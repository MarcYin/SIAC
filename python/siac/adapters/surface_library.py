"""Surface-reflectance libraries feeding the surface-driven prior.

The surface-driven retrieval predicts a scene's visible surface reflectance from
a library of clear-sky realizations of the same ground. Where that library comes
from is separable from what the predictor does with it, and the two have very
different costs and failure modes:

* building it live from a composite service is the pipeline's slowest step, and
* the radiative transfer used to atmospherically correct it decides whether the
  library is usable at all — a library corrected in one RT space and solved
  against in another injects a systematic surface offset that the solver turns
  into biased AOT.

This module makes the source pluggable behind :class:`SurfaceLibrary` and makes
the second point explicit: a library declares the :class:`~siac.domain.rt_space.RTSpace`
it was corrected in, which the pipeline checks against the solver's RT model.

A realization is deliberately the same payload shape the live composite path
already produced, so the validated reduction and prediction downstream are
untouched by the choice of source.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol

import numpy as np

from siac.adapters.bestpixel import _BESTPIXEL_DN_SCALE, _QUALITY_NODATA

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.domain.rt_space import RTSpace
    from siac.runtime.models import ObservationBundle

logger = logging.getLogger(__name__)

__all__ = [
    "PreparedSurfaceLibrary",
    "SurfaceLibrary",
    "SurfaceLibraryRealization",
    "realization_to_period",
]

#: Band order the visible predictor expects from a library realization.
PREDICTOR_BAND_ORDER = ("coastal", "blue", "green", "red", "nir", "swir16", "swir22")

#: Aliases seen in prepared libraries for the predictor's canonical band names.
_BAND_ALIASES = {
    "nir08": "nir",
    "nir8a": "nir",
    "B01": "coastal",
    "B02": "blue",
    "B03": "green",
    "B04": "red",
    "B8A": "nir",
    "B11": "swir16",
    "B12": "swir22",
}


def canonical_band_name(name: str) -> str:
    """Map a library's band name onto the predictor's canonical name."""

    return _BAND_ALIASES.get(str(name), str(name))


@dataclass(frozen=True)
class SurfaceLibraryRealization:
    """One clear-sky realization of the surface over the AOI.

    ``reflectance`` is ``(band, y, x)`` in reflectance units (not DN), aligned to
    ``band_names``. ``transform`` is the six-element affine
    ``(x_step, 0, x_origin, 0, y_step, y_origin)`` of the realization's own grid.
    """

    reflectance: np.ndarray
    band_names: tuple[str, ...]
    crs: str
    transform: tuple[float, float, float, float, float, float]
    label: str | None = None
    uncertainty: np.ndarray | None = None


def realization_to_period(realization: SurfaceLibraryRealization) -> dict[str, Any]:
    """Render a realization as the composite payload the prior reduction reads.

    Reflectance is converted back to the composite DN convention so the shared
    downstream path (reprojection, realization reduction, visible prediction)
    treats prepared and live libraries identically.
    """

    reflectance = np.asarray(realization.reflectance, dtype=np.float32)
    if reflectance.ndim != 3:
        raise ValueError(
            "Surface library realization reflectance must be (band, y, x); "
            f"got shape {reflectance.shape}."
        )
    _, height, width = reflectance.shape

    # Nodata is carried by the quality plane; mark pixels missing in every band.
    finite = np.isfinite(reflectance).any(axis=0)
    quality = np.where(finite, 0, _QUALITY_NODATA).astype(np.uint16)

    dn = np.where(np.isfinite(reflectance), reflectance / _BESTPIXEL_DN_SCALE, 0.0)
    bands = {
        str(name): dn[index].astype(np.float32) for index, name in enumerate(realization.band_names)
    }
    period: dict[str, Any] = {
        "grid": {
            "transform": list(realization.transform),
            "width": int(width),
            "height": int(height),
            "crs": realization.crs,
        },
        "bands": bands,
        "quality": quality,
        "reflectance_scale": _BESTPIXEL_DN_SCALE,
    }
    if realization.uncertainty is not None:
        uncertainty = np.asarray(realization.uncertainty, dtype=np.float32)
        period["boa_unc"] = {
            str(name): uncertainty[index] for index, name in enumerate(realization.band_names)
        }
    return period


class SurfaceLibrary(Protocol):
    """A source of surface-reflectance realizations for an observation."""

    @property
    def rt_space(self) -> RTSpace | None:
        """RT space the library's reflectance was corrected in, if managed."""

    def realizations(
        self,
        observation: ObservationBundle,
        resolution: float,
        bands: Sequence[str],
    ) -> list[SurfaceLibraryRealization]:
        """Return the library's realizations covering ``observation``."""


class PreparedSurfaceLibrary:
    """Library read from a prepared per-scene store on disk.

    The store holds one ``.npz`` per scene, keyed by the scene identifier, with
    a ``comp`` array of ``(realization, band, y, x)`` reflectance plus the grid
    (``epsg``, ``transform``) and optional ``realizations`` labels. An optional
    sidecar ``library.json`` records the RT space the library was corrected in;
    it can also be supplied directly via ``rt_space``.

    Preparing the library offline decouples the expensive acquisition and
    atmospheric correction from the retrieval, and lets the correction run in
    the same RT model the solver uses.
    """

    def __init__(
        self,
        root: str | Path,
        *,
        band_names: Sequence[str] | None = None,
        rt_space: RTSpace | None = None,
        scene_key: str | None = None,
    ) -> None:
        self._root = Path(root).expanduser()
        self._band_names = tuple(band_names) if band_names is not None else None
        self._scene_key = scene_key
        self._rt_space = rt_space if rt_space is not None else self._read_sidecar_rt_space()

    def _read_sidecar_rt_space(self) -> RTSpace | None:
        sidecar = self._root / "library.json"
        if not sidecar.is_file():
            return None
        try:
            payload = json.loads(sidecar.read_text())
        except (OSError, ValueError):
            logger.warning("Unreadable surface-library sidecar %s; RT space unknown.", sidecar)
            return None
        space = payload.get("rt_space")
        if not isinstance(space, dict):
            return None
        backend, aerosol = space.get("backend"), space.get("aerosol")
        if not backend or not aerosol:
            return None
        from siac.domain.rt_space import RTSpace as _RTSpace

        return _RTSpace(backend=str(backend), aerosol=str(aerosol))

    @property
    def rt_space(self) -> RTSpace | None:
        return self._rt_space

    def _resolve_path(self, observation: ObservationBundle) -> Path:
        key = self._scene_key or self._observation_key(observation)
        candidate = self._root / f"{key}.npz"
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(
            f"Prepared surface library has no entry for scene {key!r} under {self._root}. "
            "Build the library for this scene, or point "
            "surface_prior.prepared_library_path at the correct store."
        )

    @staticmethod
    def _observation_key(observation: ObservationBundle) -> str:
        metadata = observation.metadata or {}
        for field in ("scene_key", "matchup_id", "product_id", "scene_id"):
            value = metadata.get(field)
            if value:
                return str(value)
        raise ValueError(
            "Prepared surface library needs a scene key: set observation metadata "
            "'scene_key' (or 'product_id'), or pass scene_key explicitly."
        )

    def realizations(
        self,
        observation: ObservationBundle,
        resolution: float,  # noqa: ARG002 - library is prepared at a fixed grid
        bands: Sequence[str],  # noqa: ARG002 - library declares its own bands
    ) -> list[SurfaceLibraryRealization]:
        path = self._resolve_path(observation)
        with np.load(path, allow_pickle=False) as payload:
            comp = np.asarray(payload["comp"], dtype=np.float32)
            transform = tuple(float(v) for v in np.asarray(payload["transform"]).ravel()[:6])
            crs = self._crs_from_payload(payload)
            labels = self._labels_from_payload(payload, comp.shape[0])
            stored_bands = self._band_names_from_payload(payload, comp.shape[1])

        if comp.ndim != 4:
            raise ValueError(
                f"Prepared surface library {path} must hold a (realization, band, y, x) "
                f"'comp' array; got shape {comp.shape}."
            )
        canonical = tuple(canonical_band_name(name) for name in stored_bands)
        logger.info(
            "Prepared surface library %s: %d realizations, bands %s, rt_space=%s",
            path.name,
            comp.shape[0],
            ",".join(canonical),
            self._rt_space if self._rt_space is not None else "unmanaged",
        )
        return [
            SurfaceLibraryRealization(
                reflectance=comp[index],
                band_names=canonical,
                crs=crs,
                transform=transform,  # type: ignore[arg-type]
                label=labels[index],
            )
            for index in range(comp.shape[0])
        ]

    def _crs_from_payload(self, payload: Any) -> str:
        if "crs" in payload:
            return str(np.asarray(payload["crs"]).item())
        return f"EPSG:{int(np.asarray(payload['epsg']).item())}"

    @staticmethod
    def _labels_from_payload(payload: Any, count: int) -> list[str | None]:
        if "realizations" not in payload:
            return [None] * count
        values = [str(v) for v in np.asarray(payload["realizations"]).ravel().tolist()]
        if len(values) != count:
            return [None] * count
        return list(values)

    def _band_names_from_payload(self, payload: Any, count: int) -> tuple[str, ...]:
        if self._band_names is not None:
            names = self._band_names
        elif "band_names" in payload:
            names = tuple(str(v) for v in np.asarray(payload["band_names"]).ravel().tolist())
        else:
            names = PREDICTOR_BAND_ORDER
        if len(names) != count:
            raise ValueError(
                f"Prepared surface library declares {len(names)} band names "
                f"{list(names)} but stores {count} bands. Set "
                "surface_prior.prepared_library_bands to match the store."
            )
        return tuple(names)

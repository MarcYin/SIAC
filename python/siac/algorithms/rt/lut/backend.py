"""Zarr LUT backend implementation."""

from __future__ import annotations

import logging
import os
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar, cast

import numpy as np
import xarray as xr

from siac.algorithms.rt._run_cache import resolve_run_cache_dir
from siac.algorithms.rt.lut._compact import _CompactLUTMixin
from siac.algorithms.rt.lut._spectral import _SpectralLUTMixin
from siac.algorithms.rt.lut._spectral_math import weighted_spectral_mean
from siac.algorithms.rt.lut.store import as_local_path, build_lut_store
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.domain import SensorBand

logger = logging.getLogger(__name__)
_T = TypeVar("_T")


class ZarrLUTBackend(_CompactLUTMixin, _SpectralLUTMixin):
    """
    Look-up table backend for RT coefficients using Zarr storage.

    Implements the RTModelBackend protocol using multi-dimensional
    interpolation in a pre-computed look-up table.

    Shared SIAC RT input convention:

    - geometry angles are stored upstream in radians and converted to degrees
      before LUT interpolation
    - ``aot`` is unitless aerosol optical depth at 550 nm
    - ``tcwv`` is total column water vapour in cm precipitable water
    - ``tco3`` is total column ozone in atm-cm (DU / 1000)
    - ``elevation`` is terrain altitude in km above mean sea level

    Ozone LUTs may be indexed either in atm-cm or Dobson Units. The backend
    normalizes ``tco3`` automatically when the LUT ozone axis is in DU.

    The two supported LUT representations live in path-specific mixins:
    :class:`~siac.algorithms.rt.lut._compact._CompactLUTMixin` (legacy
    xap/xbp/xcp coefficient LUTs) and
    :class:`~siac.algorithms.rt.lut._spectral._SpectralLUTMixin` (dense
    ``TOA_rho1/2`` + ``Eg_rho1/2`` spectral LUTs). This class owns loading,
    validation, dispatch, the transient-I/O retry wrapper, and the small
    helpers shared by both paths.

    Args:
        lut_path: Path to the Zarr LUT store
        interpolation_method: Interpolation method ("linear" or "nearest")
        chunk_cache_size: Size of chunk cache in bytes
        rt_setup: Effective generic RT setup used to validate or describe the
            packaged LUT preset. The current public remote LUT is a fixed
            libRadtran preset, not a generic configurable RT family.
    """

    _SPECTRAL_IO_MAX_RETRIES = 8
    _SUPPORTED_INTERPOLATION_METHODS = frozenset({"linear", "nearest"})

    def __init__(
        self,
        lut_path: str | Path,
        interpolation_method: str = "linear",
        chunk_cache_size: int = 128 * 1024 * 1024,  # 128 MB
        storage_options: dict[str, Any] | None = None,
        rt_setup: Any | None = None,
        scene_cache_enabled: bool = True,
        scene_cache_dir: str | Path | None = None,
    ):
        self.lut_path = str(lut_path)
        if interpolation_method not in self._SUPPORTED_INTERPOLATION_METHODS:
            supported = ", ".join(sorted(self._SUPPORTED_INTERPOLATION_METHODS))
            raise ValueError(
                f"Unsupported LUT interpolation_method {interpolation_method!r}; "
                f"expected one of: {supported}"
            )
        self.interpolation_method = interpolation_method
        self.chunk_cache_size = chunk_cache_size
        self.storage_options = dict(storage_options or {})
        self._rt_setup = rt_setup
        self._scene_subset_logged = False
        # Persistent on-disk cache of the materialised scene subset. The remote
        # LUT references are already cached (lut_refs/), but the chunk BYTES are
        # otherwise re-fetched over HTTP in every fresh process; persisting the
        # small grid-snapped scene subset lets a re-run of the same scene load it
        # from local disk with no network. ``None`` disables it.
        self._scene_cache_dir: Path | None = resolve_run_cache_dir(
            scene_cache_dir, subpath="lut_subsets", enabled=scene_cache_enabled
        )

        # Lazy load the LUT
        self._lut: xr.Dataset | None = None
        self._lut_coords: dict[str, np.ndarray] = {}
        self._cache_lock = threading.Lock()
        self._scene_build_lock = threading.Lock()
        self._spectral_scene_key: tuple[float, ...] | None = None
        self._spectral_scene_subset: xr.Dataset | None = None
        self._spectral_band_grid_cache: dict[
            tuple[tuple[float, ...], tuple[str, float, float]],
            tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray],
        ] = {}

    @property
    def rt_setup(self) -> Any | None:
        return self._rt_setup

    @property
    def lut(self) -> xr.Dataset:
        """Lazily load the LUT dataset."""
        if self._lut is None:
            self._load_lut()
            self._validate_lut()
        return self._lut

    def _validate_lut(self) -> None:
        """Check that the loaded LUT contains a usable RT representation."""
        if self._lut is None:
            return
        has_compact = all(v in self._lut.data_vars for v in self._COEFFICIENT_VARS)
        has_spectral = all(v in self._lut.data_vars for v in self._SPECTRAL_VARS)
        if not has_compact and not has_spectral:
            available = sorted(self._lut.data_vars)
            raise ValueError(
                "LUT does not contain a supported RT representation. "
                f"Expected compact vars {self._COEFFICIENT_VARS} or spectral vars "
                f"{self._SPECTRAL_VARS}; found: {available}"
            )

    def _load_lut(self) -> None:
        """Load the Zarr LUT dataset."""
        local_path = as_local_path(self.lut_path)
        if local_path is not None and not local_path.exists():
            raise FileNotFoundError(f"LUT not found at {local_path}")

        logger.info(f"Loading LUT from {self.lut_path}")

        # Open with dask for lazy loading. Store may be local path, fsspec mapper, or zip mapper.
        store = build_lut_store(self.lut_path, self.storage_options)
        if self._store_contains_key(store, "zarr.json"):
            _prev = os.environ.get("ZARR_V3_EXPERIMENTAL_API")
            os.environ.setdefault("ZARR_V3_EXPERIMENTAL_API", "1")
            try:
                self._lut = xr.open_zarr(
                    store=store,
                    consolidated=False,
                    zarr_format=3,
                )
            except (KeyError, ValueError, RuntimeError, OSError) as e:
                # Narrowed from ``except Exception`` (REVIEW.md §2.1):
                # zarr raises KeyError for missing v3 entries, ValueError
                # for malformed metadata, OSError for I/O issues, and
                # xarray wraps some of these in RuntimeError. Other
                # classes (MemoryError, KeyboardInterrupt) propagate.
                logger.warning(
                    "Failed to open LUT as zarr v3 (%s); retrying with zarr v2 paths",
                    e,
                )
            finally:
                # Restore env to avoid leaking global state.
                if _prev is None:
                    os.environ.pop("ZARR_V3_EXPERIMENTAL_API", None)
                else:
                    os.environ["ZARR_V3_EXPERIMENTAL_API"] = _prev

        if self._lut is None:
            try:
                self._lut = xr.open_zarr(store=store, consolidated=True)
            except (KeyError, ValueError, RuntimeError, OSError) as e:
                # Same narrowing rationale as the zarr-v3 path above.
                logger.warning(
                    "Failed to open LUT with consolidated metadata (%s); retrying with consolidated=False",
                    e,
                )
                self._lut = xr.open_zarr(store=store, consolidated=False)

        # Cache coordinate arrays for fast interpolation
        for dim, coord in self._lut.coords.items():
            self._lut_coords[dim] = coord.values

        logger.info(
            "LUT index loaded (available axis sizes): SZA=%d, VZA=%d, RAA=%d, AOT=%d, TCWV=%d",
            len(self._lut_coords.get("sza", [])),
            len(self._lut_coords.get("vza", [])),
            len(self._lut_coords.get("raa", [])),
            len(self._lut_coords.get("aot", [])),
            len(self._lut_coords.get("tcwv", [])),
        )

    @staticmethod
    def _store_contains_key(store: Any, key: str) -> bool:
        """Best-effort key existence check for mapping/path style stores."""
        if isinstance(store, str):
            return (Path(store) / key).exists()
        try:
            return key in store
        except (TypeError, KeyError, OSError, RuntimeError, ValueError):
            # Defensive best-effort check (REVIEW.md §2.1): a custom Mapping
            # with no ``__contains__`` raises TypeError, remote stores can
            # raise OSError on a transient connection drop, fsspec wrappers
            # may surface RuntimeError, and a misconfigured zarr v3 store
            # raises ValueError on probe. ``MemoryError``/``KeyboardInterrupt``
            # still propagate.
            return False

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        """
        Compute radiative transfer coefficients via LUT interpolation.

        Args:
            geometry: Viewing geometry (sza, vza, raa in radians)
            atmo_state: Atmospheric state in the shared SIAC RT units
                (unitless aot, tcwv in cm, tco3 in atm-cm, elevation in km)
            band: Sensor band specification
            compute_jacobian: Whether to compute d(coeff)/d(aot,tcwv)

        Returns:
            RTCoefficients with xap, xbp, xcp and optional Jacobians.
        """
        template = self._require_matching_grid_shapes(geometry, atmo_state)

        sza_deg, vza_deg, raa_deg = self._geometry_to_lut_inputs(geometry)

        # Get atmospheric parameters
        aot = atmo_state.aot.values
        tcwv = atmo_state.tcwv.values
        elevation = atmo_state.elevation.values

        # Get band wavelength for LUT selection
        wavelength = band.center_wavelength

        if self._supports_coefficient_lut():
            # Interpolate compact LUT variables, then derive RT coefficients.
            path_ref, trans_down, trans_up, sph_alb = self._interpolate_lut(
                sza_deg,
                vza_deg,
                raa_deg,
                aot,
                tcwv,
                atmo_state.tco3.values,
                elevation,
                wavelength,
            )

            # Legacy compact LUTs may not carry altitude as an explicit dimension.
            if not self._coefficient_lut_has_altitude_axis():
                path_ref, trans_down, trans_up, sph_alb = self._apply_elevation_correction(
                    path_ref,
                    trans_down,
                    trans_up,
                    sph_alb,
                    elevation,
                )

            # Convert to RT coefficients
            # The correction equation is:
            #   y = xap * toa - xbp
            #   boa = y / (1 + xcp * y)
            #
            # From RT parameters:
            #   xap = 1 / (trans_down * trans_up)
            #   xbp = path_ref / (trans_down * trans_up)
            #   xcp = sph_alb
            trans_total = trans_down * trans_up
            trans_total = np.maximum(trans_total, 1e-10)  # Avoid division by zero

            xap = 1.0 / trans_total
            xbp = path_ref / trans_total
            xcp = sph_alb
        elif self._supports_spectral_lut():
            xap, xbp, xcp = self._compute_spectral_with_retry(
                sza_deg=sza_deg,
                vza_deg=vza_deg,
                raa_deg=raa_deg,
                aot=aot,
                tcwv=tcwv,
                tco3=atmo_state.tco3.values,
                elevation=elevation,
                band=band,
            )
        else:
            raise ValueError(
                "LUT does not contain a supported RT representation. "
                "Expected compact coefficient variables or dense spectral LUT variables."
            )

        # Create DataArrays
        xap_da = xr.DataArray(xap, dims=template.dims, coords=template.coords)
        xbp_da = xr.DataArray(xbp, dims=template.dims, coords=template.coords)
        xcp_da = xr.DataArray(xcp, dims=template.dims, coords=template.coords)

        # Compute Jacobians via finite differences if requested
        d_xap = None
        d_xbp = None
        d_xcp = None

        if compute_jacobian:
            d_xap, d_xbp, d_xcp = self._compute_jacobian_numerical(
                geometry, atmo_state, band, xap, xbp, xcp
            )

        return RTCoefficients(
            xap=xap_da,
            xbp=xbp_da,
            xcp=xcp_da,
            d_xap=d_xap,
            d_xbp=d_xbp,
            d_xcp=d_xcp,
        )

    @staticmethod
    def _axis_bounds(axis: np.ndarray) -> tuple[float, float]:
        axis_values = np.asarray(axis, dtype=np.float32)
        return float(np.nanmin(axis_values)), float(np.nanmax(axis_values))

    @classmethod
    def _sanitize_point_values(
        cls,
        values: np.ndarray,
        axis: np.ndarray,
    ) -> np.ndarray:
        """Clamp finite values to ``[axis_min, axis_max]`` and preserve NaN.

        Previously NaN values were replaced with the axis midpoint, which
        caused missing inputs (e.g. cloudy pixels with no atmospheric prior)
        to silently produce fabricated interpolations at the LUT's middle.
        The interpolator (``RegularGridInterpolator(... fill_value=NaN)``)
        is configured to return NaN for NaN inputs, so we now preserve NaN
        end-to-end. Finite values are still clipped to the LUT envelope —
        the LUT cannot extend beyond its sampled range, so silent clipping
        of slightly-out-of-bounds values from observation noise is
        intentional. (REVIEW.md §1.2 #4)
        """
        axis_min, axis_max = cls._axis_bounds(axis)
        arr = np.asarray(values, dtype=np.float32)
        finite_mask = np.isfinite(arr)
        clipped = np.clip(arr, axis_min, axis_max).astype(np.float32, copy=False)
        # Preserve NaN positions: where the input was non-finite (NaN/+inf/-inf)
        # we restore the original sentinel rather than the clipped boundary.
        result = np.where(finite_mask, clipped, arr).astype(np.float32, copy=False)
        return cast("np.ndarray", result)

    @staticmethod
    def _require_finite_values(
        values: np.ndarray,
        *,
        name: str,
    ) -> np.ndarray:
        arr = np.asarray(values, dtype=np.float32)
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} must contain only finite values for LUT interpolation")
        return cast("np.ndarray", arr)

    @staticmethod
    def _grid_arrays(
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> dict[str, xr.DataArray]:
        """Return the retrieval-grid fields that must align for LUT evaluation."""
        return {
            "geometry.sza": geometry.sza,
            "geometry.vza": geometry.vza,
            "geometry.raa": geometry.raa,
            "atmo_state.aot": atmo_state.aot,
            "atmo_state.tcwv": atmo_state.tcwv,
            "atmo_state.tco3": atmo_state.tco3,
            "atmo_state.elevation": atmo_state.elevation,
        }

    @staticmethod
    def _grid_shape_mismatches(
        arrays: dict[str, xr.DataArray],
        *,
        expected_shape: tuple[int, ...],
    ) -> list[str]:
        """Describe shape mismatches against the geometry template grid."""
        return [
            f"{name}={tuple(array.shape)}"
            for name, array in arrays.items()
            if tuple(array.shape) != expected_shape
        ]

    @staticmethod
    def _grid_alignment_mismatches(
        template: xr.DataArray,
        arrays: dict[str, xr.DataArray],
    ) -> tuple[list[str], list[str]]:
        """Describe dim-order and coordinate mismatches against the template grid."""
        dim_mismatches = [
            f"{name}.dims={array.dims}"
            for name, array in arrays.items()
            if array.dims != template.dims
        ]
        coord_mismatches: list[str] = []
        expected_shape = tuple(template.shape)
        for name, array in arrays.items():
            if tuple(array.shape) != expected_shape or array.dims != template.dims:
                continue
            for axis in template.dims:
                template_has_coord = axis in template.coords
                array_has_coord = axis in array.coords
                if template_has_coord != array_has_coord:
                    coord_mismatches.append(f"{name}.{axis}=missing")
                    continue
                if not template_has_coord:
                    continue
                template_values = np.asarray(template.coords[axis].values)
                array_values = np.asarray(array.coords[axis].values)
                if not np.array_equal(template_values, array_values, equal_nan=True):
                    coord_mismatches.append(f"{name}.{axis}")
        return dim_mismatches, coord_mismatches

    @staticmethod
    def _require_matching_grid_shapes(
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> xr.DataArray:
        """Validate that geometry and atmospheric fields share the same retrieval grid."""
        template = geometry.sza
        expected_shape = tuple(template.shape)
        arrays = ZarrLUTBackend._grid_arrays(geometry, atmo_state)
        mismatches = ZarrLUTBackend._grid_shape_mismatches(arrays, expected_shape=expected_shape)
        dim_mismatches, coord_mismatches = ZarrLUTBackend._grid_alignment_mismatches(
            template, arrays
        )
        if mismatches:
            details = ", ".join(mismatches)
            raise ValueError(
                "Geometry and atmospheric fields must share the same grid shape; "
                f"expected {expected_shape} from geometry.sza but got {details}"
            )
        if dim_mismatches or coord_mismatches:
            details = ", ".join(dim_mismatches + coord_mismatches)
            raise ValueError(
                "Geometry and atmospheric fields must share the same grid shape and coordinates; "
                f"expected dims {template.dims} from geometry.sza but got {details}"
            )
        return template

    @staticmethod
    def _geometry_to_lut_inputs(
        geometry: GeometryAngles,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Convert radian geometry into LUT-ready degrees.

        ``GeometryAngles`` stores angles in radians. LUT interpolation uses
        solar/view zenith and relative azimuth in degrees, with relative azimuth
        normalized to ``[0, 180]``.
        """
        sza_deg = np.rad2deg(geometry.sza.values)
        vza_deg = np.rad2deg(geometry.vza.values)
        raa_deg = np.abs(np.rad2deg(geometry.raa.values))
        return sza_deg, vza_deg, np.where(raa_deg > 180, 360 - raa_deg, raa_deg)

    @staticmethod
    def _interpolate_variable(
        var: xr.DataArray,
        coords: dict[str, xr.DataArray],
        method: str,
    ) -> xr.DataArray:
        """Interpolate one LUT variable using only coordinates it actually depends on."""
        applicable = {name: coord for name, coord in coords.items() if name in var.dims}
        if not applicable:
            return var
        if method == "linear":
            return var.interp(**applicable)
        return var.sel(**applicable, method="nearest")

    @staticmethod
    def _interpolate_variable_fast(
        var: xr.DataArray,
        coords: dict[str, xr.DataArray],
        method: str,
    ) -> np.ndarray:
        """Fast point interpolation using scipy RegularGridInterpolator.

        Replaces xr.DataArray.interp() which has heavy per-call overhead from
        coordinate alignment and broadcasting.  For N query points on a regular
        grid this is typically 50-200x faster.
        """
        from scipy.interpolate import RegularGridInterpolator

        applicable_names = [name for name in var.dims if name in coords]
        if not applicable_names:
            return np.asarray(var.values, dtype=np.float32)  # type: ignore[no-any-return]

        grid_axes = tuple(
            np.asarray(var.coords[name].values, dtype=np.float64) for name in applicable_names
        )
        values = np.asarray(var.values, dtype=np.float64)

        if method == "linear":
            interp = RegularGridInterpolator(
                grid_axes,
                values,
                method="linear",
                bounds_error=False,
                fill_value=np.nan,
            )
        else:
            interp = RegularGridInterpolator(
                grid_axes,
                values,
                method="nearest",
                bounds_error=False,
                fill_value=np.nan,
            )

        query_points = np.column_stack(
            [np.asarray(coords[name].values, dtype=np.float64) for name in applicable_names]
        )
        result: np.ndarray = interp(query_points).astype(np.float32)
        return result

    @staticmethod
    def _weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
        """Weighted mean over wavelength with coordinate-aware integration.

        Kept on the facade (not the spectral mixin) so the underlying
        ``weighted_spectral_mean`` keeps resolving through this module's
        namespace — tests monkeypatch it here.
        """
        return weighted_spectral_mean(data, weights)

    def _run_with_transient_lut_io_retry(
        self,
        fn: Callable[[], _T],
        *,
        operation: str,
    ) -> _T:
        """Retry transient remote LUT I/O for preload/read paths."""
        max_attempts = max(1, int(self._SPECTRAL_IO_MAX_RETRIES))
        for attempt in range(1, max_attempts + 1):
            try:
                return fn()
            except Exception as exc:
                # Intentionally broad: ``_is_transient_lut_io_error`` does
                # its own classification by walking ``__cause__`` /
                # ``__context__`` and matching tokens. Anything that's not
                # transient is re-raised unchanged. (REVIEW.md §2.1)
                if attempt >= max_attempts or not self._is_transient_lut_io_error(exc):
                    raise
                logger.warning(
                    "Transient LUT %s failure (%s); reloading and retrying (%d/%d).",
                    operation,
                    type(exc).__name__,
                    attempt,
                    max_attempts,
                )
                time.sleep(min(float(attempt), 5.0))
                self._reload_lut()

        raise RuntimeError(f"Failed to {operation} LUT data after retries")

    def _reload_lut(self) -> None:
        """Drop cached handles and reopen LUT metadata/chunks lazily."""
        self._lut = None
        self._lut_coords = {}
        self._scene_subset_logged = False
        with self._cache_lock:
            self._spectral_scene_key = None
            self._spectral_scene_subset = None
            self._spectral_band_grid_cache.clear()
        _ = self.lut

    @staticmethod
    def _is_transient_lut_io_error(exc: Exception) -> bool:
        """Return True for transient remote I/O errors that should be retried."""
        probes = []
        current: Exception | None = exc
        for _ in range(10):
            if current is None:
                break
            probes.append(current)
            nxt = current.__cause__ if isinstance(current.__cause__, Exception) else None
            if nxt is None and isinstance(current.__context__, Exception):
                nxt = current.__context__
            current = nxt

        tokens = (
            "referencenotreachable",
            "server disconnected",
            "serverdisconnectederror",
            "connection reset",
            "connection aborted",
            "timed out",
            "timeout",
            "temporarily unavailable",
        )
        for err in probes:
            name = type(err).__name__.lower()
            text = str(err).lower()
            if any(tok in name or tok in text for tok in tokens):
                return True
        return False

    def _supports_coefficient_lut(self) -> bool:
        """Check whether the loaded LUT exposes compact coefficient variables."""
        return all(name in self.lut.data_vars for name in self._COEFFICIENT_VARS)

    def _supports_spectral_lut(self) -> bool:
        """Check whether the loaded LUT exposes dense spectral radiative terms."""
        return all(name in self.lut.data_vars for name in self._SPECTRAL_VARS)

    def _compute_jacobian_numerical(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        xap: np.ndarray,
        xbp: np.ndarray,
        xcp: np.ndarray,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """
        Compute Jacobians via numerical finite differences.

        Args:
            geometry: Viewing geometry
            atmo_state: Atmospheric state
            band: Sensor band
            xap, xbp, xcp: Current coefficient values

        Returns:
            Tuple of Jacobian DataArrays for (d_xap, d_xbp, d_xcp)
        """
        # Forward-difference perturbation sizes — see siac.constants.
        from siac.constants import (
            DEFAULT_JACOBIAN_DELTA_AOT as delta_aot,
        )
        from siac.constants import (
            DEFAULT_JACOBIAN_DELTA_TCWV as delta_tcwv,
        )

        template = geometry.sza

        # Perturb AOT
        aot_plus = atmo_state.with_updated_aot_tcwv(
            aot=atmo_state.aot + delta_aot,
            tcwv=atmo_state.tcwv,
        )

        coeffs_aot_plus = self.compute_coefficients(
            geometry, aot_plus, band, compute_jacobian=False
        )

        d_xap_aot = (coeffs_aot_plus.xap.values - xap) / delta_aot
        d_xbp_aot = (coeffs_aot_plus.xbp.values - xbp) / delta_aot
        d_xcp_aot = (coeffs_aot_plus.xcp.values - xcp) / delta_aot

        # Perturb TCWV
        tcwv_plus = atmo_state.with_updated_aot_tcwv(
            aot=atmo_state.aot,
            tcwv=atmo_state.tcwv + delta_tcwv,
        )

        coeffs_tcwv_plus = self.compute_coefficients(
            geometry, tcwv_plus, band, compute_jacobian=False
        )

        d_xap_tcwv = (coeffs_tcwv_plus.xap.values - xap) / delta_tcwv
        d_xbp_tcwv = (coeffs_tcwv_plus.xbp.values - xbp) / delta_tcwv
        d_xcp_tcwv = (coeffs_tcwv_plus.xcp.values - xcp) / delta_tcwv

        # Create DataArrays with param dimension
        d_xap = xr.concat(
            [
                xr.DataArray(d_xap_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xap_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        d_xbp = xr.concat(
            [
                xr.DataArray(d_xbp_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xbp_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        d_xcp = xr.concat(
            [
                xr.DataArray(d_xcp_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xcp_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        return d_xap, d_xbp, d_xcp

    def supports_jacobian(self) -> bool:
        """Check if this backend supports Jacobian computation."""
        return False

    @property
    def backend_name(self) -> str:
        """Name of the RT backend."""
        return "lut"

    def is_available_for_sensor(self, _sensor_id: str, _satellite_id: str) -> bool:
        """Check if this backend has data for the specified sensor."""
        # LUT backend works for any sensor with known wavelengths
        return True

    def get_available_wavelengths(self) -> np.ndarray:
        """Get wavelengths available in the LUT."""
        wavelengths = self._lut_coords.get("wavelength")
        if wavelengths is None:
            return np.zeros(0, dtype=np.float32)
        result: np.ndarray = np.array(wavelengths, dtype=np.float32)
        return result

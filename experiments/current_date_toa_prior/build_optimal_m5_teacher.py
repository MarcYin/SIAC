#!/usr/bin/env python3
"""Build T1(tau*) labels from M5 or an explicitly provided training AOD.

The historical L1C archive remains the seasonal regression dictionary.  The
label is the B02/B03/B04 prior predicted from current-date B8A/B11/B12 after
correction at M5's accepted AOD node.  The optional provided-AOD mode is a
separate train/development-only experiment: callers supply a scene AOD and
uncertainty explicitly, and the builder never looks up AERONET itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import time
import traceback
import warnings
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
from experiments.current_date_toa_prior.build_expanded_l1c_teacher import (
    DEFAULT_CAMS,
    DEFAULT_MIE_CACHE,
    DEFAULT_MIE_MODULE,
    DEFAULT_MODULE,
    DEFAULT_RUN_CACHE,
    _backend,
    _platform,
    _scene_time,
)
from experiments.current_date_toa_prior.snow_teacher_support import (
    classify_profile,
    profile_library_support,
)
from experiments.current_date_toa_prior.surface_model_v2 import (
    _robust_median as _robust_temporal_median,
)
from rasterio.enums import Resampling

from siac.adapters.atmo.cams import CAMSProvider
from siac.adapters.atmo.fusion import FusedAODProvider
from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider
from siac.adapters.data.water_mask import load_water_mask_subset
from siac.adapters.live_l1c_library import read_mosaic_index
from siac.adapters.rsrf import load_sensor_config_with_rsrf
from siac.algorithms.cloud.providers.omnicloudmask import OmniCloudMaskProvider
from siac.algorithms.grid.assembler import _dilate_mask, _resample_external_exclusion_mask
from siac.algorithms.solver import SurfaceDrivenSolver
from siac.algorithms.solver.surface_driven import _species_candidate_rt_models
from siac.algorithms.surface.seasonal_predictor import (
    _correct_anchor_reflectance,
    predict_visible_from_tau_payload,
    seasonal_extra_tree_prior,
)
from siac.config.algorithms import SolverAlgorithmConfig
from siac.domain.rt_space import RTSpace
from siac.runtime import AtmosphericState, GeometryAngles, ObservationBundle, SurfacePrior

if TYPE_CHECKING:
    from datetime import datetime

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
CAMPAIGN = ROOT / "current_date_toa_surface_allera_stac_20260809"
SCHEMA = "siac_optimal_m5_teacher_v8"
SNOW_POLICY_SCHEMA = "siac_optimal_m5_teacher_snow_v8"
AERONET_CONDITIONED_SCHEMA = "siac_aeronet_conditioned_surface_teacher_v1"
VISIBLE = ("B02", "B03", "B04")
SOLVER_VISIBLE = ("B01", *VISIBLE)
ANCHORS = ("B8A", "B11", "B12")
#: Every Sentinel-2 land band a seasonal library can supervise. B09 is excluded:
#: the historical correction that built it used B09-derived water vapour, so it
#: is not independent of its own label.
S2_LAND_BANDS = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12")
#: Seven-band dictionaries label their planes by role; wider ones use Sentinel-2
#: band names directly. Mirrors ``seasonal_predictor._COMPOSITE_BAND_ALIASES``.
DICTIONARY_BAND_ALIASES = {
    "coastal": "B01",
    "blue": "B02",
    "green": "B03",
    "red": "B04",
    "nir": "B8A",
    "swir16": "B11",
    "swir22": "B12",
}
logger = logging.getLogger(__name__)


#: Band order of the reduced archive's ``surface`` planes
#: (``surface_model_v2.ALL_SURFACE_BANDS``). It coincides with the seven-band
#: dictionary's order, but a wider dictionary breaks that coincidence, so the
#: two are resolved separately.
REDUCED_BAND_COLUMNS = {
    "B01": 0,
    "B02": 1,
    "B03": 2,
    "B04": 3,
    "B8A": 4,
    "B11": 5,
    "B12": 6,
}


def _dictionary_columns(band_names):
    """Canonical Sentinel-2 band name -> column for a seasonal dictionary."""
    columns = {}
    for index, name in enumerate(str(value) for value in band_names):
        columns[DICTIONARY_BAND_ALIASES.get(name, name)] = index
    return columns


def _target_base_planes(t0, t0_unc, comp, target_bands, dictionary_columns):
    """Base prior planes for every target band, on the 20 m grid.

    The predictor overwrites both the reflectance and the uncertainty of every
    band it targets, so these planes are carriers rather than estimates. Bands
    the reduced archive holds keep its committed robust median; bands only the
    wider dictionary carries take that dictionary's own temporal median and
    MAD, which is the same estimator family.
    """
    boa, unc = [], []
    for name in target_bands:
        if name in REDUCED_BAND_COLUMNS and REDUCED_BAND_COLUMNS[name] < t0.shape[-1]:
            column = REDUCED_BAND_COLUMNS[name]
            boa.append(np.asarray(t0[..., column], dtype=np.float32))
            unc.append(np.asarray(t0_unc[..., column], dtype=np.float32))
            continue
        if name not in dictionary_columns:
            raise ValueError(f"seasonal dictionary cannot supervise target band {name}")
        planes = np.asarray(comp[:, dictionary_columns[name]], dtype=np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            median = np.nanmedian(planes, axis=0)
            spread = 1.4826 * np.nanmedian(np.abs(planes - median[np.newaxis]), axis=0)
        boa.append(median.astype(np.float32))
        unc.append(spread.astype(np.float32))
    return np.stack(boa, axis=0), np.stack(unc, axis=0)


def _parse_surface_target_bands(value):
    """Bands the teacher predicts and writes, independent of the AOD solve bands.

    Widening the surface output does not widen the retrieval: M5 keeps solving
    the aerosol state on its own ``--solver-solve-bands``, and the extra bands
    are read off the same fitted trees at the same solved AOD.
    """
    bands = tuple(name.strip() for name in str(value).split(",") if name.strip())
    if not bands:
        raise ValueError("--surface-target-bands cannot be empty")
    unknown = tuple(name for name in bands if name not in S2_LAND_BANDS)
    if unknown:
        raise ValueError(f"--surface-target-bands has non-land bands {unknown}")
    if len(set(bands)) != len(bands):
        raise ValueError(f"--surface-target-bands repeats a band: {bands}")
    return bands


class TeacherSceneIneligible(RuntimeError):
    """Documented label absence, distinct from an IO/solver implementation error."""
    def __init__(self, reason, diagnostics):
        super().__init__(reason)
        self.reason = reason
        self.diagnostics = diagnostics
# OmniCloudMask native labels: 0 clear, 1 thick cloud, 2 thin cloud,
# 3 shadow.  Preserve them in the archive so that the retrieval policy is
# auditable.  Only thick cloud and shadow are excluded; thin cloud is retained
# for recovery from the NIR/SWIR anchors.
OCM_RAW_IDENTITY = {0: [0], 1: [1], 2: [2], 3: [3]}
OCM_RETRIEVAL_USABLE = (0, 2)
OCM_RETRIEVAL_CLEAR_ONLY = (0,)
AOD_PRIOR_SOURCE_MISSING = np.uint8(0)
AOD_PRIOR_SOURCE_MAIAC = np.uint8(1)
AOD_PRIOR_SOURCE_CAMS = np.uint8(2)
AOD_PRIOR_SOURCE_BOTH = np.uint8(3)

COMMITTED_C0_CONTRACT = "committed-c0-exact"
CCI25_EXPERIMENT_CONTRACT = "cci25-experimental"
TEACHER_CONTRACTS: dict[str, dict[str, Any]] = {
    COMMITTED_C0_CONTRACT: {
        "aerosol_species": "cci_climatology_exact",
        "aot_axis": "acixthree",
        "aot_axis_nodes": 68,
        "aot_max": 4.0,
        "role": "promotion baseline reproducing frozen C0 (126/149 within EE)",
    },
    CCI25_EXPERIMENT_CONTRACT: {
        "aerosol_species": "cci_climatology_25pct",
        "aot_axis": "acixthree",
        "aot_axis_nodes": 68,
        "aot_max": 4.0,
        "role": "non-promoted CCI25 ablation; must pass the locked parity gate",
    },
}
DEFAULT_EXACT_MIE_CACHE = ROOT / "mie_cache_exact_harmonized_20260804"
DEFAULT_EXACT_RUN_CACHE = ROOT / "rt_run_cache_expanded_l1c_teacher_exact/rt6s/run_cache"


class _FixedAtmosphericProvider:
    """Expose an already-loaded state through the provider protocol."""

    def __init__(self, state: AtmosphericState, source_name: str) -> None:
        self._state = state
        self.source_name = source_name

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        del bounds, crs, obs_time, resolution
        return self._state


@dataclass(frozen=True)
class _FusedAODPrior:
    state: AtmosphericState
    maiac_aot: xr.DataArray
    cams_aot: xr.DataArray
    winner: xr.DataArray


@dataclass(frozen=True)
class _DeepBluePriorResult:
    """Result and provenance for an optional scalar MODIS Deep Blue fusion."""

    state: AtmosphericState
    available: bool
    aod550: float | None
    uncertainty: float | None
    metadata: dict[str, Any]


def _scalar(value: Any) -> Any:
    return np.asarray(value).item()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _grid(values: np.ndarray, transform: np.ndarray, crs: str, *, name: str) -> xr.DataArray:
    height, width = values.shape[-2:]
    a, b, c, d, e, f = (float(value) for value in transform[:6])
    if b != 0.0 or d != 0.0:
        raise ValueError("optimal-M5 teacher requires a north-up grid")
    x = c + (np.arange(width) + 0.5) * a
    y = f + (np.arange(height) + 0.5) * e
    dims = ("y", "x") if values.ndim == 2 else ("band", "y", "x")
    coords: dict[str, Any] = {"y": y, "x": x}
    return xr.DataArray(values, dims=dims, coords=coords, name=name).rio.write_crs(crs)


def _block_mean(values: np.ndarray, factor: int = 3) -> np.ndarray:
    array = np.asarray(values, dtype=np.float32)
    height, width = array.shape[-2:]
    if height % factor or width % factor:
        raise ValueError(f"shape {(height, width)} is not divisible by {factor}")
    reshaped = array.reshape(*array.shape[:-2], height // factor, factor, width // factor, factor)
    with np.errstate(invalid="ignore"):
        return np.nanmean(reshaped, axis=(-3, -1)).astype(np.float32)


def _ocm_retrieval_policy(
    raw_classes: np.ndarray,
    input_valid: np.ndarray,
    *,
    usable_classes: tuple[int, ...] = OCM_RETRIEVAL_USABLE,
) -> tuple[np.ndarray, np.ndarray]:
    """Map native OCM labels to an explicit M5 retrieval policy."""

    raw = np.asarray(raw_classes, dtype=np.uint8)
    available = np.asarray(input_valid, dtype=bool)
    if raw.shape != available.shape:
        raise ValueError(f"OCM class/input-valid shapes differ: {raw.shape} vs {available.shape}")
    unknown = np.setdiff1d(np.unique(raw[available]), np.asarray((0, 1, 2, 3)))
    if unknown.size:
        raise ValueError(f"unknown native OCM classes: {unknown.tolist()}")
    retrieval_valid = available & np.isin(raw, np.asarray(usable_classes, dtype=np.uint8))
    standardized = np.zeros(raw.shape, dtype=np.uint8)
    standardized[retrieval_valid] = np.uint8(1)
    standardized[available & (raw == 1)] = np.uint8(2)
    standardized[available & (raw == 3)] = np.uint8(3)
    return standardized, retrieval_valid


def _scalar_component_max_aod_prior(
    fused: _FusedAODPrior,
    template: xr.DataArray,
) -> AtmosphericState:
    """Approximate the locked MAXCAMS scalar backstop on the scene grid.

    The frozen 149-case harness supplied a constant AOD field equal to the
    larger of its staged scene MAIAC and CAMS values.  The expanded builder
    originally supplied a pixelwise maximum instead.  This helper makes that
    distinction testable while retaining the same source fields for audit.
    """

    component_medians: list[float] = []
    for field in (fused.maiac_aot, fused.cams_aot):
        values = np.asarray(field, dtype=np.float64)
        finite = values[np.isfinite(values)]
        if finite.size:
            component_medians.append(float(np.median(finite)))
    if not component_medians:
        raise RuntimeError("neither MAIAC nor CAMS has a finite scene AOD")
    value = np.float32(max(component_medians))
    aot = xr.full_like(template, value).astype(np.float32)
    aot_unc = xr.full_like(template, np.float32(max(0.5 * float(value), 0.05)))
    return replace(fused.state, aot=aot, aot_unc=aot_unc)


def _scene_constant(field: xr.DataArray, value: float | None = None) -> xr.DataArray:
    values = np.asarray(field, dtype=np.float64)
    finite = values[np.isfinite(values)]
    fill = (
        float(value)
        if value is not None
        else (float(np.median(finite)) if finite.size else float("nan"))
    )
    return xr.full_like(field, np.float32(fill)).astype(np.float32)


def _apply_modis_deep_blue_prior(
    base: AtmosphericState,
    path: Path,
    *,
    quality: str,
    fusion: str,
) -> _DeepBluePriorResult:
    """Fuse a staged truth-free Deep Blue scalar with the current AOD field.

    Missing Deep Blue is an explicit baseline fallback.  The staged file was
    selected using MODIS QA plus spatial/temporal proximity only; this function
    never reads AERONET or any validation label.
    """

    payload = json.loads(path.read_text())
    if payload.get("schema_version") != "siac_modis_deep_blue_prior_v1":
        raise ValueError(f"unsupported Deep Blue prior schema in {path}")
    if quality not in {"best", "qa2", "qa1"}:
        raise ValueError(f"unsupported Deep Blue quality policy {quality!r}")
    selected = payload.get("selected", {}).get(quality)
    metadata = {
        "quality": quality,
        "fusion": fusion,
        "staged_path": str(path),
        "selected": selected,
    }
    if selected is None:
        return _DeepBluePriorResult(base, False, None, None, metadata)
    aod = float(selected["aod550"])
    if not np.isfinite(aod) or not 0.0 <= aod <= 5.0:
        raise ValueError(f"invalid staged Deep Blue AOD {aod!r} in {path}")
    raw_uncertainty = selected.get("uncertainty")
    uncertainty = (
        float(raw_uncertainty)
        if raw_uncertainty is not None and np.isfinite(float(raw_uncertainty))
        else max(0.1 * aod, 0.05)
    )
    uncertainty = max(uncertainty, 0.05)
    base_aod = np.asarray(base.aot, dtype=np.float32)
    base_unc = np.asarray(base.aot_unc, dtype=np.float32)
    if fusion == "replace":
        result_aod = np.full(base_aod.shape, aod, dtype=np.float32)
        result_unc = np.full(base_unc.shape, uncertainty, dtype=np.float32)
    elif fusion == "max":
        use_deep_blue = np.float32(aod) > base_aod
        result_aod = np.maximum(base_aod, np.float32(aod))
        result_unc = np.where(use_deep_blue, np.float32(uncertainty), base_unc)
    else:
        weights = {
            "mean": 0.5,
            "db_weight_0p25": 0.25,
            "db_weight_0p75": 0.75,
        }
        try:
            weight = np.float32(weights[fusion])
        except KeyError as exc:
            raise ValueError(f"unsupported Deep Blue fusion {fusion!r}") from exc
        result_aod = (np.float32(1.0) - weight) * base_aod + weight * np.float32(aod)
        result_unc = np.sqrt(
            ((np.float32(1.0) - weight) * base_unc) ** 2 + (weight * np.float32(uncertainty)) ** 2
        )
    aot = base.aot.copy(data=np.asarray(result_aod, dtype=np.float32))
    aot_unc = base.aot_unc.copy(data=np.asarray(result_unc, dtype=np.float32))
    return _DeepBluePriorResult(
        replace(base, aot=aot, aot_unc=aot_unc),
        True,
        aod,
        uncertainty,
        metadata,
    )


def _water_mask_on_grid(
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    template: xr.DataArray,
    source: str | Path | None,
    cache_dir: Path,
    buffer_pixels: int,
) -> xr.DataArray:
    """Load and remap the same landWater2020 exclusion used by committed C0."""

    native = load_water_mask_subset(
        bounds,
        crs,
        source=source,
        cache_dir=cache_dir,
    )
    if int(buffer_pixels) > 0:
        native = native.copy(data=_dilate_mask(np.asarray(native, dtype=bool), int(buffer_pixels)))
    return _resample_external_exclusion_mask(
        native,
        tuple(int(value) for value in template.shape),
        template=template,
    ).astype(bool)


def _constant_geometry(payload: Any, template: xr.DataArray) -> GeometryAngles:
    def field(name: str) -> xr.DataArray:
        return xr.full_like(template, np.float32(_scalar(payload[name])))

    return GeometryAngles.from_degrees(
        field("mean_sza_deg"),
        field("mean_saa_deg"),
        field("mean_vza_deg"),
        field("mean_vaa_deg"),
    )


def _align(value: xr.DataArray, template: xr.DataArray, *, nearest: bool = False) -> xr.DataArray:
    return value.rio.reproject_match(
        template, resampling=Resampling.nearest if nearest else Resampling.bilinear
    ).astype(np.float32)


def _aod_state_on_grid(state: AtmosphericState, template: xr.DataArray) -> AtmosphericState:
    """Put the two aerosol fields on the solver grid before source fusion."""

    return replace(
        state,
        aot=_align(state.aot, template),
        aot_unc=_align(state.aot_unc, template),
    )


def _fuse_aod_on_grid(
    maiac: AtmosphericState | None,
    cams: AtmosphericState,
    *,
    template: xr.DataArray,
    obs_time: datetime,
    combination: str = "max",
) -> _FusedAODPrior:
    """Fuse QA-MAIAC and CAMS on the M5 grid with auditable provenance.

    The production :class:`FusedAODProvider` is used here rather than a second
    implementation of the fusion rule.  Aligning both states first makes the
    archived component fields and winner mask an exact, auditable account of
    the values supplied to M5.  If MAIAC has no usable retrieval, CAMS is
    promoted exactly as it is in the production provider.
    """

    cams_grid = _aod_state_on_grid(cams, template)
    maiac_grid = _aod_state_on_grid(maiac, template) if maiac is not None else None
    if maiac_grid is None:
        fused = cams_grid
        maiac_aot = xr.full_like(cams_grid.aot, np.float32(np.nan))
    else:
        fused = FusedAODProvider(
            _FixedAtmosphericProvider(maiac_grid, "MCD19-best-quality-policy"),
            [_FixedAtmosphericProvider(cams_grid, "CAMS")],
            op=str(combination),
        ).get_prior(
            tuple(float(value) for value in template.rio.bounds()),
            str(template.rio.crs),
            obs_time,
            abs(float(template.x.values[1] - template.x.values[0])),
        )
        maiac_aot = maiac_grid.aot

    maiac_values = np.asarray(maiac_aot, dtype=np.float32)
    cams_values = np.asarray(cams_grid.aot, dtype=np.float32)
    fused_values = np.asarray(fused.aot, dtype=np.float32)
    if str(combination) == "max":
        expected = np.fmax(maiac_values, cams_values)
    elif str(combination) == "mean":
        count = np.isfinite(maiac_values).astype(np.uint8) + np.isfinite(cams_values).astype(
            np.uint8
        )
        expected = np.nansum(np.stack((maiac_values, cams_values)), axis=0)
        expected = np.divide(
            expected,
            count,
            out=np.full(expected.shape, np.nan, dtype=np.float32),
            where=count > 0,
        )
    elif str(combination) == "plateau":
        # Recomputed here independently of the provider, as the other branches
        # are -- that is the point of this audit. The knots are imported rather
        # than restated so the two implementations cannot drift apart.
        from siac.adapters.atmo.fusion import (
            _PLATEAU_T1,
            _PLATEAU_T2,
            _PLATEAU_T3,
            _PLATEAU_T4,
        )

        low = np.fmin(maiac_values, cams_values)
        high = np.fmax(maiac_values, cams_values)
        count = np.isfinite(maiac_values).astype(np.uint8) + np.isfinite(cams_values).astype(
            np.uint8
        )
        mid = np.nansum(np.stack((maiac_values, cams_values)), axis=0)
        mid = np.divide(
            mid, count, out=np.full(mid.shape, np.nan, dtype=np.float32), where=count > 0
        )
        rise = np.clip((high - _PLATEAU_T1) / (_PLATEAU_T2 - _PLATEAU_T1), 0.0, 1.0)
        fall = np.clip((high - _PLATEAU_T3) / (_PLATEAU_T4 - _PLATEAU_T3), 0.0, 1.0)
        expected = np.where(
            high < _PLATEAU_T3,
            (1.0 - rise) * low + rise * mid,
            (1.0 - fall) * mid + fall * high,
        ).astype(np.float32)
    else:
        raise ValueError(f"unsupported AOD prior combination {combination!r}")
    if not np.allclose(fused_values, expected, rtol=0.0, atol=1.0e-6, equal_nan=True):
        raise RuntimeError(
            f"{combination}(QA-MAIAC, CAMS) fusion does not match archived components"
        )

    finite_maiac = np.isfinite(maiac_values)
    finite_cams = np.isfinite(cams_values)
    winner_values = np.full(template.shape, AOD_PRIOR_SOURCE_MISSING, dtype=np.uint8)
    if str(combination) == "mean":
        winner_values[finite_maiac & ~finite_cams] = AOD_PRIOR_SOURCE_MAIAC
        winner_values[finite_cams & ~finite_maiac] = AOD_PRIOR_SOURCE_CAMS
        winner_values[finite_maiac & finite_cams] = AOD_PRIOR_SOURCE_BOTH
    else:
        winner_values[finite_maiac & (~finite_cams | (maiac_values >= cams_values))] = (
            AOD_PRIOR_SOURCE_MAIAC
        )
        winner_values[finite_cams & (~finite_maiac | (cams_values > maiac_values))] = (
            AOD_PRIOR_SOURCE_CAMS
        )
    winner = xr.DataArray(
        winner_values,
        dims=template.dims,
        coords=template.coords,
        name="aod_prior_source",
    ).rio.write_crs(template.rio.crs)
    return _FusedAODPrior(
        state=fused,
        maiac_aot=maiac_aot,
        cams_aot=cams_grid.aot,
        winner=winner,
    )


def _reuse_fused_aod_on_grid(
    path: Path,
    *,
    template: xr.DataArray,
    combination: str,
    layout: str,
    uncertainty_mode: str,
    temporal_window_days: int,
) -> tuple[_FusedAODPrior, bool]:
    """Load a previously audited current-scene AOD prior without refetching it.

    Historical-surface ablations repeatedly solve the same current acquisition
    while changing only the historical dictionary.  Reusing the archived
    MAIAC/CAMS fields removes duplicate Earthaccess traffic.  The complete
    fusion contract and every component field are checked before reuse.
    """

    with np.load(path, allow_pickle=False) as archive:
        required = {
            "maiac_aod_prior60",
            "cams_aod_prior60",
            "aod_prior_source60",
            "fused_aod_prior_uncertainty",
            "aod_prior_fusion",
            "maiac_temporal_window_days",
            "maiac_prior_available",
        }
        missing = sorted(required - set(archive.files))
        if missing:
            raise ValueError(f"AOD-prior reuse archive lacks fields: {missing}")
        explicit_contract = {
            "pixelwise_fused_aod_prior60",
            "aod_prior_layout",
            "aod_prior_uncertainty_mode",
        }.issubset(archive.files)
        legacy_committed = bool(
            not explicit_contract
            and "schema_version" in archive.files
            and "solver_contract" in archive.files
            and str(_scalar(archive["schema_version"])) == "siac_optimal_m5_teacher_v7"
            and str(_scalar(archive["solver_contract"])) == COMMITTED_C0_CONTRACT
            and "fused_aod_prior60" in archive.files
            and not (
                {"pixelwise_fused_aod_prior60", "aod_prior_layout", "aod_prior_uncertainty_mode"}
                & set(archive.files)
            )
        )
        if not explicit_contract and not legacy_committed:
            missing_explicit = sorted(
                {
                    "pixelwise_fused_aod_prior60",
                    "aod_prior_layout",
                    "aod_prior_uncertainty_mode",
                }
                - set(archive.files)
            )
            raise ValueError(
                "AOD-prior reuse archive lacks explicit contract fields and is not "
                f"the audited legacy committed-v7 schema: {missing_explicit}"
            )
        actual_layout = (
            str(_scalar(archive["aod_prior_layout"])) if explicit_contract else "spatial"
        )
        archived_combination = str(_scalar(archive["aod_prior_fusion"]))
        # The frozen committed-C0 archive used ``pixelwise_max`` while the
        # current CLI uses ``max`` for the same spatial, per-pixel operation.
        # Canonicalise only that exact alias and only for a spatial layout;
        # scalar-component max remains a different RT prior contract.
        actual_combination = (
            "max"
            if archived_combination == "pixelwise_max" and actual_layout == "spatial"
            else archived_combination
        )
        actual_contract = {
            "combination": actual_combination,
            "layout": actual_layout,
            "uncertainty_mode": (
                str(_scalar(archive["aod_prior_uncertainty_mode"]))
                if explicit_contract
                else "provider_primary"
            ),
            "temporal_window_days": int(_scalar(archive["maiac_temporal_window_days"])),
        }
        expected_contract = {
            "combination": str(combination),
            "layout": str(layout),
            "uncertainty_mode": str(uncertainty_mode),
            "temporal_window_days": int(temporal_window_days),
        }
        if actual_contract != expected_contract:
            raise ValueError(
                f"AOD-prior reuse contract {actual_contract} "
                f"(archived combination={archived_combination!r}) != requested "
                f"{expected_contract}"
            )

        maiac_values = np.asarray(archive["maiac_aod_prior60"], dtype=np.float32)
        cams_values = np.asarray(archive["cams_aod_prior60"], dtype=np.float32)
        fused_values = np.asarray(
            archive["pixelwise_fused_aod_prior60" if explicit_contract else "fused_aod_prior60"],
            dtype=np.float32,
        )
        winner_values = np.asarray(archive["aod_prior_source60"], dtype=np.uint8)
        if "pixelwise_fused_aod_prior_uncertainty60" in archive.files:
            uncertainty_values = np.asarray(
                archive["pixelwise_fused_aod_prior_uncertainty60"], dtype=np.float32
            )
        else:
            # The legacy archive retained the exact 60->20 bilinear field.  A
            # source-grid centre occurs at every [1::3, 1::3] target pixel.
            uncertainty20 = np.asarray(archive["fused_aod_prior_uncertainty"], dtype=np.float32)
            uncertainty_values = uncertainty20[1::3, 1::3]
        maiac_available = bool(_scalar(archive["maiac_prior_available"]))

    for name, values in {
        "MAIAC": maiac_values,
        "CAMS": cams_values,
        "pixelwise fused AOD": fused_values,
        "source code": winner_values,
        "AOD uncertainty": uncertainty_values,
    }.items():
        if values.shape != template.shape:
            raise ValueError(
                f"AOD-prior reuse {name} shape {values.shape} != solver grid {template.shape}"
            )
    if np.any(~np.isin(winner_values, np.asarray((0, 1, 2, 3), dtype=np.uint8))):
        raise ValueError("AOD-prior reuse archive contains an unknown source code")

    if str(combination) == "max":
        expected = np.fmax(maiac_values, cams_values)
    elif str(combination) == "mean":
        count = np.isfinite(maiac_values).astype(np.uint8) + np.isfinite(cams_values).astype(
            np.uint8
        )
        expected = np.divide(
            np.nansum(np.stack((maiac_values, cams_values)), axis=0),
            count,
            out=np.full(template.shape, np.nan, dtype=np.float32),
            where=count > 0,
        )
    elif str(combination) == "plateau":
        # Recomputed here independently of the provider, as the other branches
        # are -- that is the point of this audit. The knots are imported rather
        # than restated so the two implementations cannot drift apart.
        from siac.adapters.atmo.fusion import (
            _PLATEAU_T1,
            _PLATEAU_T2,
            _PLATEAU_T3,
            _PLATEAU_T4,
        )

        low = np.fmin(maiac_values, cams_values)
        high = np.fmax(maiac_values, cams_values)
        count = np.isfinite(maiac_values).astype(np.uint8) + np.isfinite(cams_values).astype(
            np.uint8
        )
        mid = np.nansum(np.stack((maiac_values, cams_values)), axis=0)
        mid = np.divide(
            mid, count, out=np.full(mid.shape, np.nan, dtype=np.float32), where=count > 0
        )
        rise = np.clip((high - _PLATEAU_T1) / (_PLATEAU_T2 - _PLATEAU_T1), 0.0, 1.0)
        fall = np.clip((high - _PLATEAU_T3) / (_PLATEAU_T4 - _PLATEAU_T3), 0.0, 1.0)
        expected = np.where(
            high < _PLATEAU_T3,
            (1.0 - rise) * low + rise * mid,
            (1.0 - fall) * mid + fall * high,
        ).astype(np.float32)
    else:
        raise ValueError(f"unsupported AOD prior combination {combination!r}")
    if not np.allclose(fused_values, expected, rtol=0.0, atol=1.0e-6, equal_nan=True):
        raise ValueError("AOD-prior reuse fused field does not match its archived components")
    if not np.isfinite(uncertainty_values).any() or np.any(
        uncertainty_values[np.isfinite(uncertainty_values)] <= 0.0
    ):
        raise ValueError("AOD-prior reuse uncertainty has no valid positive support")

    def field(values: np.ndarray, name: str) -> xr.DataArray:
        return xr.DataArray(
            values,
            dims=template.dims,
            coords=template.coords,
            name=name,
        ).rio.write_crs(template.rio.crs)

    aot = field(fused_values, "aot")
    uncertainty = field(uncertainty_values, "aot_unc")
    placeholder = xr.zeros_like(template, dtype=np.float32)
    state = AtmosphericState(
        aot=aot,
        tcwv=placeholder,
        tco3=placeholder,
        aot_unc=uncertainty,
        tcwv_unc=placeholder,
        tco3_unc=placeholder,
        elevation=placeholder,
    )
    return (
        _FusedAODPrior(
            state=state,
            maiac_aot=field(maiac_values, "maiac_aot"),
            cams_aot=field(cams_values, "cams_aot"),
            winner=field(winner_values, "aod_prior_source"),
        ),
        maiac_available,
    )


def _atmo_on_grid(
    aerosol_prior: AtmosphericState,
    *,
    template: xr.DataArray,
    elevation: xr.DataArray,
    tcwv: xr.DataArray,
    tco3: xr.DataArray,
) -> AtmosphericState:
    aot = _align(aerosol_prior.aot, template)
    aot_unc = _align(aerosol_prior.aot_unc, template)
    tcwv = _align(tcwv, template)
    tco3 = _align(tco3, template)
    elevation = _align(elevation, template)
    return AtmosphericState(
        aot=aot,
        tcwv=tcwv,
        tco3=tco3,
        aot_unc=aot_unc,
        tcwv_unc=xr.apply_ufunc(np.maximum, tcwv * 0.15, np.float32(0.2)),
        tco3_unc=xr.apply_ufunc(np.maximum, tco3 * 0.1, np.float32(0.01)),
        elevation=elevation,
    )


def _teacher_contract(name: str) -> dict[str, Any]:
    try:
        return TEACHER_CONTRACTS[str(name)]
    except KeyError as exc:
        raise ValueError(f"unknown teacher solver contract {name!r}") from exc


def _parse_solver_solve_bands(value: str) -> tuple[str, ...]:
    bands = tuple(part.strip().upper() for part in str(value).split(",") if part.strip())
    if len(set(bands)) != len(bands):
        raise argparse.ArgumentTypeError("--solver-solve-bands cannot contain duplicates")
    unknown = tuple(band for band in bands if band not in SOLVER_VISIBLE)
    if unknown:
        raise argparse.ArgumentTypeError(
            f"--solver-solve-bands contains unsupported bands {unknown}; "
            f"allowed bands are {SOLVER_VISIBLE}"
        )
    missing_reference = tuple(band for band in VISIBLE if band not in bands)
    if missing_reference:
        raise argparse.ArgumentTypeError(
            "--solver-solve-bands must retain B02,B03,B04 so the teacher output "
            f"contract remains unchanged; missing {missing_reference}"
        )
    return bands


def _parse_solver_band_cost_weights(value: str) -> tuple[float, ...]:
    try:
        weights = tuple(float(part.strip()) for part in str(value).split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--solver-band-cost-weights must be comma-separated numbers"
        ) from exc
    if not weights or not np.all(np.isfinite(weights)) or any(weight < 0 for weight in weights):
        raise argparse.ArgumentTypeError(
            "--solver-band-cost-weights must be finite and non-negative"
        )
    if not any(weight > 0 for weight in weights):
        raise argparse.ArgumentTypeError(
            "--solver-band-cost-weights requires at least one positive value"
        )
    return weights


def _solver_config(
    contract: dict[str, Any],
    *,
    allow_toa_above_one: bool = False,
    solve_bands: tuple[str, ...] = VISIBLE,
    cost_aggregation: str = "sum",
    backstop_uncertainty_scale: float = 1.0,
    band_cost_weights: tuple[float, ...] | None = None,
    quadratic_refine: bool = False,
    observation_likelihood: str = "inverse_chi2",
    student_t_dof: float = 4.0,
    toa_uncertainty_floor: float = 0.003,
    rt_uncertainty_aod_scale: float = 0.0,
    unresolved_floor_aot: float = 0.05,
    unresolved_high_prior_aot: float = 0.40,
    unresolved_prior_conflict_sigma: float = 1.5,
    unresolved_band_spread: float = 0.50,
    unresolved_forward_z: float = 3.0,
) -> SolverAlgorithmConfig:
    """Construct the frozen C0 solver contract without relying on defaults."""

    return SolverAlgorithmConfig(
        method="surface_driven",
        aerosol_resolution=60.0,
        surface_driven_pool_radius_m=600.0,
        surface_driven_pool_min_count=1,
        surface_driven_backstop_calibrated=True,
        surface_driven_backstop_uncertainty_scale=backstop_uncertainty_scale,
        surface_driven_cost_aggregation=cost_aggregation,
        surface_driven_band_cost_weights=band_cost_weights,
        surface_driven_quadratic_refine=quadratic_refine,
        surface_driven_observation_likelihood=observation_likelihood,
        surface_driven_student_t_dof=student_t_dof,
        surface_driven_toa_uncertainty_floor=toa_uncertainty_floor,
        surface_driven_rt_uncertainty_aod_scale=rt_uncertainty_aod_scale,
        surface_driven_unresolved_floor_aot=unresolved_floor_aot,
        surface_driven_unresolved_high_prior_aot=unresolved_high_prior_aot,
        surface_driven_unresolved_prior_conflict_sigma=unresolved_prior_conflict_sigma,
        surface_driven_unresolved_band_spread=unresolved_band_spread,
        surface_driven_unresolved_forward_z=unresolved_forward_z,
        surface_driven_reference_tcwv=None,
        surface_driven_solve_bands=solve_bands,
        surface_driven_tau_dependent_prior=True,
        surface_driven_aerosol_species=str(contract["aerosol_species"]),
        surface_driven_scene_mean_geometry=True,
        surface_driven_aot_axis=str(contract["aot_axis"]),
        surface_driven_allow_cloud_retrieval=False,
        surface_driven_toa_upper_bound=None if allow_toa_above_one else 1.0,
        surface_driven_ignore_cloud_water=False,
    )


def _scene_mean_aod(aod: Any) -> tuple[float, int]:
    values = np.asarray(aod, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return (float(np.mean(finite)), int(finite.size)) if finite.size else (float("nan"), 0)


def _winner_aod_for_realizations(
    index_path: Path,
    realizations: tuple[str, ...],
    expected_shape: tuple[int, int],
) -> np.ndarray:
    """Map every corrected monthly realization back to its winning-day AOD.

    A monthly mosaic can contain pixels from several acquisitions, so one AOD
    scalar per realization is not faithful.  The winner index records the
    acquisition selected at every pixel; this function expands its per-day AOD
    table onto the same ``(realization, y, x)`` grid as the corrected library.
    """

    index = read_mosaic_index(index_path)
    if tuple(index.winners.shape[1:]) != tuple(expected_shape):
        raise ValueError(
            f"historical AOD-screen index grid {index.winners.shape[1:]} "
            f"!= dictionary grid {expected_shape}"
        )
    month_positions = {month: position for position, month in enumerate(index.months)}
    if len(month_positions) != len(index.months):
        raise ValueError(f"historical AOD-screen index has duplicate months: {index_path}")

    output = np.full((len(realizations), *expected_shape), np.nan, dtype=np.float32)
    for realization_index, month in enumerate(realizations):
        if month not in month_positions:
            raise ValueError(
                f"historical AOD-screen index {index_path} lacks dictionary realization {month}"
            )
        winner_plane = np.asarray(index.winners[month_positions[month]], dtype=np.int32)
        images = index.images_by_month[month]
        used = np.unique(winner_plane[winner_plane >= 0])
        if used.size and int(used.max()) >= len(images):
            raise ValueError(
                f"historical AOD-screen winner {int(used.max())} is outside the "
                f"{len(images)} images for {month}"
            )
        for image_index in used.tolist():
            day = images[int(image_index)].day
            aod = index.day_aod.get(day)
            if aod is not None and np.isfinite(aod) and float(aod) >= 0.0:
                output[realization_index][winner_plane == int(image_index)] = np.float32(aod)
    return output


def _screen_historical_comp(
    comp: np.ndarray,
    winner_aod: np.ndarray,
    *,
    max_aod: float,
    min_realizations: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Prefer low-AOD history only where enough complete samples remain.

    Pixels with at least ``min_realizations`` complete, low-AOD histories use
    only those histories.  Everywhere else the original stack is retained.
    This fallback makes the experiment coverage preserving: a sparse archive
    cannot turn into a missing teacher merely because it lacks clean days.
    """

    values = np.asarray(comp, dtype=np.float32)
    aod = np.asarray(winner_aod, dtype=np.float32)
    if values.ndim != 4:
        raise ValueError(f"historical surface stack must be 4-D, got {values.shape}")
    if aod.shape != (values.shape[0], *values.shape[2:]):
        raise ValueError(
            f"winner-AOD shape {aod.shape} does not match surface stack {values.shape}"
        )
    if not np.isfinite(max_aod) or float(max_aod) <= 0.0:
        raise ValueError("historical AOD-screen threshold must be finite and > 0")
    if int(min_realizations) < 1:
        raise ValueError("historical AOD-screen minimum realizations must be >= 1")

    complete = np.all(np.isfinite(values), axis=1)
    low_aod_complete = complete & np.isfinite(aod) & (aod <= float(max_aod))
    low_aod_count = np.count_nonzero(low_aod_complete, axis=0)
    screen_applied = low_aod_count >= int(min_realizations)
    retain = (~screen_applied)[np.newaxis] | low_aod_complete
    screened = np.where(retain[:, np.newaxis], values, np.nan).astype(np.float32)

    before_complete = int(np.count_nonzero(complete))
    after_complete = int(np.count_nonzero(np.all(np.isfinite(screened), axis=1)))
    active_counts = low_aod_count[screen_applied]
    summary = {
        "mode": "coverage_preserving_low_aod",
        "max_aod": float(max_aod),
        "min_realizations": int(min_realizations),
        "pixels": int(screen_applied.size),
        "active_pixels": int(np.count_nonzero(screen_applied)),
        "active_fraction": float(np.mean(screen_applied)),
        "complete_samples_before": before_complete,
        "complete_samples_after": after_complete,
        "retained_complete_sample_fraction": (
            float(after_complete / before_complete) if before_complete else 0.0
        ),
        "median_low_aod_realizations_when_active": (
            float(np.median(active_counts)) if active_counts.size else 0.0
        ),
        "minimum_low_aod_realizations_when_active": (
            int(np.min(active_counts)) if active_counts.size else 0
        ),
    }
    return screened, summary


def _existing_ok(
    path: Path,
    *,
    schema: str = SCHEMA,
    contract_name: str = COMMITTED_C0_CONTRACT,
    expected_policy: dict[str, Any] | None = None,
    expected_aod_extraction: str = "scene_mean",
) -> bool:
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as value:
            required = {
                "schema_version",
                "surface",
                "surface_uncertainty",
                "aod",
                "aod_uncertainty",
                "valid",
                "cloud_classes",
                "ocm_raw_classes",
                "retrieval_valid_mask",
                "surface_bands",
                "aod_axis_nodes",
                "aod_scene_mean",
                "aod_scene_valid_count",
                "aod_extraction",
                "solver_contract",
                "aerosol_species",
                "maiac_aod_prior60",
                "cams_aod_prior60",
                "fused_aod_prior60",
                "aod_prior_source60",
                "maiac_prior_available",
            }
            base_ok = bool(
                required.issubset(value.files)
                and str(_scalar(value["schema_version"])) == str(schema)
                and np.asarray(value["surface"]).shape[-1] == 3
                and str(_scalar(value["solver_contract"])) == str(contract_name)
                and int(_scalar(value["aod_axis_nodes"]))
                == int(_teacher_contract(contract_name)["aot_axis_nodes"])
                and str(_scalar(value["aerosol_species"]))
                == str(_teacher_contract(contract_name)["aerosol_species"])
                and str(_scalar(value["aod_extraction"])) == expected_aod_extraction
            )
            if not base_ok:
                return False
            for key, expected in (expected_policy or {}).items():
                if key not in value.files:
                    return False
                actual = _scalar(value[key])
                if isinstance(expected, (float, int, np.floating, np.integer)):
                    actual_number = float(actual)
                    expected_number = float(expected)
                    if np.isnan(expected_number):
                        if not np.isnan(actual_number):
                            return False
                    elif not np.isclose(actual_number, expected_number, rtol=0.0, atol=1e-9):
                        return False
                elif str(actual) != str(expected):
                    return False
            return True
    except (OSError, ValueError, KeyError, EOFError):
        return False


def build_one(matchup_id: str, args: argparse.Namespace) -> dict[str, Any]:
    contract = _teacher_contract(args.solver_contract)
    solve_bands = tuple(args.solver_solve_bands)
    target_bands = tuple(getattr(args, "surface_target_bands", VISIBLE))
    band_cost_weights = args.solver_band_cost_weights
    if band_cost_weights is not None and len(band_cost_weights) != len(solve_bands):
        raise ValueError(
            "--solver-band-cost-weights must have one value per solve band: "
            f"{len(band_cost_weights)} != {len(solve_bands)}"
        )
    current_path = args.current_root / f"{matchup_id}.npz"
    fine_path = args.fine20_root / f"{matchup_id}.npz"
    dictionary_path = args.dictionary_root / f"{matchup_id}.npz"
    reduced_path = args.reduced_root / f"{matchup_id}.npz"
    output_path = args.output_root / f"{matchup_id}.npz"
    receipt_path = args.output_root / f"{matchup_id}.receipt.json"
    source_quality_root = getattr(args, 'source_quality_root', None)
    source_quality_path = None if source_quality_root is None else Path(source_quality_root) / f'{matchup_id}.npz'
    source_quality_hash = None if source_quality_path is None else _sha256(source_quality_path)
    aod_prior_reuse_path = (
        None
        if args.aod_prior_reuse_root is None
        else args.aod_prior_reuse_root / f"{matchup_id}.npz"
    )
    deep_blue_prior_path = (
        None
        if args.deep_blue_prior_root is None
        else args.deep_blue_prior_root / f"{matchup_id}.json"
    )
    aod_prior_reuse_sha256: str | None = None
    if aod_prior_reuse_path is not None:
        if not aod_prior_reuse_path.is_file():
            raise FileNotFoundError(aod_prior_reuse_path)
        aod_prior_reuse_sha256 = _sha256(aod_prior_reuse_path)
    deep_blue_prior_sha256: str | None = None
    if deep_blue_prior_path is not None:
        if not deep_blue_prior_path.is_file():
            raise FileNotFoundError(deep_blue_prior_path)
        deep_blue_prior_sha256 = _sha256(deep_blue_prior_path)
    snow_policy_on = str(args.snow_support_policy) == "recurrent-library"
    provided_surface_aod = args.surface_target_aod550 is not None
    output_schema = (
        AERONET_CONDITIONED_SCHEMA
        if provided_surface_aod
        else SNOW_POLICY_SCHEMA
        if snow_policy_on
        else SCHEMA
    )
    aod_extraction = "provided_scene_aod550" if provided_surface_aod else "scene_mean"
    historical_screen_on = args.historical_aod_screen_max is not None
    historical_index_path: Path | None = None
    historical_index_sha256: str | None = None
    expected_policy: dict[str, Any] = {
        "aod_prior_fusion": str(args.aod_prior_combination),
        "aod_prior_layout": str(args.aod_prior_layout),
        "aod_prior_uncertainty_mode": str(args.aod_prior_uncertainty),
        "maiac_temporal_window_days": int(args.maiac_temporal_window_days),
        "tcwv_mode": str(args.tcwv_mode),
        "tcwv_missing_policy": str(args.tcwv_missing_policy),
        "tco3_mode": str(args.tco3_mode),
        "anchor_geometry": str(args.anchor_geometry),
        "seasonal_robust_clip": float(args.seasonal_robust_clip),
        # Run-level contract: the requested model. The per-scene archive records
        # the one actually used, which differs on recurrent-snow scenes.
        "predictor_model": str(args.predictor_model),
        "ocm_thin_policy": str(args.ocm_thin_policy),
        "water_mask_mode": str(args.water_mask_mode),
        "solver_toa_upper_bound": np.nan if args.allow_toa_above_one else 1.0,
        "snow_support_policy": str(args.snow_support_policy),
        "solver_quadratic_refine": bool(args.solver_quadratic_refine),
        "solver_observation_likelihood": str(args.solver_observation_likelihood),
        "solver_student_t_dof": float(args.solver_student_t_dof),
        "solver_toa_uncertainty_floor": float(args.solver_toa_uncertainty_floor),
        "solver_rt_uncertainty_aod_scale": float(args.solver_rt_uncertainty_aod_scale),
        "solver_unresolved_floor_aot": float(args.solver_unresolved_floor_aot),
        "solver_unresolved_high_prior_aot": float(args.solver_unresolved_high_prior_aot),
        "solver_unresolved_prior_conflict_sigma": float(
            args.solver_unresolved_prior_conflict_sigma
        ),
        "solver_unresolved_band_spread": float(args.solver_unresolved_band_spread),
        "solver_unresolved_forward_z": float(args.solver_unresolved_forward_z),
    }
    if provided_surface_aod:
        expected_policy.update(
            {
                "surface_target_aod_source": str(args.surface_target_aod_source),
                "surface_target_aod550": float(np.float32(args.surface_target_aod550)),
                "surface_target_aod550_uncertainty": float(
                    np.float32(args.surface_target_aod550_uncertainty)
                ),
            }
        )
    if source_quality_hash is not None:
        expected_policy['source_quality_sha256'] = source_quality_hash
    # Old committed-C0 archives predate this explicit field. Preserve their
    # cache validity, while requiring exact provenance for any experimental
    # solve-band extension.
    if solve_bands != VISIBLE:
        expected_policy["solver_solve_bands_csv"] = ",".join(solve_bands)
    if str(args.solver_cost_aggregation) != "sum":
        expected_policy["solver_cost_aggregation"] = str(args.solver_cost_aggregation)
    if band_cost_weights is not None:
        expected_policy["solver_band_cost_weights_csv"] = ",".join(
            f"{value:.12g}" for value in band_cost_weights
        )
    if not np.isclose(float(args.solver_backstop_uncertainty_scale), 1.0):
        expected_policy["solver_backstop_uncertainty_scale"] = float(
            args.solver_backstop_uncertainty_scale
        )
    if aod_prior_reuse_sha256 is not None:
        expected_policy["aod_prior_reuse_sha256"] = aod_prior_reuse_sha256
    if deep_blue_prior_sha256 is not None:
        expected_policy.update(
            {
                "deep_blue_prior_sha256": deep_blue_prior_sha256,
                "deep_blue_quality": str(args.deep_blue_quality),
                "deep_blue_fusion": str(args.deep_blue_fusion),
            }
        )
    if historical_screen_on:
        if args.historical_aod_screen_index_root is None:
            raise ValueError("historical AOD screening requires --historical-aod-screen-index-root")
        historical_index_path = args.historical_aod_screen_index_root / f"{matchup_id}.npz"
        if not historical_index_path.is_file():
            raise FileNotFoundError(historical_index_path)
        historical_index_sha256 = _sha256(historical_index_path)
        expected_policy.update(
            {
                "historical_aod_screen_mode": "coverage_preserving_low_aod",
                # The archive stores this compact scalar as float32; compare
                # against the same exact representation on cache validation.
                "historical_aod_screen_max": float(np.float32(args.historical_aod_screen_max)),
                "historical_aod_screen_min_realizations": int(
                    args.historical_aod_screen_min_realizations
                ),
                "historical_aod_screen_index_sha256": historical_index_sha256,
            }
        )
    if not args.force and _existing_ok(
        output_path,
        schema=output_schema,
        contract_name=args.solver_contract,
        expected_policy=expected_policy,
        expected_aod_extraction=aod_extraction,
    ):
        return {"matchup_id": matchup_id, "status": "exists", "output": str(output_path)}
    for path in (current_path, fine_path, dictionary_path, reduced_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    started = time.perf_counter()
    with (
        np.load(current_path, allow_pickle=False) as current,
        np.load(fine_path, allow_pickle=False) as fine,
        np.load(dictionary_path, allow_pickle=False) as dictionary,
        np.load(reduced_path, allow_pickle=False) as reduced,
    ):
        current_schema = str(_scalar(current["schema_version"]))
        fine_schema = str(_scalar(fine["schema_version"]))
        if current_schema != "siac_l1c_mixed_resolution_v3":
            raise ValueError(f"{matchup_id}: unsupported current TOA schema {current_schema!r}")
        if fine_schema != "siac_l1c_fine20_visible_v3":
            raise ValueError(f"{matchup_id}: unsupported fine TOA schema {fine_schema!r}")
        product_id = str(_scalar(current["product_id"]))
        platform = _platform(product_id)
        when = _scene_time(matchup_id)
        crs = str(_scalar(current["crs"]))
        transform20 = np.asarray(current["detail20_transform"], dtype=np.float64)[:6]
        transform60 = np.asarray(current["local60_transform"], dtype=np.float64)[:6]
        comp = np.asarray(dictionary["comp"], dtype=np.float32)
        dictionary_realizations = tuple(
            str(value) for value in np.asarray(dictionary["realizations"]).tolist()
        )
        dictionary_band_names = [
            str(value) for value in np.asarray(dictionary["band_names"]).tolist()
        ]
        dictionary_columns = _dictionary_columns(dictionary_band_names)
        dictionary_species = (
            str(_scalar(dictionary["aerosol_species"]))
            if "aerosol_species" in dictionary.files
            else None
        )
        reduced_species = (
            str(_scalar(reduced["aerosol_species"])) if "aerosol_species" in reduced.files else None
        )
        expected_species = str(contract["aerosol_species"])
        if dictionary_species != expected_species or reduced_species != expected_species:
            raise ValueError(
                f"{matchup_id}: teacher RT-space mismatch: expected {expected_species!r}, "
                f"dictionary={dictionary_species!r}, reduced={reduced_species!r}"
            )
        dictionary_epsg = int(_scalar(dictionary["epsg"]))
        t0 = np.asarray(reduced["surface"], dtype=np.float32)
        t0_unc = np.asarray(reduced["surface_uncertainty"], dtype=np.float32)
        context = np.asarray(reduced["current_context"], dtype=np.float32)
        historical_screen: dict[str, Any] = {
            "mode": "off",
            "max_aod": None,
            "min_realizations": None,
            "pixels": int(comp.shape[-2] * comp.shape[-1]),
            "active_pixels": 0,
            "active_fraction": 0.0,
            "complete_samples_before": int(np.count_nonzero(np.all(np.isfinite(comp), axis=1))),
            "complete_samples_after": int(np.count_nonzero(np.all(np.isfinite(comp), axis=1))),
            "retained_complete_sample_fraction": 1.0,
            "median_low_aod_realizations_when_active": 0.0,
            "minimum_low_aod_realizations_when_active": 0,
        }
        if historical_screen_on:
            if historical_index_path is None or historical_index_sha256 is None:
                raise RuntimeError("historical AOD-screen contract was not initialized")
            dictionary_receipt_path = dictionary_path.with_suffix(".receipt.json")
            dictionary_receipt = json.loads(dictionary_receipt_path.read_text(encoding="utf-8"))
            expected_index_sha256 = str(
                dictionary_receipt.get("winner_index", {}).get("sha256", "")
            )
            if expected_index_sha256 != historical_index_sha256:
                raise ValueError(
                    f"{matchup_id}: historical AOD-screen index does not match the "
                    "index used to construct the corrected surface dictionary"
                )
            winner_aod = _winner_aod_for_realizations(
                historical_index_path,
                dictionary_realizations,
                tuple(int(value) for value in comp.shape[-2:]),
            )
            comp, historical_screen = _screen_historical_comp(
                comp,
                winner_aod,
                max_aod=float(args.historical_aod_screen_max),
                min_realizations=int(args.historical_aod_screen_min_realizations),
            )
            screened_median, screened_uncertainty = _robust_temporal_median(
                comp,
                clip=1.5,
                floor=0.006,
            )
            t0 = np.moveaxis(screened_median, 0, -1)
            t0_unc = np.moveaxis(screened_uncertainty, 0, -1)
        toa20_values = {
            **{band: np.asarray(fine[f"fine20_{band}"], dtype=np.float32) for band in solve_bands},
            **{band: np.asarray(current[f"detail20_{band}"], dtype=np.float32) for band in ANCHORS},
        }
        toa60_values = {
            band: np.asarray(current[f"local60_{band}"], dtype=np.float32)
            for band in (*solve_bands, *ANCHORS)
        }
        template20 = _grid(t0[..., 0], transform20, crs, name="template20")
        template60 = _grid(toa60_values["B02"], transform60, crs, name="template60")
        geometry20 = _constant_geometry(current, template20)
        geometry60 = _constant_geometry(current, template60)

    source_quality_record = None
    source_retrieval_valid20 = np.ones(template20.shape, bool)
    if source_quality_path is not None:
        from experiments.current_date_toa_prior.teacher_source_validity import load_masked_observations
        toa20_values, toa60_values, source_retrieval_valid20, source_quality_record = load_masked_observations(
            source_quality_path, toa20_values, toa60_values,
            matchup_id=matchup_id, transform=transform20, crs=crs)
        if not source_retrieval_valid20.any():
            raise TeacherSceneIneligible('no_instrument_valid_visible_and_anchor_pixels', source_quality_record)

    sensor = load_sensor_config_with_rsrf("MSI", platform)
    toa20 = xr.Dataset(
        {name: _grid(value, transform20, crs, name=name) for name, value in toa20_values.items()}
    )
    toa60 = xr.Dataset(
        {name: _grid(value, transform60, crs, name=name) for name, value in toa60_values.items()}
    )
    raw_classes20 = OmniCloudMaskProvider(inference_device="cpu").predict(
        toa20["B04"],
        toa20["B03"],
        toa20["B8A"],
        class_mapping=OCM_RAW_IDENTITY,
        unmapped_to_missing=False,
    )
    ocm_inputs_valid = np.all(
        np.stack([np.isfinite(np.asarray(toa20[name])) for name in ("B04", "B03", "B8A")]),
        axis=0,
    )
    raw_class_values = np.asarray(raw_classes20, dtype=np.uint8)
    # Standard SIAC classes for downstream diagnostics: class 1 means usable
    # by this retrieval (OCM clear OR thin cloud), class 2 is thick cloud, and
    # class 3 is shadow. The separate raw array preserves the distinction.
    ocm_usable_classes = (
        OCM_RETRIEVAL_USABLE if str(args.ocm_thin_policy) == "retain" else OCM_RETRIEVAL_CLEAR_ONLY
    )
    retrieval_class_values, retrieval_valid20 = _ocm_retrieval_policy(
        raw_class_values,
        ocm_inputs_valid,
        usable_classes=ocm_usable_classes,
    )
    retrieval_valid20 &= source_retrieval_valid20
    if source_quality_path is not None and not retrieval_valid20.any():
        raise TeacherSceneIneligible('no_clear_or_thin_instrument_valid_pixels', dict(
            source_quality=source_quality_record,
            ocm_counts={str(int(k)): int(v) for k, v in zip(*np.unique(raw_class_values, return_counts=True))}))
    snow_decision: dict[str, Any] = {
        "eligible": True,
        "reason": "",
        "classification": "snow_policy_off",
        "current_snow_fraction": None,
        "supporting_realizations": None,
        "eligible_realizations": None,
        "support_ratio": None,
        "support_wilson_lower": None,
        "support_wilson_upper": None,
        "history_snow_sample_fraction": None,
    }
    if snow_policy_on:
        try:
            green_column = dictionary_columns["B03"]
            swir16_column = dictionary_columns["B11"]
        except KeyError as exc:
            raise ValueError(
                f"{matchup_id}: spectral dictionary lacks green/SWIR16 snow bands"
            ) from exc
        snow_support = profile_library_support(
            toa20_values["B03"],
            toa20_values["B11"],
            retrieval_valid20,
            comp[:, green_column],
            comp[:, swir16_column],
            ndsi_threshold=float(args.snow_ndsi_threshold),
            green_threshold=float(args.snow_green_threshold),
            realization_coverages=(float(args.snow_realization_coverage),),
        )
        snow_key = (
            f"ndsi_{float(args.snow_ndsi_threshold):.2f}_green_"
            f"{float(args.snow_green_threshold):.2f}"
        ).replace(".", "p")
        snow_decision = classify_profile(
            {"spectral_profiles": {snow_key: snow_support}},
            ndsi_threshold=float(args.snow_ndsi_threshold),
            green_threshold=float(args.snow_green_threshold),
            current_scene_snow_fraction=float(args.current_scene_snow_fraction),
            realization_coverage=float(args.snow_realization_coverage),
            recurrent_lower_bound=float(args.snow_recurrent_lower_bound),
            minimum_history_realizations=int(args.snow_minimum_history_realizations),
        )
        if not bool(snow_decision["eligible"]):
            raise TeacherSceneIneligible('current_snow_without_recurrent_library_support', snow_decision)
    classes20 = xr.DataArray(
        retrieval_class_values,
        dims=template20.dims,
        coords=template20.coords,
        name="cloud_classes",
    ).rio.write_crs(crs)
    retrieval_valid60 = retrieval_valid20.reshape(
        template60.shape[0], 3, template60.shape[1], 3
    ).all(axis=(1, 3))
    if source_quality_path is not None and not retrieval_valid60.any():
        raise TeacherSceneIneligible('no_fully_usable_60m_solver_cells',
                                     dict(usable_20m_pixels=int(retrieval_valid20.sum())))
    cloud60 = xr.DataArray(
        ~retrieval_valid60, dims=template60.dims, coords=template60.coords
    ).rio.write_crs(crs)

    bounds = (
        float(template60.x.min() - 30.0),
        float(template60.y.min() - 30.0),
        float(template60.x.max() + 30.0),
        float(template60.y.max() + 30.0),
    )
    maiac_error: dict[str, str] | None = None
    if provided_surface_aod:
        target_aod60 = xr.full_like(
            template60, np.float32(args.surface_target_aod550)
        ).rio.write_crs(crs)
        target_uncertainty60 = xr.full_like(
            template60, np.float32(args.surface_target_aod550_uncertainty)
        ).rio.write_crs(crs)
        placeholder60 = xr.zeros_like(template60, dtype=np.float32).rio.write_crs(crs)
        missing60 = xr.full_like(template60, np.float32(np.nan)).rio.write_crs(crs)
        source60 = xr.full_like(template60, AOD_PRIOR_SOURCE_MISSING, dtype=np.uint8).rio.write_crs(
            crs
        )
        fused_prior60 = _FusedAODPrior(
            state=AtmosphericState(
                aot=target_aod60,
                tcwv=placeholder60,
                tco3=placeholder60,
                aot_unc=target_uncertainty60,
                tcwv_unc=placeholder60,
                tco3_unc=placeholder60,
                elevation=placeholder60,
            ),
            maiac_aot=missing60,
            cams_aot=missing60,
            winner=source60,
        )
        maiac_available = False
        maiac_error = {
            "type": "intentionally_skipped",
            "message": "provided-AOD surface target does not read an aerosol prior",
        }
    elif aod_prior_reuse_path is not None:
        fused_prior60, maiac_available = _reuse_fused_aod_on_grid(
            aod_prior_reuse_path,
            template=template60,
            combination=str(args.aod_prior_combination),
            layout=str(args.aod_prior_layout),
            uncertainty_mode=str(args.aod_prior_uncertainty),
            temporal_window_days=int(args.maiac_temporal_window_days),
        )
    else:
        cams_provider = CAMSProvider(
            args.cams, temporal_interp=True, download_missing=False, cache_dir=args.cams_cache
        )
        cams_probe = cams_provider._load_cams_data(when)
        if cams_probe is None:
            raise RuntimeError(
                f"{matchup_id}: no CAMS current-date AOD; refusing default AOD label"
            )
        close_probe = getattr(cams_probe, "close", None)
        if close_probe is not None:
            close_probe()
        cams = cams_provider.get_prior(bounds, crs, when, 60.0)
        maiac_provider = MCD19AODProvider(
            cache_dir=args.maiac_cache,
            temporal_window_days=args.maiac_temporal_window_days,
            max_granules=args.maiac_max_granules,
            best_quality_qa=True,
            allow_default_prior=False,
        )
        maiac: AtmosphericState | None = None
        try:
            maiac = maiac_provider.get_prior(bounds, crs, when, 60.0)
        except Exception as exc:  # noqa: BLE001 - production fusion treats this source as absent
            maiac_error = {"type": type(exc).__name__, "message": str(exc)}
            logger.warning(
                "%s: no usable QA-best MCD19 prior; production fallback is CAMS-only (%s: %s)",
                matchup_id,
                type(exc).__name__,
                exc,
            )
        fused_prior60 = _fuse_aod_on_grid(
            maiac,
            cams,
            template=template60,
            obs_time=when,
            combination=str(args.aod_prior_combination),
        )
        maiac_available = maiac is not None
    aerosol_prior60 = fused_prior60.state
    if not provided_surface_aod and str(args.aod_prior_layout) == "scalar_component_max":
        aerosol_prior60 = _scalar_component_max_aod_prior(fused_prior60, template60)
    elif not provided_surface_aod and str(args.aod_prior_uncertainty) == "half_aod_floor_0p05":
        aerosol_prior60 = replace(
            aerosol_prior60,
            aot_unc=xr.apply_ufunc(
                np.maximum,
                aerosol_prior60.aot * np.float32(0.5),
                np.float32(0.05),
            ).astype(np.float32),
        )
    deep_blue_result = _DeepBluePriorResult(
        state=aerosol_prior60,
        available=False,
        aod550=None,
        uncertainty=None,
        metadata={
            "quality": str(args.deep_blue_quality),
            "fusion": str(args.deep_blue_fusion),
            "staged_path": None,
            "selected": None,
        },
    )
    if deep_blue_prior_path is not None and not provided_surface_aod:
        deep_blue_result = _apply_modis_deep_blue_prior(
            aerosol_prior60,
            deep_blue_prior_path,
            quality=str(args.deep_blue_quality),
            fusion=str(args.deep_blue_fusion),
        )
        aerosol_prior60 = deep_blue_result.state
    context20 = [
        _grid(context[..., index], transform20, crs, name=name)
        for index, name in enumerate(("elevation", "tcwv", "tco3"))
    ]
    tcwv20 = context20[1]
    tcwv_source = "current reduced-teacher context (L1C B8A/B09 CIBR)"
    tcwv_path: Path | None = None
    if str(args.tcwv_mode) == "l2a_wvp_gcs":
        if args.tcwv_root is None:
            raise ValueError("--tcwv-mode l2a_wvp_gcs requires --tcwv-root")
        requested_tcwv_path = args.tcwv_root / f"{matchup_id}.npz"
        if requested_tcwv_path.exists():
            tcwv_path = requested_tcwv_path
            with np.load(tcwv_path, allow_pickle=False) as tcwv_archive:
                if str(_scalar(tcwv_archive["schema_version"])) != "siac_l2a_wvp_gcs_v1":
                    raise ValueError(f"{matchup_id}: unsupported L2A WVP archive schema")
                tcwv_native = _grid(
                    np.asarray(tcwv_archive["tcwv60"], dtype=np.float32),
                    np.asarray(tcwv_archive["transform"], dtype=np.float64)[:6],
                    str(_scalar(tcwv_archive["crs"])),
                    name="tcwv",
                )
            tcwv20 = _align(tcwv_native, template20)
            tcwv_source = "matching public-GCS L2A WVP_60m"
        elif str(args.tcwv_missing_policy) == "error":
            raise FileNotFoundError(
                f"{matchup_id}: no captured matching public-GCS L2A WVP archive"
            )
        else:
            tcwv_source = (
                "current reduced-teacher context (L1C B8A/B09 CIBR); public-GCS L2A WVP unavailable"
            )
    elif str(args.tcwv_mode) == "scene_median":
        tcwv20 = _scene_constant(tcwv20)
    tco3_20 = context20[2]
    if str(args.tco3_mode) == "fixed_0p30":
        tco3_20 = _scene_constant(tco3_20, 0.30)
    atmo20 = _atmo_on_grid(
        aerosol_prior60,
        template=template20,
        elevation=context20[0],
        tcwv=tcwv20,
        tco3=tco3_20,
    )
    atmo60 = _atmo_on_grid(
        aerosol_prior60,
        template=template60,
        elevation=context20[0],
        tcwv=tcwv20,
        tco3=tco3_20,
    )
    observation20 = ObservationBundle(
        toa=toa20,
        geometry=geometry20,
        cloud_mask=~xr.DataArray(
            retrieval_valid20, dims=template20.dims, coords=template20.coords
        ).rio.write_crs(crs),
        sensor_config=sensor,
        metadata={"scene_key": matchup_id, "observation_time": when},
        crs=crs,
        bounds=bounds,
    )

    observation60 = ObservationBundle(
        toa=toa60,
        geometry=geometry60,
        cloud_mask=cloud60,
        sensor_config=sensor,
        metadata={"scene_key": matchup_id, "observation_time": when},
        crs=crs,
        bounds=bounds,
    )
    solve_band_columns = tuple(REDUCED_BAND_COLUMNS[name] for name in solve_bands)
    prior60_base = SurfacePrior(
        boa=_grid(
            _block_mean(np.moveaxis(t0[..., solve_band_columns], -1, 0)),
            transform60,
            crs,
            name="boa",
        ).assign_coords(band=list(solve_bands)),
        boa_unc=_grid(
            _block_mean(np.moveaxis(t0_unc[..., solve_band_columns], -1, 0)),
            transform60,
            crs,
            name="boa_unc",
        ).assign_coords(band=list(solve_bands)),
        kernels=None,
        mask=xr.DataArray(
            _block_mean(
                np.all(np.isfinite(t0[..., solve_band_columns]), axis=-1).astype(np.float32)
            )
            >= 0.999,
            dims=template60.dims,
            coords=template60.coords,
        ).rio.write_crs(crs),
        rt_space=RTSpace(backend="sixs", aerosol=str(contract["aerosol_species"])),
    )
    backend = _backend(args, platform=platform, when=when)
    finite_prior_aot = np.asarray(atmo60.aot.values)[np.isfinite(atmo60.aot.values)]
    if finite_prior_aot.size == 0:
        raise RuntimeError(f"{matchup_id}: max(QA-MAIAC, CAMS) prior contains no finite AOD")
    anchor_aot = float(np.median(finite_prior_aot))
    comp60 = _block_mean(comp)
    recurrent_snow = snow_decision["classification"] == "current_snow_recurrent_library"
    # A recurrent-snow scene aggregates with per-composite anchor weights, which a
    # pooled fit has no axis for. Fall back for that scene alone rather than
    # refusing the whole run, and record what was actually used.
    effective_predictor_model = "extra_trees_20" if recurrent_snow else str(args.predictor_model)
    predicted60 = seasonal_extra_tree_prior(
        prior60_base,
        observation60,
        seasonal_composites=comp60,
        epsg=dictionary_epsg,
        transform=transform60,
        anchor_aot=anchor_aot,
        anchor_aot_field=atmo60.aot,
        atmo_prior=atmo60,
        rt_model=backend,
        composite_band_names=dictionary_band_names,
        target_band_columns={name: dictionary_columns[name] for name in solve_bands},
        uncertainty_floor=float(args.uncertainty_floor),
        relative_uncertainty_floor=float(args.uncertainty_relative_floor),
        predictor_model=effective_predictor_model,
        robust_clip=float(args.seasonal_robust_clip),
        attach_tau_predictor=True,
        scene_mean_geometry=str(args.anchor_geometry) == "scene_mean",
        preserve_recurrent_snow_fraction=(
            float(args.snow_preserve_recurrence_fraction) if recurrent_snow else 0.0
        ),
        snow_ndsi_threshold=float(args.snow_ndsi_threshold),
        snow_green_threshold=float(args.snow_green_threshold),
        ensemble_aggregation="anchor_weighted" if recurrent_snow else "median",
    )
    if predicted60.tau_predictor is None:
        raise RuntimeError(f"{matchup_id}: seasonal predictor did not produce a tau payload")

    aggregation_weights = predicted60.tau_predictor.get("aggregation_weights")
    tau60 = {
        **predicted60.tau_predictor,
        "localizer_grid": np.asarray(predicted60.tau_predictor["localizer"]),
        "anchor_toa_grid": np.stack([toa60_values[name] for name in ANCHORS]),
        "anchor_sensor_bands": [sensor.get_band(name) for name in ANCHORS],
    }
    if aggregation_weights is not None:
        tau60["aggregation_weights_grid"] = np.asarray(aggregation_weights)
    prior60 = replace(predicted60, tau_predictor=tau60)
    water60: xr.DataArray | None = None
    if str(args.water_mask_mode) == "landwater2020":
        water60 = _water_mask_on_grid(
            bounds=bounds,
            crs=crs,
            template=template60,
            source=args.water_mask_source,
            cache_dir=args.water_mask_cache,
            buffer_pixels=int(args.water_mask_buffer_pixels),
        )
    solver_config = _solver_config(
        contract,
        allow_toa_above_one=bool(args.allow_toa_above_one),
        solve_bands=solve_bands,
        cost_aggregation=str(args.solver_cost_aggregation),
        backstop_uncertainty_scale=float(args.solver_backstop_uncertainty_scale),
        band_cost_weights=band_cost_weights,
        quadratic_refine=bool(args.solver_quadratic_refine),
        observation_likelihood=str(args.solver_observation_likelihood),
        student_t_dof=float(args.solver_student_t_dof),
        toa_uncertainty_floor=float(args.solver_toa_uncertainty_floor),
        rt_uncertainty_aod_scale=float(args.solver_rt_uncertainty_aod_scale),
        unresolved_floor_aot=float(args.solver_unresolved_floor_aot),
        unresolved_high_prior_aot=float(args.solver_unresolved_high_prior_aot),
        unresolved_prior_conflict_sigma=float(args.solver_unresolved_prior_conflict_sigma),
        unresolved_band_spread=float(args.solver_unresolved_band_spread),
        unresolved_forward_z=float(args.solver_unresolved_forward_z),
    )
    if provided_surface_aod:
        target_aod = float(args.surface_target_aod550)
        target_aod_uncertainty = float(args.surface_target_aod550_uncertainty)
        aod20 = xr.full_like(template20, np.float32(target_aod)).rio.write_crs(crs)
        aod_unc20 = xr.full_like(template20, np.float32(target_aod_uncertainty)).rio.write_crs(crs)
        aod_scene_mean = target_aod
        aod_scene_valid_count = int(np.count_nonzero(retrieval_valid60))
        actual_axis_nodes = int(contract["aot_axis_nodes"])
        diagnostics = {
            "surface_aot_axis_nodes": actual_axis_nodes,
            "surface_target_mode": "provided_scene_aod550",
            "surface_target_aod550": target_aod,
            "surface_target_aod550_uncertainty": target_aod_uncertainty,
            "m5_solver_executed": False,
            "surface_solution_unresolved": False,
            "surface_solution_unresolved_reasons": [],
        }
        teacher_success = True
    else:
        result = SurfaceDrivenSolver(solver_config).solve(
            xr.concat([toa60[name] for name in solve_bands], dim="band").assign_coords(
                band=list(solve_bands)
            ),
            prior60,
            geometry60,
            atmo60,
            backend,
            cloud60,
            [sensor.get_band(name) for name in solve_bands],
            water_mask=water60,
        )
        aod20 = _align(result.aot.rio.write_crs(crs), template20, nearest=True)
        aod_unc20 = _align(result.aot_unc.rio.write_crs(crs), template20, nearest=True)
        aod_scene_mean, aod_scene_valid_count = _scene_mean_aod(result.aot)
        diagnostics = dict(result.diagnostics)
        actual_axis_nodes = int(diagnostics.get("surface_aot_axis_nodes", -1))
        if actual_axis_nodes != int(contract["aot_axis_nodes"]):
            raise RuntimeError(
                f"{matchup_id}: solver returned {actual_axis_nodes} AOD nodes; "
                f"contract requires {contract['aot_axis_nodes']}"
            )
        teacher_success = bool(result.success)

    selected_models = _species_candidate_rt_models(
        rt_model=backend, config=solver_config, template=atmo60.aot
    )
    selected_backend = selected_models[0] if selected_models else backend
    anchors20 = {name: toa20[name] for name in ANCHORS}
    flat_anchor = np.stack([np.asarray(anchors20[name]) for name in ANCHORS], axis=-1).reshape(
        -1, 3
    )
    anchor_valid = np.all(np.isfinite(flat_anchor) & (flat_anchor > 0.0), axis=1)
    if str(args.surface_grid) == "20m":
        # The committed teacher predicts the surface on the 60 m aerosol grid
        # and bilinearly resamples it to 20 m, so only ~22% of the target's
        # spatial variability is genuinely sub-60 m. Predicting natively at
        # 20 m keeps the aerosol solve on its own 60 m grid (``predicted60``
        # above is untouched and still feeds M5) while giving the teacher real
        # 20 m structure. The atmospheric fields stay physically smooth, so
        # ``atmo20`` is the legitimate upsample of the same state.
        target_base, target_base_unc = _target_base_planes(
            t0, t0_unc, comp, target_bands, dictionary_columns
        )
        prior20_base = SurfacePrior(
            boa=_grid(target_base, transform20, crs, name="boa").assign_coords(
                band=list(target_bands)
            ),
            boa_unc=_grid(target_base_unc, transform20, crs, name="boa_unc").assign_coords(
                band=list(target_bands)
            ),
            kernels=None,
            # Validity stays on the solve bands: widening the surface output must
            # not shrink the teacher's footprint by demanding every extra band be
            # finite at a pixel the committed three-band contract accepted.
            mask=xr.DataArray(
                np.all(np.isfinite(t0[..., solve_band_columns]), axis=-1),
                dims=template20.dims,
                coords=template20.coords,
            ).rio.write_crs(crs),
            rt_space=RTSpace(backend="sixs", aerosol=str(contract["aerosol_species"])),
        )
        predicted20 = seasonal_extra_tree_prior(
            prior20_base,
            observation20,
            seasonal_composites=comp,
            epsg=dictionary_epsg,
            transform=transform20,
            anchor_aot=anchor_aot,
            anchor_aot_field=atmo20.aot,
            atmo_prior=atmo20,
            rt_model=backend,
            composite_band_names=dictionary_band_names,
            target_band_columns={name: dictionary_columns[name] for name in target_bands},
            uncertainty_floor=float(args.uncertainty_floor),
            relative_uncertainty_floor=float(args.uncertainty_relative_floor),
            predictor_model=effective_predictor_model,
            robust_clip=float(args.seasonal_robust_clip),
            attach_tau_predictor=True,
            scene_mean_geometry=str(args.anchor_geometry) == "scene_mean",
            preserve_recurrent_snow_fraction=(
                float(args.snow_preserve_recurrence_fraction) if recurrent_snow else 0.0
            ),
            snow_ndsi_threshold=float(args.snow_ndsi_threshold),
            snow_green_threshold=float(args.snow_green_threshold),
            ensemble_aggregation="anchor_weighted" if recurrent_snow else "median",
        )
        if predicted20.tau_predictor is None:
            raise RuntimeError(
                f"{matchup_id}: native 20 m seasonal predictor produced no tau payload"
            )
        base20 = predicted20.boa
        uncertainty20 = predicted20.boa_unc
        predictor_mask20 = np.asarray(predicted20.mask, dtype=bool)
        tau20 = {
            **predicted20.tau_predictor,
            "localizer_grid": np.asarray(predicted20.tau_predictor["localizer"]),
        }
        weights20 = predicted20.tau_predictor.get("aggregation_weights")
        if weights20 is not None:
            tau20["aggregation_weights_grid"] = np.asarray(weights20)
    else:
        base20 = _align(predicted60.boa, template20)
        uncertainty20 = _align(predicted60.boa_unc, template20)
        predictor_mask20 = (
            _align(predicted60.mask.astype(np.float32), template20, nearest=True) >= 0.999
        )
        tau20 = {
            **predicted60.tau_predictor,
            "localizer_grid": np.asarray(
                _align(predicted60.tau_predictor["localizer"], template20)
            ),
        }
        if aggregation_weights is not None:
            tau20["aggregation_weights_grid"] = np.stack(
                [
                    np.asarray(_align(aggregation_weights.isel(realization=index), template20))
                    for index in range(int(aggregation_weights.sizes["realization"]))
                ],
                axis=0,
            )

    def surface_and_anchor_at_aod(
        aod_field: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray]:
        corrected = _correct_anchor_reflectance(
            observation20,
            atmo_prior=atmo20,
            rt_model=selected_backend,
            template=template20,
            anchor_grids=anchors20,
            valid=anchor_valid,
            anchor_aot=anchor_aot,
            anchor_aot_field=aod_field,
            scene_mean_geometry=str(args.anchor_geometry) == "scene_mean",
        )
        anchor_boa = np.full((3, *template20.shape), np.nan, dtype=np.float64)
        anchor_boa.reshape(3, -1)[:, anchor_valid] = corrected.T
        return (
            predict_visible_from_tau_payload(
                np.asarray(base20),
                band_names=target_bands,
                tau_payload=tau20,
                anchor_boa=anchor_boa,
                aot=np.asarray(aod_field),
            ),
            anchor_boa,
        )

    def surface_at_aod(aod_field: xr.DataArray) -> np.ndarray:
        return surface_and_anchor_at_aod(aod_field)[0]

    optimal_surface_all, optimal_anchor_boa = surface_and_anchor_at_aod(aod20)
    output_indices = list(range(len(target_bands)))
    optimal_surface = np.moveaxis(np.asarray(optimal_surface_all)[output_indices], 0, -1).astype(
        np.float32
    )
    uncertainty = np.moveaxis(np.asarray(uncertainty20, dtype=np.float32)[output_indices], 0, -1)
    if provided_surface_aod and float(args.surface_target_aod550_uncertainty) > 0.0:
        target = float(args.surface_target_aod550)
        sigma = float(args.surface_target_aod550_uncertainty)
        low = max(0.0, target - sigma)
        high = min(float(contract["aot_max"]), target + sigma)
        if high > low:
            low_field = xr.full_like(template20, np.float32(low)).rio.write_crs(crs)
            high_field = xr.full_like(template20, np.float32(high)).rio.write_crs(crs)
            low_surface = np.asarray(surface_at_aod(low_field))[output_indices]
            high_surface = np.asarray(surface_at_aod(high_field))[output_indices]
            sensitivity_uncertainty = np.moveaxis(
                0.5 * np.abs(high_surface - low_surface), 0, -1
            ).astype(np.float32)
            uncertainty = np.sqrt(
                np.square(uncertainty) + np.square(sensitivity_uncertainty)
            ).astype(np.float32)
    valid = (
        retrieval_valid20
        & np.asarray(predictor_mask20, dtype=bool)
        & np.isfinite(aod20)
        & (np.asarray(aod20) >= 0.0)
        & (np.asarray(aod20) <= float(contract["aot_max"]))
        & np.all(np.isfinite(optimal_surface), axis=-1)
        & np.all(optimal_surface >= 0.0, axis=-1)
        & np.all(np.isfinite(uncertainty), axis=-1)
        & np.all(uncertainty > 0.0, axis=-1)
    )
    optimal_surface[~valid] = np.nan
    uncertainty[~valid] = np.nan
    optimal_anchor_boa_output = np.moveaxis(np.asarray(optimal_anchor_boa, dtype=np.float32), 0, -1)
    optimal_anchor_boa_output[~valid] = np.nan
    aod_values = np.asarray(aod20, dtype=np.float32)
    aod_uncertainty = np.asarray(aod_unc20, dtype=np.float32)
    aod_values[~valid] = np.nan
    aod_uncertainty[~valid] = np.nan
    prior_source_values = np.asarray(fused_prior60.winner, dtype=np.uint8)
    prior_source_support = prior_source_values != AOD_PRIOR_SOURCE_MISSING
    prior_source_pixels = int(np.count_nonzero(prior_source_support))
    prior_source_fraction = {
        "maiac": (
            float(np.count_nonzero(prior_source_values == AOD_PRIOR_SOURCE_MAIAC))
            / prior_source_pixels
            if prior_source_pixels
            else 0.0
        ),
        "cams": (
            float(np.count_nonzero(prior_source_values == AOD_PRIOR_SOURCE_CAMS))
            / prior_source_pixels
            if prior_source_pixels
            else 0.0
        ),
        "both": (
            float(np.count_nonzero(prior_source_values == AOD_PRIOR_SOURCE_BOTH))
            / prior_source_pixels
            if prior_source_pixels
            else 0.0
        ),
    }

    args.output_root.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.stem}.{os.getpid()}.partial.npz")
    np.savez_compressed(
        temporary,
        schema_version=np.asarray(output_schema),
        source_quality_sha256=np.asarray(source_quality_hash or ''),
        matchup_id=np.asarray(matchup_id),
        surface=optimal_surface,
        surface_uncertainty=uncertainty,
        surface_bands=np.asarray(target_bands),
        anchor_boa_at_solution=optimal_anchor_boa_output,
        anchor_boa_bands=np.asarray(ANCHORS),
        anchor_boa_aod_source=np.asarray(
            "caller_provided_scene_aod550" if provided_surface_aod else "m5_solution"
        ),
        solver_solve_bands_csv=np.asarray(",".join(solve_bands)),
        aod=aod_values,
        aod_uncertainty=aod_uncertainty,
        valid=valid.astype(np.uint8),
        cloud_classes=np.asarray(classes20, dtype=np.uint8),
        ocm_raw_classes=raw_class_values,
        retrieval_valid_mask=retrieval_valid20.astype(np.uint8),
        aod_axis_nodes=np.asarray(actual_axis_nodes),
        aod_scene_mean=np.asarray(aod_scene_mean, dtype=np.float32),
        aod_scene_valid_count=np.asarray(aod_scene_valid_count, dtype=np.int64),
        aod_extraction=np.asarray(aod_extraction),
        surface_target_aod_source=np.asarray(
            str(args.surface_target_aod_source) if provided_surface_aod else "m5_solution"
        ),
        surface_target_aod550=np.asarray(
            np.float32(args.surface_target_aod550) if provided_surface_aod else np.float32(np.nan)
        ),
        surface_target_aod550_uncertainty=np.asarray(
            np.float32(args.surface_target_aod550_uncertainty)
            if provided_surface_aod
            else np.float32(np.nan)
        ),
        solver_contract=np.asarray(args.solver_contract),
        solver_diagnostics_json=np.asarray(json.dumps(diagnostics, sort_keys=True)),
        transform=transform20,
        crs=np.asarray(crs),
        aerosol_species=np.asarray(str(contract["aerosol_species"])),
        fused_aod_prior=np.asarray(atmo20.aot, dtype=np.float32),
        fused_aod_prior_uncertainty=np.asarray(atmo20.aot_unc, dtype=np.float32),
        maiac_aod_prior60=np.asarray(fused_prior60.maiac_aot, dtype=np.float32),
        cams_aod_prior60=np.asarray(fused_prior60.cams_aot, dtype=np.float32),
        pixelwise_fused_aod_prior60=np.asarray(fused_prior60.state.aot, dtype=np.float32),
        pixelwise_fused_aod_prior_uncertainty60=np.asarray(
            fused_prior60.state.aot_unc, dtype=np.float32
        ),
        fused_aod_prior60=np.asarray(aerosol_prior60.aot, dtype=np.float32),
        fused_aod_prior_uncertainty60=np.asarray(aerosol_prior60.aot_unc, dtype=np.float32),
        aod_prior_source60=prior_source_values,
        aod_prior_fusion=np.asarray(str(args.aod_prior_combination)),
        aod_prior_layout=np.asarray(str(args.aod_prior_layout)),
        aod_prior_uncertainty_mode=np.asarray(str(args.aod_prior_uncertainty)),
        aod_prior_source_codes=np.asarray(
            ("0=missing", "1=MCD19-best-quality-policy", "2=CAMS", "3=both-source mean")
        ),
        maiac_prior_available=np.asarray(maiac_available),
        maiac_best_quality_qa=np.asarray(True),
        maiac_temporal_window_days=np.asarray(args.maiac_temporal_window_days),
        aod_prior_reuse_path=np.asarray(
            "" if aod_prior_reuse_path is None else str(aod_prior_reuse_path)
        ),
        aod_prior_reuse_sha256=np.asarray(aod_prior_reuse_sha256 or ""),
        deep_blue_prior_available=np.asarray(deep_blue_result.available),
        deep_blue_aod550=np.asarray(
            np.nan if deep_blue_result.aod550 is None else deep_blue_result.aod550,
            dtype=np.float32,
        ),
        deep_blue_uncertainty=np.asarray(
            np.nan if deep_blue_result.uncertainty is None else deep_blue_result.uncertainty,
            dtype=np.float32,
        ),
        deep_blue_quality=np.asarray(str(args.deep_blue_quality)),
        deep_blue_fusion=np.asarray(str(args.deep_blue_fusion)),
        deep_blue_prior_sha256=np.asarray(deep_blue_prior_sha256 or ""),
        tcwv_mode=np.asarray(str(args.tcwv_mode)),
        tcwv_missing_policy=np.asarray(str(args.tcwv_missing_policy)),
        tco3_mode=np.asarray(str(args.tco3_mode)),
        anchor_geometry=np.asarray(str(args.anchor_geometry)),
        seasonal_robust_clip=np.asarray(float(args.seasonal_robust_clip)),
        predictor_model=np.asarray(effective_predictor_model),
        ocm_thin_policy=np.asarray(str(args.ocm_thin_policy)),
        water_mask_mode=np.asarray(str(args.water_mask_mode)),
        solver_toa_upper_bound=np.asarray(
            np.nan if args.allow_toa_above_one else 1.0, dtype=np.float32
        ),
        solver_cost_aggregation=np.asarray(str(args.solver_cost_aggregation)),
        solver_band_cost_weights_csv=np.asarray(
            ""
            if band_cost_weights is None
            else ",".join(f"{value:.12g}" for value in band_cost_weights)
        ),
        solver_backstop_uncertainty_scale=np.asarray(
            float(args.solver_backstop_uncertainty_scale), dtype=np.float32
        ),
        solver_quadratic_refine=np.asarray(bool(args.solver_quadratic_refine)),
        solver_observation_likelihood=np.asarray(str(args.solver_observation_likelihood)),
        solver_student_t_dof=np.asarray(float(args.solver_student_t_dof)),
        solver_toa_uncertainty_floor=np.asarray(float(args.solver_toa_uncertainty_floor)),
        solver_rt_uncertainty_aod_scale=np.asarray(float(args.solver_rt_uncertainty_aod_scale)),
        solver_unresolved_floor_aot=np.asarray(float(args.solver_unresolved_floor_aot)),
        solver_unresolved_high_prior_aot=np.asarray(float(args.solver_unresolved_high_prior_aot)),
        solver_unresolved_prior_conflict_sigma=np.asarray(
            float(args.solver_unresolved_prior_conflict_sigma)
        ),
        solver_unresolved_band_spread=np.asarray(float(args.solver_unresolved_band_spread)),
        solver_unresolved_forward_z=np.asarray(float(args.solver_unresolved_forward_z)),
        solver_solution_unresolved=np.asarray(
            bool(diagnostics.get("surface_solution_unresolved", False))
        ),
        solver_aod_target_eligible=np.asarray(
            not bool(diagnostics.get("surface_solution_unresolved", False))
        ),
        water_mask_fraction60=np.asarray(
            0.0 if water60 is None else float(np.asarray(water60, dtype=bool).mean()),
            dtype=np.float32,
        ),
        snow_support_policy=np.asarray(str(args.snow_support_policy)),
        historical_aod_screen_mode=np.asarray(str(historical_screen["mode"])),
        historical_aod_screen_max=np.asarray(
            np.nan if historical_screen["max_aod"] is None else historical_screen["max_aod"],
            dtype=np.float32,
        ),
        historical_aod_screen_min_realizations=np.asarray(
            0
            if historical_screen["min_realizations"] is None
            else historical_screen["min_realizations"],
            dtype=np.int32,
        ),
        historical_aod_screen_index_sha256=np.asarray(historical_index_sha256 or ""),
        historical_aod_screen_active_fraction=np.asarray(
            historical_screen["active_fraction"], dtype=np.float32
        ),
        snow_policy_classification=np.asarray(str(snow_decision["classification"])),
        current_snow_fraction=np.asarray(
            np.nan
            if snow_decision.get("current_snow_fraction") is None
            else float(snow_decision["current_snow_fraction"])
        ),
        library_snow_support_ratio=np.asarray(
            np.nan
            if snow_decision.get("support_ratio") is None
            else float(snow_decision["support_ratio"])
        ),
        library_snow_support_wilson_lower=np.asarray(
            np.nan
            if snow_decision.get("support_wilson_lower") is None
            else float(snow_decision["support_wilson_lower"])
        ),
    )
    temporary.replace(output_path)
    receipt = {
        "schema_version": output_schema,
        "matchup_id": matchup_id,
        "status": "ok",
        "output": {"path": str(output_path), "sha256": _sha256(output_path)},
        "teacher": (
            "T1(tau_AERONET) from current B8A/B11/B12 corrected at an explicitly "
            "provided train/development scene AOD"
            if provided_surface_aod
            else "T1(tau*) from current B8A/B11/B12 and final cloud-masked M5 AOD"
        ),
        "solver_contract": {
            "name": str(args.solver_contract),
            **contract,
            "aerosol_resolution_m": 60.0,
            "solve_bands": list(solve_bands),
            "pool_radius_m": 600.0,
            "cost_aggregation": str(args.solver_cost_aggregation),
            "band_cost_weights": (None if band_cost_weights is None else list(band_cost_weights)),
            "backstop_uncertainty_scale": float(args.solver_backstop_uncertainty_scale),
            "quadratic_refine": bool(args.solver_quadratic_refine),
            "observation_likelihood": str(args.solver_observation_likelihood),
            "student_t_dof": float(args.solver_student_t_dof),
            "toa_uncertainty_floor": float(args.solver_toa_uncertainty_floor),
            "rt_uncertainty_aod_scale": float(args.solver_rt_uncertainty_aod_scale),
            "unresolved_thresholds": {
                "floor_aot": float(args.solver_unresolved_floor_aot),
                "high_prior_aot": float(args.solver_unresolved_high_prior_aot),
                "prior_conflict_sigma": float(args.solver_unresolved_prior_conflict_sigma),
                "band_spread": float(args.solver_unresolved_band_spread),
                "forward_z": float(args.solver_unresolved_forward_z),
            },
            "aod_target_eligible": not bool(diagnostics.get("surface_solution_unresolved", False)),
            "retrieval_extraction": (
                "caller-provided scene AOD550; M5 not executed"
                if provided_surface_aod
                else "mean of finite solved 60 m AOD pixels"
            ),
            "scene_mean_aod": aod_scene_mean,
            "scene_valid_count": aod_scene_valid_count,
            "diagnostics": diagnostics,
        },
        "snow_support_policy": {
            "mode": str(args.snow_support_policy),
            "decision": snow_decision,
            "recurrent_snow_predictor": {
                "active": recurrent_snow,
                "temporal_clip": float(args.seasonal_robust_clip),
                "preserve_recurrence_fraction": (
                    float(args.snow_preserve_recurrence_fraction) if recurrent_snow else None
                ),
                "aggregation": "anchor_weighted" if recurrent_snow else "median",
            },
        },
        "historical_aod_screen": {
            **historical_screen,
            "index_path": None if historical_index_path is None else str(historical_index_path),
            "index_sha256": historical_index_sha256,
            "background_reduction": {
                "temporal_clip_mad": 1.5,
                "uncertainty_floor": 0.006,
            },
            "seasonal_predictor_temporal_clip_mad": float(args.seasonal_robust_clip),
        },
        "surface_uncertainty_model": {
            "absolute_floor": float(args.uncertainty_floor),
            "relative_floor": float(args.uncertainty_relative_floor),
            "formula": "max(MAD*1.4826, absolute_floor, relative_floor * |prediction|)",
        },
        "surface_grid": str(args.surface_grid),
        "aeronet_role": (
            "provided explicitly by the caller as a train/development-only label; "
            "the builder performs no AERONET lookup"
            if provided_surface_aod
            else "none; not read by builder"
        ),
        "surface_target_aod": {
            "source": (
                str(args.surface_target_aod_source) if provided_surface_aod else "m5_solution"
            ),
            "aod550": (float(args.surface_target_aod550) if provided_surface_aod else None),
            "uncertainty": (
                float(args.surface_target_aod550_uncertainty) if provided_surface_aod else None
            ),
            "m5_solver_executed": not provided_surface_aod,
        },
        "aod_prior": {
            "operation": (
                "not_used_for_provided_surface_target"
                if provided_surface_aod
                else str(args.aod_prior_layout)
            ),
            "source_fusion": (
                "none; explicit train/development target AOD"
                if provided_surface_aod
                else f"pixelwise {args.aod_prior_combination}(QA-best MCD19A2, CAMS)"
            ),
            "solver_field": (
                "constant explicitly provided scene AOD; M5 not executed"
                if provided_surface_aod
                else "pixelwise fused field with native uncertainty"
                if str(args.aod_prior_layout) == "spatial"
                else (
                    "constant max(scene-median MAIAC, scene-median CAMS), with "
                    "sigma=max(0.5*AOD,0.05)"
                )
            ),
            "maiac": {
                "available": maiac_available,
                "product": "MCD19A2",
                "best_quality_qa": True,
                "granule_local_loose_qa_fallback_when_no_best_pixels": True,
                "temporal_window_days": args.maiac_temporal_window_days,
                "nearest_valid_orbit_per_pixel": True,
                "allow_default_prior": False,
                "cache_dir": str(args.maiac_cache),
                "error": maiac_error,
            },
            "cams": {
                "available": not provided_surface_aod,
                "current_date": True,
                "temporal_interpolation": True,
                "allow_default_prior": False,
                "path": str(args.cams),
            },
            "winner_fraction_on_60m_solver_grid": prior_source_fraction,
            "reused_from": (
                None
                if aod_prior_reuse_path is None
                else {
                    "path": str(aod_prior_reuse_path),
                    "sha256": aod_prior_reuse_sha256,
                }
            ),
            "source_codes": {
                "0": "missing",
                "1": "MCD19A2 best-quality policy",
                "2": "CAMS",
                "3": "both-source mean",
            },
            "modis_deep_blue": {
                "active": deep_blue_prior_path is not None,
                "available": deep_blue_result.available,
                "quality": str(args.deep_blue_quality),
                "fusion": str(args.deep_blue_fusion),
                "aod550": deep_blue_result.aod550,
                "uncertainty": deep_blue_result.uncertainty,
                "selection": deep_blue_result.metadata,
                "product": "MOD04_L2/MYD04_L2 Collection 6.1 dedicated Deep Blue land SDS",
                "aeronet_role": "none",
            },
        },
        "toa_reflectance_domain": {
            "minimum": 0.0,
            "maximum": None,
            "upper_clipping": False,
            "solver_validity_upper_bound": (None if args.allow_toa_above_one else 1.0),
        },
        "input_contract": {
            "source_quality": source_quality_record,
            "aod_prior_layout": str(args.aod_prior_layout),
            "aod_prior_combination": str(args.aod_prior_combination),
            "aod_prior_uncertainty": str(args.aod_prior_uncertainty),
            "deep_blue_quality": str(args.deep_blue_quality),
            "deep_blue_fusion": str(args.deep_blue_fusion),
            "maiac_temporal_window_days": int(args.maiac_temporal_window_days),
            "tcwv_mode": str(args.tcwv_mode),
            "tcwv_missing_policy": str(args.tcwv_missing_policy),
            "tcwv_source": tcwv_source,
            "tcwv_path": None if tcwv_path is None else str(tcwv_path),
            "tcwv_sha256": None if tcwv_path is None else _sha256(tcwv_path),
            "tco3_mode": str(args.tco3_mode),
            "anchor_geometry": str(args.anchor_geometry),
            "solver_geometry": "scene_mean",
            "seasonal_robust_clip": float(args.seasonal_robust_clip),
            "predictor_model": effective_predictor_model,
            "ocm_thin_policy": str(args.ocm_thin_policy),
            "ocm_grid_resolution_m": 20,
            "water_mask_mode": str(args.water_mask_mode),
            "water_mask_buffer_pixels_native": int(args.water_mask_buffer_pixels),
            "allow_toa_above_one": bool(args.allow_toa_above_one),
            "solver_solve_bands": list(solve_bands),
            "solver_cost_aggregation": str(args.solver_cost_aggregation),
            "solver_band_cost_weights": (
                None if band_cost_weights is None else list(band_cost_weights)
            ),
            "solver_backstop_uncertainty_scale": float(args.solver_backstop_uncertainty_scale),
            "solver_quadratic_refine": bool(args.solver_quadratic_refine),
            "solver_observation_likelihood": str(args.solver_observation_likelihood),
            "solver_student_t_dof": float(args.solver_student_t_dof),
            "solver_toa_uncertainty_floor": float(args.solver_toa_uncertainty_floor),
            "solver_rt_uncertainty_aod_scale": float(args.solver_rt_uncertainty_aod_scale),
            "solver_unresolved_floor_aot": float(args.solver_unresolved_floor_aot),
            "solver_unresolved_high_prior_aot": float(args.solver_unresolved_high_prior_aot),
            "solver_unresolved_prior_conflict_sigma": float(
                args.solver_unresolved_prior_conflict_sigma
            ),
            "solver_unresolved_band_spread": float(args.solver_unresolved_band_spread),
            "solver_unresolved_forward_z": float(args.solver_unresolved_forward_z),
            "historical_aod_screen_mode": str(historical_screen["mode"]),
            "historical_aod_screen_max": historical_screen["max_aod"],
            "historical_aod_screen_min_realizations": historical_screen["min_realizations"],
        },
        "cloud_contract": (
            "OCM native clear+thin-cloud pixels are retrieval-usable; only thick cloud, "
            "shadow, and missing OCM inputs are masked; conservative 20m-to-60m "
            "all-usable aggregation"
            if str(args.ocm_thin_policy) == "retain"
            else (
                "Only OCM native clear pixels are retrieval-usable; thin/thick cloud, "
                "shadow, and missing inputs are masked; conservative 20m-to-60m "
                "all-usable aggregation"
            )
        ),
        "water_mask_fraction60": (
            0.0 if water60 is None else float(np.asarray(water60, dtype=bool).mean())
        ),
        "valid_pixels": int(valid.sum()),
        "valid_fraction": float(valid.mean()),
        "retrieval_usable_fraction": float(retrieval_valid20.mean()),
        "ocm_clear_fraction": float((ocm_inputs_valid & (raw_class_values == 0)).mean()),
        "ocm_thin_cloud_fraction": float((ocm_inputs_valid & (raw_class_values == 2)).mean()),
        "ocm_thick_cloud_fraction": float((ocm_inputs_valid & (raw_class_values == 1)).mean()),
        "ocm_shadow_fraction": float((ocm_inputs_valid & (raw_class_values == 3)).mean()),
        "ocm_missing_fraction": float((~ocm_inputs_valid).mean()),
        "m5_executed": not provided_surface_aod,
        "m5_success": None if provided_surface_aod else teacher_success,
        "teacher_success": teacher_success,
        "elapsed_s": round(time.perf_counter() - started, 3),
        "inputs": {
            name: {"path": str(path), "sha256": _sha256(path)}
            for name, path in {
                "current": current_path,
                "fine20": fine_path,
                "dictionary": dictionary_path,
                "reduced_t0": reduced_path,
                **(
                    {"historical_aod_screen_index": historical_index_path}
                    if historical_index_path is not None
                    else {}
                ),
                **(
                    {"current_aod_prior_reuse": aod_prior_reuse_path}
                    if aod_prior_reuse_path is not None
                    else {}
                ),
                **(
                    {"modis_deep_blue_prior": deep_blue_prior_path}
                    if deep_blue_prior_path is not None
                    else {}
                ),
            }.items()
        },
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "matchup_id": matchup_id,
        "status": "ok",
        "valid_fraction": receipt["valid_fraction"],
        "elapsed_s": receipt["elapsed_s"],
        "output": str(output_path),
    }


def parser() -> argparse.ArgumentParser:
    """Build the teacher CLI.

    Exposed separately from ``main`` so the committed defaults are testable
    without executing a build.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_ids", nargs="+")
    parser.add_argument("--current-root", type=Path, default=CAMPAIGN / "current_l1c_multiscale")
    parser.add_argument("--fine20-root", type=Path, default=CAMPAIGN / "current_l1c_fine20_visible")
    parser.add_argument('--source-quality-root', type=Path,
                        help='Optional exact-grid instrument QA, applied before OCM/M5 and recorded by hash; no TOA clipping.')
    parser.add_argument(
        "--solver-contract",
        choices=tuple(TEACHER_CONTRACTS),
        default=COMMITTED_C0_CONTRACT,
    )
    parser.add_argument(
        "--solver-solve-bands",
        type=_parse_solver_solve_bands,
        default=VISIBLE,
        help=(
            "Comma-separated surface-driven solve bands. B02,B03,B04 is the frozen "
            "committed contract; B01,B02,B03,B04 is an experimental deep-blue arm."
        ),
    )
    parser.add_argument(
        "--surface-target-bands",
        type=_parse_surface_target_bands,
        default=VISIBLE,
        help=(
            "Comma-separated bands the teacher predicts and writes. B02,B03,B04 is the "
            "committed three-band output; the full land set needs a dictionary that "
            "carries those bands. The aerosol solve is unaffected."
        ),
    )
    parser.add_argument(
        "--solver-cost-aggregation",
        choices=("sum", "median_scaled"),
        default="sum",
        help=(
            "Per-node visible-band cost aggregation. sum is the frozen path; "
            "median_scaled is a development-only robust-band experiment."
        ),
    )
    parser.add_argument(
        "--solver-band-cost-weights",
        type=_parse_solver_band_cost_weights,
        help=(
            "Optional comma-separated relative cost weights aligned with "
            "--solver-solve-bands. They are normalized to unit mean inside the "
            "solver; omission preserves equal weights."
        ),
    )
    parser.add_argument(
        "--solver-backstop-uncertainty-scale",
        type=float,
        default=1.0,
        help=(
            "Experimental multiplier for the aerosol-prior backstop sigma; 1.0 "
            "preserves the frozen solver."
        ),
    )
    parser.add_argument(
        "--solver-quadratic-refine",
        action="store_true",
        help=(
            "Refine only well-conditioned interior AOD-node minima with a bounded "
            "three-point parabola; boundary solutions are never extrapolated."
        ),
    )
    parser.add_argument(
        "--solver-observation-likelihood",
        choices=("inverse_chi2", "forward_student_t"),
        default="inverse_chi2",
    )
    parser.add_argument("--solver-student-t-dof", type=float, default=4.0)
    parser.add_argument("--solver-toa-uncertainty-floor", type=float, default=0.003)
    parser.add_argument(
        "--solver-rt-uncertainty-aod-scale",
        type=float,
        default=0.0,
        help=(
            "BOA uncertainty inflation per unit AOD at 550 nm; scales inversely "
            "with wavelength and never clips signed/negative BOA."
        ),
    )
    parser.add_argument("--solver-unresolved-floor-aot", type=float, default=0.05)
    parser.add_argument("--solver-unresolved-high-prior-aot", type=float, default=0.40)
    parser.add_argument("--solver-unresolved-prior-conflict-sigma", type=float, default=1.5)
    parser.add_argument("--solver-unresolved-band-spread", type=float, default=0.50)
    parser.add_argument("--solver-unresolved-forward-z", type=float, default=3.0)
    parser.add_argument(
        "--dictionary-root", type=Path, default=CAMPAIGN / "l1c_teacher20_cci_exact"
    )
    parser.add_argument(
        "--reduced-root", type=Path, default=CAMPAIGN / "l1c_teacher20_cci_exact_reduced"
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=(CAMPAIGN / "optimal_m5_teacher20_cci_exact_committed_c0"),
    )
    parser.add_argument(
        "--surface-target-aod550",
        type=float,
        help=(
            "Build the current-anchor surface at this explicitly provided scene AOD550 "
            "instead of running M5. This mode is for site-disjoint train/development "
            "labels only and must never be used for a locked holdout target."
        ),
    )
    parser.add_argument(
        "--surface-target-aod550-uncertainty",
        type=float,
        default=0.02,
        help=(
            "One-sigma uncertainty of the provided scene AOD. The target surface "
            "uncertainty includes the surface change at AOD +/- this value."
        ),
    )
    parser.add_argument(
        "--surface-target-aod-source",
        default="caller_provided_train_development_aeronet",
        help="Auditable provenance label for --surface-target-aod550.",
    )
    parser.add_argument("--cams", type=Path, default=DEFAULT_CAMS)
    parser.add_argument("--cams-cache", type=Path, default=ROOT / "cache/cams_optimal_m5_teacher")
    parser.add_argument("--maiac-cache", type=Path, default=ROOT / "cache/maiac_day_aod")
    parser.add_argument("--maiac-temporal-window-days", type=int, default=2)
    parser.add_argument("--maiac-max-granules", type=int, default=8)
    parser.add_argument(
        "--aod-prior-reuse-root",
        type=Path,
        help=(
            "Reuse audited current-scene MAIAC/CAMS prior fields from this teacher root. "
            "The full fusion/layout/uncertainty/time-window contract must match."
        ),
    )
    parser.add_argument(
        "--deep-blue-prior-root",
        type=Path,
        help=(
            "Optional staged MOD04_L2/MYD04_L2 Collection 6.1 Deep Blue JSON root. "
            "Selection in those files is QA/time/distance based and does not use AERONET."
        ),
    )
    parser.add_argument(
        "--deep-blue-quality",
        choices=("best", "qa2", "qa1"),
        default="best",
    )
    parser.add_argument(
        "--deep-blue-fusion",
        choices=("replace", "mean", "max", "db_weight_0p25", "db_weight_0p75"),
        default="mean",
        help="Fusion of an available scalar Deep Blue AOD with the current AOD prior field.",
    )
    parser.add_argument(
        "--aod-prior-layout",
        choices=("spatial", "scalar_component_max"),
        default="spatial",
        help="AOD field supplied to M5; scalar_component_max approximates frozen MAXCAMS.",
    )
    parser.add_argument(
        "--aod-prior-combination",
        choices=("max", "mean", "plateau"),
        default="max",
        help=(
            "Pixelwise MAIAC/CAMS fusion; max is the committed C0 rule. "
            "plateau switches by regime (min when clean, mean when moderate, "
            "max when heavy), which scores 81.2% within EE against max's 59.8% "
            "on 95,619 AERONET matchups."
        ),
    )
    parser.add_argument(
        "--aod-prior-uncertainty",
        choices=("provider_primary", "half_aod_floor_0p05"),
        default="provider_primary",
        help=(
            "Uncertainty for a spatially fused AOD field. provider_primary preserves "
            "the current behaviour; half_aod_floor_0p05 isolates the locked C0 rule."
        ),
    )
    parser.add_argument(
        "--tcwv-mode",
        choices=("current_context", "scene_median", "l2a_wvp_gcs"),
        default="current_context",
    )
    parser.add_argument("--tcwv-root", type=Path)
    parser.add_argument(
        "--tcwv-missing-policy",
        choices=("error", "current_context"),
        default="error",
        help="Policy when an exact public-GCS L2A WVP archive is unavailable.",
    )
    parser.add_argument(
        "--tco3-mode",
        choices=("current_context", "fixed_0p30"),
        default="current_context",
    )
    parser.add_argument(
        "--anchor-geometry",
        choices=("scene_mean", "native"),
        default="scene_mean",
    )
    parser.add_argument(
        "--predictor-model",
        choices=("extra_trees_20", "extra_trees_20_pooled", "extra_tree"),
        default="extra_trees_20",
        help=(
            "How the seasonal library is fitted. The default trains one model per "
            "composite and takes a median across them, which bounds each member by "
            "one composite's values and inherits the seasonal middle's dark bias. "
            "extra_trees_20_pooled trains a single model over every composite at "
            "once; it cannot be combined with a recurrent-snow policy, whose "
            "anchor weighting needs a per-composite axis."
        ),
    )
    parser.add_argument(
        "--seasonal-robust-clip",
        type=float,
        default=0.0,
        help=(
            "MAD clip applied to the ET20 realization-training stack. The committed "
            "C0 predictor uses 0 (unclipped); the separately prepared background "
            "prior may still use its validated 1.5-MAD reduction."
        ),
    )
    parser.add_argument(
        "--historical-aod-screen-index-root",
        type=Path,
        help=(
            "Winner-index root used to map each corrected historical pixel to its "
            "acquisition-day AOD. Must be the exact index recorded by the dictionary."
        ),
    )
    parser.add_argument(
        "--historical-aod-screen-max",
        type=float,
        help=(
            "Use only complete historical realizations at or below this AOD where "
            "enough are available; otherwise retain the original history."
        ),
    )
    parser.add_argument(
        "--historical-aod-screen-min-realizations",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--ocm-thin-policy",
        choices=("retain", "mask"),
        default="retain",
    )
    parser.add_argument(
        "--water-mask-mode",
        choices=("none", "landwater2020"),
        default="none",
    )
    parser.add_argument("--water-mask-source")
    parser.add_argument(
        "--water-mask-cache",
        type=Path,
        default=ROOT / "cache/water-mask",
    )
    parser.add_argument("--water-mask-buffer-pixels", type=int, default=32)
    parser.add_argument(
        "--allow-toa-above-one",
        action="store_true",
        help=(
            "Retain finite positive surface-driven solve pixels when scaled L1C "
            "TOA exceeds 1.0; the captured TOA is never clipped."
        ),
    )
    parser.add_argument("--sixs-module", type=Path, default=DEFAULT_MODULE)
    parser.add_argument("--sixs-threads", type=int, default=4)
    parser.add_argument("--mie-cache", type=Path)
    parser.add_argument("--mie-module", type=Path, default=DEFAULT_MIE_MODULE)
    parser.add_argument("--run-cache", type=Path)
    parser.add_argument(
        "--snow-support-policy",
        choices=("off", "recurrent-library"),
        default="off",
        help=(
            "Reject unsupported current-snow labels and use a recurrent-snow-aware "
            "predictor for supported scenes. Off preserves the frozen v5 builder."
        ),
    )
    parser.add_argument("--snow-ndsi-threshold", type=float, default=0.40)
    parser.add_argument("--snow-green-threshold", type=float, default=0.20)
    parser.add_argument("--current-scene-snow-fraction", type=float, default=0.10)
    parser.add_argument("--snow-realization-coverage", type=float, default=0.50)
    parser.add_argument("--snow-recurrent-lower-bound", type=float, default=0.20)
    parser.add_argument("--snow-minimum-history-realizations", type=int, default=3)
    parser.add_argument("--snow-preserve-recurrence-fraction", type=float, default=0.20)
    parser.add_argument(
        "--uncertainty-floor",
        type=float,
        default=0.006,
        help=(
            "Absolute floor on the seasonal-predictor surface uncertainty. The "
            "committed teacher uses 0.006; lower it only together with "
            "--uncertainty-relative-floor, which supplies the brightness scale."
        ),
    )
    parser.add_argument(
        "--uncertainty-relative-floor",
        type=float,
        default=0.0,
        help=(
            "Fraction of predicted reflectance used as an additional, "
            "brightness-proportional uncertainty floor. Zero reproduces the "
            "committed absolute-floor teacher exactly."
        ),
    )
    parser.add_argument(
        "--surface-grid",
        choices=("60m", "20m"),
        default="60m",
        help=(
            "Grid on which the seasonal surface prior is predicted. '60m' is "
            "the committed behaviour: the prior is block-meaned to the aerosol "
            "grid and the result bilinearly resampled back to 20 m. '20m' "
            "predicts natively at 20 m against the upsampled atmospheric "
            "fields, so the target carries genuine sub-60 m structure."
        ),
    )
    parser.add_argument("--force", action="store_true")
    return parser


def validate_args(args: argparse.Namespace, cli: argparse.ArgumentParser | None = None) -> None:
    """Reject CLI combinations the builder cannot honour.

    Split out of ``main`` so the guards are exercised directly by tests.
    """

    fail = cli.error if cli is not None else _fail
    if args.uncertainty_floor < 0.0:
        fail("--uncertainty-floor must be >= 0")
    if args.uncertainty_relative_floor < 0.0:
        fail("--uncertainty-relative-floor must be >= 0")


def _fail(message: str) -> None:
    raise SystemExit(message)


def main() -> None:
    cli = parser()
    args = cli.parse_args()
    validate_args(args, cli)
    if args.seasonal_robust_clip < 0.0:
        cli.error("--seasonal-robust-clip must be >= 0")
    if args.surface_target_aod550 is not None and (
        not np.isfinite(args.surface_target_aod550)
        or args.surface_target_aod550 < 0.0
        or args.surface_target_aod550 > float(_teacher_contract(args.solver_contract)["aot_max"])
    ):
        cli.error("--surface-target-aod550 must be finite and inside the RT AOD axis")
    if (
        not np.isfinite(args.surface_target_aod550_uncertainty)
        or args.surface_target_aod550_uncertainty < 0.0
    ):
        cli.error("--surface-target-aod550-uncertainty must be finite and >= 0")
    if (
        not np.isfinite(args.solver_backstop_uncertainty_scale)
        or args.solver_backstop_uncertainty_scale <= 0.0
    ):
        cli.error("--solver-backstop-uncertainty-scale must be finite and > 0")
    positive_solver_values = {
        "--solver-student-t-dof": args.solver_student_t_dof,
        "--solver-unresolved-prior-conflict-sigma": (args.solver_unresolved_prior_conflict_sigma),
        "--solver-unresolved-forward-z": args.solver_unresolved_forward_z,
    }
    for option, value in positive_solver_values.items():
        if not np.isfinite(value) or value <= 0.0:
            cli.error(f"{option} must be finite and > 0")
    nonnegative_solver_values = {
        "--solver-toa-uncertainty-floor": args.solver_toa_uncertainty_floor,
        "--solver-rt-uncertainty-aod-scale": args.solver_rt_uncertainty_aod_scale,
        "--solver-unresolved-floor-aot": args.solver_unresolved_floor_aot,
        "--solver-unresolved-high-prior-aot": args.solver_unresolved_high_prior_aot,
        "--solver-unresolved-band-spread": args.solver_unresolved_band_spread,
    }
    for option, value in nonnegative_solver_values.items():
        if not np.isfinite(value) or value < 0.0:
            cli.error(f"{option} must be finite and >= 0")
    if args.historical_aod_screen_max is not None:
        if not np.isfinite(args.historical_aod_screen_max) or args.historical_aod_screen_max <= 0:
            cli.error("--historical-aod-screen-max must be finite and > 0")
        if args.historical_aod_screen_index_root is None:
            cli.error("--historical-aod-screen-max requires --historical-aod-screen-index-root")
    elif args.historical_aod_screen_index_root is not None:
        cli.error("--historical-aod-screen-index-root requires --historical-aod-screen-max")
    if args.historical_aod_screen_min_realizations < 1:
        cli.error("--historical-aod-screen-min-realizations must be >= 1")
    if args.water_mask_buffer_pixels < 0:
        cli.error("--water-mask-buffer-pixels must be >= 0")
    if args.aod_prior_layout == "scalar_component_max" and args.aod_prior_combination != "max":
        cli.error(
            f"--aod-prior-combination={args.aod_prior_combination} "
            "requires --aod-prior-layout=spatial"
        )
    if args.mie_cache is None:
        args.mie_cache = (
            DEFAULT_EXACT_MIE_CACHE
            if args.solver_contract == COMMITTED_C0_CONTRACT
            else DEFAULT_MIE_CACHE
        )
    if args.run_cache is None:
        args.run_cache = (
            DEFAULT_EXACT_RUN_CACHE
            if args.solver_contract == COMMITTED_C0_CONTRACT
            else DEFAULT_RUN_CACHE
        )
    args.cams_cache.mkdir(parents=True, exist_ok=True)
    args.maiac_cache.mkdir(parents=True, exist_ok=True)
    args.water_mask_cache.mkdir(parents=True, exist_ok=True)
    args.run_cache.mkdir(parents=True, exist_ok=True)
    failures = 0
    for matchup_id in args.matchup_ids:
        # One bad scene must not abort the rest of the shard's list: a shard
        # carries many scenes and losing the tail costs hours of recompute.
        try:
            record = build_one(matchup_id, args)
        except Exception as exc:  # noqa: BLE001 - reported per scene, not raised
            failures += 1
            traceback.print_exc()
            record = {
                "matchup_id": matchup_id,
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        print(json.dumps(record, sort_keys=True), flush=True)
    # Exit non-zero so a shard that lost scenes does not report COMPLETED.
    if failures:
        raise SystemExit(f"{failures} of {len(args.matchup_ids)} scenes failed")


if __name__ == "__main__":
    main()

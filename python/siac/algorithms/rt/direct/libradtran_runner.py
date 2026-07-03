"""Scene-scoped grid runner driving libRadtran ``uvspec`` for RT coefficients.

``uvspec`` is far too slow to invoke per pixel, so this runner mirrors 6S's
``scene_lut`` mode: it builds a small grid over the scene's atmosphere ranges
(AOT, TCWV) at the scene-mean geometry, runs ``uvspec`` once per grid node per
surface albedo (the two-albedo method), and assembles the four dense-spectral
terms ``TOA_rho1/2`` + ``Eg_rho1/2`` into an in-memory dataset in the *exact*
schema of the remote libRadtran LUT. That dataset is then handed to a thin
:class:`~siac.algorithms.rt.lut.backend.ZarrLUTBackend` subclass, so the entire
RSRF-convolution + ``derive_standard_rt_coefficients`` path is reused verbatim
and the output matches the remote LUT's code path.

Coefficient algebra note: ``derive_standard_rt_coefficients`` is *scale
invariant* in the ``Eg`` terms at each wavelength (a per-λ constant cancels in
``S``, ``path_ref`` and ``T_total``), so ``Eg_rho`` only needs to be
proportional to the surface global irradiance. ``TOA_rho`` must be a true
reflectance — obtained directly via uvspec ``output_quantity reflectivity``.
"""

from __future__ import annotations

import hashlib
import logging
import math
import os
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.algorithms.rt._run_cache import (
    load_cache_entry,
    resolve_run_cache_dir,
    store_cache_entry,
)
from siac.algorithms.rt.direct._build_common import prepend_conda_lib_path
from siac.algorithms.rt.direct.libradtran_build import LibRadtranPaths, ensure_libradtran
from siac.algorithms.rt.lut.backend import ZarrLUTBackend
from siac.algorithms.rt.lut.constants import TCO3_ATMCM_TO_DU, TCWV_CM_TO_LUT_MM
from siac.runtime import AtmosphericState, GeometryAngles

if TYPE_CHECKING:
    from datetime import datetime

    from siac.config.algorithms import LibRadtranAlgorithmConfig
    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime import RTCoefficients

logger = logging.getLogger(__name__)

#: LUT dimension order shared with the remote spectral LUT.
_LUT_DIMS = ("sza", "vza", "raa", "ozone", "altitude", "aot", "tcwv", "wavelength")
#: Relative paths inside the libRadtran ``data`` dir for the fixed preset.
_US_STANDARD_ATMOSPHERE = "atmmod/afglus.dat"
_KURUDZ_SOLAR = "solar_flux/kurudz_1.0nm.dat"
#: Hard ceiling on concurrent ``uvspec`` processes when ``native_threads`` is
#: unset. Each ``reptran fine`` uvspec over the full window is multi-GB, so the
#: default must NOT scale with core count (a 12-core box would otherwise launch
#: 12 multi-GB processes and exhaust RAM). The reference LUT generation ran ONE
#: uvspec per cluster task with 63 GB reserved, parallelising across nodes — so
#: a single box should run very few at once. Concurrency is further bounded by
#: available memory; see ``_resolve_worker_count``.
_DEFAULT_MAX_WORKERS = 3
#: Fraction of *available* RAM the worker pool may target.
_MEMORY_BUDGET_FRACTION = 0.7
#: Per-``uvspec`` peak-RAM model by reptran resolution: a base for the
#: band-model tables plus a strong per-window slope (the 1 nm output grid means
#: DISORT holds large per-wavelength arrays, and ``fine`` has far more internal
#: representative wavelengths). Calibrated from measured peaks at 16 streams over
#: a 1900 nm window: coarse 3.9 GB, medium 11.5 GB, **fine 35-45 GB** (fine was
#: 6.7 GB over 100 nm -> slope ~22 GB/1000 nm). Deliberately conservative
#: (estimates >= measured) so the memory cap errs toward fewer workers. NOTE:
#: a single fine run over the full S2 range needs ~45-55 GB, so it fits at most
#: ONE-at-a-time on a 64 GB box and not at all on <=32 GB.
_UVSPEC_BASE_GB = {"fine": 5.0, "medium": 5.0, "coarse": 2.0}
_UVSPEC_GB_PER_1000NM = {"fine": 22.0, "medium": 4.0, "coarse": 1.5}

#: Recommended deep H2O absorption bands (lo_nm, hi_nm) to run at a finer
#: resolution via ``mol_abs_regions``. Boundaries sit in strong absorption away
#: from S2 surface-band centres, so a per-region scheme like
#: ``[(lo, hi, "reptran fine") for lo, hi in DEEP_WATER_H2O_BANDS_NM]`` over a
#: ``"reptran"`` base does not split any S2 band's response support.
DEEP_WATER_H2O_BANDS_NM: tuple[tuple[float, float], ...] = (
    (915.0, 985.0),  # 0.94 um (lower bound in the B08 end ~906 / B09 start ~932 gap)
    (1085.0, 1180.0),  # 1.13 um
    (1330.0, 1460.0),  # 1.38 um (contains the unused B10 cirrus band)
    (1790.0, 1970.0),  # 1.88 um
)


def _cgroup_mem_limit_gb() -> float | None:
    """Memory limit (GB) from a cgroup v2/v1 controller, or None if unlimited/absent.

    ``psutil.virtual_memory().available`` reports HOST physical RAM and ignores
    container/SLURM cgroup limits, so on a constrained task it would over-size the
    worker pool (a cause of the >100 GB blow-up on a memory-limited node). Reading
    the actual cgroup ceiling lets the budget respect it.
    """
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            text = Path(path).read_text().strip()
        except OSError:
            continue
        if not text or text == "max":
            continue
        try:
            value = int(text)
        except ValueError:
            continue
        # cgroup v1 uses a huge sentinel (~2**63) for "unlimited"; ignore it.
        if 0 < value < (1 << 62):
            return value / 1.0e9
    return None


def _slurm_mem_gb() -> float | None:
    """Memory allocation (GB) from SLURM env (``--mem`` / ``--mem-per-cpu``), or None."""
    per_node = os.environ.get("SLURM_MEM_PER_NODE")
    if per_node and per_node.isdigit():
        return float(per_node) / 1024.0  # SLURM reports MB
    per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    cpus = os.environ.get("SLURM_CPUS_ON_NODE")
    if per_cpu and per_cpu.isdigit() and cpus and cpus.isdigit():
        return float(per_cpu) * float(cpus) / 1024.0
    return None


class _InMemorySpectralLUTBackend(ZarrLUTBackend):
    """A :class:`ZarrLUTBackend` fed an already-in-memory spectral dataset.

    Overriding :meth:`_load_lut` (the only method that touches ``lut_path``)
    means the whole spectral interpolation + coefficient-derivation path is
    inherited with no zarr/network IO.
    """

    def __init__(
        self,
        dataset: xr.Dataset,
        *,
        interpolation_method: str = "linear",
        rt_setup: Any | None = None,
    ) -> None:
        super().__init__(
            "<in-memory-libradtran>",
            interpolation_method=interpolation_method,
            rt_setup=rt_setup,
            # The on-disk scene-subset cache keys on ``lut_path`` (here a constant
            # placeholder) + grid-snapped scene coords. It exists only to avoid
            # re-fetching the REMOTE LUT's HTTP chunk bytes - an in-memory dataset
            # has no network cost, and the placeholder is NOT content-identifying,
            # so two different in-memory scene LUTs (e.g. coarse vs fine, or
            # different scenes that snap to the same key) would COLLIDE and serve
            # each other's stale subset - silently freezing the aot interpolation.
            # Disable it: the within-process ``_spectral_scene_subset`` cache still
            # serves every band of this scene, and the expensive uvspec runs are
            # already cached at the run-cache layer.
            scene_cache_enabled=False,
        )
        self._dataset = dataset

    def _load_lut(self) -> None:
        self._lut = self._dataset
        for dim, coord in self._lut.coords.items():
            self._lut_coords[dim] = coord.values


def build_uvspec_deck(
    *,
    data_dir: str,
    aot550: float,
    tcwv_cm: float,
    tco3_atmcm: float,
    elevation_km: float,
    sza_deg: float,
    vza_deg: float,
    raa_deg: float,
    albedo: float,
    wavelength_min_nm: float,
    wavelength_max_nm: float,
    mol_abs_param: str,
    number_of_streams: int,
    aerosol_species: str = "continental_average",
) -> str:
    """Build a uvspec input deck for one grid node + surface albedo.

    Units are converted to libRadtran conventions: TCWV cm -> mm (``MM``), TCO3
    atm-cm -> Dobson (``DU``), VZA -> ``umu = cos(VZA)`` (>0 = upward/satellite
    view). ``output_quantity reflectivity`` makes ``uu`` come out as TOA
    reflectance; ``eglo`` at the surface level is the global-irradiance term.
    """
    umu = math.cos(math.radians(vza_deg))
    # umu must be strictly within (0, 1]; clamp away from the grazing limit.
    umu = min(1.0, max(1.0e-6, umu))
    lines = [
        f"data_files_path {data_dir}",
        f"atmosphere_file {Path(data_dir) / _US_STANDARD_ATMOSPHERE}",
        f"source solar {Path(data_dir) / _KURUDZ_SOLAR}",
        f"mol_abs_param {mol_abs_param}",
        "rte_solver disort",
        f"number_of_streams {int(number_of_streams)}",
        "aerosol_default",
        f"aerosol_species_file {aerosol_species}",
        f"aerosol_set_tau_at_wvl 550 {aot550:.6f}",
        f"mol_modify H2O {tcwv_cm * TCWV_CM_TO_LUT_MM:.6f} MM",
        f"mol_modify O3 {tco3_atmcm * TCO3_ATMCM_TO_DU:.6f} DU",
        f"sza {sza_deg:.4f}",
        "phi0 0",
        f"phi {raa_deg:.4f}",
        f"umu {umu:.6f}",
        f"altitude {max(0.0, elevation_km):.4f}",
        "zout sur toa",
        f"albedo {albedo:.4f}",
        f"wavelength {wavelength_min_nm:.1f} {wavelength_max_nm:.1f}",
        "output_quantity reflectivity",
        "output_user lambda eglo uu",
        "quiet",
    ]
    return "\n".join(lines) + "\n"


def parse_uvspec_table(stdout: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse ``output_user lambda eglo uu`` with ``zout sur toa``.

    uvspec emits two rows per wavelength (surface then TOA, in ``zout`` order),
    each ``lambda eglo uu``. Returns ``(wavelength_nm, eg_surface, toa_uu)``
    where ``eg_surface`` is the surface global-irradiance term (``Eg_rho``) and
    ``toa_uu`` is the TOA reflectance (``TOA_rho``).
    """
    rows: list[list[float]] = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
        except ValueError:
            continue
    if not rows or len(rows) % 2 != 0:
        raise ValueError(
            f"Unexpected uvspec output: parsed {len(rows)} numeric rows "
            "(expected an even count = 2 per wavelength)."
        )
    arr = np.asarray(rows, dtype=np.float64).reshape(-1, 2, 3)
    # Each pair shares a wavelength; row 0 = surface (zout sur), row 1 = toa.
    wavelength = arr[:, 0, 0]
    eg_surface = arr[:, 0, 1]
    toa_uu = arr[:, 1, 2]
    return wavelength, eg_surface, toa_uu


def _fold_raa_deg(value: float) -> float:
    raa = abs(value) % 360.0
    return 360.0 - raa if raa > 180.0 else raa


def _scene_mean_deg(values: xr.DataArray) -> float:
    return float(np.degrees(np.nanmean(np.asarray(values.values, dtype=np.float64))))


def _scene_mean(values: xr.DataArray) -> float:
    return float(np.nanmean(np.asarray(values.values, dtype=np.float64)))


def _axis_from_scene(values: xr.DataArray, n_nodes: int) -> np.ndarray:
    """Build a strictly-increasing axis spanning the scene's finite range."""
    finite = np.asarray(values.values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot build a libRadtran grid axis from an all-NaN scene field.")
    lo = float(np.min(finite))
    hi = float(np.max(finite))
    if not hi > lo:
        pad = max(abs(lo) * 0.1, 1.0e-3)
        lo, hi = lo - pad, hi + pad
    # AOT and TCWV are physically non-negative; never emit a negative grid node
    # (uvspec rejects a negative aerosol optical depth / column water amount).
    if lo < 0.0:
        lo, hi = 0.0, max(hi, 1.0e-3)
    axis = np.linspace(lo, hi, max(2, int(n_nodes)), dtype=np.float64)
    # Quantize for run-cache key stability (the scene min/max inherit any
    # cross-process jitter in upstream resampling; see _DECK_AXIS_DECIMALS).
    axis = _round_axis(axis)
    if axis.size == 1:
        # A scene range narrower than the quantum collapsed to one node; keep
        # this function's >=2-node contract with a one-quantum-wide axis.
        axis = np.array([axis[0], axis[0] + 10.0**-_DECK_AXIS_DECIMALS], dtype=np.float64)
    return cast("np.ndarray", axis)


#: Node ceilings for the joint grid-search LUT. The solver's aot axis is
#: ``grid_search_aot_points`` nodes (~40); each node is one uvspec run per albedo,
#: so we thin to ~16 (off-node candidates are linearly interpolated - the same
#: approximation the scene-LUT path already makes in geometry). tcwv collapses to
#: one node when the stage fixes tcwv, so a small ceiling is plenty.
_JOINT_LUT_MAX_AOT_NODES = 16
_JOINT_LUT_MAX_TCWV_NODES = 4

#: Quantization (decimal places) for every continuous scalar embedded in a
#: uvspec deck. The persistent run cache keys on a SHA-256 of the deck text, so
#: any input that is not bit-reproducible across processes busts all keys at
#: once. The scene means are NOT bit-reproducible: on T33KWP the solver-grid
#: sza mean was observed to flip between 37.575276... and 37.575361... deg
#: across identical re-runs (~9e-5 deg apart, straddling the deck's 4-dp
#: format boundary), recomputing all 32 joint-LUT uvspec runs (~14 min) on
#: every retrieval. Rounding to these quanta is far above that jitter and far
#: below RT sensitivity (1e-3 deg in sza shifts sec(sza) by ~2e-5 relative;
#: 1e-4 in aot shifts TOA reflectance by ~1e-5), so deck text - and therefore
#: the cache key - is stable across re-runs of the same scene.
_DECK_ANGLE_DECIMALS = 3
_DECK_AXIS_DECIMALS = 4
_DECK_TCO3_ATMCM_DECIMALS = 5
_DECK_ELEVATION_KM_DECIMALS = 4


def _round_axis(axis: np.ndarray) -> np.ndarray:
    """Quantize axis values for deck/cache-key stability (see _DECK_*_DECIMALS).

    ``np.unique`` collapses any nodes that the rounding made identical, so the
    returned axis is strictly increasing (duplicate coords would break the
    in-memory LUT's interpolation).
    """
    return cast(
        "np.ndarray",
        np.unique(np.round(np.asarray(axis, dtype=np.float64), _DECK_AXIS_DECIMALS)),
    )


def _subsample_axis(axis: np.ndarray, max_nodes: int) -> np.ndarray:
    """Unique, sorted axis thinned to <= ``max_nodes`` nodes (endpoints kept).

    Axis values are quantized first (``_round_axis``) so the downstream deck
    text - and the persistent run-cache key derived from it - cannot vary with
    sub-1e-4 jitter in the solver's candidate derivation.
    """
    axis = _round_axis(axis)
    if axis.size <= max_nodes or axis.size == 0:
        return axis
    idx = np.unique(np.linspace(0, axis.size - 1, max_nodes).round().astype(int))
    return cast("np.ndarray", axis[idx])


def _scalar_cell(value: float) -> xr.DataArray:
    """A 1x1 ``(y, x)`` field holding one value (single-point scene inputs)."""
    return xr.DataArray(np.full((1, 1), float(value), dtype=np.float64), dims=("y", "x"))


def _nanmean_value(values: xr.DataArray) -> float:
    """Scene-mean of a field as a plain float (NaN when nothing is finite)."""
    return float(np.nanmean(np.asarray(values.values, dtype=np.float64)))


class _LibRadtranJointGridSearchLUT:
    """Joint (aot x tcwv) grid-search LUT backed by one in-memory spectral LUT.

    Built once by :meth:`LibRadtranRunner.build_joint_grid_search_lut`; each
    block-grid-search candidate is served by :meth:`evaluate`, which interpolates
    the prebuilt scene LUT instead of launching a fresh ``uvspec`` build. This is
    the libRadtran analogue of the 6S ``JointGridSearchLUT`` amortisation - it is
    what makes a full per-pixel AOD retrieval tractable (otherwise every one of
    the ~40 aot candidates rebuilds the multi-minute scene LUT).
    """

    def __init__(
        self,
        *,
        backend: _InMemorySpectralLUTBackend,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        node_count: int,
    ) -> None:
        self._backend = backend
        self._bands = bands
        finite = int(np.count_nonzero(np.isfinite(np.asarray(atmo_state.aot.values))))
        self.plan = SimpleNamespace(lut_case_count=int(node_count), direct_case_count=finite)
        # The in-memory scene LUT collapses geometry/ozone/altitude to the
        # scene mean (singleton axes), so for a spatially-uniform (aot, tcwv)
        # candidate EVERY pixel maps to the same LUT query. Evaluate once on a
        # 1x1 scene at the scene means and broadcast the three scalars, instead
        # of paying the full per-pixel interpolation pipeline for
        # ~grid_pixels x bands x candidates identical points.
        self._template_shape = tuple(np.asarray(atmo_state.aot.values).shape)
        self._template_dims = tuple(atmo_state.aot.dims)
        # saa=0 / vaa=mean(raa) reproduces the scene-mean relative azimuth the
        # full-field path would compute (raa = |vaa - saa| % 2pi, and the mean
        # of per-pixel raa already lies in [0, 2pi)).
        self._point_geometry = GeometryAngles(
            sza=_scalar_cell(_nanmean_value(geometry.sza)),
            saa=_scalar_cell(0.0),
            vza=_scalar_cell(_nanmean_value(geometry.vza)),
            vaa=_scalar_cell(_nanmean_value(geometry.raa)),
        )
        zero = _scalar_cell(0.0)
        self._point_state = {
            "tco3": _scalar_cell(_nanmean_value(atmo_state.tco3)),
            "elevation": _scalar_cell(_nanmean_value(atmo_state.elevation)),
            "aot_unc": zero,
            "tcwv_unc": zero,
            "tco3_unc": zero,
        }

    def evaluate(self, aot_val: float, tcwv_val: float) -> list[dict[str, xr.DataArray]]:
        """Per-band xap/xbp/xcp for one (aot, tcwv) candidate via interpolation."""
        atmo = AtmosphericState(
            aot=_scalar_cell(float(aot_val)),
            tcwv=_scalar_cell(float(tcwv_val)),
            **self._point_state,
        )
        outputs: list[dict[str, xr.DataArray]] = []
        for band in self._bands:
            coeffs = self._backend.compute_coefficients(self._point_geometry, atmo, band)
            outputs.append(
                {
                    name: xr.DataArray(
                        np.full(
                            self._template_shape,
                            float(np.asarray(values.values).reshape(-1)[0]),
                            dtype=np.float32,
                        ),
                        dims=self._template_dims,
                    )
                    for name, values in (
                        ("xap", coeffs.xap),
                        ("xbp", coeffs.xbp),
                        ("xcp", coeffs.xcp),
                    )
                }
            )
        return outputs


class LibRadtranRunner:
    """Build a scene-scoped libRadtran spectral LUT and derive RT coefficients."""

    def __init__(
        self,
        *,
        libradtran_config: LibRadtranAlgorithmConfig,
        sensor_config: SensorConfig | None = None,
        rt_setup: Any | None = None,
    ) -> None:
        self._config = libradtran_config
        self._sensor_config = sensor_config
        self._rt_setup = rt_setup
        self._paths: LibRadtranPaths | None = None
        self._observation_time: datetime | None = None
        self._scene_key: tuple[float, ...] | None = None
        self._scene_backend: _InMemorySpectralLUTBackend | None = None
        self._straddle_warned: set[str] = set()
        # Serializes scene-LUT construction. The solver evaluates bands on worker
        # threads that can all miss the empty cache at once and each launch a full
        # (~5-9 GB) uvspec build of the SAME scene (signature ignores the band) -
        # 3-4 concurrent builds blew past RAM and were jetsam-killed. The lock
        # collapses that race to a single build; the rest reuse the warm cache.
        self._build_lock = threading.Lock()
        # Persistent on-disk cache of uvspec run outputs (deck text -> spectra).
        # A hit skips the uvspec subprocess entirely, so re-running a scene - or
        # reusing nodes across the retrieval/diagnostics/correction builds - costs
        # no radiative-transfer time. ``None`` disables it.
        self._run_cache_dir: Path | None = resolve_run_cache_dir(
            getattr(libradtran_config, "run_cache_dir", None),
            subpath="libradtran/uvspec_run_cache",
            enabled=bool(getattr(libradtran_config, "run_cache_enabled", True)),
        )
        self._cache_stats_lock = threading.Lock()
        self._run_cache_hits = 0
        self._run_cache_misses = 0

    def set_observation_time(self, observation_time: datetime | None) -> None:
        self._observation_time = observation_time

    def _segments(self) -> list[tuple[float, float, str]]:
        """Partition the output window into contiguous (lo, hi, model) segments.

        ``mol_abs_param`` is the base; each region overrides it over its range.
        Regions come from explicit ``mol_abs_regions`` or, when those are unset
        and ``adaptive_deep_water_fine`` is on (default), from the deep H2O bands
        (``DEEP_WATER_H2O_BANDS_NM``) run at ``"reptran fine"`` - the adaptive
        scheme that keeps each ``uvspec`` process small. Returns segments
        covering [wmin, wmax] with no gaps/overlap.
        """
        cfg = self._config
        wmin, wmax = float(cfg.wavelength_min_nm), float(cfg.wavelength_max_nm)
        base = str(cfg.mol_abs_param)
        regions: tuple[tuple[float, float, str], ...] | None = cfg.mol_abs_regions
        if not regions and getattr(cfg, "adaptive_deep_water_fine", False):
            regions = tuple(
                (lo, hi, "reptran fine")
                for lo, hi in DEEP_WATER_H2O_BANDS_NM
                if min(hi, wmax) > max(lo, wmin)
            )
        if not regions:
            return [(wmin, wmax, base)]
        segs: list[tuple[float, float, str]] = []
        cur = wmin
        for lo_r, hi_r, model in sorted(regions, key=lambda r: float(r[0])):
            lo, hi = max(float(lo_r), wmin), min(float(hi_r), wmax)
            if hi <= cur:  # region outside the window or already covered
                continue
            if lo > cur:
                segs.append((cur, lo, base))
            segs.append((lo, hi, str(model)))
            cur = hi
        if cur < wmax:
            segs.append((cur, wmax, base))
        return segs

    def _warn_band_straddle(self, band: SensorBand) -> None:
        """Warn once if a band's response support crosses a band-model seam."""
        segments = self._segments()
        if len(segments) < 2:  # single segment -> no interior seams to straddle
            return
        name = str(getattr(band, "name", "?"))
        if name in self._straddle_warned:
            return
        self._straddle_warned.add(name)
        wl = getattr(band, "rsrf_wavelengths_nm", None)
        resp = getattr(band, "rsrf_response", None)
        if wl is not None and resp is not None and len(wl):
            wl_arr = np.asarray(wl, dtype=np.float64)
            resp_arr = np.asarray(resp, dtype=np.float64)
            # Material support (>1% of peak) — negligible RSRF tails crossing a
            # seam contribute <1% of the band integral, so they don't warrant a
            # warning.
            support = wl_arr[resp_arr > 0.01 * float(np.max(resp_arr))]
            lo, hi = float(np.min(support)), float(np.max(support))
        else:
            c, bw = float(band.center_wavelength), float(band.bandwidth or 0.0)
            lo, hi = c - bw, c + bw
        for seg_lo, _seg_hi, _m in segments[1:]:  # interior boundaries = each segment start
            if lo < seg_lo < hi:
                logger.warning(
                    "libradtran: band %s response support [%.0f, %.0f] nm straddles a "
                    "mol_abs_regions boundary at %.0f nm; its band integral mixes two band "
                    "models. Move the region boundary outside this band's support.",
                    name,
                    lo,
                    hi,
                    seg_lo,
                )
                return

    def _ensure_paths(self) -> LibRadtranPaths:
        if self._paths is None:
            self._paths = ensure_libradtran(self._config)
        return self._paths

    def compute_coefficients(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
    ) -> RTCoefficients:
        self._warn_band_straddle(band)
        backend = self._scene_backend_for(geometry, atmo_state)
        return backend.compute_coefficients(geometry, atmo_state, band)

    def build_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
    ) -> _LibRadtranJointGridSearchLUT:
        """Build one scene LUT spanning the grid-search range, served by interp.

        The block-grid-search calls the RT model with many ``(aot, tcwv)``
        candidates. Each ``compute_coefficients`` call would otherwise rebuild
        the multi-minute ``uvspec`` scene LUT (its cache signature varies with
        the per-candidate aot/tcwv range). Here we build a single in-memory
        spectral LUT whose aot/tcwv nodes ARE the grid-search axes - so lookups
        are exact at the nodes and within-range elsewhere - then serve every
        candidate by interpolation. Mirrors the 6S joint-LUT amortisation.
        """
        for band in bands:
            self._warn_band_straddle(band)
        aot_nodes = _subsample_axis(aot_axis, _JOINT_LUT_MAX_AOT_NODES)
        tcwv_nodes = _subsample_axis(tcwv_axis, _JOINT_LUT_MAX_TCWV_NODES)
        # Full-precision trace of the node axes: any change here (vs a prior
        # run of the same scene) flips every uvspec deck's run-cache key, so a
        # mysterious "0 from cache" re-run is diagnosable from DEBUG logs.
        logger.debug(
            "libradtran: joint grid-search nodes aot=%s tcwv=%s "
            "(thinned from %d aot / %d tcwv solver candidates)",
            np.array2string(aot_nodes, precision=17, max_line_width=10_000),
            np.array2string(tcwv_nodes, precision=17, max_line_width=10_000),
            int(np.asarray(aot_axis).size),
            int(np.asarray(tcwv_axis).size),
        )
        with self._build_lock:
            # Cache as the scene backend so a same-scene compute_coefficients call
            # reuses it (the joint LUT spans the full search range, so prior-range
            # values fall inside it) rather than launching another build.
            backend = self._build_and_cache_backend(
                geometry,
                atmo_state,
                aot_axis_override=aot_nodes,
                tcwv_axis_override=tcwv_nodes,
            )
        return _LibRadtranJointGridSearchLUT(
            backend=backend,
            geometry=geometry,
            atmo_state=atmo_state,
            bands=list(bands),
            node_count=int(aot_nodes.size * max(1, tcwv_nodes.size)),
        )

    def _scene_backend_for(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> _InMemorySpectralLUTBackend:
        key = self._scene_signature(geometry, atmo_state)
        if self._scene_backend is not None and key == self._scene_key:
            return self._scene_backend
        with self._build_lock:
            # Re-check under the lock: a concurrent caller may have just built
            # this exact scene while we waited, so reuse it rather than launch a
            # second redundant (and memory-doubling) uvspec build.
            if self._scene_backend is not None and key == self._scene_key:
                return self._scene_backend
            return self._build_and_cache_backend(geometry, atmo_state)

    def _build_and_cache_backend(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        *,
        aot_axis_override: np.ndarray | None = None,
        tcwv_axis_override: np.ndarray | None = None,
    ) -> _InMemorySpectralLUTBackend:
        """Build the scene LUT + backend and publish it as the cached scene.

        Callers must hold ``_build_lock`` (this is the single site that
        constructs the in-memory backend and assigns the scene-cache slot, so
        the two cannot drift apart).
        """
        dataset = self._build_scene_lut(
            geometry,
            atmo_state,
            aot_axis_override=aot_axis_override,
            tcwv_axis_override=tcwv_axis_override,
        )
        backend = _InMemorySpectralLUTBackend(
            dataset,
            interpolation_method=str(self._config.interpolation_method),
            rt_setup=self._rt_setup,
        )
        self._scene_backend = backend
        self._scene_key = self._scene_signature(geometry, atmo_state)
        return backend

    def _scene_signature(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> tuple[float, ...]:
        return (
            round(_scene_mean_deg(geometry.sza), 3),
            round(_scene_mean_deg(geometry.vza), 3),
            round(_fold_raa_deg(_scene_mean_deg(geometry.raa)), 3),
            round(_scene_mean(atmo_state.tco3), 6),
            round(_scene_mean(atmo_state.elevation), 4),
            round(float(np.nanmin(atmo_state.aot.values)), 4),
            round(float(np.nanmax(atmo_state.aot.values)), 4),
            round(float(np.nanmin(atmo_state.tcwv.values)), 4),
            round(float(np.nanmax(atmo_state.tcwv.values)), 4),
        )

    def _build_scene_lut(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        *,
        aot_axis_override: np.ndarray | None = None,
        tcwv_axis_override: np.ndarray | None = None,
    ) -> xr.Dataset:
        paths = self._ensure_paths()
        cfg = self._config
        data_dir = str(paths.data_dir)

        # Geometry + slow-varying gases collapse to the scene mean (one node):
        # the spectral-LUT consumer selects geometry by nearest-to-scene-mean
        # and averages ozone/altitude over their range anyway. Every scalar is
        # quantized (see _DECK_*_DECIMALS) because it is embedded in the deck
        # text that keys the persistent run cache: upstream resampling is not
        # bit-reproducible across processes (observed ~9e-5 deg flips in the
        # T33KWP solver-grid sza mean), and an unrounded mean recomputes every
        # uvspec run on every retrieval of the same scene.
        sza_deg = round(_scene_mean_deg(geometry.sza), _DECK_ANGLE_DECIMALS)
        vza_deg = round(_scene_mean_deg(geometry.vza), _DECK_ANGLE_DECIMALS)
        raa_deg = round(_fold_raa_deg(_scene_mean_deg(geometry.raa)), _DECK_ANGLE_DECIMALS)
        tco3_atmcm = round(_scene_mean(atmo_state.tco3), _DECK_TCO3_ATMCM_DECIMALS)
        ozone_du = tco3_atmcm * TCO3_ATMCM_TO_DU
        altitude_km = round(
            max(0.0, _scene_mean(atmo_state.elevation)), _DECK_ELEVATION_KM_DECIMALS
        )

        n_nodes = min(
            int(cfg.scene_lut_max_nodes_per_axis),
            max(2, int(math.isqrt(int(cfg.scene_lut_max_cases) // 2))),
        )
        # Explicit overrides let the joint grid-search LUT span the solver's full
        # (aot, tcwv) search range with asymmetric node counts (many aot nodes,
        # one fixed-tcwv node) instead of the per-scene prior range. Overrides
        # are re-quantized here (idempotent for the joint path, which rounds in
        # _subsample_axis) so no caller can reintroduce jittery deck values.
        aot_axis = (
            _round_axis(aot_axis_override)
            if aot_axis_override is not None
            else _axis_from_scene(atmo_state.aot, n_nodes)
        )
        tcwv_axis = (
            _round_axis(tcwv_axis_override)
            if tcwv_axis_override is not None
            else _axis_from_scene(atmo_state.tcwv, n_nodes)
        )
        # Full-precision trace of every deck-embedded scalar (the run-cache key
        # inputs); compare across runs when diagnosing unexpected cache misses.
        logger.debug(
            "libradtran: deck scalars sza=%r vza=%r raa=%r tco3_atmcm=%r altitude_km=%r "
            "aot_axis=%s tcwv_axis=%s",
            sza_deg,
            vza_deg,
            raa_deg,
            tco3_atmcm,
            altitude_km,
            np.array2string(aot_axis, precision=17, max_line_width=10_000),
            np.array2string(tcwv_axis, precision=17, max_line_width=10_000),
        )

        # Endpoint-safe 1 nm (or wavelength_step_nm) grid: linspace pins both ends
        # exactly, so a non-integer step can't overshoot wavelength_max_nm (which
        # np.arange(stop=max+step) could, leaving an out-of-window cell).
        wmin_grid = float(cfg.wavelength_min_nm)
        wmax_grid = float(cfg.wavelength_max_nm)
        step_grid = float(cfg.wavelength_step_nm)
        n_pts = max(2, int(round((wmax_grid - wmin_grid) / step_grid)) + 1)
        wl_grid: np.ndarray = np.linspace(wmin_grid, wmax_grid, n_pts, dtype=np.float64)

        # Partition the window into (lo, hi, model) segments and run one uvspec
        # per (aot, tcwv, albedo, segment), stitching the 1 nm pieces. With no
        # mol_abs_regions this is a single full-window segment (= prior behaviour).
        segments = self._segments()
        seg_masks: list[np.ndarray] = []
        for j, (lo, hi, _m) in enumerate(segments):
            last = j == len(segments) - 1
            seg_masks.append((wl_grid >= lo) & (wl_grid <= hi if last else wl_grid < hi))

        jobs: list[tuple[int, int, int, float, int]] = []
        for i_aot in range(len(aot_axis)):
            for i_tcwv in range(len(tcwv_axis)):
                for i_rho, rho in enumerate((cfg.rho1, cfg.rho2)):
                    for seg_idx in range(len(segments)):
                        jobs.append((i_aot, i_tcwv, i_rho, float(rho), seg_idx))

        def _run_job(
            job: tuple[int, int, int, float, int],
        ) -> tuple[int, int, int, int, np.ndarray, np.ndarray]:
            i_aot, i_tcwv, i_rho, rho, seg_idx = job
            seg_lo, seg_hi, seg_model = segments[seg_idx]
            deck = build_uvspec_deck(
                data_dir=data_dir,
                aot550=float(aot_axis[i_aot]),
                tcwv_cm=float(tcwv_axis[i_tcwv]),
                tco3_atmcm=tco3_atmcm,
                elevation_km=altitude_km,
                sza_deg=sza_deg,
                vza_deg=vza_deg,
                raa_deg=raa_deg,
                albedo=rho,
                wavelength_min_nm=seg_lo,
                wavelength_max_nm=seg_hi,
                mol_abs_param=seg_model,
                number_of_streams=int(cfg.number_of_streams),
                aerosol_species=str(getattr(cfg, "aerosol_species", "continental_average")),
            )
            wl, eg, toa = self._run_uvspec(paths, deck)
            target = wl_grid[seg_masks[seg_idx]]
            # uvspec/reptran emit at internal representative wavelengths that can
            # start a hair above seg_lo (and end a hair below seg_hi), so clamp at
            # the edges (np.interp default) rather than leaving NaN holes: a sub-nm
            # extrapolation at a band-model-chosen seam is far safer than a NaN that
            # would poison any band whose RSRF support reaches that seam (and with
            # the adaptive default every scene is segmented).
            eg_seg = np.interp(target, wl, eg)
            toa_seg = np.interp(target, wl, toa)
            return i_aot, i_tcwv, i_rho, seg_idx, toa_seg, eg_seg

        n_aot, n_tcwv, n_wl = len(aot_axis), len(tcwv_axis), len(wl_grid)
        toa: np.ndarray = np.full((2, n_aot, n_tcwv, n_wl), np.nan, dtype=np.float32)
        eg: np.ndarray = np.full((2, n_aot, n_tcwv, n_wl), np.nan, dtype=np.float32)
        # Concurrency is bounded by the HEAVIEST segment so a coarse+fine mix is
        # safe (a narrow fine segment is ~7 GB, not the ~45 GB of full-range fine).
        heaviest_gb = max(self._estimate_uvspec_gb_for(m, hi - lo) for lo, hi, m in segments)
        workers = self._resolve_worker_count(len(jobs), heaviest_gb)
        seg_desc = ", ".join(f"{m.split()[-1]}[{int(lo)}-{int(hi)}]" for lo, hi, m in segments)
        logger.info(
            "libradtran: scene LUT = %d uvspec run(s) [%d node(s) x %d segment(s): %s] "
            "across %d worker(s) (<=%.1f GB/process)",
            len(jobs),
            len(jobs) // max(1, len(segments)),
            len(segments),
            seg_desc,
            workers,
            heaviest_gb,
        )
        cache_hits0, cache_misses0 = self._run_cache_hits, self._run_cache_misses
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for i_aot, i_tcwv, i_rho, seg_idx, toa_vals, eg_vals in pool.map(_run_job, jobs):
                toa[i_rho, i_aot, i_tcwv][seg_masks[seg_idx]] = toa_vals
                eg[i_rho, i_aot, i_tcwv][seg_masks[seg_idx]] = eg_vals
        if self._run_cache_dir is not None:
            logger.info(
                "libradtran: scene LUT done - %d uvspec run(s): %d from cache, %d computed (cache=%s)",
                len(jobs),
                self._run_cache_hits - cache_hits0,
                self._run_cache_misses - cache_misses0,
                self._run_cache_dir,
            )

        return self._assemble_dataset(
            toa=toa,
            eg=eg,
            sza_deg=sza_deg,
            vza_deg=vza_deg,
            raa_deg=raa_deg,
            ozone_du=ozone_du,
            altitude_km=altitude_km,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            wl_grid=wl_grid,
        )

    def _assemble_dataset(
        self,
        *,
        toa: np.ndarray,
        eg: np.ndarray,
        sza_deg: float,
        vza_deg: float,
        raa_deg: float,
        ozone_du: float,
        altitude_km: float,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        wl_grid: np.ndarray,
    ) -> xr.Dataset:
        # Promote (rho, aot, tcwv, wl) -> full (sza,vza,raa,ozone,altitude,aot,tcwv,wl)
        # with singleton geometry/gas axes (the LUT consumer collapses them).
        def _cube(values: np.ndarray) -> np.ndarray:
            return values[np.newaxis, np.newaxis, np.newaxis, np.newaxis, np.newaxis, :, :, :]

        # Axis units follow the remote spectral-LUT schema (ozone in DU, tcwv in
        # MM) — the inherited ZarrLUTBackend evaluation converts the SIAC state
        # (tco3 atm-cm, tcwv cm) into these units, so the in-memory axes must
        # match or every per-pixel tcwv lookup lands ~10x too dry. ``tcwv_axis``
        # arrives in cm (state units, also used for the uvspec decks).
        coords = {
            "sza": np.array([sza_deg], dtype=np.float32),
            "vza": np.array([vza_deg], dtype=np.float32),
            "raa": np.array([raa_deg], dtype=np.float32),
            "ozone": np.array([ozone_du], dtype=np.float32),
            "altitude": np.array([altitude_km], dtype=np.float32),
            "aot": aot_axis.astype(np.float32),
            "tcwv": (tcwv_axis * TCWV_CM_TO_LUT_MM).astype(np.float32),
            "wavelength": wl_grid.astype(np.float32),
        }
        # No ``solar_irradiance`` variable: the remote libRadtran LUT does not
        # carry one, so its RSRF convolution is bandpass-only. Matching that
        # keeps this backend's band integration identical to the LUT route
        # (uvspec ``output_quantity reflectivity`` already cancels the solar
        # spectrum out of TOA_rho/Eg_rho, so no solar weighting is needed).
        data_vars: dict[str, Any] = {
            "TOA_rho1": (_LUT_DIMS, _cube(toa[0])),
            "TOA_rho2": (_LUT_DIMS, _cube(toa[1])),
            "Eg_rho1": (_LUT_DIMS, _cube(eg[0])),
            "Eg_rho2": (_LUT_DIMS, _cube(eg[1])),
        }
        return xr.Dataset(
            data_vars=data_vars,
            coords=coords,
            attrs={
                "rho1": float(self._config.rho1),
                "rho2": float(self._config.rho2),
                "generator": "libRadtran-uvspec",
                "aerosol_profile": "continental_average",
                "atmospheric_profile_shape": "us_standard",
            },
        )

    @staticmethod
    def _estimate_uvspec_gb_for(mol_abs_param: str, span_nm: float) -> float:
        """Conservative peak-RAM estimate (GB) for one ``uvspec`` run."""
        span = max(1.0, float(span_nm))
        res = mol_abs_param.lower()
        base = next(
            (gb for t, gb in _UVSPEC_BASE_GB.items() if t in res), _UVSPEC_BASE_GB["coarse"]
        )
        slope = next(
            (gb for t, gb in _UVSPEC_GB_PER_1000NM.items() if t in res),
            _UVSPEC_GB_PER_1000NM["coarse"],
        )
        return max(0.5, base + span / 1000.0 * slope)

    def _estimate_uvspec_gb(self) -> float:
        """Peak-RAM estimate for the global (single-segment) ``mol_abs_param``."""
        cfg = self._config
        return self._estimate_uvspec_gb_for(
            str(cfg.mol_abs_param),
            float(cfg.wavelength_max_nm) - float(cfg.wavelength_min_nm),
        )

    def _effective_available_gb(self) -> float | None:
        """Smallest of host-available RAM and any cgroup/SLURM memory ceiling."""
        candidates: list[float] = []
        try:
            import psutil  # type: ignore[import-untyped]

            candidates.append(float(psutil.virtual_memory().available) / 1.0e9)
        except Exception:  # noqa: BLE001 - psutil import/probe failure must not be fatal
            pass
        for limit in (_cgroup_mem_limit_gb(), _slurm_mem_gb()):
            if limit is not None and limit > 0:
                candidates.append(limit)
        return min(candidates) if candidates else None

    def _memory_worker_cap(self, per_proc_gb: float | None = None) -> int:
        """Max concurrent uvspec processes within the memory budget.

        Budget = ``min(_MEMORY_BUDGET_FRACTION x effective-available-RAM,
        memory_budget_gb)``; effective-available respects cgroup/SLURM limits. A
        single process larger than the whole budget still runs (one at a time)
        with a warning - we cannot launch fewer than one.
        """
        per_proc = self._estimate_uvspec_gb() if per_proc_gb is None else per_proc_gb
        per_proc = max(0.5, float(per_proc))
        available_gb = self._effective_available_gb()
        budget: float | None = None
        if available_gb is not None:
            budget = _MEMORY_BUDGET_FRACTION * available_gb
        cap_gb = getattr(self._config, "memory_budget_gb", None)
        if cap_gb is not None:
            budget = float(cap_gb) if budget is None else min(budget, float(cap_gb))
        if budget is None:  # no probe and no explicit cap
            return _DEFAULT_MAX_WORKERS
        if per_proc > budget:
            logger.warning(
                "libradtran: one '%s' uvspec is ~%.0f GB but the memory budget is ~%.0f GB "
                "(available ~%s GB, cap %s GB); it will run one-at-a-time and may still swap/OOM. "
                "Use a narrower window or a coarser base (adaptive deep-water-fine keeps it small).",
                self._config.mol_abs_param,
                per_proc,
                budget,
                f"{available_gb:.0f}" if available_gb is not None else "?",
                f"{float(cap_gb):.0f}" if cap_gb is not None else "unset",
            )
        return max(1, int(budget / per_proc))

    def _resolve_worker_count(self, n_jobs: int, per_proc_gb: float | None = None) -> int:
        """Bound concurrent ``uvspec`` processes by available RAM.

        Each ``reptran fine`` uvspec over the full window is multi-GB, so the
        worker count must NOT scale with core count (12 cores x multi-GB was the
        original >100 GB out-of-memory blow-up). When ``native_threads`` is unset
        the default is small (``_DEFAULT_MAX_WORKERS``) and further capped by
        memory (``per_proc_gb`` = the heaviest segment when mol_abs_regions are
        used); an explicit ``native_threads`` is honoured but still memory-capped
        (with a warning) so it cannot exhaust RAM.
        """
        cfg = self._config
        cpu = os.cpu_count() or 4
        requested = (
            int(cfg.native_threads) if cfg.native_threads else min(cpu, _DEFAULT_MAX_WORKERS)
        )
        workers = max(1, min(requested, self._memory_worker_cap(per_proc_gb), n_jobs))
        if cfg.native_threads and workers < int(cfg.native_threads):
            logger.warning(
                "libradtran: capping uvspec concurrency at %d (native_threads=%d requested); the "
                "heaviest '%s' run is ~%.1f GB and only ~%.0f%% of available RAM may be used.",
                workers,
                int(cfg.native_threads),
                cfg.mol_abs_param,
                per_proc_gb if per_proc_gb is not None else self._estimate_uvspec_gb(),
                _MEMORY_BUDGET_FRACTION * 100.0,
            )
        return workers

    def _run_uvspec(
        self,
        paths: LibRadtranPaths,
        deck: str,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Persistent run cache: the deck text fully determines uvspec's output
        # (disort is deterministic; the deck embeds the versioned data path), so a
        # SHA-256 of the deck is an exact key. A hit skips the subprocess entirely.
        cache_file: Path | None = None
        if self._run_cache_dir is not None:
            digest = hashlib.sha256(deck.encode("utf-8")).hexdigest()
            cache_file = self._run_cache_dir / f"{digest}.npz"
            cached = self._load_cached_run(cache_file)
            if cached is not None:
                self._record_cache(hit=True)
                return cached
        result = self._invoke_uvspec(paths, deck)
        if cache_file is not None:
            self._record_cache(hit=False)
            self._store_cached_run(cache_file, result)
        return result

    def _invoke_uvspec(
        self,
        paths: LibRadtranPaths,
        deck: str,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Run in a private scratch cwd: uvspec writes a ``randomseed`` artifact
        # to its working directory, and the thread pool runs several at once, so
        # each needs its own dir to avoid littering the project / racing.
        timeout_s = float(getattr(self._config, "uvspec_timeout_s", 1800.0))
        with tempfile.TemporaryDirectory(prefix="siac-uvspec-") as scratch:
            try:
                result = subprocess.run(  # noqa: S603 - fixed argv (uvspec), deck on stdin
                    [str(paths.uvspec)],
                    input=deck,
                    capture_output=True,
                    text=True,
                    check=False,
                    env=self._uvspec_env(),
                    cwd=scratch,
                    timeout=timeout_s,
                )
            except subprocess.TimeoutExpired as exc:
                # A wedged uvspec must not block a pool worker forever.
                raise RuntimeError(
                    f"uvspec timed out after {timeout_s:.0f}s (uvspec_timeout_s). The run may be "
                    "pathological; raise uvspec_timeout_s or narrow the wavelength window."
                ) from exc
        if result.returncode != 0:
            tail = result.stderr.strip()[-600:]
            raise RuntimeError(f"uvspec failed (exit {result.returncode}): {tail}")
        return parse_uvspec_table(result.stdout)

    @staticmethod
    def _load_cached_run(
        cache_file: Path,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        def _read(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            with np.load(path) as data:
                return (
                    np.asarray(data["wl"], dtype=np.float64),
                    np.asarray(data["eg"], dtype=np.float64),
                    np.asarray(data["toa"], dtype=np.float64),
                )

        cached: tuple[np.ndarray, np.ndarray, np.ndarray] | None = load_cache_entry(
            cache_file, _read
        )
        return cached

    @staticmethod
    def _store_cached_run(
        cache_file: Path,
        result: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        wl, eg, toa = result

        def _write(tmp: Path) -> None:
            # Passing a file handle to np.savez avoids its ".npz" auto-suffix
            # on the temp file.
            with tmp.open("wb") as handle:
                np.savez(handle, wl=wl, eg=eg, toa=toa)

        store_cache_entry(cache_file, _write)

    def _record_cache(self, *, hit: bool) -> None:
        with self._cache_stats_lock:
            if hit:
                self._run_cache_hits += 1
            else:
                self._run_cache_misses += 1

    def _uvspec_env(self) -> dict[str, str]:
        return prepend_conda_lib_path(dict(os.environ))

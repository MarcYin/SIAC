"""Dense spectral-LUT path of :class:`~siac.algorithms.rt.lut.backend.ZarrLUTBackend`.

Implements the SPECTRAL representation (``TOA_rho1/2`` + ``Eg_rho1/2`` terms
on a wavelength axis): grid-snapped scene subsetting with its in-memory and
on-disk caches, per-band RSRF convolution into aot/tcwv grids, point
interpolation of the band-compressed terms, and recovery of the standard
xap/xbp/xcp coefficients via the two-albedo algebra. The mixin is composed
into ``ZarrLUTBackend`` (which owns loading, dispatch, the retry wrapper, and
the shared sanitizer/interpolation helpers declared in the ``TYPE_CHECKING``
block below); the cache attributes and locks are initialised in
``ZarrLUTBackend.__init__``.
"""

from __future__ import annotations

import hashlib
import logging
import time
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.adapters.rsrf import coerce_band_rsrf
from siac.algorithms.rt._run_cache import load_cache_entry, store_cache_entry
from siac.algorithms.rt.lut._spectral_math import (
    build_point_interpolation_coords,
    build_spectral_integration_weights,
    derive_standard_rt_coefficients,
    finite_mean,
    finite_range,
    spectral_scene_cache_key,
    summarize_spectral_scene,
)
from siac.algorithms.rt.lut.constants import state_tco3_to_lut
from siac.algorithms.rt.lut.rsrf_kernel import build_aligned_rsrf_kernel
from siac.domain.spectral import RelativeSpectralResponse

if TYPE_CHECKING:
    import threading
    from collections.abc import Callable
    from typing import TypeVar

    from siac.domain import SensorBand
    from siac.runtime import AtmosphericState, GeometryAngles

    _T = TypeVar("_T")

logger = logging.getLogger(__name__)


class _SpectralLUTMixin:
    """Dense spectral TOA/Eg LUT methods for ``ZarrLUTBackend``."""

    _SPECTRAL_VARS = (
        "TOA_rho1",
        "TOA_rho2",
        "Eg_rho1",
        "Eg_rho2",
    )
    _SOLAR_IRRADIANCE_NAMES = (
        "extraterrestrial_solar_irradiance",
        "solar_irradiance",
    )
    _DEFAULT_SURFACE_RHO1 = 0.15
    _DEFAULT_SURFACE_RHO2 = 0.5

    if TYPE_CHECKING:
        # Attributes (initialised in ZarrLUTBackend.__init__) and shared
        # helpers provided by the ZarrLUTBackend facade.
        lut_path: str
        interpolation_method: str
        _lut_coords: dict[str, np.ndarray]
        _scene_subset_logged: bool
        _scene_cache_dir: Path | None
        _cache_lock: threading.Lock
        _scene_build_lock: threading.Lock
        _spectral_scene_key: tuple[float, ...] | None
        _spectral_scene_subset: xr.Dataset | None
        _spectral_band_grid_cache: dict[
            tuple[tuple[float, ...], tuple[str, float, float]],
            tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray],
        ]

        @property
        def lut(self) -> xr.Dataset: ...

        def _run_with_transient_lut_io_retry(
            self, fn: Callable[[], _T], *, operation: str
        ) -> _T: ...

        def _supports_spectral_lut(self) -> bool: ...

        @staticmethod
        def _axis_bounds(axis: np.ndarray) -> tuple[float, float]: ...

        @staticmethod
        def _require_finite_values(values: np.ndarray, *, name: str) -> np.ndarray: ...

        @staticmethod
        def _interpolate_variable_fast(
            var: xr.DataArray,
            coords: dict[str, xr.DataArray],
            method: str,
        ) -> np.ndarray: ...

        @staticmethod
        def _require_matching_grid_shapes(
            geometry: GeometryAngles,
            atmo_state: AtmosphericState,
        ) -> xr.DataArray: ...

        @staticmethod
        def _geometry_to_lut_inputs(
            geometry: GeometryAngles,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...

        @staticmethod
        def _weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray: ...

    @classmethod
    def _axis_midpoint(cls, axis: np.ndarray) -> float:
        lo, hi = cls._axis_bounds(axis)
        return (lo + hi) * 0.5

    @staticmethod
    def _finite_range(
        values: np.ndarray,
        *,
        fallback: tuple[float, float],
    ) -> tuple[float, float]:
        return cast("tuple[float, float]", finite_range(values, fallback=fallback))

    @staticmethod
    def _finite_mean(values: np.ndarray, *, fallback: float) -> float:
        return cast("float", finite_mean(values, fallback=fallback))

    def _compute_spectral_with_retry(
        self,
        *,
        sza_deg: np.ndarray,
        vza_deg: np.ndarray,
        raa_deg: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        band: SensorBand,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute spectral coefficients with retry for transient remote LUT read failures."""
        return self._run_with_transient_lut_io_retry(
            lambda: self._compute_coefficients_from_spectral_lut(
                sza_deg,
                vza_deg,
                raa_deg,
                aot,
                tcwv,
                tco3,
                elevation,
                band,
            ),
            operation="read",
        )

    def _compute_coefficients_from_spectral_lut(
        self,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        band: SensorBand,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Derive standard RT coefficients from dense spectral LUT terms."""
        if "wavelength" not in self._lut_coords:
            raise ValueError("Spectral LUT must define a wavelength coordinate")

        scene_key, lut_scene = self._resolve_spectral_scene_subset(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )

        target_shape = aot.shape
        toa_rho1, toa_rho2, eg_rho1, eg_rho2 = self._interpolate_spectral_band_terms(
            scene_key,
            lut_scene,
            band,
            aot=aot.ravel(),
            tcwv=tcwv.ravel(),
        )
        xap, xbp, xcp = self._derive_coefficients_from_spectral_terms(
            toa_rho1=toa_rho1,
            toa_rho2=toa_rho2,
            eg_rho1=eg_rho1,
            eg_rho2=eg_rho2,
        )

        return (
            xap.reshape(target_shape).astype(np.float32),
            xbp.reshape(target_shape).astype(np.float32),
            xcp.reshape(target_shape).astype(np.float32),
        )

    def _resolve_spectral_scene_subset(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[tuple[float, ...], xr.Dataset]:
        """Return the cached spectral scene subset and emit the one-time summary log."""
        scene_key, lut_scene = self._get_or_build_spectral_scene_subset(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        if self._scene_subset_logged:
            return scene_key, lut_scene

        sza_mean, vza_mean, raa_mean, tco3_bounds, elevation_bounds = summarize_spectral_scene(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        logger.info(
            (
                "Scene LUT subset: sza=%.3f deg, vza=%.3f deg, raa=%.3f deg, "
                "ozone=[%.3f, %.3f], altitude=[%.3f, %.3f]"
            ),
            sza_mean,
            vza_mean,
            raa_mean,
            tco3_bounds[0],
            tco3_bounds[1],
            elevation_bounds[0],
            elevation_bounds[1],
        )
        self._scene_subset_logged = True
        return scene_key, lut_scene

    def _interpolate_spectral_band_terms(
        self,
        scene_key: tuple[float, ...],
        lut_scene: xr.Dataset,
        band: SensorBand,
        *,
        aot: np.ndarray,
        tcwv: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Interpolate band-compressed spectral LUT terms for each retrieval point."""
        coords = self._build_aot_tcwv_point_coords(
            lut_scene,
            aot=aot,
            tcwv=tcwv,
        )
        toa_rho1_grid, toa_rho2_grid, eg_rho1_grid, eg_rho2_grid = (
            self._get_or_build_spectral_band_grids(
                scene_key,
                lut_scene,
                band,
            )
        )
        return (
            self._interpolate_variable_fast(toa_rho1_grid, coords, self.interpolation_method),
            self._interpolate_variable_fast(toa_rho2_grid, coords, self.interpolation_method),
            self._interpolate_variable_fast(eg_rho1_grid, coords, self.interpolation_method),
            self._interpolate_variable_fast(eg_rho2_grid, coords, self.interpolation_method),
        )

    def _derive_coefficients_from_spectral_terms(
        self,
        *,
        toa_rho1: np.ndarray,
        toa_rho2: np.ndarray,
        eg_rho1: np.ndarray,
        eg_rho2: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Convert dense spectral radiative terms into standard RT coefficients."""
        rho1, rho2 = self._spectral_reference_reflectances()
        return cast(
            "tuple[np.ndarray, np.ndarray, np.ndarray]",
            derive_standard_rt_coefficients(
                toa_rho1=toa_rho1,
                toa_rho2=toa_rho2,
                eg_rho1=eg_rho1,
                eg_rho2=eg_rho2,
                rho1=rho1,
                rho2=rho2,
            ),
        )

    def _spectral_scene_cache_key(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[float, ...]:
        """Build a stable scene cache key snapped to the LUT grid.

        Angle means are snapped to the nearest LUT coordinate (mirroring
        ``_subset_spectral_lut_for_scene``) so that different spatial
        resolutions of the same scene (e.g. the atmospheric grid used by
        preload vs. the coarse grid used by M3) always produce the *same*
        key and share a single download.
        """
        return cast(
            "tuple[float, ...]",
            spectral_scene_cache_key(
                sza=sza,
                vza=vza,
                raa=raa,
                tco3=tco3,
                elevation=elevation,
                snap_axes=self._lut_coords,
            ),
        )

    def _spectral_scene_summary(
        cls,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[float, float, float, tuple[float, float], tuple[float, float]]:
        """Summarize the scene geometry and optional axes for cache/logging reuse."""
        tco3_arr = np.asarray(tco3, dtype=np.float32)
        elevation_arr = np.asarray(elevation, dtype=np.float32)
        return (
            round(cls._finite_mean(sza, fallback=0.0), 3),
            round(cls._finite_mean(vza, fallback=0.0), 3),
            round(cls._finite_mean(raa, fallback=0.0), 3),
            (
                round(cls._finite_range(tco3_arr, fallback=(0.0, 0.0))[0], 3),
                round(cls._finite_range(tco3_arr, fallback=(0.0, 0.0))[1], 3),
            ),
            (
                round(cls._finite_range(elevation_arr, fallback=(0.0, 0.0))[0], 3),
                round(cls._finite_range(elevation_arr, fallback=(0.0, 0.0))[1], 3),
            ),
        )

    @staticmethod
    def _spectral_band_cache_key(band: SensorBand) -> tuple[str, float, float]:
        """Build a stable per-band cache key for spectral compression grids."""
        return (band.name, float(band.center_wavelength), float(band.bandwidth))

    def _get_or_build_spectral_scene_subset(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[tuple[float, ...], xr.Dataset]:
        """Return cached scene LUT subset or build and cache it once."""
        scene_key = self._spectral_scene_cache_key(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )

        # Fast path: already cached.
        with self._cache_lock:
            if self._spectral_scene_key == scene_key and self._spectral_scene_subset is not None:
                return scene_key, self._spectral_scene_subset

        # Slow path: hold the build lock so only one thread downloads.
        # A second thread arriving here will block until the first finishes,
        # then hit the fast-path cache check above on re-entry.
        with self._scene_build_lock:
            # Re-check after acquiring the build lock — another thread may
            # have completed the download while we waited.
            with self._cache_lock:
                if (
                    self._spectral_scene_key == scene_key
                    and self._spectral_scene_subset is not None
                ):
                    return scene_key, self._spectral_scene_subset

            # Persistent disk cache: a prior run materialised this scene subset,
            # so load it from local disk instead of re-fetching chunks over HTTP.
            subset = self._load_scene_subset_from_disk(scene_key)
            loaded_from_disk = subset is not None

            if subset is None:
                _t0 = time.perf_counter()
                subset = self._subset_spectral_lut_for_scene(
                    self.lut,
                    sza=sza,
                    vza=vza,
                    raa=raa,
                    tco3=tco3,
                    elevation=elevation,
                )
                logger.info(
                    "_subset_spectral_lut_for_scene (lazy graph) %.3f s", time.perf_counter() - _t0
                )

                # Eagerly materialise the scene subset into memory.  Without this,
                # downstream band-grid operations (.integrate, .values) each trigger
                # independent HTTP range-request round-trips to the remote Zarr store.
                # Using scheduler="threads" parallelises chunk fetches over HTTP,
                # cutting wall-clock time from ~33 s (sequential) to a few seconds.
                _t0 = time.perf_counter()
                if hasattr(subset, "compute"):
                    subset = subset.compute(scheduler="threads", num_workers=8)
                logger.info("subset.compute() (materialise) %.3f s", time.perf_counter() - _t0)

            # Single publication site (we hold the build lock, so no other
            # thread can be populating concurrently).
            with self._cache_lock:
                self._spectral_scene_key = scene_key
                self._spectral_scene_subset = subset
                self._spectral_band_grid_cache.clear()
                self._scene_subset_logged = False

        # Persist OUTSIDE the build lock: waiters only need the in-memory copy,
        # so the netCDF write must not extend their blockage (writes are
        # atomic + best-effort, see store_cache_entry).
        if not loaded_from_disk:
            self._store_scene_subset_to_disk(scene_key, subset)
        return scene_key, subset

    #: Disk-format version of the scene-subset cache. Bump whenever
    #: ``_subset_spectral_lut_for_scene`` semantics change (axis selection,
    #: averaging, unit handling): the key carries no semantic inputs beyond the
    #: scene coordinates, so without a bump old caches would silently serve
    #: subsets built by the old logic. v2: ozone range converted atm-cm -> DU
    #: before slicing (pre-v2 entries were pinned to the 200 DU edge node).
    _SCENE_SUBSET_CACHE_VERSION = "v2"

    def _scene_subset_cache_path(self, scene_key: tuple[float, ...]) -> Path | None:
        """On-disk path for this LUT + grid-snapped scene subset (or None)."""
        if self._scene_cache_dir is None:
            return None
        stem = Path(self.lut_path).name or "lut"
        # Key on LUT identity + the grid-snapped scene key, so different LUTs or
        # scenes never collide and a re-run of the same scene hits.
        raw = (
            self._SCENE_SUBSET_CACHE_VERSION
            + "|"
            + self.lut_path
            + "|"
            + ",".join(f"{v:.6g}" for v in scene_key)
        )
        digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]
        return self._scene_cache_dir / f"{stem}.{digest}.subset.nc"

    def _load_scene_subset_from_disk(self, scene_key: tuple[float, ...]) -> xr.Dataset | None:
        path = self._scene_subset_cache_path(scene_key)
        if path is None:
            return None

        def _read(entry: Path) -> xr.Dataset:
            subset = xr.open_dataset(entry).load()
            logger.info("LUT scene subset loaded from disk cache %s", entry)
            return subset

        return load_cache_entry(path, _read)

    def _store_scene_subset_to_disk(self, scene_key: tuple[float, ...], subset: xr.Dataset) -> None:
        path = self._scene_subset_cache_path(scene_key)
        if path is not None:
            store_cache_entry(path, subset.to_netcdf)

    def _get_or_build_spectral_band_grids(
        self,
        scene_key: tuple[float, ...],
        lut_scene: xr.Dataset,
        band: SensorBand,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
        """Return cached per-band spectral-compressed grids or build once."""
        band_key = self._spectral_band_cache_key(band)
        cache_key = (scene_key, band_key)

        with self._cache_lock:
            cached = self._spectral_band_grid_cache.get(cache_key)
        if cached is not None:
            return cached

        _t0 = time.perf_counter()
        lut_band = self._subset_wavelength_for_band(lut_scene, band)
        weights = self._spectral_integration_weights(band, lut_band)
        grids = (
            self._weighted_spectral_mean(lut_band["TOA_rho1"], weights),
            self._weighted_spectral_mean(lut_band["TOA_rho2"], weights),
            self._weighted_spectral_mean(lut_band["Eg_rho1"], weights),
            self._weighted_spectral_mean(lut_band["Eg_rho2"], weights),
        )
        logger.info(
            "_get_or_build_spectral_band_grids(%s) %.3f s", band.name, time.perf_counter() - _t0
        )

        with self._cache_lock:
            self._spectral_band_grid_cache[cache_key] = grids
        return grids

    def preload_scene_subset(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand] | None = None,
    ) -> None:
        """
        Preload LUT metadata and cache scene/band spectral subsets.

        Intended for background execution while other modules are running.
        """
        self._run_with_transient_lut_io_retry(
            lambda: self._preload_scene_subset_once(
                geometry=geometry,
                atmo_state=atmo_state,
                bands=bands,
            ),
            operation="preload",
        )

    def _preload_scene_subset_once(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand] | None,
    ) -> None:
        """Populate scene-level spectral LUT caches once for the requested bands."""
        # Ensure LUT metadata index is loaded exactly once.
        _ = self.lut
        if not self._supports_spectral_lut():
            return
        self._require_matching_grid_shapes(geometry, atmo_state)

        sza, vza, raa = self._geometry_to_lut_inputs(geometry)
        tco3 = atmo_state.tco3.values
        elevation = atmo_state.elevation.values

        scene_key, lut_scene = self._get_or_build_spectral_scene_subset(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        for band in bands or []:
            self._get_or_build_spectral_band_grids(scene_key, lut_scene, band)

    def _subset_spectral_lut_for_scene(
        self,
        lut: xr.Dataset,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> xr.Dataset:
        """Subset spectral LUT using scene-angle means and ozone/elevation ranges."""
        out = lut

        for dim, scene_values in {"sza": sza, "vza": vza, "raa": raa}.items():
            if dim in out.coords:
                axis = np.asarray(out.coords[dim].values, dtype=np.float32)
                finite_scene_values = self._require_finite_values(scene_values, name=dim)
                value = self._finite_mean(finite_scene_values, fallback=self._axis_midpoint(axis))
                out = out.sel({dim: value}, method="nearest")
            if dim in out.dims and out.sizes.get(dim, 0) == 1:
                out = out.squeeze(dim=dim, drop=True)

        # The spectral-LUT schema stores ozone in Dobson units while the SIAC
        # state carries tco3 in atm-cm (DU / 1000). Convert BEFORE slicing:
        # slicing the DU axis [200..600] with atm-cm values (~0.3) selected an
        # empty range, and the nearest-value fallback then silently pinned every
        # scene to the 200 DU edge node (visibly wrong Chappuis-band absorption,
        # AOT biased low). Altitude is km on both sides.
        ozone_du = state_tco3_to_lut(np.asarray(tco3, dtype=np.float64))
        for dim, scene_values in {"ozone": ozone_du, "altitude": elevation}.items():
            if dim not in out.coords:
                continue
            axis = np.asarray(out.coords[dim].values, dtype=np.float32)
            if axis.size == 0:
                continue
            axis_min, axis_max = self._axis_bounds(axis)
            vmin, vmax = self._finite_range(scene_values, fallback=(axis_min, axis_max))
            lo = float(np.clip(min(vmin, vmax), axis_min, axis_max))
            hi = float(np.clip(max(vmin, vmax), axis_min, axis_max))
            subset = out.sel({dim: slice(lo, hi)})
            if subset.sizes.get(dim, 0) == 0:
                subset = out.sel({dim: (lo + hi) * 0.5}, method="nearest")
            if dim in subset.dims and subset.sizes.get(dim, 1) > 1:
                subset = subset.mean(dim=dim)
            if dim in subset.dims and subset.sizes.get(dim, 0) == 1:
                subset = subset.squeeze(dim=dim, drop=True)
            out = subset

        return out

    @staticmethod
    def _build_aot_tcwv_point_coords(
        lut: xr.Dataset,
        *,
        aot: np.ndarray,
        tcwv: np.ndarray,
    ) -> dict[str, xr.DataArray]:
        """Build per-pixel LUT coordinates (aot unitless; tcwv cm -> LUT mm axis)."""
        # Late import: the shared sanitizer helpers live on the facade class,
        # whose module imports this mixin at load time.
        from siac.algorithms.rt.lut.backend import ZarrLUTBackend

        return cast(
            "dict[str, xr.DataArray]",
            build_point_interpolation_coords(
                lut,
                aot=aot,
                tcwv=tcwv,
                require_finite_values=ZarrLUTBackend._require_finite_values,
                sanitize_point_values=ZarrLUTBackend._sanitize_point_values,
            ),
        )

    def _subset_wavelength_for_band(
        self,
        lut: xr.Dataset,
        band: SensorBand,
        sigma_factor: float = 4.0,
    ) -> xr.Dataset:
        """Trim spectral LUT to a narrow wavelength window around the target band."""
        if "wavelength" not in lut.coords:
            return lut

        wl_axis = np.asarray(lut.coords["wavelength"].values, dtype=np.float32)
        if wl_axis.size == 0:
            return lut

        del sigma_factor
        rsrf_band = (
            self._band_rsrf(band)
            if band.has_rsrf
            else coerce_band_rsrf(
                band,
                sensor_id="UNKNOWN",
                satellite_id="UNKNOWN",
            )
        )
        kernel = build_aligned_rsrf_kernel(
            rsrf_band,
            lut_wavelengths_nm=wl_axis,
            lut_id=self.lut_path,
            support_padding=1,
        )
        return lut.isel(wavelength=slice(kernel.start_index, kernel.end_index))

    def _spectral_integration_weights(
        self, band: SensorBand, lut: xr.Dataset | None = None
    ) -> xr.DataArray:
        """Build wavelength weights for spectral convolution (bandpass * optional solar spectrum)."""
        source = lut if lut is not None else self.lut
        return build_spectral_integration_weights(
            source,
            band,
            lut_path=self.lut_path,
            solar_irradiance_names=self._SOLAR_IRRADIANCE_NAMES,
            band_rsrf=self._band_rsrf,
        )

    @staticmethod
    def _band_rsrf(band: SensorBand) -> RelativeSpectralResponse:
        """Build a canonical RSRF object from a band-carried tabulated response."""
        if not band.has_rsrf:
            raise ValueError(f"Band {band.name!r} does not define a tabulated RSRF")
        return RelativeSpectralResponse.from_tabulated(
            sensor_id="UNKNOWN",
            satellite_id="UNKNOWN",
            band_name=band.name,
            wavelengths_nm=np.asarray(band.rsrf_wavelengths_nm, dtype=np.float32),
            response=np.asarray(band.rsrf_response, dtype=np.float32),
        )

    def _spectral_reference_reflectances(self) -> tuple[float, float]:
        """Return reference surface reflectances used by dense spectral LUT variables."""
        attrs = self.lut.attrs
        rho1 = float(
            attrs.get("rho1", attrs.get("reference_reflectance_1", self._DEFAULT_SURFACE_RHO1))
        )
        rho2 = float(
            attrs.get("rho2", attrs.get("reference_reflectance_2", self._DEFAULT_SURFACE_RHO2))
        )
        return rho1, rho2

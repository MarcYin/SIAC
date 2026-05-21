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

from siac.adapters.rsrf import coerce_band_rsrf
from siac.algorithms.rt.lut._spectral_math import (
    build_point_interpolation_coords,
    build_spectral_integration_weights,
    derive_standard_rt_coefficients,
    finite_mean,
    finite_range,
    weighted_spectral_mean,
)
from siac.algorithms.rt.lut.rsrf_kernel import build_aligned_rsrf_kernel
from siac.algorithms.rt.lut.store import as_local_path, build_lut_store
from siac.domain.spectral import RelativeSpectralResponse
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.domain import SensorBand

logger = logging.getLogger(__name__)
_T = TypeVar("_T")


class ZarrLUTBackend:
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

    Args:
        lut_path: Path to the Zarr LUT store
        interpolation_method: Interpolation method ("linear" or "nearest")
        chunk_cache_size: Size of chunk cache in bytes
        rt_setup: Effective generic RT setup used to validate or describe the
            packaged LUT preset. The current public remote LUT is a fixed
            libRadtran preset, not a generic configurable RT family.
    """

    _COEFFICIENT_VARS = (
        "path_reflectance",
        "transmittance_down",
        "transmittance_up",
        "spherical_albedo",
    )
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
    _SPECTRAL_IO_MAX_RETRIES = 8
    _SUPPORTED_INTERPOLATION_METHODS = frozenset({"linear", "nearest"})

    def __init__(
        self,
        lut_path: str | Path,
        interpolation_method: str = "linear",
        chunk_cache_size: int = 128 * 1024 * 1024,  # 128 MB
        storage_options: dict[str, Any] | None = None,
        rt_setup: Any | None = None,
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

    def _interpolate_lut(
        self,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        wavelength: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Interpolate LUT to get RT parameters.

        Returns:
            Tuple of (path_reflectance, trans_down, trans_up, sph_albedo)
        """
        # Ensure LUT and coordinate cache are initialized.
        if not self._lut_coords:
            _ = self.lut
        wavelength_axis = np.asarray(
            self._lut_coords.get("wavelength", np.array([], dtype=np.float32)),
            dtype=np.float32,
        )
        if wavelength_axis.size == 0:
            raise ValueError("Coefficient LUT must define a non-empty wavelength coordinate")

        original_shape = sza.shape

        # Flatten inputs
        sza_flat = sza.ravel()
        vza_flat = vza.ravel()
        raa_flat = raa.ravel()
        aot_flat = aot.ravel()
        tcwv_flat = tcwv.ravel()

        # Find nearest wavelength in LUT
        wl_idx = np.argmin(np.abs(wavelength_axis - wavelength))
        wl_value = wavelength_axis[wl_idx]

        logger.debug(f"Using LUT wavelength {wl_value:.1f}nm for band at {wavelength:.1f}nm")

        # Select wavelength slice
        lut_wl = self.lut.sel(wavelength=wl_value, method="nearest")
        coords = self._build_point_coords(
            sza=sza_flat,
            vza=vza_flat,
            raa=raa_flat,
            aot=aot_flat,
            tcwv=tcwv_flat,
            tco3=tco3.ravel(),
            elevation=elevation.ravel(),
        )

        # Perform interpolation for each pixel
        # For efficiency, we batch pixels with similar coordinates
        if self.interpolation_method == "linear":
            path_ref, trans_down, trans_up, sph_alb = self._linear_interpolate(
                lut_wl,
                coords,
            )
        else:
            path_ref, trans_down, trans_up, sph_alb = self._nearest_interpolate(
                lut_wl,
                coords,
            )

        return (
            path_ref.reshape(original_shape),
            trans_down.reshape(original_shape),
            trans_up.reshape(original_shape),
            sph_alb.reshape(original_shape),
        )

    def _linear_interpolate(
        self,
        lut: xr.Dataset,
        coords: dict[str, xr.DataArray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Multi-dimensional linear interpolation."""
        path_ref = self._interpolate_variable_fast(lut["path_reflectance"], coords, "linear")
        trans_down = self._interpolate_variable_fast(lut["transmittance_down"], coords, "linear")
        trans_up = self._interpolate_variable_fast(lut["transmittance_up"], coords, "linear")
        sph_alb = self._interpolate_variable_fast(lut["spherical_albedo"], coords, "linear")

        return (path_ref, trans_down, trans_up, sph_alb)

    def _nearest_interpolate(
        self,
        lut: xr.Dataset,
        coords: dict[str, xr.DataArray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Nearest-neighbor interpolation."""
        path_ref = self._interpolate_variable_fast(lut["path_reflectance"], coords, "nearest")
        trans_down = self._interpolate_variable_fast(lut["transmittance_down"], coords, "nearest")
        trans_up = self._interpolate_variable_fast(lut["transmittance_up"], coords, "nearest")
        sph_alb = self._interpolate_variable_fast(lut["spherical_albedo"], coords, "nearest")

        return (path_ref, trans_down, trans_up, sph_alb)

    @staticmethod
    def _axis_bounds(axis: np.ndarray) -> tuple[float, float]:
        axis_values = np.asarray(axis, dtype=np.float32)
        return float(np.nanmin(axis_values)), float(np.nanmax(axis_values))

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
    def _point_indexer(values: np.ndarray) -> xr.DataArray:
        """Wrap flattened interpolation values in the LUT point dimension."""
        return xr.DataArray(np.asarray(values, dtype=np.float32), dims=["point"])

    @classmethod
    def _required_point_indexer(
        cls,
        values: np.ndarray,
        *,
        name: str,
    ) -> xr.DataArray:
        """Build a point indexer for a required LUT interpolation axis."""
        return cls._point_indexer(cls._require_finite_values(values, name=name))

    def _normalize_ozone_point_values(self, tco3: np.ndarray) -> np.ndarray:
        """Normalize ozone values to the LUT unit convention when needed.

        SIAC carries ozone in atm-cm. Some LUTs index ozone in DU, so values are
        upscaled by 1000 when the LUT axis is clearly in the DU regime.
        """
        ozone = np.asarray(tco3, dtype=np.float32)
        ozone_axis = np.asarray(self._lut_coords["ozone"], dtype=np.float32)
        # Atmospheric ozone often arrives in atm-cm (~0.3), while LUTs may use DU (~300).
        finite_ozone = ozone[np.isfinite(ozone)]
        if (
            ozone_axis.size
            and finite_ozone.size
            and np.nanmax(np.abs(ozone_axis)) > 20
            and np.nanmax(np.abs(finite_ozone)) < 10
        ):
            return cast("np.ndarray", ozone * 1000.0)
        return cast("np.ndarray", ozone)

    def _sanitize_point_indexers(self, coords: dict[str, xr.DataArray]) -> dict[str, xr.DataArray]:
        """Clamp available point indexers to the finite domain of each LUT axis."""
        sanitized = dict(coords)
        for name, coord in list(sanitized.items()):
            if name not in self._lut_coords:
                continue
            axis = np.asarray(self._lut_coords[name], dtype=np.float32)
            if axis.size == 0:
                continue
            sanitized[name] = self._point_indexer(self._sanitize_point_values(coord.values, axis))
        return sanitized

    def _build_point_coords(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> dict[str, xr.DataArray]:
        """Build point-indexer coordinates for interpolation.

        Expected units:

        - ``sza``, ``vza``, ``raa`` in degrees
        - ``aot`` unitless
        - ``tcwv`` in cm precipitable water
        - ``tco3`` in atm-cm; converted to DU automatically when the LUT ozone
          axis uses DU
        - ``elevation`` in km
        """
        coords: dict[str, xr.DataArray] = {
            "sza": self._required_point_indexer(sza, name="sza"),
            "vza": self._required_point_indexer(vza, name="vza"),
            "raa": self._required_point_indexer(raa, name="raa"),
            "aot": self._required_point_indexer(aot, name="aot"),
            "tcwv": self._required_point_indexer(tcwv, name="tcwv"),
        }

        if "ozone" in self._lut_coords:
            coords["ozone"] = self._point_indexer(self._normalize_ozone_point_values(tco3))

        if "altitude" in self._lut_coords:
            coords["altitude"] = self._point_indexer(elevation)

        return self._sanitize_point_indexers(coords)

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

        sza_mean, vza_mean, raa_mean, tco3_bounds, elevation_bounds = self._spectral_scene_summary(
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

        Angle means are snapped to the nearest LUT coordinate so that
        different spatial resolutions of the same scene (e.g. the
        atmospheric grid used by preload vs. the coarse grid used by M3)
        always produce the *same* key and share a single download.
        """
        sza_mean = self._finite_mean(sza, fallback=0.0)
        vza_mean = self._finite_mean(vza, fallback=0.0)
        raa_mean = self._finite_mean(raa, fallback=0.0)

        # Snap to nearest LUT grid point (mirrors _subset_spectral_lut_for_scene)
        coords = self._lut_coords
        if "sza" in coords and coords["sza"].size:
            sza_mean = float(coords["sza"][np.argmin(np.abs(coords["sza"] - sza_mean))])
        if "vza" in coords and coords["vza"].size:
            vza_mean = float(coords["vza"][np.argmin(np.abs(coords["vza"] - vza_mean))])
        if "raa" in coords and coords["raa"].size:
            raa_mean = float(coords["raa"][np.argmin(np.abs(coords["raa"] - raa_mean))])

        tco3_arr = np.asarray(tco3, dtype=np.float32)
        elevation_arr = np.asarray(elevation, dtype=np.float32)
        tco3_bounds = self._finite_range(tco3_arr, fallback=(0.0, 0.0))
        elevation_bounds = self._finite_range(elevation_arr, fallback=(0.0, 0.0))
        return (
            round(sza_mean, 3),
            round(vza_mean, 3),
            round(raa_mean, 3),
            round(tco3_bounds[0], 3),
            round(tco3_bounds[1], 3),
            round(elevation_bounds[0], 3),
            round(elevation_bounds[1], 3),
        )

    @classmethod
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

            with self._cache_lock:
                if self._spectral_scene_key != scene_key:
                    self._spectral_scene_key = scene_key
                    self._spectral_scene_subset = subset
                    self._spectral_band_grid_cache.clear()
                    self._scene_subset_logged = False
                if self._spectral_scene_subset is None:
                    self._spectral_scene_subset = subset
                return scene_key, self._spectral_scene_subset

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

        for dim, scene_values in {"ozone": tco3, "altitude": elevation}.items():
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
        """Build per-pixel LUT coordinates for unitless AOT and TCWV in cm."""
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

    @staticmethod
    def _weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
        """Weighted mean over wavelength with coordinate-aware integration."""
        return weighted_spectral_mean(data, weights)

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

    def _coefficient_lut_has_altitude_axis(self) -> bool:
        """Return True when coefficient LUT variables already model altitude explicitly."""
        return all(
            "altitude" in self.lut[name].dims for name in self._COEFFICIENT_VARS if name in self.lut
        )

    def _apply_elevation_correction(
        self,
        path_ref: np.ndarray,
        trans_down: np.ndarray,
        trans_up: np.ndarray,
        sph_alb: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Apply elevation correction to RT parameters.

        Uses Beer-Lambert law approximation for pressure-altitude relationship.

        Args:
            path_ref: Path reflectance
            trans_down: Downward transmittance
            trans_up: Upward transmittance
            sph_alb: Spherical albedo
            elevation: Surface elevation in km

        Returns:
            Corrected RT parameters
        """
        # Atmospheric Rayleigh scale height — see siac.constants.
        from siac.constants import ATMOSPHERIC_SCALE_HEIGHT_KM

        # Elevation correction factor
        correction = np.exp(-elevation / ATMOSPHERIC_SCALE_HEIGHT_KM)

        # Apply correction (reduce path effects at altitude)
        path_ref_corr = path_ref * correction
        sph_alb_corr = sph_alb * correction

        # Transmittance increases with altitude
        trans_correction = np.power(correction, 0.5)  # Approximate
        trans_down_corr = 1.0 - (1.0 - trans_down) * trans_correction
        trans_up_corr = 1.0 - (1.0 - trans_up) * trans_correction

        return path_ref_corr, trans_down_corr, trans_up_corr, sph_alb_corr

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

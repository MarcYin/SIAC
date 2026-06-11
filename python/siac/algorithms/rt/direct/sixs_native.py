"""Native 6SV2.1 Python-extension bridge and batched xarray runner."""

from __future__ import annotations

import ctypes
import hashlib
import importlib.machinery
import importlib.util
import logging
import os
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

from siac.algorithms.rt._run_cache import (
    load_cache_entry,
    resolve_run_cache_dir,
    store_cache_entry,
)
from siac.algorithms.rt.direct.sixs_build import ensure_native_sixs_module
from siac.geo.resample import resample_field_to_template, shares_template_grid
from siac.runtime.models import copy_spatial_metadata_like
from siac.sixs_outputs import (
    SIXS_BASE_OUTPUTS,
    SIXS_NATIVE_OUTPUT_NAMES,
    SIXS_OUTPUT_VARIABLE_CHOICES,
)

if TYPE_CHECKING:
    from datetime import datetime

    from siac.config.algorithms import RTSetupConfig, SixSAlgorithmConfig
    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime import AtmosphericState, GeometryAngles

logger = logging.getLogger(__name__)

_SIXS_GRID_START_UM = 0.25
_SIXS_GRID_STEP_UM = 0.0025
_SIXS_GRID_SIZE = 1501
_SIXS_GRID_UM = _SIXS_GRID_START_UM + _SIXS_GRID_STEP_UM * np.arange(
    _SIXS_GRID_SIZE, dtype=np.float64
)
_BASE_OUTPUTS = SIXS_BASE_OUTPUTS
_ALL_OUTPUTS = SIXS_OUTPUT_VARIABLE_CHOICES
_ATMOSPHERIC_PROFILE_CODES: dict[str, int] = {
    "no_gas": 0,
    "tropical": 1,
    "midlatitude_summer": 2,
    "midlatitude_winter": 3,
    "subarctic_summer": 4,
    "subarctic_winter": 5,
    "us_standard_62": 6,
    "user_profile": 7,
    "user_water_ozone": 8,
}
_AEROSOL_PROFILE_CODES: dict[str, int] = {
    "layered_profile": -1,
    "none": 0,
    "continental": 1,
    "maritime": 2,
    "urban": 3,
    "user_mixture": 4,
    "desert": 5,
    "biomass_burning": 6,
    "stratospheric": 7,
    "multimodal_log_normal": 8,
    "modified_gamma": 9,
    "junge_power_law": 10,
    "sun_photometer": 11,
    "user_model": 12,
}
_AEROSOL_LAYER_TYPE_CODES: dict[str, int] = {
    "none": 0,
    "continental": 1,
    "maritime": 2,
    "urban": 3,
    "desert": 5,
    "biomass_burning": 6,
    "stratospheric": 7,
}
_BUILTIN_REFLECTANCE_CODES: dict[str, int] = {
    "green_vegetation": 1,
    "clear_water": 2,
    "sand": 3,
    "lake_water": 4,
}
_BRDF_MODEL_CODES: dict[str, int] = {
    "user_defined": 0,
    "hapke": 1,
    "verstraete": 2,
    "roujean": 3,
    "walthall": 4,
    "minnaert": 5,
    "ocean": 6,
    "iaquinta_pinty": 7,
    "rahman": 8,
    "kuusk": 9,
    "modis": 10,
    "ross_li_maignan": 11,
}
_MODULE_CACHE: dict[tuple[Path, int, int], Any] = {}
_OPENMP_RUNTIME_PRELOADED = False


def _preload_openmp_runtime() -> None:
    """Load OpenMP symbols globally for native 6S extensions built on macOS."""
    global _OPENMP_RUNTIME_PRELOADED
    if _OPENMP_RUNTIME_PRELOADED:
        return

    candidates: list[Path] = []
    env_path = os.getenv("SIAC_SIXS_OPENMP_RUNTIME")
    if env_path:
        candidates.append(Path(env_path).expanduser())

    lib_dir = Path(sys.prefix) / "lib"
    for name in (
        "libgomp.dylib",
        "libgomp.1.dylib",
        "libgomp.so.1",
        "libgomp.so",
        "libomp.dylib",
        "libomp.so",
    ):
        candidates.append(lib_dir / name)

    for candidate in candidates:
        if not candidate.exists():
            continue
        try:
            ctypes.CDLL(os.fspath(candidate), mode=ctypes.RTLD_GLOBAL)
        except OSError as exc:
            logger.debug("Failed to preload OpenMP runtime %s: %s", candidate, exc)
            continue
        _OPENMP_RUNTIME_PRELOADED = True
        return


_CASE_ARRAY_NAMES: tuple[str, ...] = (
    "sza_deg",
    "saa_deg",
    "vza_deg",
    "vaa_deg",
    "aot550",
    "tcwv_cm",
    "tco3_atmcm",
    "elevation_km",
)


def _module_isolation_root() -> Path:
    root = Path.home().expanduser() / ".cache" / "siac" / "rt6s" / "module_copies"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _preferred_module_isolation_root(source: Path) -> Path:
    parent = source.expanduser().resolve().parent
    if os.access(parent, os.W_OK | os.X_OK):
        return parent
    return _module_isolation_root()


def _copy_isolated_module(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
        return
    except OSError:
        pass
    last_error: PermissionError | None = None
    for attempt in range(3):
        try:
            shutil.copyfile(source, destination)
            shutil.copystat(source, destination)
            return
        except PermissionError as exc:
            last_error = exc
            if attempt == 2:
                break
            time.sleep(0.05 * (attempt + 1))
    if last_error is not None:
        raise last_error


def _default_native_threads() -> int:
    cpu_total = os.cpu_count() or 1
    return max(1, cpu_total)


def _resample_geometry_to_template(
    geometry: GeometryAngles, template: xr.DataArray
) -> GeometryAngles:
    if all(
        shares_template_grid(field, template)
        for field in (geometry.sza, geometry.saa, geometry.vza, geometry.vaa)
    ):
        return geometry
    return type(geometry)(
        sza=resample_field_to_template(geometry.sza, template),
        saa=resample_field_to_template(geometry.saa, template),
        vza=resample_field_to_template(geometry.vza, template),
        vaa=resample_field_to_template(geometry.vaa, template),
    )


def _wrap_azimuth_degrees(angle: xr.DataArray, name: str) -> xr.DataArray:
    return xr.DataArray(
        np.mod(np.degrees(angle.values), 360.0),
        dims=angle.dims,
        coords=angle.coords,
        name=name,
    )


def _to_native_output_names(
    output_variables: tuple[str, ...] | list[str] | None,
) -> tuple[str, ...]:
    if output_variables is None:
        return _ALL_OUTPUTS
    selected = [name for name in output_variables if name not in {"xbp", "xcp"}]
    normalized: list[str] = ["xap", "xbp", "xcp"]
    for name in selected:
        if name not in normalized:
            normalized.append(name)
    return tuple(normalized)


def _build_spectral_response(band: SensorBand) -> tuple[np.ndarray, float, float]:
    if band.has_rsrf and band.rsrf_wavelengths_nm is not None and band.rsrf_response is not None:
        wavelengths_um = np.asarray(band.rsrf_wavelengths_nm, dtype=np.float64) / 1000.0
        response = np.asarray(band.rsrf_response, dtype=np.float64)
    else:
        half_width_um = max(float(band.bandwidth) / 2000.0, _SIXS_GRID_STEP_UM)
        center_um = float(band.center_wavelength) / 1000.0
        wavelengths_um = np.array(
            [center_um - half_width_um, center_um + half_width_um],
            dtype=np.float64,
        )
        response = np.array([1.0, 1.0], dtype=np.float64)

    order = np.argsort(wavelengths_um)
    wavelengths_um = wavelengths_um[order]
    response = np.clip(response[order], 0.0, None)

    native = np.interp(_SIXS_GRID_UM, wavelengths_um, response, left=0.0, right=0.0)
    positive = np.flatnonzero(native > 0.0)
    if positive.size == 0:
        center_index = int(
            np.clip(
                round((band.wavelength_um - _SIXS_GRID_START_UM) / _SIXS_GRID_STEP_UM),
                0,
                _SIXS_GRID_SIZE - 1,
            )
        )
        native[center_index] = 1.0
        positive = np.array([center_index], dtype=np.int64)

    area = float(np.trapezoid(native, _SIXS_GRID_UM))
    if area > 0.0:
        native /= area

    wlinf = float(_SIXS_GRID_UM[int(positive[0])])
    wlsup = float(_SIXS_GRID_UM[int(positive[-1])])
    return np.ascontiguousarray(native, dtype=np.float64), wlinf, wlsup


def _empty_radiosonde_profile() -> dict[str, np.ndarray]:
    return {
        "altitude_km": np.zeros(34, dtype=np.float64),
        "pressure_mb": np.zeros(34, dtype=np.float64),
        "temperature_k": np.zeros(34, dtype=np.float64),
        "water_g_m3": np.zeros(34, dtype=np.float64),
        "ozone_g_m3": np.zeros(34, dtype=np.float64),
    }


def _empty_aerosol_inputs() -> dict[str, np.ndarray | float | int]:
    return {
        "mixture": np.zeros(4, dtype=np.float64),
        "dist_rmin": 0.0,
        "dist_rmax": 0.0,
        "dist_component_count": 0,
        "dist_x1": np.zeros(4, dtype=np.float64),
        "dist_x2": np.zeros(4, dtype=np.float64),
        "dist_x3": np.zeros(4, dtype=np.float64),
        "dist_cij": np.zeros(4, dtype=np.float64),
        "dist_rn": np.zeros((20, 4), dtype=np.float64),
        "dist_ri": np.zeros((20, 4), dtype=np.float64),
        "sun_count": 0,
        "sun_radius": np.zeros(50, dtype=np.float64),
        "sun_dvlogr": np.zeros(50, dtype=np.float64),
        "layer_count": 0,
        "layer_height": np.zeros(50, dtype=np.float64),
        "layer_aot": np.zeros(50, dtype=np.float64),
        "layer_type": np.zeros(50, dtype=np.int32),
    }


def _empty_surface_inputs() -> dict[str, np.ndarray | float | int]:
    return {
        "inhomo": 0,
        "idirec": 0,
        "target_mode": 0,
        "target_constant": 0.0,
        "target_spectrum": np.zeros(_SIXS_GRID_SIZE, dtype=np.float64),
        "env_mode": 0,
        "env_constant": 0.0,
        "env_spectrum": np.zeros(_SIXS_GRID_SIZE, dtype=np.float64),
        "radius_km": 1.0,
        "brdf_model": 0,
        "brdf_params": np.zeros(12, dtype=np.float64),
        "brdf_options": np.zeros(5, dtype=np.int32),
        "brdf_struct": np.zeros(4, dtype=np.float64),
        "brdf_optics": np.zeros(3, dtype=np.float64),
        "brdf_table_solar": np.zeros((10, 13), dtype=np.float64),
        "brdf_table_view": np.zeros((10, 13), dtype=np.float64),
        "brdf_spherical_albedo": 0.0,
        "brdf_directional_reflectance": 0.0,
    }


def _auto_atmospheric_profile_from_latitude_and_month(latitude: float, month: int) -> str:
    rounded_lat = int(round(latitude, -1))
    rounded_lat = max(-80, min(80, rounded_lat))
    profiles = {
        1: {
            80: "subarctic_winter",
            70: "subarctic_winter",
            60: "midlatitude_winter",
            50: "midlatitude_winter",
            40: "subarctic_summer",
            30: "midlatitude_summer",
            20: "tropical",
            10: "tropical",
            0: "tropical",
            -10: "tropical",
            -20: "tropical",
            -30: "midlatitude_summer",
            -40: "subarctic_summer",
            -50: "subarctic_summer",
            -60: "midlatitude_winter",
            -70: "midlatitude_winter",
            -80: "midlatitude_winter",
        },
        5: {
            80: "subarctic_winter",
            70: "midlatitude_winter",
            60: "midlatitude_winter",
            50: "subarctic_summer",
            40: "subarctic_summer",
            30: "midlatitude_summer",
            20: "tropical",
            10: "tropical",
            0: "tropical",
            -10: "tropical",
            -20: "tropical",
            -30: "midlatitude_summer",
            -40: "subarctic_summer",
            -50: "subarctic_summer",
            -60: "midlatitude_winter",
            -70: "midlatitude_winter",
            -80: "midlatitude_winter",
        },
        7: {
            80: "midlatitude_winter",
            70: "midlatitude_winter",
            60: "subarctic_summer",
            50: "subarctic_summer",
            40: "midlatitude_summer",
            30: "tropical",
            20: "tropical",
            10: "tropical",
            0: "tropical",
            -10: "tropical",
            -20: "midlatitude_summer",
            -30: "midlatitude_summer",
            -40: "subarctic_summer",
            -50: "midlatitude_winter",
            -60: "midlatitude_winter",
            -70: "midlatitude_winter",
            -80: "subarctic_winter",
        },
        9: {
            80: "midlatitude_winter",
            70: "midlatitude_winter",
            60: "subarctic_summer",
            50: "subarctic_summer",
            40: "midlatitude_summer",
            30: "tropical",
            20: "tropical",
            10: "tropical",
            0: "tropical",
            -10: "tropical",
            -20: "midlatitude_summer",
            -30: "midlatitude_summer",
            -40: "subarctic_summer",
            -50: "midlatitude_winter",
            -60: "midlatitude_winter",
            -70: "midlatitude_winter",
            -80: "midlatitude_winter",
        },
        11: {
            80: "subarctic_winter",
            70: "subarctic_winter",
            60: "midlatitude_winter",
            50: "subarctic_summer",
            40: "subarctic_summer",
            30: "midlatitude_summer",
            20: "tropical",
            10: "tropical",
            0: "tropical",
            -10: "tropical",
            -20: "tropical",
            -30: "midlatitude_summer",
            -40: "subarctic_summer",
            -50: "subarctic_summer",
            -60: "midlatitude_winter",
            -70: "midlatitude_winter",
            -80: "midlatitude_winter",
        },
    }
    if month in {1, 2, 3, 4}:
        lookup = profiles[1]
    elif month in {5, 6}:
        lookup = profiles[5]
    elif month in {7, 8}:
        lookup = profiles[7]
    elif month in {9, 10}:
        lookup = profiles[9]
    else:
        lookup = profiles[11]
    return lookup[rounded_lat]


def _resolve_atmospheric_inputs(
    rt_setup: RTSetupConfig, *, month: int
) -> tuple[int, dict[str, np.ndarray]]:
    atmosphere = rt_setup.atmosphere
    if atmosphere is None or atmosphere.profile is None:
        raise ValueError(
            "Native 6S requires rt_setup.atmosphere.profile to be resolved before execution."
        )
    profile = atmosphere.profile
    radiosonde = _empty_radiosonde_profile()
    if profile == "auto_latitude_date":
        if atmosphere.profile_latitude is None:
            raise ValueError(
                "rt.setup.atmosphere.profile_latitude is required when profile='auto_latitude_date'."
            )
        profile = _auto_atmospheric_profile_from_latitude_and_month(
            float(atmosphere.profile_latitude),
            int(month),
        )
    if profile == "user_profile":
        if atmosphere.radiosonde_profile is None:
            raise ValueError(
                "rt.setup.atmosphere.radiosonde_profile is required for profile='user_profile'."
            )
        radiosonde = {
            "altitude_km": np.asarray(atmosphere.radiosonde_profile.altitude_km, dtype=np.float64),
            "pressure_mb": np.asarray(atmosphere.radiosonde_profile.pressure_mb, dtype=np.float64),
            "temperature_k": np.asarray(
                atmosphere.radiosonde_profile.temperature_k, dtype=np.float64
            ),
            "water_g_m3": np.asarray(atmosphere.radiosonde_profile.water_g_m3, dtype=np.float64),
            "ozone_g_m3": np.asarray(atmosphere.radiosonde_profile.ozone_g_m3, dtype=np.float64),
        }
    return _ATMOSPHERIC_PROFILE_CODES[profile], radiosonde


def _resolve_atmospheric_mode(rt_setup: RTSetupConfig, *, month: int) -> int:
    mode, _ = _resolve_atmospheric_inputs(rt_setup, month=month)
    return mode


def _resolve_atmospheric_columns_mode(rt_setup: RTSetupConfig) -> int:
    atmosphere = rt_setup.atmosphere
    if atmosphere is None or atmosphere.profile is None:
        raise ValueError("Native 6S requires rt_setup.atmosphere to be resolved before execution.")
    if atmosphere.columns_mode != "input_columns":
        return 0
    if atmosphere.profile == "no_gas":
        return 0
    return 1


def _resolve_aerosol_mode(rt_setup: RTSetupConfig) -> int:
    aerosol = rt_setup.aerosol
    if aerosol is None or aerosol.profile is None:
        raise ValueError(
            "Native 6S requires rt_setup.aerosol.profile to be resolved before execution."
        )
    return _AEROSOL_PROFILE_CODES[aerosol.profile]


def _resolve_aerosol_inputs(
    rt_setup: RTSetupConfig,
) -> tuple[int, dict[str, np.ndarray | float | int]]:
    aerosol = rt_setup.aerosol
    if aerosol is None or aerosol.profile is None:
        raise ValueError("Native 6S requires rt_setup.aerosol to be resolved before execution.")
    mode = _resolve_aerosol_mode(rt_setup)
    payload = _empty_aerosol_inputs()
    if aerosol.profile == "user_mixture":
        payload["mixture"] = np.asarray(aerosol.mixture or (0.0, 0.0, 0.0, 0.0), dtype=np.float64)
    elif aerosol.profile in {"multimodal_log_normal", "modified_gamma", "junge_power_law"}:
        if aerosol.distribution is None:
            raise ValueError("rt.setup.aerosol.distribution is required for aerosol distributions.")
        payload["dist_rmin"] = float(aerosol.distribution.rmin)
        payload["dist_rmax"] = float(aerosol.distribution.rmax)
        payload["dist_component_count"] = len(aerosol.distribution.components)
        x1 = np.zeros(4, dtype=np.float64)
        x2 = np.zeros(4, dtype=np.float64)
        x3 = np.zeros(4, dtype=np.float64)
        cij = np.zeros(4, dtype=np.float64)
        rn = np.zeros((20, 4), dtype=np.float64)
        ri = np.zeros((20, 4), dtype=np.float64)
        for idx, component in enumerate(aerosol.distribution.components):
            x1[idx] = float(component.rmean)
            x2[idx] = float(component.sigma)
            rn[:, idx] = np.asarray(component.refr_real, dtype=np.float64)
            ri[:, idx] = np.asarray(component.refr_imag, dtype=np.float64)
        if aerosol.profile == "multimodal_log_normal":
            for idx, component in enumerate(aerosol.distribution.components):
                cij[idx] = float(component.percentage_density)
        elif aerosol.profile == "modified_gamma":
            x3[0] = float(aerosol.distribution.components[0].percentage_density)
            cij[0] = 1.0
        else:
            cij[0] = 1.0
        payload["dist_x1"] = x1
        payload["dist_x2"] = x2
        payload["dist_x3"] = x3
        payload["dist_cij"] = cij
        payload["dist_rn"] = rn
        payload["dist_ri"] = ri
    elif aerosol.profile == "sun_photometer":
        if aerosol.sun_photometer_aerosol is None:
            raise ValueError(
                "rt.setup.aerosol.sun_photometer_aerosol is required for sun-photometer aerosols."
            )
        count = len(aerosol.sun_photometer_aerosol.radii_um)
        radius = np.zeros(50, dtype=np.float64)
        dvlogr = np.zeros(50, dtype=np.float64)
        radius[:count] = np.asarray(aerosol.sun_photometer_aerosol.radii_um, dtype=np.float64)
        dvlogr[:count] = np.asarray(aerosol.sun_photometer_aerosol.dv_dlogr, dtype=np.float64)
        rn = np.zeros((20, 4), dtype=np.float64)
        ri = np.zeros((20, 4), dtype=np.float64)
        rn[:, 0] = np.asarray(aerosol.sun_photometer_aerosol.refr_real, dtype=np.float64)
        ri[:, 0] = np.asarray(aerosol.sun_photometer_aerosol.refr_imag, dtype=np.float64)
        payload["sun_count"] = count
        payload["sun_radius"] = radius
        payload["sun_dvlogr"] = dvlogr
        payload["dist_rn"] = rn
        payload["dist_ri"] = ri
    elif aerosol.profile == "layered_profile":
        if not aerosol.layer_profile:
            raise ValueError(
                "rt.setup.aerosol.layer_profile is required for layered aerosol profiles."
            )
        count = len(aerosol.layer_profile)
        height = np.zeros(50, dtype=np.float64)
        aot = np.zeros(50, dtype=np.float64)
        layer_type = np.zeros(50, dtype=np.int32)
        for idx, layer in enumerate(aerosol.layer_profile):
            height[idx] = float(layer.height_km)
            aot[idx] = float(layer.optical_thickness)
            layer_type[idx] = _AEROSOL_LAYER_TYPE_CODES[layer.aerosol_type]
        payload["layer_count"] = count
        payload["layer_height"] = height
        payload["layer_aot"] = aot
        payload["layer_type"] = layer_type
    return mode, payload


def _build_surface_spectrum(
    values: tuple[float, ...] | None, wavelengths_um: tuple[float, ...] | None
) -> np.ndarray:
    if values is None:
        return np.zeros(_SIXS_GRID_SIZE, dtype=np.float64)
    data = np.asarray(values, dtype=np.float64)
    if wavelengths_um is None:
        if data.size != _SIXS_GRID_SIZE:
            raise ValueError(
                "Surface spectra without explicit wavelengths must contain 1501 samples at 2.5 nm spacing."
            )
        return np.ascontiguousarray(np.clip(data, 0.0, None), dtype=np.float64)
    wavelength_values = np.asarray(wavelengths_um, dtype=np.float64)
    order = np.argsort(wavelength_values)
    interpolated = np.interp(
        _SIXS_GRID_UM,
        wavelength_values[order],
        data[order],
        left=0.0,
        right=0.0,
    )
    return np.ascontiguousarray(np.clip(interpolated, 0.0, None), dtype=np.float64)


def _resolve_surface_reflectance(
    reflectance: Any,
) -> tuple[int, float, np.ndarray]:
    if reflectance.kind == "constant":
        return 0, float(reflectance.constant or 0.0), np.zeros(_SIXS_GRID_SIZE, dtype=np.float64)
    if reflectance.kind == "built_in":
        return (
            _BUILTIN_REFLECTANCE_CODES[reflectance.built_in],
            0.0,
            np.zeros(_SIXS_GRID_SIZE, dtype=np.float64),
        )
    return -1, 0.0, _build_surface_spectrum(reflectance.values, reflectance.wavelengths_um)


def _prepare_brdf_table(table: tuple[tuple[float, ...], ...] | None) -> np.ndarray:
    if table is None:
        return np.zeros((10, 13), dtype=np.float64)
    values = np.asarray(table, dtype=np.float64)
    if values.shape != (13, 10):
        raise ValueError("BRDF tables must be provided as 13 azimuth rows by 10 zenith columns.")
    prepared = np.zeros((10, 13), dtype=np.float64)
    for k in range(13):
        prepared[:, k] = values[k, ::-1]
    return np.ascontiguousarray(prepared, dtype=np.float64)


def _resolve_surface_inputs(rt_setup: RTSetupConfig) -> dict[str, np.ndarray | float | int]:
    payload = _empty_surface_inputs()
    surface = rt_setup.surface
    if surface is None:
        raise ValueError("Native 6S requires rt_setup.surface to be resolved before execution.")
    target_mode, target_constant, target_spectrum = _resolve_surface_reflectance(surface.target)
    payload["target_mode"] = target_mode
    payload["target_constant"] = target_constant
    payload["target_spectrum"] = target_spectrum

    if surface.mode == "homogeneous_lambertian":
        payload["inhomo"] = 0
        payload["idirec"] = 0
        return payload

    if surface.mode == "heterogeneous_lambertian":
        if surface.environment is None:
            raise ValueError(
                "Heterogeneous Lambertian surfaces require an environment reflectance."
            )
        env_mode, env_constant, env_spectrum = _resolve_surface_reflectance(surface.environment)
        payload["inhomo"] = 1
        payload["idirec"] = 0
        payload["env_mode"] = env_mode
        payload["env_constant"] = env_constant
        payload["env_spectrum"] = env_spectrum
        payload["radius_km"] = float(surface.radius_km)
        return payload

    if surface.brdf is None:
        raise ValueError("BRDF surfaces require a brdf configuration.")
    brdf = surface.brdf
    payload["inhomo"] = 0
    payload["idirec"] = 1
    payload["brdf_model"] = _BRDF_MODEL_CODES[brdf.model]
    params = np.zeros(12, dtype=np.float64)
    options = np.zeros(5, dtype=np.int32)
    struct = np.zeros(4, dtype=np.float64)
    optics = np.zeros(3, dtype=np.float64)

    if brdf.model == "hapke":
        params[:4] = [
            brdf.parameters["albedo"],
            brdf.parameters["asymmetry"],
            brdf.parameters["amplitude"],
            brdf.parameters["width"],
        ]
    elif brdf.model == "verstraete":
        options[2] = int(brdf.options["kappa_parameterisation"])
        options[3] = int(brdf.options["phase_function"])
        options[4] = int(brdf.options["multiple_scattering"])
        struct[:] = [
            brdf.parameters["leaf_area_density"],
            brdf.parameters["sun_flecks_radius"],
            brdf.parameters.get("kappa1", brdf.parameters.get("chil", 0.0)),
            brdf.parameters.get("kappa2", 0.0),
        ]
        optics[:] = [
            brdf.parameters["single_scattering_albedo"],
            brdf.parameters.get("phase_parameter_1", 0.0),
            brdf.parameters.get("phase_parameter_2", 0.0),
        ]
    elif brdf.model == "roujean":
        params[:3] = [brdf.parameters["albedo"], brdf.parameters["k1"], brdf.parameters["k2"]]
    elif brdf.model == "walthall":
        params[:4] = [
            brdf.parameters["param1"],
            brdf.parameters["param2"],
            brdf.parameters["param3"],
            brdf.parameters["albedo"],
        ]
    elif brdf.model == "minnaert":
        params[:2] = [brdf.parameters["k"], brdf.parameters["albedo"]]
    elif brdf.model == "ocean":
        params[:4] = [
            brdf.parameters["wind_speed"],
            brdf.parameters["wind_azimuth"],
            brdf.parameters["salinity"],
            brdf.parameters["pigment_concentration"],
        ]
    elif brdf.model == "iaquinta_pinty":
        options[0] = int(brdf.options["leaf_distribution"])
        options[1] = int(brdf.options["hot_spot"])
        params[:5] = [
            brdf.parameters["lai"],
            brdf.parameters["hot_spot_parameter"],
            brdf.parameters["leaf_reflectance"],
            brdf.parameters["leaf_transmittance"],
            brdf.parameters["soil_albedo"],
        ]
    elif brdf.model == "rahman":
        params[:3] = [
            brdf.parameters["intensity"],
            brdf.parameters["asymmetry_factor"],
            brdf.parameters["structural_parameter"],
        ]
    elif brdf.model == "kuusk":
        params[:9] = [
            brdf.parameters["lai"],
            brdf.parameters["lad_eps"],
            brdf.parameters["lad_thm"],
            brdf.parameters["relative_leaf_size"],
            brdf.parameters["chlorophyll_content"],
            brdf.parameters["leaf_water_equiv_thickness"],
            brdf.parameters["effective_num_layers"],
            brdf.parameters["ratio_refractive_indices"],
            brdf.parameters["weight_first_price_function"],
        ]
    elif brdf.model in {"modis", "ross_li_maignan"}:
        params[:3] = [brdf.parameters["k_iso"], brdf.parameters["k_vol"], brdf.parameters["k_geo"]]
    elif brdf.model == "user_defined":
        payload["brdf_table_solar"] = _prepare_brdf_table(brdf.table_solar_zenith)
        payload["brdf_table_view"] = _prepare_brdf_table(brdf.table_view_zenith)
        payload["brdf_spherical_albedo"] = float(brdf.spherical_albedo or 0.0)
        payload["brdf_directional_reflectance"] = float(brdf.directional_reflectance or 0.0)

    payload["brdf_params"] = params
    payload["brdf_options"] = options
    payload["brdf_struct"] = struct
    payload["brdf_optics"] = optics
    return payload


def _resolve_atmospheric_correction_inputs(rt_setup: RTSetupConfig) -> tuple[int, float]:
    correction = rt_setup.atmospheric_correction
    correction_mode = (
        correction.mode if correction is not None and correction.mode is not None else "none"
    )
    reference_reflectance = (
        rt_setup.reference_reflectance if rt_setup.reference_reflectance is not None else 0.1
    )
    value = float(
        correction.value
        if correction is not None and correction.value is not None
        else reference_reflectance
    )
    mapping = {
        "none": (-1, 0.0),
        "lambertian_reflectance": (0, -value),
        "lambertian_radiance": (0, value),
        "brdf_reflectance": (1, -value),
        "brdf_radiance": (1, value),
    }
    return mapping[correction_mode]


def _encode_aerosol_model_path(path: Path | None) -> str:
    if path is None:
        return ""
    return os.fspath(path)


def _empty_output_bundle(length: int) -> dict[str, np.ndarray]:
    return {name: np.empty(length, dtype=np.float64) for name in _ALL_OUTPUTS}


def _as_float64(array: np.ndarray | xr.DataArray) -> np.ndarray:
    values = np.asarray(array, dtype=np.float64)
    return np.ascontiguousarray(values.reshape(-1), dtype=np.float64)


def _as_output_array(values: np.ndarray, template: xr.DataArray, name: str) -> xr.DataArray:
    data = xr.DataArray(
        values.reshape(template.shape),
        dims=template.dims,
        coords=template.coords,
        name=name,
    )
    return copy_spatial_metadata_like(data, template)


def _nan_output_arrays(
    template: xr.DataArray, selected_names: tuple[str, ...]
) -> dict[str, xr.DataArray]:
    return {
        name: _as_output_array(np.full(template.size, np.nan, dtype=np.float64), template, name)
        for name in selected_names
    }


def _slice_case_kwargs(kwargs: dict[str, Any], start: int, stop: int) -> dict[str, Any]:
    sliced = dict(kwargs)
    for name in _CASE_ARRAY_NAMES:
        sliced[name] = kwargs[name][start:stop]
    return sliced


def _build_scene_lut_axis(values: np.ndarray, target_size: int) -> np.ndarray:
    finite = np.asarray(values[np.isfinite(values)], dtype=np.float64)
    if finite.size == 0:
        return np.zeros(1, dtype=np.float64)
    unique = np.unique(finite)
    if unique.size <= max(1, target_size):
        return np.ascontiguousarray(unique, dtype=np.float64)
    if target_size <= 1:
        return np.ascontiguousarray(unique[:1], dtype=np.float64)
    quantiles = np.linspace(0.0, 1.0, target_size, dtype=np.float64)
    axis = np.quantile(finite, quantiles, method="linear")
    axis[0] = float(finite.min())
    axis[-1] = float(finite.max())
    axis = np.unique(np.asarray(axis, dtype=np.float64))
    if axis.size == 1 and unique.size > 1:
        axis = np.array([unique[0], unique[-1]], dtype=np.float64)
    return np.ascontiguousarray(axis, dtype=np.float64)


def _scene_lut_case_count(axes: dict[str, np.ndarray]) -> int:
    count = 1
    for axis in axes.values():
        count *= max(1, int(axis.size))
    return count


def _build_scene_lut_plan(
    case_arrays: dict[str, np.ndarray],
    *,
    max_nodes_per_axis: int,
    max_cases: int,
) -> _SceneLUTPlan:
    axes = {
        name: _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64), max_nodes_per_axis
        )
        for name in _CASE_ARRAY_NAMES
    }
    while _scene_lut_case_count(axes) > max_cases:
        reducible = [name for name, axis in axes.items() if axis.size > 1]
        if not reducible:
            break
        name = max(reducible, key=lambda item: axes[item].size)
        axes[name] = _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64), axes[name].size - 1
        )

    mesh = np.meshgrid(*(axes[name] for name in _CASE_ARRAY_NAMES), indexing="ij")
    grid_case_arrays = {
        name: np.ascontiguousarray(mesh[idx].reshape(-1), dtype=np.float64)
        for idx, name in enumerate(_CASE_ARRAY_NAMES)
    }
    return _SceneLUTPlan(
        axes=axes,
        grid_case_arrays=grid_case_arrays,
        direct_case_count=int(np.asarray(case_arrays[_CASE_ARRAY_NAMES[0]]).size),
        lut_case_count=int(grid_case_arrays[_CASE_ARRAY_NAMES[0]].size),
    )


def _interpolate_scene_lut_outputs(
    plan: _SceneLUTPlan,
    native_outputs: _NativeBatchResult,
    case_arrays: dict[str, np.ndarray],
    selected_names: tuple[str, ...],
) -> dict[str, np.ndarray]:
    from scipy.interpolate import RegularGridInterpolator

    axes_order = _CASE_ARRAY_NAMES
    varying = [name for name in axes_order if plan.axes[name].size > 1]
    n_cases = int(np.asarray(case_arrays[axes_order[0]]).size)
    result: dict[str, np.ndarray] = {}
    full_shape = tuple(int(plan.axes[name].size) for name in axes_order)
    if not varying:
        for name in selected_names:
            value = float(np.asarray(native_outputs.outputs[name], dtype=np.float64).reshape(-1)[0])
            result[name] = np.full(n_cases, value, dtype=np.float64)
        return result

    sample_points = np.column_stack(
        [np.asarray(case_arrays[name], dtype=np.float64) for name in varying]
    )
    for name in selected_names:
        values = np.asarray(native_outputs.outputs[name], dtype=np.float64).reshape(full_shape)
        reduced = values
        for axis_index in reversed(range(len(axes_order))):
            if plan.axes[axes_order[axis_index]].size == 1:
                reduced = np.take(reduced, 0, axis=axis_index)
        interpolator = RegularGridInterpolator(
            tuple(plan.axes[axis_name] for axis_name in varying),
            reduced,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        result[name] = np.ascontiguousarray(interpolator(sample_points), dtype=np.float64)
    return result


def _should_use_scene_lut(
    *,
    mode: str,
    direct_case_count: int,
    lut_case_count: int,
    min_pixels: int,
    required_speedup: float,
) -> bool:
    if mode == "direct":
        return False
    if mode == "scene_lut":
        return True
    if direct_case_count < min_pixels or lut_case_count <= 0 or lut_case_count >= direct_case_count:
        return False
    return (float(direct_case_count) / float(lut_case_count)) >= required_speedup


def _build_joint_grid_search_lut_plan(
    case_arrays: dict[str, np.ndarray],
    *,
    aot_axis: np.ndarray,
    tcwv_axis: np.ndarray,
    max_nodes_per_axis: int,
    max_cases: int,
) -> _SceneLUTPlan:
    """Build a scene-LUT plan with explicit aot/tcwv axes for joint grid-search reuse.

    Unlike :func:`_build_scene_lut_plan`, this builder takes the aot550 and
    tcwv_cm axes as inputs (rather than deriving them from per-pixel
    candidate values). The remaining six geometric/atmospheric axes are
    derived from the per-pixel ``case_arrays`` as usual. The trimming step
    that reduces total case count to fit ``max_cases`` only shrinks the
    geometric axes — the explicit aot/tcwv axes are preserved because their
    nodes must coincide with the grid-search candidate values for the
    block-grid-search reuse to be numerically exact at the grid points.
    """
    aot_axis_arr = np.ascontiguousarray(np.unique(np.asarray(aot_axis, dtype=np.float64)))
    tcwv_axis_arr = np.ascontiguousarray(np.unique(np.asarray(tcwv_axis, dtype=np.float64)))
    if aot_axis_arr.size == 0:
        aot_axis_arr = np.zeros(1, dtype=np.float64)
    if tcwv_axis_arr.size == 0:
        tcwv_axis_arr = np.zeros(1, dtype=np.float64)
    fixed_axes = {"aot550", "tcwv_cm"}
    axes: dict[str, np.ndarray] = {}
    for name in _CASE_ARRAY_NAMES:
        if name == "aot550":
            axes[name] = aot_axis_arr
        elif name == "tcwv_cm":
            axes[name] = tcwv_axis_arr
        else:
            axes[name] = _build_scene_lut_axis(
                np.asarray(case_arrays[name], dtype=np.float64), max_nodes_per_axis
            )

    # Shrink only the geometric axes to fit the case budget.
    while _scene_lut_case_count(axes) > max_cases:
        reducible = [
            name for name in _CASE_ARRAY_NAMES if name not in fixed_axes and axes[name].size > 1
        ]
        if not reducible:
            break
        name = max(reducible, key=lambda item: axes[item].size)
        axes[name] = _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64), axes[name].size - 1
        )

    mesh = np.meshgrid(*(axes[name] for name in _CASE_ARRAY_NAMES), indexing="ij")
    grid_case_arrays = {
        name: np.ascontiguousarray(mesh[idx].reshape(-1), dtype=np.float64)
        for idx, name in enumerate(_CASE_ARRAY_NAMES)
    }
    return _SceneLUTPlan(
        axes=axes,
        grid_case_arrays=grid_case_arrays,
        direct_case_count=int(np.asarray(case_arrays[_CASE_ARRAY_NAMES[0]]).size),
        lut_case_count=int(grid_case_arrays[_CASE_ARRAY_NAMES[0]].size),
    )


@dataclass(frozen=True)
class _NativeBatchResult:
    outputs: dict[str, np.ndarray]
    status: np.ndarray


@dataclass(frozen=True)
class _PreparedSceneInputs:
    template: xr.DataArray
    valid_mask: np.ndarray
    flat_valid_mask: np.ndarray
    common_kwargs: dict[str, Any]
    case_arrays: dict[str, np.ndarray]


@dataclass(frozen=True)
class _SceneLUTPlan:
    axes: dict[str, np.ndarray]
    grid_case_arrays: dict[str, np.ndarray]
    direct_case_count: int
    lut_case_count: int


@dataclass(frozen=True)
class PreparedCorrectionScene:
    """Scene-level prep + LUT plan shared by all bands in M6 correction.

    Wave 18e: in the previous design every per-band ``compute_coefficients``
    call inside :class:`siac.algorithms.correction.AtmosphericCorrector`
    re-ran ``_prepare_scene_inputs`` and ``_build_scene_lut_plan`` from
    scratch — work that doesn't depend on the spectral band and is
    redundant across 13 bands. This object captures that prep once so the
    per-band executor can dispatch only the band-specific 6S work.

    The object is intentionally **read-only** and safe to pass to multiple
    threads concurrently — the 6S call reads ``prepared.case_arrays`` and
    ``plan.grid_case_arrays`` but never mutates them.
    """

    prepared: _PreparedSceneInputs
    plan: _SceneLUTPlan
    use_scene_lut: bool


class JointGridSearchLUT:
    """Precomputed RT-coefficient LUT spanning a 2-D (aot, tcwv) grid.

    Built once via :meth:`SixSNativeRunner.build_joint_grid_search_lut`, then
    queried via :meth:`evaluate` to retrieve per-band ``xap``/``xbp``/``xcp``
    DataArrays matching the scene template — without invoking the native 6S
    runner again.

    The LUT spans the cross-product of the caller-supplied ``aot_axis`` and
    ``tcwv_axis`` with the (coarsened) per-pixel geometric/atmospheric axes
    (sza, saa, vza, vaa, tco3, elevation). At evaluate time the per-pixel
    geometric/atmospheric coordinates are looked up directly while the
    aot/tcwv coordinates are broadcast scalars — typically equal to one of
    the LUT node values, in which case the lookup is exact (no aot/tcwv
    interpolation error). This is the key invariant the block-grid-search
    relies on for numerical equivalence with the per-candidate scene-LUT
    path.
    """

    def __init__(
        self,
        *,
        prepared: _PreparedSceneInputs,
        plan: _SceneLUTPlan,
        band_native_outputs: list[_NativeBatchResult],
        selected_names: tuple[str, ...],
    ) -> None:
        from scipy.interpolate import RegularGridInterpolator

        self._prepared = prepared
        self._plan = plan
        self._selected_names = selected_names
        self._n_pixels = int(np.count_nonzero(prepared.valid_mask))
        # Names of axes that genuinely vary in the LUT (some may collapse to
        # a single node, e.g. when the scene is uniform in that dimension).
        self._varying_axes: tuple[str, ...] = tuple(
            name for name in _CASE_ARRAY_NAMES if plan.axes[name].size > 1
        )
        self._varying_is_atmo: tuple[bool, ...] = tuple(
            name in {"aot550", "tcwv_cm"} for name in self._varying_axes
        )

        # Per-pixel sample points for the varying geometric axes are constant
        # across all evaluate() calls — precompute the columns once.
        self._geom_columns: list[np.ndarray] = []
        for name in self._varying_axes:
            if name in {"aot550", "tcwv_cm"}:
                self._geom_columns.append(np.empty(0, dtype=np.float64))  # placeholder
            else:
                self._geom_columns.append(
                    np.ascontiguousarray(np.asarray(prepared.case_arrays[name], dtype=np.float64))
                )

        # Build a RegularGridInterpolator per (band, output_name). These reuse
        # the underlying LUT values; only the per-pixel sample points change
        # between evaluate() calls.
        full_shape = tuple(int(plan.axes[name].size) for name in _CASE_ARRAY_NAMES)
        varying_axes_values = tuple(plan.axes[name] for name in self._varying_axes)
        self._band_interpolators: list[dict[str, Any]] = []
        for native_outputs in band_native_outputs:
            band_interp: dict[str, Any] = {}
            for name in selected_names:
                values = np.asarray(native_outputs.outputs[name], dtype=np.float64).reshape(
                    full_shape
                )
                reduced = values
                for axis_index in reversed(range(len(_CASE_ARRAY_NAMES))):
                    if plan.axes[_CASE_ARRAY_NAMES[axis_index]].size == 1:
                        reduced = np.take(reduced, 0, axis=axis_index)
                if not self._varying_axes:
                    band_interp[name] = _ScalarInterpolator(float(reduced))
                else:
                    band_interp[name] = RegularGridInterpolator(
                        varying_axes_values,
                        np.ascontiguousarray(reduced, dtype=np.float64),
                        method="linear",
                        bounds_error=False,
                        fill_value=np.nan,
                    )
            self._band_interpolators.append(band_interp)

        # If the scene has zero valid pixels we can short-circuit evaluate().
        self._empty = self._n_pixels == 0

    @property
    def band_count(self) -> int:
        return len(self._band_interpolators)

    @property
    def template(self) -> xr.DataArray:
        return self._prepared.template

    @property
    def selected_names(self) -> tuple[str, ...]:
        return self._selected_names

    @property
    def plan(self) -> _SceneLUTPlan:
        return self._plan

    def evaluate(self, aot_val: float, tcwv_val: float) -> list[dict[str, xr.DataArray]]:
        """Return per-band xap/xbp/xcp DataArrays for the given (aot, tcwv)."""
        prepared = self._prepared
        if self._empty:
            return [
                _nan_output_arrays(prepared.template, self._selected_names)
                for _ in self._band_interpolators
            ]

        if self._varying_axes:
            n_pixels = self._n_pixels
            sample_columns: list[np.ndarray] = []
            for idx, name in enumerate(self._varying_axes):
                if name == "aot550":
                    sample_columns.append(np.full(n_pixels, float(aot_val), dtype=np.float64))
                elif name == "tcwv_cm":
                    sample_columns.append(np.full(n_pixels, float(tcwv_val), dtype=np.float64))
                else:
                    sample_columns.append(self._geom_columns[idx])
            sample_points = np.column_stack(sample_columns)
        else:
            sample_points = np.empty((0, 0), dtype=np.float64)

        results: list[dict[str, xr.DataArray]] = []
        flat_valid_mask = prepared.flat_valid_mask
        template = prepared.template
        for band_interp in self._band_interpolators:
            outputs_by_name: dict[str, xr.DataArray] = {}
            for name in self._selected_names:
                interpolator = band_interp[name]
                if self._varying_axes:
                    interpolated = np.ascontiguousarray(
                        interpolator(sample_points), dtype=np.float64
                    )
                else:
                    interpolated = np.full(self._n_pixels, float(interpolator()), dtype=np.float64)
                full = np.full(template.size, np.nan, dtype=np.float64)
                full[flat_valid_mask] = interpolated
                outputs_by_name[name] = _as_output_array(full, template, name)
            results.append(outputs_by_name)
        return results


class _ScalarInterpolator:
    """Stand-in for RegularGridInterpolator when the LUT is fully degenerate."""

    __slots__ = ("_value",)

    def __init__(self, value: float) -> None:
        self._value = float(value)

    def __call__(self, points: np.ndarray | None = None) -> float:  # noqa: ARG002
        # ``points`` is accepted but unused — kept for signature parity
        # with ``scipy.interpolate.RegularGridInterpolator.__call__``,
        # which the JointGridSearchLUT calls polymorphically.
        return self._value


def _load_extension_module(module_path: Path) -> Any:
    resolved_path = Path(module_path).resolve()
    stat = resolved_path.stat()
    cache_key = (resolved_path, stat.st_mtime_ns, stat.st_size)
    cached = _MODULE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    for stale_key in [key for key in _MODULE_CACHE if key[0] == resolved_path and key != cache_key]:
        _MODULE_CACHE.pop(stale_key, None)

    module_name = resolved_path.name.split(".", 1)[0]
    _preload_openmp_runtime()

    loader = importlib.machinery.ExtensionFileLoader(
        module_name,
        os.fspath(resolved_path),
    )
    spec = importlib.util.spec_from_file_location(
        module_name,
        os.fspath(resolved_path),
        loader=loader,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not create an import spec for native 6S module: {resolved_path}")

    sys.modules.pop(module_name, None)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    _MODULE_CACHE[cache_key] = module
    return module


class _SixSExtensionModule:
    """Thin loader around the compiled F2PY extension."""

    def __init__(self, module_path: Path, *, isolate: bool = False) -> None:
        self.path = Path(module_path)
        self._temp_dir: Path | None = None
        load_path = self.path
        if isolate:
            self._temp_dir = Path(
                tempfile.mkdtemp(
                    prefix="siac_rt6s_module_",
                    dir=os.fspath(_preferred_module_isolation_root(self.path)),
                )
            )
            load_path = self._temp_dir / self.path.name
            _copy_isolated_module(self.path, load_path)
        self._module = _load_extension_module(load_path)
        self._run_batch = self._module.sixs_f2py_run_batch

    def close(self) -> None:
        if self._temp_dir is not None:
            shutil.rmtree(self._temp_dir, ignore_errors=True)
            self._temp_dir = None

    def __del__(self) -> None:
        self.close()

    def run_batch(
        self,
        *,
        n_threads: int,
        month: int,
        day: int,
        atmospheric_mode: int,
        atmospheric_columns_mode: int,
        radiosonde_altitude_km: np.ndarray,
        radiosonde_pressure_mb: np.ndarray,
        radiosonde_temperature_k: np.ndarray,
        radiosonde_water_g_m3: np.ndarray,
        radiosonde_ozone_g_m3: np.ndarray,
        aerosol_mode: int,
        aerosol_mixture: np.ndarray,
        aerosol_distribution_rmin: float,
        aerosol_distribution_rmax: float,
        aerosol_distribution_component_count: int,
        aerosol_distribution_x1: np.ndarray,
        aerosol_distribution_x2: np.ndarray,
        aerosol_distribution_x3: np.ndarray,
        aerosol_distribution_cij: np.ndarray,
        aerosol_distribution_rn: np.ndarray,
        aerosol_distribution_ri: np.ndarray,
        aerosol_sun_count: int,
        aerosol_sun_radius: np.ndarray,
        aerosol_sun_dvlogr: np.ndarray,
        aerosol_layer_count: int,
        aerosol_layer_height: np.ndarray,
        aerosol_layer_aot: np.ndarray,
        aerosol_layer_type: np.ndarray,
        reference_reflectance: float,
        spectral_wlinf: float,
        spectral_wlsup: float,
        spectral_response: np.ndarray,
        aerosol_model_path: str,
        surface_inhomo: int,
        surface_idirec: int,
        surface_target_mode: int,
        surface_target_constant: float,
        surface_target_spectrum: np.ndarray,
        surface_env_mode: int,
        surface_env_constant: float,
        surface_env_spectrum: np.ndarray,
        surface_radius_km: float,
        surface_brdf_model: int,
        surface_brdf_params: np.ndarray,
        surface_brdf_options: np.ndarray,
        surface_brdf_struct: np.ndarray,
        surface_brdf_optics: np.ndarray,
        surface_brdf_table_solar: np.ndarray,
        surface_brdf_table_view: np.ndarray,
        surface_brdf_spherical_albedo: float,
        surface_brdf_directional_reflectance: float,
        atmospheric_correction_mode: int,
        atmospheric_correction_value: float,
        sza_deg: np.ndarray,
        saa_deg: np.ndarray,
        vza_deg: np.ndarray,
        vaa_deg: np.ndarray,
        aot550: np.ndarray,
        tcwv_cm: np.ndarray,
        tco3_atmcm: np.ndarray,
        elevation_km: np.ndarray,
    ) -> _NativeBatchResult:
        n_cases = int(sza_deg.size)
        if n_cases == 0:
            return _NativeBatchResult(
                outputs=_empty_output_bundle(0), status=np.zeros(0, dtype=np.int32)
            )

        native_output_matrix = np.empty(
            (len(SIXS_NATIVE_OUTPUT_NAMES), n_cases),
            dtype=np.float64,
            order="F",
        )
        status = np.zeros(n_cases, dtype=np.int32)

        self._run_batch(
            int(month),
            int(day),
            int(atmospheric_mode),
            int(atmospheric_columns_mode),
            np.asarray(radiosonde_altitude_km, dtype=np.float64),
            np.asarray(radiosonde_pressure_mb, dtype=np.float64),
            np.asarray(radiosonde_temperature_k, dtype=np.float64),
            np.asarray(radiosonde_water_g_m3, dtype=np.float64),
            np.asarray(radiosonde_ozone_g_m3, dtype=np.float64),
            int(aerosol_mode),
            np.asarray(aerosol_mixture, dtype=np.float64),
            float(aerosol_distribution_rmin),
            float(aerosol_distribution_rmax),
            int(aerosol_distribution_component_count),
            np.asarray(aerosol_distribution_x1, dtype=np.float64),
            np.asarray(aerosol_distribution_x2, dtype=np.float64),
            np.asarray(aerosol_distribution_x3, dtype=np.float64),
            np.asarray(aerosol_distribution_cij, dtype=np.float64),
            np.asarray(aerosol_distribution_rn, dtype=np.float64),
            np.asarray(aerosol_distribution_ri, dtype=np.float64),
            int(aerosol_sun_count),
            np.asarray(aerosol_sun_radius, dtype=np.float64),
            np.asarray(aerosol_sun_dvlogr, dtype=np.float64),
            int(aerosol_layer_count),
            np.asarray(aerosol_layer_height, dtype=np.float64),
            np.asarray(aerosol_layer_aot, dtype=np.float64),
            np.asarray(aerosol_layer_type, dtype=np.int32),
            float(reference_reflectance),
            float(spectral_wlinf),
            float(spectral_wlsup),
            np.asarray(spectral_response, dtype=np.float64),
            str(aerosol_model_path),
            int(surface_inhomo),
            int(surface_idirec),
            int(surface_target_mode),
            float(surface_target_constant),
            np.asarray(surface_target_spectrum, dtype=np.float64),
            int(surface_env_mode),
            float(surface_env_constant),
            np.asarray(surface_env_spectrum, dtype=np.float64),
            float(surface_radius_km),
            int(surface_brdf_model),
            np.asarray(surface_brdf_params, dtype=np.float64),
            np.asarray(surface_brdf_options, dtype=np.int32),
            np.asarray(surface_brdf_struct, dtype=np.float64),
            np.asarray(surface_brdf_optics, dtype=np.float64),
            np.asarray(surface_brdf_table_solar, dtype=np.float64),
            np.asarray(surface_brdf_table_view, dtype=np.float64),
            float(surface_brdf_spherical_albedo),
            float(surface_brdf_directional_reflectance),
            int(atmospheric_correction_mode),
            float(atmospheric_correction_value),
            np.asarray(sza_deg, dtype=np.float64),
            np.asarray(saa_deg, dtype=np.float64),
            np.asarray(vza_deg, dtype=np.float64),
            np.asarray(vaa_deg, dtype=np.float64),
            np.asarray(aot550, dtype=np.float64),
            np.asarray(tcwv_cm, dtype=np.float64),
            np.asarray(tco3_atmcm, dtype=np.float64),
            np.asarray(elevation_km, dtype=np.float64),
            int(max(1, n_threads)),
            native_output_matrix,
            status,
            n_cases=int(n_cases),
        )
        outputs = _empty_output_bundle(n_cases)
        for row_index, name in enumerate(SIXS_NATIVE_OUTPUT_NAMES):
            outputs[name] = np.ascontiguousarray(
                native_output_matrix[row_index, :], dtype=np.float64
            )
        outputs["xbp"] = np.ascontiguousarray(outputs["xb"], dtype=np.float64)
        outputs["xcp"] = np.ascontiguousarray(outputs["xc"], dtype=np.float64)
        return _NativeBatchResult(
            outputs={
                name: np.ascontiguousarray(values, dtype=np.float64)
                for name, values in outputs.items()
            },
            status=np.ascontiguousarray(status, dtype=np.int32),
        )


class SixSNativeRunner:
    """High-level xarray runner around the native 6SV2.1 Python extension."""

    def __init__(
        self,
        *,
        sixs_config: SixSAlgorithmConfig,
        rt_setup: RTSetupConfig,
        sensor_config: SensorConfig | None = None,
    ) -> None:
        self._config = sixs_config
        self._rt_setup = rt_setup
        self._sensor_config = sensor_config
        self._module_path = ensure_native_sixs_module(sixs_config)
        self._native_threads = int(sixs_config.native_threads or _default_native_threads())
        self._observation_time: datetime | None = None
        self._openmp_session: _SixSExtensionModule | None = None
        self._worker_sessions: list[_SixSExtensionModule] = []
        self._worker_sessions_available: bool | None = None
        #: Single-shot cache for a prebuilt joint grid-search LUT (wave 18).
        #: ``preload_joint_grid_search_lut`` writes here; the next call to
        #: ``build_joint_grid_search_lut`` reads + clears so subsequent
        #: requests (e.g. across solver stages) recompute fresh inputs.
        self._cached_joint_lut: JointGridSearchLUT | None = None
        self._cached_joint_lut_signature: tuple | None = None
        # Persistent on-disk cache of the native scene-LUT GRID batch (the
        # Fortran RT compute). Keyed by the exact batch inputs + module identity,
        # so a re-run of the same scene reuses the small grid coefficients
        # instead of recomputing. Only the bounded grid batch is cached; the
        # cheap per-pixel interpolation is not. ``None`` disables it.
        self._run_cache_dir: Path | None = resolve_run_cache_dir(
            getattr(sixs_config, "run_cache_dir", None),
            subpath="rt6s/run_cache",
            enabled=bool(getattr(sixs_config, "run_cache_enabled", True)),
        )

    def set_observation_time(self, observation_time: datetime | None) -> None:
        self._observation_time = observation_time

    def close(self) -> None:
        if self._openmp_session is not None:
            self._openmp_session.close()
            self._openmp_session = None
        for session in self._worker_sessions:
            session.close()
        self._worker_sessions = []
        self._worker_sessions_available = None

    def __del__(self) -> None:
        self.close()

    def compute_coefficients(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        output_variables: tuple[str, ...] | list[str] | None = None,
    ) -> dict[str, xr.DataArray]:
        return self._compute_band_outputs_multi(
            geometry=geometry,
            atmo_state=atmo_state,
            bands=[band],
            output_variables=output_variables,
        )[0]

    def compute_coefficients_multi(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        compute_jacobian: bool = False,
        output_variables: tuple[str, ...] | list[str] | None = None,
    ) -> list[dict[str, xr.DataArray]]:
        if compute_jacobian:
            raise NotImplementedError("The native 6S backend does not expose Jacobians.")
        return self._compute_band_outputs_multi(
            geometry=geometry,
            atmo_state=atmo_state,
            bands=bands,
            output_variables=output_variables,
        )

    def preload_scene_subset(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        output_variables: tuple[str, ...] | list[str] | None = None,
    ) -> None:
        if not bands:
            return None
        prepared = self._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
        if not np.any(prepared.valid_mask):
            return None
        plan = _build_scene_lut_plan(
            prepared.case_arrays,
            max_nodes_per_axis=int(self._config.scene_lut_max_nodes_per_axis),
            max_cases=int(self._config.scene_lut_max_cases),
        )
        if not _should_use_scene_lut(
            mode=str(self._config.mode),
            direct_case_count=plan.direct_case_count,
            lut_case_count=plan.lut_case_count,
            min_pixels=int(self._config.scene_lut_min_pixels),
            required_speedup=float(self._config.scene_lut_required_speedup),
        ):
            return None
        selected_names = _to_native_output_names(
            output_variables or getattr(self._config, "output_variables", None)
        )
        for band in bands:
            self._run_scene_lut_batch(
                prepared=prepared, band=band, selected_names=selected_names, plan=plan
            )
        return None

    @staticmethod
    def _joint_lut_signature(
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
        output_variables: tuple[str, ...] | list[str] | None,
    ) -> tuple:
        """Return a hashable identity for a joint-LUT request.

        Wave 18 (opt 4): used to match a precomputed (preload-stage) joint
        LUT against the build call inside the solver. The signature folds
        in the call's shape-relevant inputs so a precompute that doesn't
        match the eventual solver request can't be reused incorrectly.
        """
        return (
            tuple(geometry.sza.shape),
            tuple(getattr(atmo_state.tco3, "shape", ())),
            tuple(np.asarray(aot_axis, dtype=np.float64).tolist()),
            tuple(np.asarray(tcwv_axis, dtype=np.float64).tolist()),
            tuple(band.name for band in bands),
            tuple(output_variables) if output_variables is not None else None,
        )

    def preload_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
        output_variables: tuple[str, ...] | list[str] | None = None,
    ) -> JointGridSearchLUT | None:
        """Eagerly build a joint LUT and cache it for a single later read.

        Designed to be invoked during the prior-fetch stage so the joint
        LUT's 6S work overlaps with M2/M3 wall-clock instead of running
        sequentially after them. The very next ``build_joint_grid_search_lut``
        call with matching inputs picks it up, then the cache is cleared
        (so staged solvers that rebuild with different inputs aren't
        served a stale entry).
        """
        joint = self.build_joint_grid_search_lut(
            geometry=geometry,
            atmo_state=atmo_state,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            bands=bands,
            output_variables=output_variables,
            _bypass_cache=True,
        )
        if joint is None:
            return None
        self._cached_joint_lut = joint
        self._cached_joint_lut_signature = self._joint_lut_signature(
            geometry=geometry,
            atmo_state=atmo_state,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            bands=bands,
            output_variables=output_variables,
        )
        logger.info(
            "Preloaded joint grid-search LUT in parallel with prior fetch (%d bands, %d cases).",
            len(bands),
            joint.plan.lut_case_count,
        )
        return joint

    def build_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
        output_variables: tuple[str, ...] | list[str] | None = None,
        _bypass_cache: bool = False,
    ) -> JointGridSearchLUT | None:
        """Build a joint scene-LUT covering an explicit (aot, tcwv) grid.

        Returns ``None`` when the optimization is disabled by config, the
        scene has no valid pixels, or the runner is configured to bypass the
        scene-LUT entirely (``sixs.mode == "direct"``). In all of those
        cases the caller should fall back to the per-candidate
        :meth:`compute_coefficients` path.
        """
        if not bands:
            return None
        if not bool(getattr(self._config, "joint_grid_search_lut_enabled", True)):
            return None
        if str(self._config.mode) == "direct":
            return None

        # Wave 18: serve a precomputed LUT if one matches this call exactly.
        # Single-shot — the cache is cleared on hit so staged solvers that
        # rebuild with different inputs aren't served a stale entry.
        if not _bypass_cache and self._cached_joint_lut is not None:
            sig = self._joint_lut_signature(
                geometry=geometry,
                atmo_state=atmo_state,
                aot_axis=aot_axis,
                tcwv_axis=tcwv_axis,
                bands=bands,
                output_variables=output_variables,
            )
            if sig == self._cached_joint_lut_signature:
                logger.info(
                    "Using preloaded joint grid-search LUT (%d cases, %d bands).",
                    self._cached_joint_lut.plan.lut_case_count,
                    self._cached_joint_lut.band_count,
                )
                cached = self._cached_joint_lut
                self._cached_joint_lut = None
                self._cached_joint_lut_signature = None
                return cached

        prepared = self._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
        if not np.any(prepared.valid_mask):
            return None
        max_nodes_per_axis = int(
            getattr(
                self._config,
                "joint_grid_search_lut_max_nodes_per_axis",
                self._config.scene_lut_max_nodes_per_axis,
            )
        )
        max_cases = int(
            getattr(
                self._config,
                "joint_grid_search_lut_max_cases",
                max(self._config.scene_lut_max_cases, 16384),
            )
        )
        plan = _build_joint_grid_search_lut_plan(
            prepared.case_arrays,
            aot_axis=np.asarray(aot_axis, dtype=np.float64),
            tcwv_axis=np.asarray(tcwv_axis, dtype=np.float64),
            max_nodes_per_axis=max_nodes_per_axis,
            max_cases=max_cases,
        )
        if plan.lut_case_count <= 0:
            return None
        selected_names = _to_native_output_names(
            output_variables or getattr(self._config, "output_variables", None)
        )
        logger.debug(
            "Building joint grid-search LUT: %d pixels, %d cases per band, %d bands.",
            plan.direct_case_count,
            plan.lut_case_count,
            len(bands),
        )
        band_outputs = self._run_joint_lut_bands(prepared=prepared, plan=plan, bands=bands)
        return JointGridSearchLUT(
            prepared=prepared,
            plan=plan,
            band_native_outputs=band_outputs,
            selected_names=selected_names,
        )

    def _run_joint_lut_bands(
        self,
        *,
        prepared: _PreparedSceneInputs,
        plan: _SceneLUTPlan,
        bands: list[SensorBand],
    ) -> list[_NativeBatchResult]:
        """Run the joint LUT's per-band 6S batches.

        When ``parallel_backend == "worker_libraries"`` we dispatch the
        per-band batches concurrently across isolated library sessions —
        each 6S call has a substantial serial portion that OpenMP can't
        cover, so band-level parallelism gives a meaningful wall-clock win
        on top of the per-call OpenMP scaling. With the OpenMP backend the
        runner already saturates cores for a single batch, so we keep the
        sequential band loop to avoid oversubscription.
        """
        backend = str(getattr(self._config, "parallel_backend", "openmp"))
        n_bands = len(bands)
        if backend != "worker_libraries" or n_bands <= 1:
            band_outputs: list[_NativeBatchResult] = []
            for band in bands:
                native_kwargs = self._band_native_kwargs(prepared, band)
                native_kwargs.update(plan.grid_case_arrays)
                band_outputs.append(self._run_native_batch(**native_kwargs))
            return band_outputs

        worker_count = min(n_bands, self._worker_library_count())
        sessions = self._ensure_worker_sessions(worker_count)
        if not sessions:
            band_outputs = []
            for band in bands:
                native_kwargs = self._band_native_kwargs(prepared, band)
                native_kwargs.update(plan.grid_case_arrays)
                band_outputs.append(self._run_native_batch(**native_kwargs))
            return band_outputs

        # Each session uses a fraction of the available OpenMP threads so
        # the sum across sessions stays close to (but does not exceed) the
        # physical core budget. The Fortran 6S kernel is largely serial
        # once you go past 4 threads per call anyway, so this allocation
        # captures most of the OpenMP win plus all of the band-parallelism.
        per_worker_threads = max(1, int(self._native_threads) // worker_count)
        results_by_band: dict[int, _NativeBatchResult] = {}

        # Pre-build per-band kwargs outside the executor to keep the worker
        # closures small and to surface band-spec errors deterministically.
        band_kwargs: list[dict[str, Any]] = []
        for band in bands:
            kwargs = self._band_native_kwargs(prepared, band)
            kwargs.update(plan.grid_case_arrays)
            band_kwargs.append(kwargs)

        def _run_single_band(session: _SixSExtensionModule, band_index: int) -> _NativeBatchResult:
            kwargs = band_kwargs[band_index]
            result = session.run_batch(n_threads=per_worker_threads, **kwargs)
            failed = result.status != 0
            if np.any(failed):
                for _name, values in result.outputs.items():
                    values[failed] = np.nan
                logger.warning(
                    "6S native runner returned %d non-zero status values for joint-LUT band #%d.",
                    int(np.count_nonzero(failed)),
                    band_index,
                )
            return result

        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            future_to_band = {
                executor.submit(_run_single_band, sessions[idx % worker_count], idx): idx
                for idx in range(n_bands)
            }
            for future in future_to_band:
                band_index = future_to_band[future]
                results_by_band[band_index] = future.result()

        return [results_by_band[idx] for idx in range(n_bands)]

    def _compute_band_outputs_multi(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        output_variables: tuple[str, ...] | list[str] | None,
    ) -> list[dict[str, xr.DataArray]]:
        if not bands:
            return []
        prepared = self._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
        selected_names = _to_native_output_names(output_variables)
        if not np.any(prepared.valid_mask):
            return [_nan_output_arrays(prepared.template, selected_names) for _ in bands]
        plan = _build_scene_lut_plan(
            prepared.case_arrays,
            max_nodes_per_axis=int(self._config.scene_lut_max_nodes_per_axis),
            max_cases=int(self._config.scene_lut_max_cases),
        )
        use_scene_lut = _should_use_scene_lut(
            mode=str(self._config.mode),
            direct_case_count=plan.direct_case_count,
            lut_case_count=plan.lut_case_count,
            min_pixels=int(self._config.scene_lut_min_pixels),
            required_speedup=float(self._config.scene_lut_required_speedup),
        )
        outputs: list[dict[str, xr.DataArray]] = []
        for band in bands:
            if use_scene_lut:
                native_outputs = self._run_scene_lut_batch(
                    prepared=prepared,
                    band=band,
                    selected_names=selected_names,
                    plan=plan,
                )
            else:
                native_outputs = self._run_native_batch(**self._band_native_kwargs(prepared, band))
            outputs.append(self._render_native_outputs(prepared, selected_names, native_outputs))
        return outputs

    def prepare_correction_scene(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> PreparedCorrectionScene:
        """Build a shared scene prep + LUT plan once for all bands in M6.

        Wave 18e: returns a frozen :class:`PreparedCorrectionScene` that
        :meth:`compute_coefficients_with_prepared` can consume for each
        band without redoing the per-band-redundant prep work (scene
        inputs assembly, valid-mask construction, scene-LUT plan).

        The object is safe to share across threads — every consumer reads
        from it but never mutates it.
        """
        prepared = self._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
        plan = _build_scene_lut_plan(
            prepared.case_arrays,
            max_nodes_per_axis=int(self._config.scene_lut_max_nodes_per_axis),
            max_cases=int(self._config.scene_lut_max_cases),
        )
        use_scene_lut = _should_use_scene_lut(
            mode=str(self._config.mode),
            direct_case_count=plan.direct_case_count,
            lut_case_count=plan.lut_case_count,
            min_pixels=int(self._config.scene_lut_min_pixels),
            required_speedup=float(self._config.scene_lut_required_speedup),
        )
        return PreparedCorrectionScene(
            prepared=prepared,
            plan=plan,
            use_scene_lut=use_scene_lut,
        )

    def compute_coefficients_with_prepared(
        self,
        *,
        prepared_scene: PreparedCorrectionScene,
        band: SensorBand,
        output_variables: tuple[str, ...] | list[str] | None = None,
    ) -> dict[str, xr.DataArray]:
        """Run RT for a single band using a pre-computed scene prep.

        Wave 18e: in M6 correction this is called per band inside the
        per-band executor — the 6S kernel work is unavoidably band-specific
        (each band has its own spectral response), but the scene
        preparation that wraps it is now done once via
        :meth:`prepare_correction_scene` and shared.
        """
        selected_names = _to_native_output_names(output_variables)
        if not np.any(prepared_scene.prepared.valid_mask):
            return _nan_output_arrays(prepared_scene.prepared.template, selected_names)
        if prepared_scene.use_scene_lut:
            native_outputs = self._run_scene_lut_batch(
                prepared=prepared_scene.prepared,
                band=band,
                selected_names=selected_names,
                plan=prepared_scene.plan,
            )
        else:
            native_outputs = self._run_native_batch(
                **self._band_native_kwargs(prepared_scene.prepared, band)
            )
        return self._render_native_outputs(prepared_scene.prepared, selected_names, native_outputs)

    def _prepare_scene_inputs(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> _PreparedSceneInputs:
        template = atmo_state.aot
        geometry_on_grid = _resample_geometry_to_template(geometry, template)

        sza_deg = xr.DataArray(
            np.degrees(geometry_on_grid.sza.values), dims=template.dims, coords=template.coords
        )
        saa_deg = _wrap_azimuth_degrees(geometry_on_grid.saa, name="saa_deg")
        vza_deg = xr.DataArray(
            np.degrees(geometry_on_grid.vza.values), dims=template.dims, coords=template.coords
        )
        vaa_deg = _wrap_azimuth_degrees(geometry_on_grid.vaa, name="vaa_deg")
        valid = (
            np.isfinite(sza_deg.values)
            & np.isfinite(saa_deg.values)
            & np.isfinite(vza_deg.values)
            & np.isfinite(vaa_deg.values)
            & np.isfinite(atmo_state.aot.values)
            & np.isfinite(atmo_state.tcwv.values)
            & np.isfinite(atmo_state.tco3.values)
            & np.isfinite(atmo_state.elevation.values)
        )
        month = (
            self._observation_time.month
            if self._observation_time is not None
            else self._config.month
        )
        day = self._observation_time.day if self._observation_time is not None else self._config.day
        atmospheric_mode, radiosonde = _resolve_atmospheric_inputs(self._rt_setup, month=month)
        aerosol_mode, aerosol_inputs = _resolve_aerosol_inputs(self._rt_setup)
        surface_inputs = _resolve_surface_inputs(self._rt_setup)
        atmospheric_correction_mode, atmospheric_correction_value = (
            _resolve_atmospheric_correction_inputs(self._rt_setup)
        )
        reference_reflectance = (
            float(self._rt_setup.reference_reflectance)
            if self._rt_setup.reference_reflectance is not None
            else 0.1
        )
        aerosol_model_path = ""
        if self._rt_setup.aerosol is not None and self._rt_setup.aerosol.model_path is not None:
            aerosol_model_path = _encode_aerosol_model_path(self._rt_setup.aerosol.model_path)
        return _PreparedSceneInputs(
            template=template,
            valid_mask=valid,
            flat_valid_mask=np.ravel(valid),
            common_kwargs={
                "month": month,
                "day": day,
                "atmospheric_mode": atmospheric_mode,
                "atmospheric_columns_mode": _resolve_atmospheric_columns_mode(self._rt_setup),
                "radiosonde_altitude_km": np.ascontiguousarray(
                    radiosonde["altitude_km"], dtype=np.float64
                ),
                "radiosonde_pressure_mb": np.ascontiguousarray(
                    radiosonde["pressure_mb"], dtype=np.float64
                ),
                "radiosonde_temperature_k": np.ascontiguousarray(
                    radiosonde["temperature_k"], dtype=np.float64
                ),
                "radiosonde_water_g_m3": np.ascontiguousarray(
                    radiosonde["water_g_m3"], dtype=np.float64
                ),
                "radiosonde_ozone_g_m3": np.ascontiguousarray(
                    radiosonde["ozone_g_m3"], dtype=np.float64
                ),
                "aerosol_mode": aerosol_mode,
                "aerosol_mixture": np.ascontiguousarray(
                    aerosol_inputs["mixture"], dtype=np.float64
                ),
                "aerosol_distribution_rmin": float(aerosol_inputs["dist_rmin"]),
                "aerosol_distribution_rmax": float(aerosol_inputs["dist_rmax"]),
                "aerosol_distribution_component_count": int(aerosol_inputs["dist_component_count"]),
                "aerosol_distribution_x1": np.ascontiguousarray(
                    aerosol_inputs["dist_x1"], dtype=np.float64
                ),
                "aerosol_distribution_x2": np.ascontiguousarray(
                    aerosol_inputs["dist_x2"], dtype=np.float64
                ),
                "aerosol_distribution_x3": np.ascontiguousarray(
                    aerosol_inputs["dist_x3"], dtype=np.float64
                ),
                "aerosol_distribution_cij": np.ascontiguousarray(
                    aerosol_inputs["dist_cij"], dtype=np.float64
                ),
                "aerosol_distribution_rn": np.ascontiguousarray(
                    aerosol_inputs["dist_rn"], dtype=np.float64
                ),
                "aerosol_distribution_ri": np.ascontiguousarray(
                    aerosol_inputs["dist_ri"], dtype=np.float64
                ),
                "aerosol_sun_count": int(aerosol_inputs["sun_count"]),
                "aerosol_sun_radius": np.ascontiguousarray(
                    aerosol_inputs["sun_radius"], dtype=np.float64
                ),
                "aerosol_sun_dvlogr": np.ascontiguousarray(
                    aerosol_inputs["sun_dvlogr"], dtype=np.float64
                ),
                "aerosol_layer_count": int(aerosol_inputs["layer_count"]),
                "aerosol_layer_height": np.ascontiguousarray(
                    aerosol_inputs["layer_height"], dtype=np.float64
                ),
                "aerosol_layer_aot": np.ascontiguousarray(
                    aerosol_inputs["layer_aot"], dtype=np.float64
                ),
                "aerosol_layer_type": np.ascontiguousarray(
                    aerosol_inputs["layer_type"], dtype=np.int32
                ),
                "reference_reflectance": reference_reflectance,
                "aerosol_model_path": aerosol_model_path,
                "surface_inhomo": int(surface_inputs["inhomo"]),
                "surface_idirec": int(surface_inputs["idirec"]),
                "surface_target_mode": int(surface_inputs["target_mode"]),
                "surface_target_constant": float(surface_inputs["target_constant"]),
                "surface_target_spectrum": np.ascontiguousarray(
                    surface_inputs["target_spectrum"], dtype=np.float64
                ),
                "surface_env_mode": int(surface_inputs["env_mode"]),
                "surface_env_constant": float(surface_inputs["env_constant"]),
                "surface_env_spectrum": np.ascontiguousarray(
                    surface_inputs["env_spectrum"], dtype=np.float64
                ),
                "surface_radius_km": float(surface_inputs["radius_km"]),
                "surface_brdf_model": int(surface_inputs["brdf_model"]),
                "surface_brdf_params": np.ascontiguousarray(
                    surface_inputs["brdf_params"], dtype=np.float64
                ),
                "surface_brdf_options": np.ascontiguousarray(
                    surface_inputs["brdf_options"], dtype=np.int32
                ),
                "surface_brdf_struct": np.ascontiguousarray(
                    surface_inputs["brdf_struct"], dtype=np.float64
                ),
                "surface_brdf_optics": np.ascontiguousarray(
                    surface_inputs["brdf_optics"], dtype=np.float64
                ),
                "surface_brdf_table_solar": np.ascontiguousarray(
                    surface_inputs["brdf_table_solar"], dtype=np.float64
                ),
                "surface_brdf_table_view": np.ascontiguousarray(
                    surface_inputs["brdf_table_view"], dtype=np.float64
                ),
                "surface_brdf_spherical_albedo": float(surface_inputs["brdf_spherical_albedo"]),
                "surface_brdf_directional_reflectance": float(
                    surface_inputs["brdf_directional_reflectance"]
                ),
                "atmospheric_correction_mode": int(atmospheric_correction_mode),
                "atmospheric_correction_value": float(atmospheric_correction_value),
            },
            case_arrays={
                "sza_deg": _as_float64(sza_deg.values[valid]),
                "saa_deg": _as_float64(saa_deg.values[valid]),
                "vza_deg": _as_float64(vza_deg.values[valid]),
                "vaa_deg": _as_float64(vaa_deg.values[valid]),
                "aot550": _as_float64(atmo_state.aot.values[valid]),
                "tcwv_cm": _as_float64(atmo_state.tcwv.values[valid]),
                "tco3_atmcm": _as_float64(atmo_state.tco3.values[valid]),
                "elevation_km": _as_float64(atmo_state.elevation.values[valid]),
            },
        )

    def _band_native_kwargs(
        self, prepared: _PreparedSceneInputs, band: SensorBand
    ) -> dict[str, Any]:
        response, wlinf, wlsup = _build_spectral_response(band)
        kwargs = dict(prepared.common_kwargs)
        kwargs.update(
            {
                "spectral_wlinf": wlinf,
                "spectral_wlsup": wlsup,
                "spectral_response": response,
            }
        )
        kwargs.update(prepared.case_arrays)
        return kwargs

    def _render_native_outputs(
        self,
        prepared: _PreparedSceneInputs,
        selected_names: tuple[str, ...],
        native_outputs: _NativeBatchResult,
    ) -> dict[str, xr.DataArray]:
        outputs_by_name: dict[str, xr.DataArray] = {}
        for name in selected_names:
            full = np.full(prepared.template.size, np.nan, dtype=np.float64)
            full[prepared.flat_valid_mask] = native_outputs.outputs[name]
            outputs_by_name[name] = _as_output_array(full, prepared.template, name)
        return outputs_by_name

    def _run_scene_lut_batch(
        self,
        *,
        prepared: _PreparedSceneInputs,
        band: SensorBand,
        selected_names: tuple[str, ...],
        plan: _SceneLUTPlan,
    ) -> _NativeBatchResult:
        native_kwargs = self._band_native_kwargs(prepared, band)
        native_kwargs.update(plan.grid_case_arrays)
        lut_outputs = self._run_native_batch_cached(native_kwargs)
        interpolated_outputs = _interpolate_scene_lut_outputs(
            plan,
            lut_outputs,
            prepared.case_arrays,
            selected_names,
        )
        status = np.zeros(plan.direct_case_count, dtype=np.int32)
        return _NativeBatchResult(outputs=interpolated_outputs, status=status)

    def _run_native_batch_cached(self, native_kwargs: dict[str, Any]) -> _NativeBatchResult:
        """Run the native grid batch, served from disk when the inputs repeat.

        The grid batch is a pure, deterministic function of ``native_kwargs`` and
        the compiled module, so a SHA-256 of (all kwargs + module identity) is an
        exact key. A hit returns the cached grid coefficients without invoking
        the Fortran kernel.
        """
        key = self._native_grid_cache_key(native_kwargs)
        if key is not None:
            cached = self._load_native_grid(key)
            if cached is not None:
                return cached
        result = self._run_native_batch(**native_kwargs)
        if key is not None:
            self._store_native_grid(key, result)
        return result

    def _native_grid_cache_key(self, native_kwargs: dict[str, Any]) -> str | None:
        if self._run_cache_dir is None:
            return None
        h = hashlib.sha256()
        h.update(b"sixs-grid-batch-cache-v1\x00")
        # Module identity: a rebuilt 6S binary must invalidate the cache.
        h.update(str(self._module_path).encode("utf-8"))
        try:
            stat = Path(self._module_path).stat()
            h.update(f"{stat.st_mtime_ns}:{stat.st_size}".encode())
        except OSError:
            pass
        for name in sorted(native_kwargs):
            h.update(name.encode("utf-8"))
            h.update(b"\x00")
            arr = np.asarray(native_kwargs[name])
            if arr.dtype == object:
                h.update(repr(native_kwargs[name]).encode("utf-8"))
            else:
                h.update(str(arr.dtype).encode("utf-8"))
                h.update(str(arr.shape).encode("utf-8"))
                h.update(np.ascontiguousarray(arr).tobytes())
            h.update(b"\x00")
        return h.hexdigest()

    def _load_native_grid(self, key: str) -> _NativeBatchResult | None:
        if self._run_cache_dir is None:
            return None

        def _read(path: Path) -> _NativeBatchResult:
            with np.load(path) as data:
                status = np.asarray(data["__status__"])
                outputs = {
                    str(name)[4:]: np.asarray(data[name])
                    for name in data.files
                    if str(name).startswith("out_")
                }
            return _NativeBatchResult(outputs=outputs, status=status)

        return load_cache_entry(self._run_cache_dir / f"{key}.npz", _read)

    def _store_native_grid(self, key: str, result: _NativeBatchResult) -> None:
        if self._run_cache_dir is None:
            return
        payload = {f"out_{name}": np.asarray(arr) for name, arr in result.outputs.items()}
        payload["__status__"] = np.asarray(result.status)

        def _write(tmp: Path) -> None:
            # File handle avoids np.savez's ".npz" auto-suffix on the temp.
            with tmp.open("wb") as handle:
                np.savez(handle, **payload)

        store_cache_entry(self._run_cache_dir / f"{key}.npz", _write)

    def _ensure_openmp_session(self) -> _SixSExtensionModule:
        if self._openmp_session is None:
            try:
                self._openmp_session = _SixSExtensionModule(self._module_path, isolate=True)
            except PermissionError:
                logger.warning(
                    "Could not create an isolated native 6S OpenMP session; loading the compiled module in place."
                )
                self._openmp_session = _SixSExtensionModule(self._module_path, isolate=False)
        return self._openmp_session

    def _ensure_worker_sessions(self, worker_count: int) -> list[_SixSExtensionModule]:
        if self._worker_sessions_available is False:
            return []
        if len(self._worker_sessions) == worker_count and worker_count > 0:
            return self._worker_sessions
        for session in self._worker_sessions:
            session.close()
        try:
            self._worker_sessions = [
                _SixSExtensionModule(self._module_path, isolate=True)
                for _ in range(max(1, worker_count))
            ]
        except PermissionError:
            self._worker_sessions = []
            self._worker_sessions_available = False
            logger.warning(
                "Could not create isolated native 6S worker sessions; falling back to OpenMP execution."
            )
            return []
        self._worker_sessions_available = True
        return self._worker_sessions

    def _worker_library_count(self) -> int:
        configured = getattr(self._config, "worker_libraries", None)
        if configured is not None:
            return max(1, int(configured))
        return max(1, int(self._native_threads))

    def _run_native_batch(self, **kwargs: Any) -> _NativeBatchResult:
        backend = str(getattr(self._config, "parallel_backend", "openmp"))
        if backend == "worker_libraries":
            result = self._run_native_batch_worker_libraries(**kwargs)
        else:
            extension = self._ensure_openmp_session()
            result = extension.run_batch(n_threads=self._native_threads, **kwargs)
        failed = result.status != 0
        if np.any(failed):
            for _name, values in result.outputs.items():
                values[failed] = np.nan
            logger.warning(
                "6S native runner returned %d non-zero status values.",
                int(np.count_nonzero(failed)),
            )
        return result

    def _run_native_batch_worker_libraries(self, **kwargs: Any) -> _NativeBatchResult:
        n_cases = int(np.asarray(kwargs["sza_deg"]).size)
        if n_cases == 0:
            return _NativeBatchResult(
                outputs=_empty_output_bundle(0), status=np.zeros(0, dtype=np.int32)
            )
        chunk_size = max(1, int(getattr(self._config, "chunk_size", 4096)))
        chunks = [
            (start, min(start + chunk_size, n_cases)) for start in range(0, n_cases, chunk_size)
        ]
        worker_count = min(len(chunks), self._worker_library_count())
        sessions = self._ensure_worker_sessions(worker_count)
        if not sessions:
            extension = self._ensure_openmp_session()
            return extension.run_batch(n_threads=self._native_threads, **kwargs)
        outputs = _empty_output_bundle(n_cases)
        status = np.zeros(n_cases, dtype=np.int32)
        assignments: list[list[tuple[int, int]]] = [[] for _ in range(worker_count)]
        for idx, chunk in enumerate(chunks):
            assignments[idx % worker_count].append(chunk)

        def _run_assigned(
            session: _SixSExtensionModule, assigned: list[tuple[int, int]]
        ) -> list[tuple[int, int, _NativeBatchResult]]:
            completed: list[tuple[int, int, _NativeBatchResult]] = []
            for start, stop in assigned:
                result = session.run_batch(n_threads=1, **_slice_case_kwargs(kwargs, start, stop))
                completed.append((start, stop, result))
            return completed

        if worker_count == 1:
            completed_runs = [_run_assigned(sessions[0], assignments[0])]
        else:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                completed_runs = list(
                    executor.map(
                        _run_assigned,
                        sessions,
                        assignments,
                    )
                )

        for worker_runs in completed_runs:
            for start, stop, result in worker_runs:
                status[start:stop] = result.status
                for name, values in result.outputs.items():
                    outputs[name][start:stop] = values
        return _NativeBatchResult(outputs=outputs, status=status)

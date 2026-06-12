"""Compiled-module loading and isolated sessions for the native 6S extension."""

from __future__ import annotations

import ctypes
import importlib.machinery
import importlib.util
import logging
import os
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from siac.sixs_outputs import SIXS_NATIVE_OUTPUT_NAMES, SIXS_OUTPUT_VARIABLE_CHOICES

logger = logging.getLogger(__name__)

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


def _empty_output_bundle(length: int) -> dict[str, np.ndarray]:
    return {name: np.empty(length, dtype=np.float64) for name in SIXS_OUTPUT_VARIABLE_CHOICES}


@dataclass(frozen=True)
class _NativeBatchResult:
    outputs: dict[str, np.ndarray]
    status: np.ndarray


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

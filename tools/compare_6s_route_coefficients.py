"""Compare xap/xbp/xcp from direct 6S, native 6S scene-LUT, and remote ZIP/Zarr LUT."""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import build_native_sixs_module
from siac.algorithms.rt.direct.sixs_native import _build_scene_lut_plan
from siac.algorithms.rt.lut import DEFAULT_LUT_URL, ZarrLUTBackend
from siac.config import SixSAlgorithmConfig
from siac.config.schema import (
    SixSAtmosphericCorrectionConfig,
    SixSAtmosphericProfile,
    SixSSpectralReflectanceConfig,
    SixSSurfaceConfig,
)
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("tmp/6s_upstream"), help="Path to unpacked 6SV2.1 sources.")
    parser.add_argument("--build-dir", type=Path, default=Path("tmp/6s_route_coefficients"), help="Build root for the native 6S module.")
    parser.add_argument("--module-path", type=Path, default=None, help="Optional prebuilt native 6S Python extension.")
    parser.add_argument("--compiler", default="gfortran", help="Fortran compiler executable or path.")
    parser.add_argument("--native-threads", type=int, default=4, help="OpenMP thread count for native 6S.")
    parser.add_argument("--size", type=int, default=32, help="Square scene size in pixels.")
    parser.add_argument(
        "--atmospheric-profile",
        default="us_standard_62",
        help="Native 6S atmospheric profile shape for the local direct/scene-LUT runs.",
    )
    parser.add_argument(
        "--atmospheric-columns-mode",
        choices=("input_columns", "profile_default"),
        default="input_columns",
        help="Whether local 6S scales the selected atmospheric profile to the scene tcwv/tco3 inputs.",
    )
    parser.add_argument("--tcwv", type=float, default=2.0, help="Scene total column water vapour (g/cm^2).")
    parser.add_argument("--tco3", type=float, default=0.3, help="Scene total column ozone (cm-atm).")
    parser.add_argument(
        "--scenario",
        choices=("smooth_geometry", "mixed_geometry"),
        default="smooth_geometry",
        help="Synthetic scene pattern for geometry/AOT/elevation variation.",
    )
    parser.add_argument("--lut-path", default=DEFAULT_LUT_URL, help="Path or URL to the remote/local Zarr LUT.")
    parser.add_argument(
        "--bands",
        default="B02,B03,B04,B08",
        help="Comma-separated synthetic band names to compare.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tmp/6s_route_coefficients/report.json"),
        help="Where to write the JSON report.",
    )
    return parser.parse_args()


def _band_catalog() -> dict[str, SensorBand]:
    return {
        "B02": SensorBand(name="B02", center_wavelength=490.0, bandwidth=20.0, resolution=10.0, band_index=1),
        "B03": SensorBand(name="B03", center_wavelength=560.0, bandwidth=35.0, resolution=10.0, band_index=2),
        "B04": SensorBand(name="B04", center_wavelength=665.0, bandwidth=30.0, resolution=10.0, band_index=3),
        "B08": SensorBand(name="B08", center_wavelength=842.0, bandwidth=115.0, resolution=10.0, band_index=7),
    }


def _select_bands(request: str) -> list[SensorBand]:
    catalog = _band_catalog()
    names = [item.strip().upper() for item in request.split(",") if item.strip()]
    bands: list[SensorBand] = []
    for name in names:
        if name not in catalog:
            raise SystemExit(f"Unknown band {name!r}; choose from {', '.join(sorted(catalog))}.")
        bands.append(catalog[name])
    if not bands:
        raise SystemExit("No bands selected.")
    return bands


def _linspace_grid(base: float, delta_x: float, delta_y: float, shape: tuple[int, int]) -> np.ndarray:
    y = np.linspace(0.0, 1.0, shape[0], dtype=np.float32)
    x = np.linspace(0.0, 1.0, shape[1], dtype=np.float32)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return np.asarray(base + delta_x * xx + delta_y * yy, dtype=np.float32)


def _build_scene(
    shape: tuple[int, int],
    scenario: str,
    *,
    tcwv_cm: float,
    tco3_atmcm: float,
) -> tuple[GeometryAngles, AtmosphericState]:
    if scenario == "smooth_geometry":
        sza = _linspace_grid(28.0, 2.0, 1.0, shape)
        saa = _linspace_grid(150.0, 4.0, 2.0, shape)
        vza = _linspace_grid(5.0, 0.8, 0.6, shape)
        vaa = _linspace_grid(110.0, 3.0, 2.0, shape)
        aot = _linspace_grid(0.12, 0.08, 0.04, shape)
        elevation = _linspace_grid(0.15, 0.05, 0.03, shape)
    else:
        sza = _linspace_grid(28.0, 6.0, 3.0, shape)
        saa = _linspace_grid(150.0, 18.0, 12.0, shape)
        vza = _linspace_grid(5.0, 2.5, 1.5, shape)
        vaa = _linspace_grid(110.0, 15.0, 10.0, shape)
        aot = _linspace_grid(0.12, 0.20, 0.12, shape)
        elevation = _linspace_grid(0.15, 0.25, 0.12, shape)

    tcwv = np.full(shape, float(tcwv_cm), dtype=np.float32)
    tco3 = np.full(shape, float(tco3_atmcm), dtype=np.float32)
    geometry = GeometryAngles.from_degrees(
        xr.DataArray(sza, dims=("y", "x")),
        xr.DataArray(saa, dims=("y", "x")),
        xr.DataArray(vza, dims=("y", "x")),
        xr.DataArray(vaa, dims=("y", "x")),
    )
    atmo_state = AtmosphericState(
        aot=xr.DataArray(aot, dims=("y", "x")),
        tcwv=xr.DataArray(tcwv, dims=("y", "x")),
        tco3=xr.DataArray(tco3, dims=("y", "x")),
        aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.005, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(elevation, dims=("y", "x")),
    )
    return geometry, atmo_state


def _sixs_config(
    *,
    source_dir: Path,
    build_dir: Path,
    module_path: Path | None,
    compiler: str,
    native_threads: int,
    mode: str,
    atmospheric_profile: SixSAtmosphericProfile,
    atmospheric_columns_mode: str,
) -> SixSAlgorithmConfig:
    return SixSAlgorithmConfig(
        source_dir=source_dir,
        build_dir=build_dir,
        module_path=module_path,
        auto_build=module_path is None,
        compiler=compiler,
        build_profile="release",
        mode=mode,
        parallel_backend="openmp",
        native_threads=native_threads,
        month=7,
        day=12,
        atmospheric_profile=atmospheric_profile,
        atmospheric_columns_mode=atmospheric_columns_mode,
        aerosol_profile="continental",
        reference_reflectance=0.1,
        surface=SixSSurfaceConfig(
            mode="homogeneous_lambertian",
            target=SixSSpectralReflectanceConfig(kind="constant", constant=0.1),
        ),
        atmospheric_correction=SixSAtmosphericCorrectionConfig(
            mode="lambertian_reflectance",
            value=0.1,
        ),
        output_variables=("xap", "xbp", "xcp"),
    )


def _coeff_summary(coeffs: RTCoefficients) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    for name in ("xap", "xbp", "xcp"):
        values = np.asarray(coeffs.get_output(name).values, dtype=np.float64)
        finite = values[np.isfinite(values)]
        summary[name] = {
            "min": float(np.min(finite)) if finite.size else float("nan"),
            "mean": float(np.mean(finite)) if finite.size else float("nan"),
            "max": float(np.max(finite)) if finite.size else float("nan"),
        }
    return summary


def _comparison_stats(reference: xr.DataArray, candidate: xr.DataArray) -> dict[str, float]:
    ref = np.asarray(reference.values, dtype=np.float64)
    cand = np.asarray(candidate.values, dtype=np.float64)
    mask = np.isfinite(ref) & np.isfinite(cand)
    if not np.any(mask):
        return {
            "valid_fraction": 0.0,
            "mean_abs": float("nan"),
            "median_abs": float("nan"),
            "max_abs": float("nan"),
            "rmse": float("nan"),
            "mean_rel": float("nan"),
            "max_rel": float("nan"),
        }

    diff = cand[mask] - ref[mask]
    abs_diff = np.abs(diff)
    denom = np.maximum(np.abs(ref[mask]), 1.0e-12)
    rel = abs_diff / denom
    return {
        "valid_fraction": float(mask.mean()),
        "mean_abs": float(np.mean(abs_diff)),
        "median_abs": float(np.median(abs_diff)),
        "max_abs": float(np.max(abs_diff)),
        "rmse": float(np.sqrt(np.mean(diff * diff))),
        "mean_rel": float(np.mean(rel)),
        "max_rel": float(np.max(rel)),
    }


def _run_native_backend(
    backend: SixSBackend,
    *,
    geometry: GeometryAngles,
    atmo_state: AtmosphericState,
    bands: list[SensorBand],
) -> tuple[dict[str, RTCoefficients], float]:
    backend.compute_coefficients_multi(geometry, atmo_state, bands)
    start = time.perf_counter()
    outputs = backend.compute_coefficients_multi(geometry, atmo_state, bands)
    elapsed = time.perf_counter() - start
    return {band.name: coeffs for band, coeffs in zip(bands, outputs, strict=True)}, elapsed


def _run_remote_lut_backend(
    backend: ZarrLUTBackend,
    *,
    geometry: GeometryAngles,
    atmo_state: AtmosphericState,
    bands: list[SensorBand],
) -> tuple[dict[str, RTCoefficients], float]:
    for band in bands:
        backend.compute_coefficients(geometry, atmo_state, band)
    start = time.perf_counter()
    outputs = {band.name: backend.compute_coefficients(geometry, atmo_state, band) for band in bands}
    elapsed = time.perf_counter() - start
    return outputs, elapsed


def _scene_lut_counts(backend: SixSBackend, geometry: GeometryAngles, atmo_state: AtmosphericState) -> dict[str, int]:
    prepared = backend._runner._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
    plan = _build_scene_lut_plan(
        prepared.case_arrays,
        max_nodes_per_axis=int(backend._runner._config.scene_lut_max_nodes_per_axis),
        max_cases=int(backend._runner._config.scene_lut_max_cases),
    )
    return {
        "direct_case_count": int(plan.direct_case_count),
        "scene_lut_case_count": int(plan.lut_case_count),
    }


def _scene_summary(geometry: GeometryAngles, atmo_state: AtmosphericState) -> dict[str, dict[str, float]]:
    def _range(values: xr.DataArray) -> dict[str, float]:
        array = np.asarray(values.values, dtype=np.float64)
        return {
            "min": float(np.nanmin(array)),
            "mean": float(np.nanmean(array)),
            "max": float(np.nanmax(array)),
        }

    return {
        "sza_deg": _range(geometry.sza * (180.0 / np.pi)),
        "saa_deg": _range(geometry.saa * (180.0 / np.pi)),
        "vza_deg": _range(geometry.vza * (180.0 / np.pi)),
        "vaa_deg": _range(geometry.vaa * (180.0 / np.pi)),
        "aot550": _range(atmo_state.aot),
        "tcwv_cm": _range(atmo_state.tcwv),
        "tco3_atmcm": _range(atmo_state.tco3),
        "elevation_km": _range(atmo_state.elevation),
    }


def main() -> int:
    args = _parse_args()
    source_dir = args.source_dir.expanduser().resolve()
    build_dir = args.build_dir.expanduser().resolve()
    module_path = args.module_path.expanduser().resolve() if args.module_path is not None else None
    bands = _select_bands(str(args.bands))
    shape = (int(args.size), int(args.size))
    geometry, atmo_state = _build_scene(
        shape,
        str(args.scenario),
        tcwv_cm=float(args.tcwv),
        tco3_atmcm=float(args.tco3),
    )
    observation_time = datetime(2025, 7, 12, 10, 30)

    built_module = module_path
    if built_module is None:
        built_module = build_native_sixs_module(
            _sixs_config(
                source_dir=source_dir,
                build_dir=build_dir / "native_build",
                module_path=None,
                compiler=str(args.compiler),
                native_threads=int(args.native_threads),
                mode="direct",
                atmospheric_profile=str(args.atmospheric_profile),
                atmospheric_columns_mode=str(args.atmospheric_columns_mode),
            )
        )

    direct_backend = SixSBackend(
        sixs_config=_sixs_config(
            source_dir=source_dir,
            build_dir=build_dir / "direct",
            module_path=built_module,
            compiler=str(args.compiler),
            native_threads=int(args.native_threads),
            mode="direct",
            atmospheric_profile=str(args.atmospheric_profile),
            atmospheric_columns_mode=str(args.atmospheric_columns_mode),
        )
    )
    scene_lut_backend = SixSBackend(
        sixs_config=_sixs_config(
            source_dir=source_dir,
            build_dir=build_dir / "scene_lut",
            module_path=built_module,
            compiler=str(args.compiler),
            native_threads=int(args.native_threads),
            mode="scene_lut",
            atmospheric_profile=str(args.atmospheric_profile),
            atmospheric_columns_mode=str(args.atmospheric_columns_mode),
        )
    )
    remote_backend = ZarrLUTBackend(str(args.lut_path), interpolation_method="linear", storage_options={"timeout": 20.0})

    direct_backend.set_observation_time(observation_time)
    scene_lut_backend.set_observation_time(observation_time)

    direct_outputs, direct_seconds = _run_native_backend(
        direct_backend,
        geometry=geometry,
        atmo_state=atmo_state,
        bands=bands,
    )
    scene_lut_outputs, scene_lut_seconds = _run_native_backend(
        scene_lut_backend,
        geometry=geometry,
        atmo_state=atmo_state,
        bands=bands,
    )
    remote_outputs, remote_seconds = _run_remote_lut_backend(
        remote_backend,
        geometry=geometry,
        atmo_state=atmo_state,
        bands=bands,
    )
    plan_counts = _scene_lut_counts(scene_lut_backend, geometry, atmo_state)

    report: dict[str, Any] = {
        "experiment": "compare_6s_route_coefficients",
        "source_dir": str(source_dir),
        "module_path": str(built_module),
        "remote_lut_path": str(args.lut_path),
        "grid_shape": list(shape),
        "scenario": str(args.scenario),
        "native_threads": int(args.native_threads),
        "matching_assumptions": {
            "native_atmospheric_profile": str(args.atmospheric_profile),
            "native_atmospheric_columns_mode": str(args.atmospheric_columns_mode),
            "native_aerosol_profile": "continental",
            "surface_mode": "homogeneous_lambertian",
            "reference_reflectance": 0.1,
            "scene_columns": {
                "tcwv_cm": float(args.tcwv),
                "tco3_atmcm": float(args.tco3),
            },
            "note": (
                "The remote ZIP/Zarr LUT is queried at the same tcwv/tco3 scene inputs as the "
                "native 6S routes. The local direct and scene-LUT 6S runs use the selected "
                "atmospheric profile shape together with atmospheric_columns_mode to determine "
                "whether those scene columns override the profile defaults."
            ),
        },
        "scene_summary": _scene_summary(geometry, atmo_state),
        "route_timings_seconds": {
            "direct_6s": direct_seconds,
            "scene_lut_6s": scene_lut_seconds,
            "remote_zip_zarr_lut": remote_seconds,
        },
        "scene_lut_plan": plan_counts,
        "bands": {},
    }

    for band in bands:
        direct_coeffs = direct_outputs[band.name]
        scene_coeffs = scene_lut_outputs[band.name]
        remote_coeffs = remote_outputs[band.name]
        report["bands"][band.name] = {
            "band": {
                "center_wavelength_nm": float(band.center_wavelength),
                "bandwidth_nm": float(band.bandwidth),
            },
            "route_summaries": {
                "direct_6s": _coeff_summary(direct_coeffs),
                "scene_lut_6s": _coeff_summary(scene_coeffs),
                "remote_zip_zarr_lut": _coeff_summary(remote_coeffs),
            },
            "comparisons_vs_direct_6s": {
                "scene_lut_6s": {
                    name: _comparison_stats(direct_coeffs.get_output(name), scene_coeffs.get_output(name))
                    for name in ("xap", "xbp", "xcp")
                },
                "remote_zip_zarr_lut": {
                    name: _comparison_stats(direct_coeffs.get_output(name), remote_coeffs.get_output(name))
                    for name in ("xap", "xbp", "xcp")
                },
            },
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(args.output)

    direct_backend._runner.close()
    scene_lut_backend._runner.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Benchmark native direct 6S against the scene-generated LUT route."""

from __future__ import annotations

import argparse
import json
import statistics
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_native import _build_scene_lut_plan, _SixSExtensionModule
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles
from siac.sixs_upstream_parity import SixSParityCase, default_parity_cases

if TYPE_CHECKING:
    from collections.abc import Callable


@dataclass(frozen=True)
class BenchmarkResult:
    case_name: str
    scenario: str
    grid_shape: tuple[int, int]
    route: str
    native_threads: int
    repeats: int
    median_seconds: float
    mean_seconds: float
    min_seconds: float
    max_seconds: float
    direct_case_count: int
    scene_lut_case_count: int
    outputs: tuple[str, ...]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("tmp/6s_upstream"), help="Path to unpacked 6S sources.")
    parser.add_argument("--build-dir", type=Path, default=Path("tmp/6s_benchmark_build"), help="Build root for compiled modules.")
    parser.add_argument("--module-path", type=Path, default=None, help="Optional prebuilt native module path.")
    parser.add_argument("--compiler", default="gfortran", help="Fortran compiler executable or path.")
    parser.add_argument("--native-threads", type=int, default=4, help="Native threads for OpenMP/direct runs.")
    parser.add_argument("--repeats", type=int, default=3, help="Measured repeats per route after one warmup.")
    parser.add_argument(
        "--band-count",
        type=int,
        default=4,
        help="How many benchmark bands to include from the predefined set.",
    )
    parser.add_argument(
        "--module-isolation",
        choices=("isolate", "no-isolate"),
        default="no-isolate",
        help="How benchmark runners should load the compiled Python extension.",
    )
    parser.add_argument(
        "--sizes",
        default="32,64,96",
        help="Comma-separated square atmospheric-grid sizes to benchmark.",
    )
    parser.add_argument(
        "--cases",
        default="lambertian_user_water_ozone_continental,rahman_brdf_biomass_burning",
        help="Comma-separated parity case names to benchmark.",
    )
    parser.add_argument(
        "--scenarios",
        default="smooth_atmo,mixed_geometry",
        help="Comma-separated synthetic scene scenarios to benchmark.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tmp/6s_route_benchmark/report.json"),
        help="Where to write the JSON benchmark report.",
    )
    return parser.parse_args()


def _benchmark_bands() -> list[SensorBand]:
    return [
        SensorBand(name="B02", center_wavelength=490.0, bandwidth=20.0, resolution=10.0, band_index=1),
        SensorBand(name="B03", center_wavelength=560.0, bandwidth=35.0, resolution=10.0, band_index=2),
        SensorBand(name="B04", center_wavelength=665.0, bandwidth=30.0, resolution=10.0, band_index=3),
        SensorBand(name="B08", center_wavelength=842.0, bandwidth=115.0, resolution=10.0, band_index=7),
    ]


def _linspace_grid(base: float, delta_x: float, delta_y: float, shape: tuple[int, int]) -> np.ndarray:
    y = np.linspace(0.0, 1.0, shape[0], dtype=np.float32)
    x = np.linspace(0.0, 1.0, shape[1], dtype=np.float32)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return np.asarray(base + delta_x * xx + delta_y * yy, dtype=np.float32)


def _geometry_constant(case: SixSParityCase, shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles.from_degrees(
        xr.DataArray(np.full(shape, case.sza_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.saa_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.vza_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.vaa_deg, dtype=np.float32), dims=("y", "x")),
    )


def _scene_smooth_atmo(case: SixSParityCase, shape: tuple[int, int]) -> tuple[GeometryAngles, AtmosphericState]:
    geometry = _geometry_constant(case, shape)
    aot = _linspace_grid(case.aot550, 0.15, 0.05, shape)
    tcwv = _linspace_grid(case.tcwv_cm, 0.5, 0.3, shape)
    tco3 = np.full(shape, case.tco3_atmcm, dtype=np.float32)
    elevation = np.full(shape, case.elevation_km, dtype=np.float32)
    return geometry, AtmosphericState(
        aot=xr.DataArray(aot, dims=("y", "x")),
        tcwv=xr.DataArray(tcwv, dims=("y", "x")),
        tco3=xr.DataArray(tco3, dims=("y", "x")),
        aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(elevation, dims=("y", "x")),
    )


def _scene_mixed_geometry(case: SixSParityCase, shape: tuple[int, int]) -> tuple[GeometryAngles, AtmosphericState]:
    geometry = GeometryAngles.from_degrees(
        xr.DataArray(_linspace_grid(case.sza_deg, 4.0, 2.0, shape), dims=("y", "x")),
        xr.DataArray(_linspace_grid(case.saa_deg, 20.0, 15.0, shape), dims=("y", "x")),
        xr.DataArray(_linspace_grid(case.vza_deg, 3.0, 2.0, shape), dims=("y", "x")),
        xr.DataArray(_linspace_grid(case.vaa_deg, 18.0, 12.0, shape), dims=("y", "x")),
    )
    aot = _linspace_grid(case.aot550, 0.15, 0.10, shape)
    tcwv = _linspace_grid(case.tcwv_cm, 0.6, 0.4, shape)
    tco3 = np.full(shape, case.tco3_atmcm, dtype=np.float32)
    elevation = _linspace_grid(case.elevation_km, 0.2, 0.1, shape)
    return geometry, AtmosphericState(
        aot=xr.DataArray(aot, dims=("y", "x")),
        tcwv=xr.DataArray(tcwv, dims=("y", "x")),
        tco3=xr.DataArray(tco3, dims=("y", "x")),
        aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(elevation, dims=("y", "x")),
    )


_SCENARIOS: dict[str, Callable[[SixSParityCase, tuple[int, int]], tuple[GeometryAngles, AtmosphericState]]] = {
    "smooth_atmo": _scene_smooth_atmo,
    "mixed_geometry": _scene_mixed_geometry,
}


def _prepare_backend(
    case: SixSParityCase,
    *,
    route: str,
    source_dir: Path,
    build_dir: Path,
    module_path: Path | None,
    compiler: str,
    native_threads: int,
    outputs: tuple[str, ...],
    module_isolation: str,
) -> SixSBackend:
    config = case.config.model_copy(
        update={
            "source_dir": source_dir,
            "build_dir": build_dir / route,
            "module_path": module_path,
            "auto_build": module_path is None,
            "compiler": compiler,
            "build_profile": "release",
            "mode": route,
            "parallel_backend": "openmp",
            "native_threads": native_threads,
            "output_variables": outputs,
            "scene_lut_max_nodes_per_axis": 4,
            "scene_lut_max_cases": 4096,
            "scene_lut_min_pixels": 512,
            "scene_lut_required_speedup": 1.5,
        }
    )
    backend = SixSBackend(sixs_config=config)
    if module_path is not None and module_isolation == "no-isolate":
        backend._runner._openmp_session = _SixSExtensionModule(module_path, isolate=False)
    backend.set_observation_time(datetime(2025, case.config.month, case.config.day, 10, 30))
    return backend


def _timed_run(
    backend: SixSBackend,
    geometry: GeometryAngles,
    atmo_state: AtmosphericState,
    bands: list[SensorBand],
    repeats: int,
) -> list[float]:
    backend.compute_coefficients_multi(geometry, atmo_state, bands)
    samples: list[float] = []
    for _ in range(repeats):
        start = time.perf_counter()
        backend.compute_coefficients_multi(geometry, atmo_state, bands)
        samples.append(time.perf_counter() - start)
    return samples


def _plan_counts(backend: SixSBackend, geometry: GeometryAngles, atmo_state: AtmosphericState) -> tuple[int, int]:
    runner = backend._runner
    prepared = runner._prepare_scene_inputs(geometry=geometry, atmo_state=atmo_state)
    plan = _build_scene_lut_plan(
        prepared.case_arrays,
        max_nodes_per_axis=int(runner._config.scene_lut_max_nodes_per_axis),
        max_cases=int(runner._config.scene_lut_max_cases),
    )
    return plan.direct_case_count, plan.lut_case_count


def main() -> int:
    args = _parse_args()
    source_dir = args.source_dir.expanduser().resolve()
    build_dir = args.build_dir.expanduser().resolve()
    module_path = args.module_path.expanduser().resolve() if args.module_path is not None else None
    sizes = [int(item.strip()) for item in str(args.sizes).split(",") if item.strip()]
    case_names = {item.strip() for item in str(args.cases).split(",") if item.strip()}
    scenario_names = [item.strip() for item in str(args.scenarios).split(",") if item.strip()]
    outputs = ("xap", "xbp", "xcp")

    cases = [case for case in default_parity_cases() if case.name in case_names]
    if not cases:
        raise SystemExit("No benchmark cases matched --cases.")

    results: list[BenchmarkResult] = []
    bands = _benchmark_bands()[: max(1, int(args.band_count))]
    for case in cases:
        for scenario_name in scenario_names:
            scene_builder = _SCENARIOS[scenario_name]
            for size in sizes:
                shape = (size, size)
                geometry, atmo_state = scene_builder(case, shape)
                direct_backend = _prepare_backend(
                    case,
                    route="direct",
                    source_dir=source_dir,
                    build_dir=build_dir,
                    module_path=module_path,
                    compiler=args.compiler,
                    native_threads=args.native_threads,
                    outputs=outputs,
                    module_isolation=args.module_isolation,
                )
                scene_lut_backend = _prepare_backend(
                    case,
                    route="scene_lut",
                    source_dir=source_dir,
                    build_dir=build_dir,
                    module_path=module_path,
                    compiler=args.compiler,
                    native_threads=args.native_threads,
                    outputs=outputs,
                    module_isolation=args.module_isolation,
                )
                direct_case_count, lut_case_count = _plan_counts(scene_lut_backend, geometry, atmo_state)
                for route, backend in (("direct", direct_backend), ("scene_lut", scene_lut_backend)):
                    samples = _timed_run(backend, geometry, atmo_state, bands, args.repeats)
                    results.append(
                        BenchmarkResult(
                            case_name=case.name,
                            scenario=scenario_name,
                            grid_shape=shape,
                            route=route,
                            native_threads=args.native_threads,
                            repeats=args.repeats,
                            median_seconds=statistics.median(samples),
                            mean_seconds=float(statistics.fmean(samples)),
                            min_seconds=min(samples),
                            max_seconds=max(samples),
                            direct_case_count=direct_case_count,
                            scene_lut_case_count=lut_case_count,
                            outputs=outputs,
                        )
                    )
                direct_backend._runner.close()
                scene_lut_backend._runner.close()

    report = {
        "native_threads": args.native_threads,
        "repeats": args.repeats,
        "module_isolation": args.module_isolation,
        "cases": sorted(case_names),
        "scenarios": scenario_names,
        "sizes": sizes,
        "results": [asdict(item) for item in results],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

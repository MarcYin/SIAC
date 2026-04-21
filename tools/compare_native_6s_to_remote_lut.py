"""Compare native direct 6S against the remote libRadtran ZIP/Zarr LUT."""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import build_native_sixs_module
from siac.algorithms.rt.lut import DEFAULT_LUT_URL, ZarrLUTBackend
from siac.config import SixSAlgorithmConfig
from siac.config.schema import (
    SixSAtmosphericCorrectionConfig,
    SixSSpectralReflectanceConfig,
    SixSSurfaceConfig,
)
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients


@dataclass(frozen=True)
class PointCase:
    name: str
    sza_deg: float
    saa_deg: float
    vza_deg: float
    vaa_deg: float
    aot550: float
    tcwv_cm: float
    tco3_atmcm: float
    elevation_km: float


@dataclass(frozen=True)
class SceneCase:
    name: str
    scenario: str
    shape: tuple[int, int]
    base: PointCase


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("tmp/6s_upstream"), help="Path to unpacked 6SV2.1 sources.")
    parser.add_argument("--build-dir", type=Path, default=Path("tmp/rt_model_family_report"), help="Build root for the native 6S module.")
    parser.add_argument("--module-path", type=Path, default=None, help="Optional prebuilt native 6S Python extension.")
    parser.add_argument("--compiler", default="gfortran", help="Fortran compiler executable or path.")
    parser.add_argument("--native-threads", type=int, default=4, help="OpenMP thread count for native 6S.")
    parser.add_argument("--lut-path", default=DEFAULT_LUT_URL, help="Path or URL to the remote ZIP/Zarr LUT.")
    parser.add_argument(
        "--bands",
        default="B02,B03,B04,B08",
        help="Comma-separated band names from the built-in benchmark catalog.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tmp/rt_model_family_report/report.json"),
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


def _point_cases() -> list[PointCase]:
    return [
        PointCase(
            name="reference_lowland",
            sza_deg=30.0,
            saa_deg=150.0,
            vza_deg=5.0,
            vaa_deg=110.0,
            aot550=0.15,
            tcwv_cm=2.0,
            tco3_atmcm=0.30,
            elevation_km=0.1,
        ),
        PointCase(
            name="humid_highland",
            sza_deg=30.0,
            saa_deg=150.0,
            vza_deg=5.0,
            vaa_deg=110.0,
            aot550=0.30,
            tcwv_cm=4.0,
            tco3_atmcm=0.33,
            elevation_km=1.0,
        ),
        PointCase(
            name="oblique_dry_geometry",
            sza_deg=45.0,
            saa_deg=170.0,
            vza_deg=18.0,
            vaa_deg=120.0,
            aot550=0.12,
            tcwv_cm=1.2,
            tco3_atmcm=0.28,
            elevation_km=0.3,
        ),
    ]


def _scene_cases() -> list[SceneCase]:
    points = {case.name: case for case in _point_cases()}
    return [
        SceneCase(
            name="reference_lowland_smooth_2x2",
            scenario="smooth_geometry",
            shape=(2, 2),
            base=points["reference_lowland"],
        ),
        SceneCase(
            name="oblique_dry_mixed_2x2",
            scenario="mixed_geometry",
            shape=(2, 2),
            base=points["oblique_dry_geometry"],
        ),
    ]


def _build_point_scene(case: PointCase) -> tuple[GeometryAngles, AtmosphericState]:
    shape = (1, 1)
    geometry = GeometryAngles.from_degrees(
        xr.DataArray(np.full(shape, case.sza_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.saa_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.vza_deg, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, case.vaa_deg, dtype=np.float32), dims=("y", "x")),
    )
    atmo_state = AtmosphericState(
        aot=xr.DataArray(np.full(shape, case.aot550, dtype=np.float32), dims=("y", "x")),
        tcwv=xr.DataArray(np.full(shape, case.tcwv_cm, dtype=np.float32), dims=("y", "x")),
        tco3=xr.DataArray(np.full(shape, case.tco3_atmcm, dtype=np.float32), dims=("y", "x")),
        aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(np.full(shape, case.elevation_km, dtype=np.float32), dims=("y", "x")),
    )
    return geometry, atmo_state


def _build_scene_case(case: SceneCase) -> tuple[GeometryAngles, AtmosphericState]:
    shape = case.shape
    base = case.base
    if case.scenario == "smooth_geometry":
        sza = _linspace_grid(base.sza_deg, 1.0, 0.5, shape)
        saa = _linspace_grid(base.saa_deg, 3.0, 2.0, shape)
        vza = _linspace_grid(base.vza_deg, 0.8, 0.5, shape)
        vaa = _linspace_grid(base.vaa_deg, 2.0, 1.5, shape)
        aot = _linspace_grid(base.aot550, 0.04, 0.02, shape)
        tcwv = _linspace_grid(base.tcwv_cm, 0.25, 0.15, shape)
        tco3 = _linspace_grid(base.tco3_atmcm, 0.01, 0.005, shape)
        elevation = _linspace_grid(base.elevation_km, 0.05, 0.02, shape)
    else:
        sza = _linspace_grid(base.sza_deg, 5.0, 3.0, shape)
        saa = _linspace_grid(base.saa_deg, 12.0, 8.0, shape)
        vza = _linspace_grid(base.vza_deg, 4.0, 3.0, shape)
        vaa = _linspace_grid(base.vaa_deg, 15.0, 10.0, shape)
        aot = _linspace_grid(base.aot550, 0.10, 0.06, shape)
        tcwv = _linspace_grid(base.tcwv_cm, 0.35, 0.20, shape)
        tco3 = _linspace_grid(base.tco3_atmcm, 0.015, 0.01, shape)
        elevation = _linspace_grid(base.elevation_km, 0.10, 0.05, shape)
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
        tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(elevation, dims=("y", "x")),
    )
    return geometry, atmo_state


def _linspace_grid(base: float, delta_x: float, delta_y: float, shape: tuple[int, int]) -> np.ndarray:
    y = np.linspace(0.0, 1.0, shape[0], dtype=np.float32)
    x = np.linspace(0.0, 1.0, shape[1], dtype=np.float32)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return np.asarray(base + delta_x * xx + delta_y * yy, dtype=np.float32)


def _sixs_config(
    *,
    source_dir: Path,
    build_dir: Path,
    module_path: Path | None,
    compiler: str,
    native_threads: int,
) -> SixSAlgorithmConfig:
    return SixSAlgorithmConfig(
        source_dir=source_dir,
        build_dir=build_dir,
        module_path=module_path,
        auto_build=module_path is None,
        compiler=compiler,
        build_profile="release",
        mode="direct",
        parallel_backend="openmp",
        native_threads=native_threads,
        month=7,
        day=12,
        atmospheric_profile="us_standard_62",
        atmospheric_columns_mode="input_columns",
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


def _coeffs_to_scalars(coeffs: RTCoefficients) -> dict[str, float]:
    return {
        name: float(coeffs.get_output(name).values[0, 0])
        for name in ("xap", "xbp", "xcp")
    }


def _scalar_differences(reference: dict[str, float], candidate: dict[str, float]) -> dict[str, dict[str, float]]:
    metrics: dict[str, dict[str, float]] = {}
    for name, ref_value in reference.items():
        cand_value = candidate[name]
        abs_diff = abs(cand_value - ref_value)
        denom = max(abs(ref_value), 1.0e-12)
        metrics[name] = {
            "direct_value": ref_value,
            "candidate_value": cand_value,
            "abs_diff": abs_diff,
            "rel_diff": abs_diff / denom,
        }
    return metrics


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
        "max_abs": float(np.max(abs_diff)),
        "rmse": float(np.sqrt(np.mean(diff * diff))),
        "mean_rel": float(np.mean(rel)),
        "max_rel": float(np.max(rel)),
    }


def _aggregate_relative_differences(entries: list[dict[str, Any]]) -> dict[str, dict[str, dict[str, float]]]:
    grouped: dict[str, dict[str, list[float]]] = {}
    for entry in entries:
        for band_name, band_report in entry["bands"].items():
            bucket = grouped.setdefault(band_name, {"xap": [], "xbp": [], "xcp": []})
            metrics = band_report["remote_vs_direct"]
            for coeff_name in ("xap", "xbp", "xcp"):
                if isinstance(metrics[coeff_name], dict) and "rel_diff" in metrics[coeff_name]:
                    bucket[coeff_name].append(float(metrics[coeff_name]["rel_diff"]))
                else:
                    bucket[coeff_name].append(float(metrics[coeff_name]["mean_rel"]))

    summary: dict[str, dict[str, dict[str, float]]] = {}
    for band_name, coeff_map in grouped.items():
        summary[band_name] = {}
        for coeff_name, values in coeff_map.items():
            summary[band_name][coeff_name] = {
                "min_rel": float(min(values)),
                "mean_rel": float(sum(values) / len(values)),
                "max_rel": float(max(values)),
            }
    return summary


def main() -> int:
    args = _parse_args()
    source_dir = args.source_dir.expanduser().resolve()
    build_dir = args.build_dir.expanduser().resolve()
    module_path = args.module_path.expanduser().resolve() if args.module_path is not None else None
    bands = _select_bands(str(args.bands))
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
            )
        )

    sixs_backend = SixSBackend(
        sixs_config=_sixs_config(
            source_dir=source_dir,
            build_dir=build_dir / "direct",
            module_path=built_module,
            compiler=str(args.compiler),
            native_threads=int(args.native_threads),
        )
    )
    sixs_backend.set_observation_time(observation_time)
    remote_backend = ZarrLUTBackend(str(args.lut_path), interpolation_method="linear", storage_options={"timeout": 20.0})

    report: dict[str, Any] = {
        "experiment": "compare_native_6s_to_remote_lut",
        "module_path": str(built_module),
        "remote_lut_path": str(args.lut_path),
        "remote_lut_identity": {
            "generator": "libRadtran",
            "aerosol_profile": "continental_average",
            "atmospheric_profile_shape": "us_standard_62",
            "note": (
                "This remote ZIP/Zarr LUT is the SIAC remote spectral LUT derived from libRadtran. "
                "It is not a 6S LUT, so non-zero coefficient differences are expected even when "
                "tcwv, tco3, geometry, altitude, and AOT inputs are matched."
            ),
        },
        "native_6s_identity": {
            "atmospheric_profile": "us_standard_62",
            "atmospheric_columns_mode": "input_columns",
            "aerosol_profile": "continental",
            "surface_mode": "homogeneous_lambertian",
        },
        "bands": [
            {
                "name": band.name,
                "center_wavelength_nm": float(band.center_wavelength),
                "bandwidth_nm": float(band.bandwidth),
            }
            for band in bands
        ],
        "point_cases": [],
        "scene_cases": [],
    }

    for case in _point_cases():
        geometry, atmo_state = _build_point_scene(case)
        case_report: dict[str, Any] = {
            "case": asdict(case),
            "bands": {},
        }
        for band in bands:
            start = time.perf_counter()
            direct_coeffs = sixs_backend.compute_coefficients(geometry, atmo_state, band)
            direct_seconds = time.perf_counter() - start

            start = time.perf_counter()
            remote_coeffs = remote_backend.compute_coefficients(geometry, atmo_state, band)
            remote_seconds = time.perf_counter() - start

            direct_values = _coeffs_to_scalars(direct_coeffs)
            remote_values = _coeffs_to_scalars(remote_coeffs)
            case_report["bands"][band.name] = {
                "timing_seconds": {
                    "native_direct_6s": direct_seconds,
                    "remote_zip_zarr_lut": remote_seconds,
                },
                "direct_6s": direct_values,
                "remote_zip_zarr_lut": remote_values,
                "remote_vs_direct": _scalar_differences(direct_values, remote_values),
            }
        report["point_cases"].append(case_report)

    for case in _scene_cases():
        geometry, atmo_state = _build_scene_case(case)
        case_report = {
            "case": {
                "name": case.name,
                "scenario": case.scenario,
                "shape": list(case.shape),
                "base": asdict(case.base),
            },
            "bands": {},
        }
        for band in bands:
            start = time.perf_counter()
            direct_coeffs = sixs_backend.compute_coefficients(geometry, atmo_state, band)
            direct_seconds = time.perf_counter() - start

            start = time.perf_counter()
            remote_coeffs = remote_backend.compute_coefficients(geometry, atmo_state, band)
            remote_seconds = time.perf_counter() - start

            case_report["bands"][band.name] = {
                "timing_seconds": {
                    "native_direct_6s": direct_seconds,
                    "remote_zip_zarr_lut": remote_seconds,
                },
                "direct_6s_summary": _coeff_summary(direct_coeffs),
                "remote_zip_zarr_lut_summary": _coeff_summary(remote_coeffs),
                "remote_vs_direct": {
                    name: _comparison_stats(direct_coeffs.get_output(name), remote_coeffs.get_output(name))
                    for name in ("xap", "xbp", "xcp")
                },
            }
        report["scene_cases"].append(case_report)

    report["summary"] = {
        "point_case_relative_differences_by_band": _aggregate_relative_differences(report["point_cases"]),
        "scene_case_relative_differences_by_band": _aggregate_relative_differences(report["scene_cases"]),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(args.output)

    sixs_backend._runner.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

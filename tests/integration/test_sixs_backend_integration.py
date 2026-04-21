from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.adapters.rt import build_rt_model
from siac.algorithms.rt.direct.sixs_build import build_native_sixs_module
from siac.config import RTSetupConfig, SixSAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients
from siac.sixs_upstream_parity import default_parity_cases, run_upstream_parity_suite


def _require_native_6s_stack() -> Path:
    source_dir = Path("tmp/6s_upstream").resolve()
    if not (source_dir / "main.f").exists():
        pytest.skip("Upstream 6S source tree is not available under tmp/6s_upstream.")
    missing = [tool for tool in ("gfortran", "meson", "ninja") if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"Native 6S integration tests require toolchain components on PATH: {', '.join(missing)}.")
    return source_dir


def _adapter_config(sixs_config: SixSAlgorithmConfig, rt_setup: RTSetupConfig) -> SimpleNamespace:
    return SimpleNamespace(
        sensor="auto",
        paths=SimpleNamespace(emulator_dir=None, lut_path=None),
        rt_model=SimpleNamespace(
            backend="sixs",
            setup=rt_setup,
            sixs=sixs_config,
            emulator_dir=None,
            lut_path=None,
            fallback_to_lut=False,
            lut_storage_options={},
            lut_interpolation="linear",
        ),
    )


def _scene_geometry(case_name: str, shape: tuple[int, int] = (2, 2)) -> GeometryAngles:
    cases = {case.name: case for case in default_parity_cases()}
    case = cases[case_name]
    sza = np.full(shape, case.sza_deg, dtype=np.float32)
    saa = np.full(shape, case.saa_deg, dtype=np.float32)
    vza = np.full(shape, case.vza_deg, dtype=np.float32)
    vaa = np.full(shape, case.vaa_deg, dtype=np.float32)
    return GeometryAngles.from_degrees(
        xr.DataArray(sza, dims=("y", "x")),
        xr.DataArray(saa, dims=("y", "x")),
        xr.DataArray(vza, dims=("y", "x")),
        xr.DataArray(vaa, dims=("y", "x")),
    )


def _scene_atmosphere(case_name: str, shape: tuple[int, int] = (2, 2)) -> AtmosphericState:
    cases = {case.name: case for case in default_parity_cases()}
    case = cases[case_name]
    aot = np.full(shape, case.aot550, dtype=np.float32)
    tcwv = np.full(shape, case.tcwv_cm, dtype=np.float32)
    tco3 = np.full(shape, case.tco3_atmcm, dtype=np.float32)
    elevation = np.full(shape, case.elevation_km, dtype=np.float32)
    aot[0, 1] = np.nan
    return AtmosphericState(
        aot=xr.DataArray(aot, dims=("y", "x")),
        tcwv=xr.DataArray(tcwv, dims=("y", "x")),
        tco3=xr.DataArray(tco3, dims=("y", "x")),
        aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(elevation, dims=("y", "x")),
    )


def _integration_bands() -> list[SensorBand]:
    return [
        SensorBand(
            name="B04",
            center_wavelength=665.0,
            bandwidth=30.0,
            resolution=10.0,
            band_index=3,
            rsrf_wavelengths_nm=np.arange(650.0, 680.0 + 2.5, 2.5, dtype=np.float64),
            rsrf_response=np.ones(13, dtype=np.float64),
        ),
        SensorBand(
            name="B03",
            center_wavelength=560.0,
            bandwidth=35.0,
            resolution=10.0,
            band_index=2,
            rsrf_wavelengths_nm=np.arange(545.0, 575.0 + 2.5, 2.5, dtype=np.float64),
            rsrf_response=np.array([0.0, 0.2, 0.6, 1.0, 0.8, 0.4, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        ),
    ]


@pytest.fixture(scope="module")
def sixs_source_dir() -> Path:
    return _require_native_6s_stack()


@pytest.fixture(scope="module")
def native_module_path(tmp_path_factory: pytest.TempPathFactory, sixs_source_dir: Path) -> Path:
    build_root = tmp_path_factory.mktemp("sixs_backend_integration_build")
    config = SixSAlgorithmConfig(
        source_dir=sixs_source_dir,
        build_dir=build_root,
        compiler="gfortran",
    )
    return build_native_sixs_module(config)


@pytest.mark.integration
class TestSixSBackendIntegration:
    def test_native_backend_matches_original_case_suite(
        self,
        tmp_path: Path,
        sixs_source_dir: Path,
    ) -> None:
        report = run_upstream_parity_suite(
            source_dir=sixs_source_dir,
            upstream_build_dir=tmp_path / "upstream_exec",
            native_build_dir=tmp_path / "native_build",
            output_dir=tmp_path / "cases",
            compiler="gfortran",
        )

        assert report["all_cases_matched"] is True
        assert len(report["cases"]) == len(default_parity_cases())
        for case_report in report["cases"]:
            assert case_report["matched"] is True
            assert case_report["matched_variable_count"] == case_report["compared_variable_count"]
            assert case_report["compared_variable_count"] >= 100
            assert len(case_report["mismatches"]) == 0

    def test_backend_scene_multiband_outputs_match_single_band_contract(
        self,
        native_module_path: Path,
        sixs_source_dir: Path,
    ) -> None:
        case = next(case for case in default_parity_cases() if case.name == "ocean_brdf_maritime")
        requested_outputs = (
            "xap",
            "xbp",
            "xcp",
            "tgasm",
            "sutott",
            "sast",
            "rocave",
            "rfoamave",
            "rwatave",
            "rglitave",
            "rpfet",
            "plumet",
            "xpol",
        )
        sixs_config = case.sixs_config.model_copy(
            update={
                "source_dir": sixs_source_dir,
                "module_path": native_module_path,
                "auto_build": False,
                "native_threads": 1,
                "output_variables": requested_outputs,
            }
        )
        backend = build_rt_model(_adapter_config(sixs_config, case.rt_setup))
        geometry = _scene_geometry(case.name)
        atmosphere = _scene_atmosphere(case.name)
        bands = _integration_bands()

        backend.set_observation_time(datetime(2025, case.sixs_config.month, case.sixs_config.day, 10, 30))
        assert backend.backend_name == "sixs"
        assert backend.preload_scene_subset(geometry, atmosphere, bands) is None

        multi = backend.compute_coefficients_multi(geometry, atmosphere, bands)
        singles = [backend.compute_coefficients(geometry, atmosphere, band) for band in bands]

        assert len(multi) == len(bands) == len(singles)
        valid_mask = np.isfinite(atmosphere.aot.values)
        invalid_mask = ~valid_mask

        for combined, single in zip(multi, singles, strict=True):
            assert isinstance(combined, RTCoefficients)
            assert combined.output_names == requested_outputs
            assert single.output_names == requested_outputs

            for output_name in requested_outputs:
                combined_arr = combined.get_output(output_name)
                single_arr = single.get_output(output_name)

                assert combined_arr.dims == ("y", "x")
                assert single_arr.dims == ("y", "x")
                assert np.all(np.isnan(combined_arr.values[invalid_mask]))
                assert np.all(np.isnan(single_arr.values[invalid_mask]))
                assert np.all(np.isfinite(combined_arr.values[valid_mask]))
                assert np.all(np.isfinite(single_arr.values[valid_mask]))
                np.testing.assert_allclose(
                    combined_arr.values,
                    single_arr.values,
                    rtol=1.0e-7,
                    atol=1.0e-9,
                    equal_nan=True,
                )

    def test_backend_reuses_live_native_session_across_repeated_calls(
        self,
        native_module_path: Path,
        sixs_source_dir: Path,
    ) -> None:
        case = next(case for case in default_parity_cases() if case.name == "rahman_brdf_biomass_burning")
        sixs_config = case.sixs_config.model_copy(
            update={
                "source_dir": sixs_source_dir,
                "module_path": native_module_path,
                "auto_build": False,
                "native_threads": 1,
                "output_variables": ("xap", "xbp", "xcp", "tgasm", "sutott", "sast"),
            }
        )
        backend = build_rt_model(_adapter_config(sixs_config, case.rt_setup))
        geometry = _scene_geometry(case.name)
        atmosphere = _scene_atmosphere(case.name)
        band = _integration_bands()[0]

        first = backend.compute_coefficients(geometry, atmosphere, band)
        second = backend.compute_coefficients(geometry, atmosphere, band)

        for name in first.output_names:
            np.testing.assert_allclose(
                first.get_output(name).values,
                second.get_output(name).values,
                rtol=1.0e-7,
                atol=1.0e-9,
                equal_nan=True,
            )

    def test_scene_lut_matches_direct_for_grid_aligned_scene(
        self,
        native_module_path: Path,
        sixs_source_dir: Path,
    ) -> None:
        case = next(case for case in default_parity_cases() if case.name == "ocean_brdf_maritime")
        geometry = GeometryAngles.from_degrees(
            xr.DataArray([[case.sza_deg, case.sza_deg], [case.sza_deg, case.sza_deg]], dims=("y", "x")),
            xr.DataArray([[case.saa_deg, case.saa_deg], [case.saa_deg, case.saa_deg]], dims=("y", "x")),
            xr.DataArray([[case.vza_deg, case.vza_deg], [case.vza_deg, case.vza_deg]], dims=("y", "x")),
            xr.DataArray([[case.vaa_deg, case.vaa_deg], [case.vaa_deg, case.vaa_deg]], dims=("y", "x")),
        )
        atmosphere = AtmosphericState(
            aot=xr.DataArray([[case.aot550, case.aot550 + 0.03], [case.aot550, case.aot550 + 0.03]], dims=("y", "x")),
            tcwv=xr.DataArray([[case.tcwv_cm, case.tcwv_cm + 0.4], [case.tcwv_cm, case.tcwv_cm + 0.4]], dims=("y", "x")),
            tco3=xr.DataArray([[case.tco3_atmcm, case.tco3_atmcm], [case.tco3_atmcm, case.tco3_atmcm]], dims=("y", "x")),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=("y", "x")),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.05, dtype=np.float32), dims=("y", "x")),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=("y", "x")),
            elevation=xr.DataArray([[case.elevation_km, case.elevation_km], [case.elevation_km, case.elevation_km]], dims=("y", "x")),
        )
        requested_outputs = ("xap", "xbp", "xcp", "tgasm", "sutott", "sast", "rocave")
        direct_backend = build_rt_model(
            _adapter_config(
                case.sixs_config.model_copy(
                    update={
                        "source_dir": sixs_source_dir,
                        "module_path": native_module_path,
                        "auto_build": False,
                        "mode": "direct",
                        "native_threads": 1,
                        "output_variables": requested_outputs,
                    }
                ),
                case.rt_setup,
            )
        )
        scene_lut_backend = build_rt_model(
            _adapter_config(
                case.sixs_config.model_copy(
                    update={
                        "source_dir": sixs_source_dir,
                        "module_path": native_module_path,
                        "auto_build": False,
                        "mode": "scene_lut",
                        "native_threads": 1,
                        "scene_lut_max_nodes_per_axis": 2,
                        "scene_lut_max_cases": 16,
                        "output_variables": requested_outputs,
                    }
                ),
                case.rt_setup,
            )
        )
        band = _integration_bands()[0]

        direct = direct_backend.compute_coefficients(geometry, atmosphere, band)
        scene_lut = scene_lut_backend.compute_coefficients(geometry, atmosphere, band)

        for name in requested_outputs:
            np.testing.assert_allclose(
                direct.get_output(name).values,
                scene_lut.get_output(name).values,
                rtol=1.0e-6,
                atol=1.0e-8,
                equal_nan=True,
            )

    def test_worker_libraries_matches_openmp_backend(
        self,
        native_module_path: Path,
        sixs_source_dir: Path,
    ) -> None:
        case = next(case for case in default_parity_cases() if case.name == "rahman_brdf_biomass_burning")
        requested_outputs = ("xap", "xbp", "xcp", "tgasm", "sutott", "sast")
        common_update = {
            "source_dir": sixs_source_dir,
            "module_path": native_module_path,
            "auto_build": False,
            "output_variables": requested_outputs,
        }
        openmp_backend = build_rt_model(
            _adapter_config(
                case.sixs_config.model_copy(
                    update={
                        **common_update,
                        "parallel_backend": "openmp",
                        "native_threads": 2,
                    }
                ),
                case.rt_setup,
            )
        )
        worker_backend = build_rt_model(
            _adapter_config(
                case.sixs_config.model_copy(
                    update={
                        **common_update,
                        "parallel_backend": "worker_libraries",
                        "native_threads": 2,
                        "worker_libraries": 2,
                        "chunk_size": 2,
                    }
                ),
                case.rt_setup,
            )
        )
        geometry = _scene_geometry(case.name, shape=(3, 3))
        atmosphere = _scene_atmosphere(case.name, shape=(3, 3))
        band = _integration_bands()[0]

        openmp = openmp_backend.compute_coefficients(geometry, atmosphere, band)
        worker = worker_backend.compute_coefficients(geometry, atmosphere, band)

        for name in requested_outputs:
            np.testing.assert_allclose(
                openmp.get_output(name).values,
                worker.get_output(name).values,
                rtol=1.0e-6,
                atol=1.0e-8,
                equal_nan=True,
            )

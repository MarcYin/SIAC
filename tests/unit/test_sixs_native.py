from __future__ import annotations

import subprocess
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.direct import sixs_build as sixs_build_module
from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import (
    _compile_f2py_extension,
    _distutils_backend_supported,
    _distutils_extra_link_args,
    _EncodedStringIO,
    _find_built_extension,
    _resolve_f2py_backends,
    patch_aeroso_source,
    patch_discom_source,
    patch_kernelpol_source,
    patch_main_source,
    patch_mie_source,
    patch_threadprivate_directives,
    resolve_build_paths,
)
from siac.algorithms.rt.direct.sixs_native import (
    _build_scene_lut_plan,
    _build_spectral_response,
    _NativeBatchResult,
    _resolve_atmospheric_mode,
    _should_use_scene_lut,
)
from siac.config import SixSAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles
from siac.sixs_outputs import SIXS_CORE_OUTPUT_SPECS


def _sample_geometry(shape: tuple[int, int] = (2, 2)) -> GeometryAngles:
    return GeometryAngles.from_degrees(
        xr.DataArray(np.full(shape, 30.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, 150.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, 5.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full(shape, 110.0, dtype=np.float32), dims=("y", "x")),
    )


def _sample_atmo(shape: tuple[int, int] = (2, 2)) -> AtmosphericState:
    ones = np.ones(shape, dtype=np.float32)
    return AtmosphericState(
        aot=xr.DataArray(0.15 * ones, dims=("y", "x")),
        tcwv=xr.DataArray(2.0 * ones, dims=("y", "x")),
        tco3=xr.DataArray(0.3 * ones, dims=("y", "x")),
        aot_unc=xr.DataArray(0.01 * ones, dims=("y", "x")),
        tcwv_unc=xr.DataArray(0.05 * ones, dims=("y", "x")),
        tco3_unc=xr.DataArray(0.01 * ones, dims=("y", "x")),
        elevation=xr.DataArray(0.1 * ones, dims=("y", "x")),
    )


def test_patch_discom_source_applies_boundary_fix() -> None:
    original = "\n".join(
        [
            "        if ((l.lt.20).and.(wldis(l).lt.wlinf).and.",
            "     a     (wldis(l+1).lt.wlinf)) goto 50",
            "        if ((l.gt.1).and.(wldis(l).gt.wlsup).and.",
            "     a      (wldis(l-1).gt.wlsup)) goto 50",
        ]
    )

    patched = patch_discom_source(original)

    assert "if (l.lt.20) then" in patched
    assert "if (l.gt.1) then" in patched


def test_patch_kernelpol_source_initializes_zero_order_terms() -> None:
    original = "\n".join(
        [
            "        psl(2,j)=xdb",
            "        psl(2,-j)=xdb",
            "        rsl(1,j)=0.0D+00",
            "        rsl(1,-j)=0.0D+00",
            "        xdb=3.D+00*(1.D+00-c*c)/2.D+00/sqrt(6.D+00)",
            "\tif (abs(xdb).lt.1.e-30) xdb =0.0D+00",
            "        rsl(2,j)=xdb",
            "        rsl(2,-j)=xdb",
            "        tsl(1,j)=0.0D+00",
            "        tsl(1,-j)=0.0D+00",
            "        tsl(2,j)=0.0D+00",
            "        tsl(2,-j)=0.0D+00",
        ]
    )

    patched = patch_kernelpol_source(original)

    assert "rsl(0,j)=0.0D+00" in patched
    assert "tsl(0,j)=0.0D+00" in patched


def test_patch_aeroso_source_widens_mie_path_buffer() -> None:
    patched = patch_aeroso_source(
        "\n".join(
            [
                "      character FILE*80",
                "      if (iaer.ge.8.and.iaer.le.11) then",
                "        open(10,file=FILE)",
            ]
        )
    )

    assert "character FILE*1024" in patched
    assert "len_trim(FILE).gt.0" in patched


def test_patch_mie_source_uses_dynamic_scratch_storage() -> None:
    original = "\n".join(
        [
            "      do i=1,icp",
            "        np(i)=0.D+00",
            "      enddo",
            "      parameter (nser=10000000)      ",
            "      double precision Ren,Imn,X,Y,Up,XnumRDnY,XnumIDnY",
            "      double precision XdenDnY,coxj,Qsca,Qext,xJonH,XdenGNX",
            "      double precision Xnum1An,Xnum2An,XdenAn,Xden1An,Xden2An,RAnb,IAnb",
            "      double precision Xnum1Bn,Xnum2Bn,XdenBn,Xden1Bn,Xden2Bn,RBnb,IBnb",
            "      double precision xmud,xpond,RS1,RS2,IS1,IS2,co_n,test",
            "      double precision xj(0:nser),xy(-1:nser),Rn(0:nser)",
            "      double precision IDnY(0:nser),RDnX(0:nser),RDnY(0:nser)",
            "      double precision IGnX(0:nser),RGnX(0:nser)",
            "      double precision RAn(0:nser),IAn(0:nser),RBn(0:nser),IBn(0:nser)",
            "      double precision TAUn(0:nser),PIn(0:nser)",
            "      if (mu.ge.nser) then",
            '         write(6,*) " Error, nser is too small, mu is equal to : ",mu',
            "         Stop",
            "         endif",
            "",
            "",
            "      if (mu.le.0) then",
            '         write(6,*) " Error, mu is too small, mu is equal to : ",mu',
            "         Stop",
            "         endif",
            "      enddo",
            "               ",
            "      return",
            "      end",
        ]
    )

    patched = patch_mie_source(original)

    assert "allocatable ::" in patched
    assert "scratch_n=max0(mu+1,2)" in patched
    assert "vi(i)=0.0" in patched
    assert "deallocate(xj,xy,Rn,IDnY,RDnX,RDnY,IGnX,RGnX," in patched


def test_patch_main_source_rewrites_program_entrypoint() -> None:
    upstream_main = Path("tmp/6s_upstream/main.f")
    if not upstream_main.exists():
        pytest.skip("Upstream 6S source tree is not available for main.f transform validation.")

    original = upstream_main.read_text()

    patched = patch_main_source(original)

    assert "subroutine sixs_case_core" in patched
    assert "ier_out=.FALSE." in patched
    assert "igeom=0" in patched
    assert "iwave=1" in patched
    assert "output_values_out(" in patched
    assert f"parameter (sixs_output_count={len(SIXS_CORE_OUTPUT_SPECS)})" in patched
    assert "return" in patched


def test_patch_threadprivate_directives_handles_common_block_styles() -> None:
    original = "\n".join(
        [
            "      common /foo/ a, b",
            "      common/bar/c, d",
            "      common /baz/ e,",
            "     &f, g",
        ]
    )

    patched = patch_threadprivate_directives(original)

    assert "c$omp threadprivate(/foo/)" in patched
    assert "c$omp threadprivate(/bar/)" in patched
    assert "c$omp threadprivate(/baz/)" in patched


def test_encoded_string_io_reports_utf8_encoding() -> None:
    stream = _EncodedStringIO()

    assert stream.encoding == "utf-8"


def test_default_native_profile_maps_to_user_water_ozone() -> None:
    assert _resolve_atmospheric_mode(SixSAlgorithmConfig()) == 8


def test_spectral_response_build_uses_band_rsrf_support() -> None:
    band = SensorBand(
        name="B04",
        center_wavelength=665.0,
        bandwidth=30.0,
        resolution=10.0,
        band_index=3,
        rsrf_wavelengths_nm=np.array([640.0, 650.0, 665.0, 680.0, 690.0], dtype=np.float32),
        rsrf_response=np.array([0.0, 0.4, 1.0, 0.5, 0.0], dtype=np.float32),
    )

    response, wlinf, wlsup = _build_spectral_response(band)

    assert response.shape == (1501,)
    assert 0.64 <= wlinf <= 0.67
    assert 0.66 <= wlsup <= 0.69
    assert np.isclose(np.trapezoid(response, 0.25 + 0.0025 * np.arange(1501)), 1.0, atol=1e-3)


def test_backend_propagates_observation_time_to_runner() -> None:
    class _FakeRunner:
        def __init__(self) -> None:
            self.observation_time: datetime | None = None

        def set_observation_time(self, observation_time: datetime | None) -> None:
            self.observation_time = observation_time

        def compute_coefficients(self, **kwargs):
            _ = kwargs
            template = _sample_atmo().aot
            return {
                "xap": xr.ones_like(template),
                "xbp": xr.zeros_like(template),
                "xcp": xr.zeros_like(template),
                "tgasm": xr.full_like(template, 0.8),
            }

        def compute_coefficients_multi(self, **kwargs):
            return [self.compute_coefficients(**kwargs) for _ in kwargs["bands"]]

        def preload_scene_subset(self, *args, **kwargs):
            _ = (args, kwargs)
            return None

    runner = _FakeRunner()
    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(output_variables=("xap", "xbp", "xcp", "tgasm")),
        runner=runner,
    )
    obs_time = datetime(2025, 7, 12, 10, 30)

    backend.set_observation_time(obs_time)
    coeffs = backend.compute_coefficients(
        _sample_geometry(),
        _sample_atmo(),
        SensorBand(
            name="B02",
            center_wavelength=490.0,
            bandwidth=20.0,
            resolution=10.0,
            band_index=1,
        ),
    )

    assert runner.observation_time == obs_time
    assert "tgasm" in coeffs.extras


def test_resolve_build_paths_separates_release_and_parity_roots() -> None:
    release_paths = resolve_build_paths(SixSAlgorithmConfig())
    parity_paths = resolve_build_paths(SixSAlgorithmConfig(build_profile="parity"))

    assert release_paths.root_dir != parity_paths.root_dir
    assert release_paths.root_dir.name == "release"
    assert parity_paths.root_dir.name == "parity"


def test_resolve_f2py_backends_skips_unavailable_distutils(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sixs_build_module,
        "_distutils_backend_supported",
        lambda: (False, "setuptools 82 is too new"),
    )
    monkeypatch.setattr(sixs_build_module.shutil, "which", lambda tool: "/usr/bin/" + tool)

    assert _resolve_f2py_backends() == ("meson",)


def test_distutils_backend_supported_rejects_new_setuptools(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("numpy.f2py.f2py2e.MESON_ONLY_VER", False)
    monkeypatch.setattr(
        sixs_build_module.importlib.metadata,
        "version",
        lambda name: "82.0.1" if name == "setuptools" else "0",
    )

    supported, reason = _distutils_backend_supported()

    assert supported is False
    assert reason is not None
    assert "setuptools 82.0.1 is too new" in reason


def test_resolve_f2py_backends_rejects_invalid_distutils_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        sixs_build_module,
        "_distutils_backend_supported",
        lambda: (False, "NumPy F2PY is Meson-only on Python 3.12+"),
    )
    monkeypatch.setenv("SIAC_SIXS_F2PY_BACKEND", "distutils")

    with pytest.raises(RuntimeError, match="distutils backend"):
        _resolve_f2py_backends()


def test_distutils_extra_link_args_include_gcc_runtime_for_openmp() -> None:
    command = _distutils_extra_link_args(["-O3", "-fopenmp"])

    assert command == ["-lgomp", "-lquadmath", "-lm"]


def test_compile_f2py_extension_accepts_built_module_after_nonzero_exit(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "main.f").write_text("      end\n", encoding="utf-8")

    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))
    built_module = build_paths.root_dir / f"{build_paths.module_name}.cpython-test.so"
    find_calls = {"count": 0}

    def _fake_find_built_extension(_paths):
        find_calls["count"] += 1
        if find_calls["count"] < 2:
            return None
        return built_module

    def _fake_run_distutils_backend(**_kwargs):
        built_module.parent.mkdir(parents=True, exist_ok=True)
        built_module.write_text("built", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["python", "-m", "numpy.f2py"],
            returncode=1,
            stdout="running build",
            stderr='buildmodule: Could not find the body of interfaced routine "sixs_case_core". Skipping.',
        )

    monkeypatch.setattr(sixs_build_module, "parse_makefile_sources", lambda _path: [source_dir / "main.f"])
    monkeypatch.setattr(sixs_build_module, "_resolve_f2py_backends", lambda: ("distutils",))
    monkeypatch.setattr(sixs_build_module, "_run_distutils_backend", _fake_run_distutils_backend)
    monkeypatch.setattr(sixs_build_module, "find_built_extension", _fake_find_built_extension)

    _compile_f2py_extension(
        source_dir=source_dir,
        build_paths=build_paths,
        compiler="gfortran",
        build_profile="release",
    )

    assert built_module.exists()


def test_find_built_extension_searches_recent_extra_roots(tmp_path: Path) -> None:
    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))
    extra_root = tmp_path / "workspace"
    extra_root.mkdir(parents=True, exist_ok=True)
    module_path = extra_root / f"{build_paths.module_name}.cpython-test.so"
    module_path.write_text("module", encoding="utf-8")

    found = _find_built_extension(
        build_paths,
        extra_roots=(extra_root,),
        min_mtime=time.time() - 5.0,
    )

    assert found == module_path


def test_scene_lut_plan_and_auto_selection_reduce_native_case_count() -> None:
    case_arrays = {
        "sza_deg": np.full(32, 25.0, dtype=np.float64),
        "saa_deg": np.full(32, 110.0, dtype=np.float64),
        "vza_deg": np.full(32, 5.0, dtype=np.float64),
        "vaa_deg": np.full(32, 95.0, dtype=np.float64),
        "aot550": np.linspace(0.05, 0.4, 32, dtype=np.float64),
        "tcwv_cm": np.linspace(0.8, 3.0, 32, dtype=np.float64),
        "tco3_atmcm": np.full(32, 0.3, dtype=np.float64),
        "elevation_km": np.full(32, 0.2, dtype=np.float64),
    }

    plan = _build_scene_lut_plan(case_arrays, max_nodes_per_axis=3, max_cases=81)

    assert plan.direct_case_count == 32
    assert plan.lut_case_count <= 81
    assert _should_use_scene_lut(
        mode="auto",
        direct_case_count=plan.direct_case_count,
        lut_case_count=plan.lut_case_count,
        min_pixels=8,
        required_speedup=1.1,
    ) is True
    assert _should_use_scene_lut(
        mode="direct",
        direct_case_count=plan.direct_case_count,
        lut_case_count=plan.lut_case_count,
        min_pixels=8,
        required_speedup=1.1,
    ) is False


def test_scene_lut_mode_interpolates_linear_native_outputs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    geometry = GeometryAngles.from_degrees(
        xr.DataArray(np.full((3, 3), 25.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((3, 3), 110.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((3, 3), 5.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((3, 3), 95.0, dtype=np.float32), dims=("y", "x")),
    )
    atmo = AtmosphericState(
        aot=xr.DataArray(
            [[0.10, 0.12, 0.14], [0.11, 0.13, 0.15], [0.12, 0.14, 0.16]],
            dims=("y", "x"),
        ),
        tcwv=xr.DataArray(
            [[1.5, 1.6, 1.7], [1.55, 1.65, 1.75], [1.6, 1.7, 1.8]],
            dims=("y", "x"),
        ),
        tco3=xr.DataArray(
            np.full((3, 3), 0.30, dtype=np.float32),
            dims=("y", "x"),
        ),
        aot_unc=xr.DataArray(np.full((3, 3), 0.01, dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(np.full((3, 3), 0.05, dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(np.full((3, 3), 0.01, dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(
            np.full((3, 3), 0.2, dtype=np.float32),
            dims=("y", "x"),
        ),
    )
    band = SensorBand(
        name="B04",
        center_wavelength=665.0,
        bandwidth=30.0,
        resolution=10.0,
        band_index=3,
    )
    direct_runner = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="direct",
            output_variables=("xap", "xbp", "xcp", "tgasm"),
        )
    )._runner
    scene_runner = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            output_variables=("xap", "xbp", "xcp", "tgasm"),
            scene_lut_max_nodes_per_axis=2,
            scene_lut_max_cases=256,
        )
    )._runner

    call_sizes: list[int] = []

    def _fake_run_native_batch(**kwargs):
        n_cases = int(np.asarray(kwargs["sza_deg"]).size)
        call_sizes.append(n_cases)
        sza = np.asarray(kwargs["sza_deg"], dtype=np.float64)
        aot = np.asarray(kwargs["aot550"], dtype=np.float64)
        tcwv = np.asarray(kwargs["tcwv_cm"], dtype=np.float64)
        outputs = {
            "xap": 0.01 * sza + 0.1 * aot + 0.05 * tcwv,
            "xbp": 0.001 * sza + 0.02 * aot,
            "xcp": 0.002 * tcwv + 0.01 * aot,
            "tgasm": 0.8 + 0.001 * sza - 0.01 * aot,
        }
        return _NativeBatchResult(
            outputs={name: np.ascontiguousarray(values, dtype=np.float64) for name, values in outputs.items()},
            status=np.zeros(n_cases, dtype=np.int32),
        )

    monkeypatch.setattr(direct_runner, "_run_native_batch", _fake_run_native_batch)
    monkeypatch.setattr(scene_runner, "_run_native_batch", _fake_run_native_batch)

    direct = direct_runner.compute_coefficients(geometry=geometry, atmo_state=atmo, band=band, output_variables=("xap", "xbp", "xcp", "tgasm"))
    direct_case_count = call_sizes[-1]
    scene = scene_runner.compute_coefficients(geometry=geometry, atmo_state=atmo, band=band, output_variables=("xap", "xbp", "xcp", "tgasm"))
    scene_case_count = call_sizes[-1]

    assert scene_case_count < direct_case_count
    for name in ("xap", "xbp", "xcp", "tgasm"):
        np.testing.assert_allclose(scene[name].values, direct[name].values, rtol=1.0e-10, atol=1.0e-10)


def test_worker_library_backend_slices_batches_and_merges_outputs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )
    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            parallel_backend="worker_libraries",
            worker_libraries=2,
            chunk_size=2,
        )
    )
    runner = backend._runner

    class _FakeSession:
        def __init__(self) -> None:
            self.chunk_sizes: list[int] = []

        def run_batch(self, *, n_threads: int, **kwargs):
            _ = n_threads
            sza = np.asarray(kwargs["sza_deg"], dtype=np.float64)
            self.chunk_sizes.append(int(sza.size))
            outputs = {
                name: np.ascontiguousarray(np.zeros(sza.size, dtype=np.float64), dtype=np.float64)
                for name in ("xap", "xbp", "xcp")
            }
            outputs["xap"] = sza.copy()
            outputs["xbp"] = sza + 1.0
            outputs["xcp"] = sza + 2.0
            return _NativeBatchResult(outputs=outputs, status=np.zeros(sza.size, dtype=np.int32))

        def close(self) -> None:
            return None

    sessions = [_FakeSession(), _FakeSession()]
    monkeypatch.setattr(runner, "_ensure_worker_sessions", lambda worker_count: sessions[:worker_count])

    result = runner._run_native_batch_worker_libraries(
        month=1,
        day=1,
        atmospheric_mode=8,
        radiosonde_altitude_km=np.zeros(34, dtype=np.float64),
        radiosonde_pressure_mb=np.zeros(34, dtype=np.float64),
        radiosonde_temperature_k=np.zeros(34, dtype=np.float64),
        radiosonde_water_g_m3=np.zeros(34, dtype=np.float64),
        radiosonde_ozone_g_m3=np.zeros(34, dtype=np.float64),
        aerosol_mode=1,
        aerosol_mixture=np.zeros(4, dtype=np.float64),
        aerosol_distribution_rmin=0.0,
        aerosol_distribution_rmax=0.0,
        aerosol_distribution_component_count=0,
        aerosol_distribution_x1=np.zeros(4, dtype=np.float64),
        aerosol_distribution_x2=np.zeros(4, dtype=np.float64),
        aerosol_distribution_x3=np.zeros(4, dtype=np.float64),
        aerosol_distribution_cij=np.zeros(4, dtype=np.float64),
        aerosol_distribution_rn=np.zeros((20, 4), dtype=np.float64),
        aerosol_distribution_ri=np.zeros((20, 4), dtype=np.float64),
        aerosol_sun_count=0,
        aerosol_sun_radius=np.zeros(50, dtype=np.float64),
        aerosol_sun_dvlogr=np.zeros(50, dtype=np.float64),
        aerosol_layer_count=0,
        aerosol_layer_height=np.zeros(50, dtype=np.float64),
        aerosol_layer_aot=np.zeros(50, dtype=np.float64),
        aerosol_layer_type=np.zeros(50, dtype=np.int32),
        reference_reflectance=0.1,
        spectral_wlinf=0.64,
        spectral_wlsup=0.68,
        spectral_response=np.zeros(1501, dtype=np.float64),
        aerosol_model_path="",
        surface_inhomo=0,
        surface_idirec=0,
        surface_target_mode=0,
        surface_target_constant=0.0,
        surface_target_spectrum=np.zeros(1501, dtype=np.float64),
        surface_env_mode=0,
        surface_env_constant=0.0,
        surface_env_spectrum=np.zeros(1501, dtype=np.float64),
        surface_radius_km=1.0,
        surface_brdf_model=0,
        surface_brdf_params=np.zeros(12, dtype=np.float64),
        surface_brdf_options=np.zeros(5, dtype=np.int32),
        surface_brdf_struct=np.zeros(4, dtype=np.float64),
        surface_brdf_optics=np.zeros(3, dtype=np.float64),
        surface_brdf_table_solar=np.zeros((10, 13), dtype=np.float64),
        surface_brdf_table_view=np.zeros((10, 13), dtype=np.float64),
        surface_brdf_spherical_albedo=0.0,
        surface_brdf_directional_reflectance=0.0,
        atmospheric_correction_mode=0,
        atmospheric_correction_value=-0.1,
        sza_deg=np.arange(5, dtype=np.float64),
        saa_deg=np.arange(5, dtype=np.float64),
        vza_deg=np.arange(5, dtype=np.float64),
        vaa_deg=np.arange(5, dtype=np.float64),
        aot550=np.arange(5, dtype=np.float64),
        tcwv_cm=np.arange(5, dtype=np.float64),
        tco3_atmcm=np.arange(5, dtype=np.float64),
        elevation_km=np.arange(5, dtype=np.float64),
    )

    np.testing.assert_allclose(result.outputs["xap"], np.arange(5, dtype=np.float64))
    np.testing.assert_allclose(result.outputs["xbp"], np.arange(5, dtype=np.float64) + 1.0)
    np.testing.assert_allclose(result.outputs["xcp"], np.arange(5, dtype=np.float64) + 2.0)
    assert sessions[0].chunk_sizes == [2, 1]
    assert sessions[1].chunk_sizes == [2]

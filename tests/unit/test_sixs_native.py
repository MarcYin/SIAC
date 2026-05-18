from __future__ import annotations

import io
import subprocess
import tarfile
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.direct import sixs_build as sixs_build_module
from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import (
    _build_f2py_command,
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
    JointGridSearchLUT,
    _build_joint_grid_search_lut_plan,
    _build_scene_lut_plan,
    _build_spectral_response,
    _NativeBatchResult,
    _resolve_atmospheric_columns_mode,
    _resolve_atmospheric_mode,
    _should_use_scene_lut,
)
from siac.config import RTSetupConfig, SixSAlgorithmConfig
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


def test_patch_main_source_supports_profile_column_scaling() -> None:
    upstream_main = Path("tmp/6s_upstream/main.f")
    if not upstream_main.exists():
        pytest.skip("Upstream 6S source tree is not available for main.f transform validation.")

    patched = patch_main_source(upstream_main.read_text())

    assert "atmospheric_columns_mode_in" in patched
    assert "idatm_abstra=idatm" in patched
    assert "uwus=max(uw,1.0e-06)" in patched
    assert "if (idatm_abstra.eq.8) then" in patched
    assert "call abstra(idatm_abstra,wl,xmus,xmuv,uw,uo3,uwus,uo3us," in patched


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


def test_default_native_profile_maps_to_us_standard_62() -> None:
    assert (
        _resolve_atmospheric_mode(
            RTSetupConfig(
                atmosphere={"profile": "us_standard_62", "columns_mode": "input_columns"}
            ),
            month=1,
        )
        == 6
    )


def test_default_native_columns_mode_uses_scene_inputs() -> None:
    assert (
        _resolve_atmospheric_columns_mode(
            RTSetupConfig(atmosphere={"profile": "us_standard_62", "columns_mode": "input_columns"})
        )
        == 1
    )
    assert (
        _resolve_atmospheric_columns_mode(
            RTSetupConfig(
                atmosphere={"profile": "us_standard_62", "columns_mode": "profile_default"}
            )
        )
        == 0
    )
    assert (
        _resolve_atmospheric_columns_mode(
            RTSetupConfig(atmosphere={"profile": "no_gas", "columns_mode": "input_columns"})
        )
        == 0
    )


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


def test_backend_merges_partial_rt_setup_with_native_defaults() -> None:
    class _FakeRunner:
        def set_observation_time(self, observation_time: datetime | None) -> None:
            _ = observation_time

        def compute_coefficients(self, **kwargs):
            _ = kwargs
            template = _sample_atmo().aot
            return {
                "xap": xr.ones_like(template),
                "xbp": xr.zeros_like(template),
                "xcp": xr.zeros_like(template),
            }

        def compute_coefficients_multi(self, **kwargs):
            return [self.compute_coefficients(**kwargs) for _ in kwargs["bands"]]

        def preload_scene_subset(self, *args, **kwargs):
            _ = (args, kwargs)
            return None

    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(),
        rt_setup=RTSetupConfig(
            surface={
                "mode": "homogeneous_brdf",
                "brdf": {
                    "model": "rahman",
                    "parameters": {
                        "intensity": 0.12,
                        "asymmetry_factor": 0.03,
                        "structural_parameter": 0.45,
                    },
                },
            },
            atmospheric_correction={"mode": "brdf_reflectance", "value": 0.2},
            reference_reflectance=0.2,
        ),
        runner=_FakeRunner(),
    )

    assert backend.rt_setup.atmosphere is not None
    assert backend.rt_setup.atmosphere.profile == "us_standard_62"
    assert backend.rt_setup.aerosol is not None
    assert backend.rt_setup.aerosol.profile == "continental"
    assert backend.rt_setup.surface is not None
    assert backend.rt_setup.surface.mode == "homogeneous_brdf"
    assert backend.rt_setup.atmospheric_correction is not None
    assert backend.rt_setup.atmospheric_correction.mode == "brdf_reflectance"


def test_resolve_build_paths_separates_release_and_parity_roots() -> None:
    release_paths = resolve_build_paths(SixSAlgorithmConfig())
    parity_paths = resolve_build_paths(SixSAlgorithmConfig(build_profile="parity"))

    assert release_paths.root_dir != parity_paths.root_dir
    assert release_paths.root_dir.name == "release"
    assert parity_paths.root_dir.name == "parity"


def test_fetch_and_unpack_source_rejects_path_traversal_archive(tmp_path: Path) -> None:
    archive_path = tmp_path / "bad-6s.tar"
    payload = b"unsafe"
    with tarfile.open(archive_path, "w") as archive:
        member = tarfile.TarInfo("../escape.txt")
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    with pytest.raises(RuntimeError, match="unsafe 6S source archive member"):
        sixs_build_module._fetch_and_unpack_source(
            "https://example.invalid/6s.tar",
            archive_path,
            tmp_path / "upstream",
        )

    assert not (tmp_path / "escape.txt").exists()


def test_resolve_f2py_backends_skips_unavailable_distutils(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sixs_build_module,
        "_distutils_backend_supported",
        lambda: (False, "setuptools 82 is too new"),
    )
    monkeypatch.setattr(sixs_build_module.shutil, "which", lambda tool: "/usr/bin/" + tool)

    assert _resolve_f2py_backends() == ("meson",)


def test_distutils_backend_supported_rejects_new_setuptools(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    monkeypatch.chdir(tmp_path)
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "main.f").write_text("      end\n", encoding="utf-8")

    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))
    built_module = build_paths.root_dir / f"{build_paths.module_name}.cpython-test.so"
    find_calls = {"count": 0}

    def _fake_find_built_extension(_paths, **_kwargs):
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

    monkeypatch.setattr(
        sixs_build_module, "parse_makefile_sources", lambda _path: [source_dir / "main.f"]
    )
    monkeypatch.setattr(
        sixs_build_module,
        "_generate_f2py_signature",
        lambda **_kwargs: source_dir / "_siac_rt6s_native.pyf",
    )
    monkeypatch.setattr(sixs_build_module, "_resolve_f2py_backends", lambda: ("distutils",))
    monkeypatch.setattr(sixs_build_module, "_run_distutils_backend", _fake_run_distutils_backend)
    monkeypatch.setattr(sixs_build_module, "_find_built_extension", _fake_find_built_extension)
    monkeypatch.setattr(
        sixs_build_module,
        "_validate_extension_import",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            args=["python", "-c", "import module"],
            returncode=0,
            stdout=str(built_module),
            stderr="",
        ),
    )

    _compile_f2py_extension(
        source_dir=source_dir,
        build_paths=build_paths,
        compiler="gfortran",
        build_profile="release",
    )

    assert built_module.exists()


def test_compile_f2py_extension_persists_backend_diagnostics_on_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "main.f").write_text("      end\n", encoding="utf-8")

    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))

    def _fake_run_distutils_backend(**_kwargs):
        return subprocess.CompletedProcess(
            args=["python", "-m", "numpy.f2py"],
            returncode=1,
            stdout="distutils stdout",
            stderr="distutils stderr",
        )

    monkeypatch.setattr(
        sixs_build_module, "parse_makefile_sources", lambda _path: [source_dir / "main.f"]
    )
    monkeypatch.setattr(
        sixs_build_module,
        "_generate_f2py_signature",
        lambda **_kwargs: source_dir / "_siac_rt6s_native.pyf",
    )
    monkeypatch.setattr(sixs_build_module, "_resolve_f2py_backends", lambda: ("distutils",))
    monkeypatch.setattr(sixs_build_module, "_run_distutils_backend", _fake_run_distutils_backend)
    monkeypatch.setattr(sixs_build_module, "_find_built_extension", lambda *_args, **_kwargs: None)

    with pytest.raises(RuntimeError, match="6S native Python extension build failed"):
        _compile_f2py_extension(
            source_dir=source_dir,
            build_paths=build_paths,
            compiler="gfortran",
            build_profile="release",
        )

    diagnostics_dir = build_paths.root_dir / "diagnostics"
    assert (diagnostics_dir / "f2py-distutils.stdout.txt").read_text(
        encoding="utf-8"
    ) == "distutils stdout"
    assert (diagnostics_dir / "f2py-distutils.stderr.txt").read_text(
        encoding="utf-8"
    ) == "distutils stderr"
    assert "status=failed" in (diagnostics_dir / "f2py-distutils.summary.txt").read_text(
        encoding="utf-8"
    )
    assert "Backend distutils" in (diagnostics_dir / "build_failure_summary.txt").read_text(
        encoding="utf-8"
    )


def test_compile_f2py_extension_accepts_module_found_in_environment_roots(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "main.f").write_text("      end\n", encoding="utf-8")

    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))
    env_site = tmp_path / "site-packages"
    built_module = env_site / f"{build_paths.module_name}.cpython-test.so"

    def _fake_run_distutils_backend(**_kwargs):
        built_module.parent.mkdir(parents=True, exist_ok=True)
        built_module.write_text("built", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["python", "-m", "numpy.f2py"],
            returncode=1,
            stdout="running build",
            stderr='buildmodule: Could not find the body of interfaced routine "sixs_case_core". Skipping.',
        )

    monkeypatch.setattr(
        sixs_build_module, "parse_makefile_sources", lambda _path: [source_dir / "main.f"]
    )
    monkeypatch.setattr(
        sixs_build_module,
        "_generate_f2py_signature",
        lambda **_kwargs: source_dir / "_siac_rt6s_native.pyf",
    )
    monkeypatch.setattr(sixs_build_module, "_resolve_f2py_backends", lambda: ("distutils",))
    monkeypatch.setattr(sixs_build_module, "_run_distutils_backend", _fake_run_distutils_backend)
    monkeypatch.setattr(sixs_build_module, "_environment_extension_roots", lambda: (env_site,))
    monkeypatch.setattr(
        sixs_build_module,
        "_validate_extension_import",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            args=["python", "-c", "import module"],
            returncode=0,
            stdout=str(built_module),
            stderr="",
        ),
    )

    _compile_f2py_extension(
        source_dir=source_dir,
        build_paths=build_paths,
        compiler="gfortran",
        build_profile="release",
    )

    assert built_module.exists()


def test_compile_f2py_extension_falls_back_after_import_validation_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "main.f").write_text("      end\n", encoding="utf-8")

    build_paths = resolve_build_paths(SixSAlgorithmConfig(build_dir=str(tmp_path / "build")))
    meson_module = build_paths.root_dir / f"{build_paths.module_name}.cpython-meson.so"
    distutils_module = build_paths.root_dir / f"{build_paths.module_name}.cpython-distutils.so"
    find_calls = {"count": 0}

    def _fake_find_built_extension(_paths, **_kwargs):
        find_calls["count"] += 1
        if find_calls["count"] == 1:
            return None
        if find_calls["count"] == 2:
            return meson_module
        return distutils_module

    def _fake_subprocess_run(*args, **kwargs):
        meson_module.parent.mkdir(parents=True, exist_ok=True)
        meson_module.write_text("meson", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["python", "-m", "numpy.f2py", "--backend", "meson"],
            returncode=0,
            stdout="meson build ok",
            stderr="",
        )

    def _fake_run_distutils_backend(**_kwargs):
        distutils_module.parent.mkdir(parents=True, exist_ok=True)
        distutils_module.write_text("distutils", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["python", "-m", "numpy.f2py", "--backend", "distutils"],
            returncode=0,
            stdout="distutils build ok",
            stderr="",
        )

    def _fake_validate_extension_import(_module_name, module_path):
        module_path = Path(module_path)
        if module_path == meson_module:
            return subprocess.CompletedProcess(
                args=["python", "-c", "import module"],
                returncode=1,
                stdout="",
                stderr="ImportError: undefined symbol: GOMP_parallel",
            )
        return subprocess.CompletedProcess(
            args=["python", "-c", "import module"],
            returncode=0,
            stdout=str(module_path),
            stderr="",
        )

    monkeypatch.setattr(
        sixs_build_module, "parse_makefile_sources", lambda _path: [source_dir / "main.f"]
    )
    monkeypatch.setattr(
        sixs_build_module,
        "_generate_f2py_signature",
        lambda **_kwargs: source_dir / "_siac_rt6s_native.pyf",
    )
    monkeypatch.setattr(sixs_build_module, "_resolve_f2py_backends", lambda: ("meson", "distutils"))
    monkeypatch.setattr(sixs_build_module.subprocess, "run", _fake_subprocess_run)
    monkeypatch.setattr(sixs_build_module, "_run_distutils_backend", _fake_run_distutils_backend)
    monkeypatch.setattr(sixs_build_module, "_find_built_extension", _fake_find_built_extension)
    monkeypatch.setattr(sixs_build_module, "_environment_extension_roots", lambda: ())
    monkeypatch.setattr(
        sixs_build_module, "_validate_extension_import", _fake_validate_extension_import
    )

    _compile_f2py_extension(
        source_dir=source_dir,
        build_paths=build_paths,
        compiler="gfortran",
        build_profile="release",
    )

    diagnostics_dir = build_paths.root_dir / "diagnostics"
    meson_summary = (diagnostics_dir / "f2py-meson.summary.txt").read_text(encoding="utf-8")
    assert "status=import-failed" in meson_summary
    assert "import_check_returncode=1" in meson_summary
    assert "GOMP_parallel" in (diagnostics_dir / "f2py-meson.import_check.stderr.txt").read_text(
        encoding="utf-8"
    )
    assert not meson_module.exists()
    assert distutils_module.exists()


def test_build_f2py_command_uses_signature_file_before_sources(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    signature_path = tmp_path / "_siac_rt6s_native.pyf"
    first_source = source_dir / "main.f"
    second_source = source_dir / "siac_rt6s_bridge.f90"
    first_source.write_text("      end\n", encoding="utf-8")
    second_source.write_text("end\n", encoding="utf-8")

    command = _build_f2py_command(
        backend="meson",
        module_name="_siac_rt6s_native",
        source_dir=source_dir,
        signature_path=signature_path,
        compile_sources=[first_source, second_source],
        flags=["-O3", "-fopenmp"],
        f2py_build_dir=tmp_path / "f2py_build",
    )

    assert str(signature_path.resolve()) in command
    assert command.index(str(signature_path.resolve())) < command.index(str(first_source.resolve()))
    assert "only:" not in command


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
    assert (
        _should_use_scene_lut(
            mode="auto",
            direct_case_count=plan.direct_case_count,
            lut_case_count=plan.lut_case_count,
            min_pixels=8,
            required_speedup=1.1,
        )
        is True
    )
    assert (
        _should_use_scene_lut(
            mode="direct",
            direct_case_count=plan.direct_case_count,
            lut_case_count=plan.lut_case_count,
            min_pixels=8,
            required_speedup=1.1,
        )
        is False
    )


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
    # parallel_backend="openmp" forces the path that goes through
    # ``_run_native_batch`` (which the test below monkey-patches). The
    # default is ``worker_libraries``, which would skip the mock and try
    # to load real .so copies — see wave 18 for the band-parallel
    # implementation that introduced that.
    direct_runner = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="direct",
            output_variables=("xap", "xbp", "xcp", "tgasm"),
            parallel_backend="openmp",
        )
    )._runner
    scene_runner = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            output_variables=("xap", "xbp", "xcp", "tgasm"),
            scene_lut_max_nodes_per_axis=2,
            scene_lut_max_cases=256,
            parallel_backend="openmp",
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
            outputs={
                name: np.ascontiguousarray(values, dtype=np.float64)
                for name, values in outputs.items()
            },
            status=np.zeros(n_cases, dtype=np.int32),
        )

    monkeypatch.setattr(direct_runner, "_run_native_batch", _fake_run_native_batch)
    monkeypatch.setattr(scene_runner, "_run_native_batch", _fake_run_native_batch)

    direct = direct_runner.compute_coefficients(
        geometry=geometry,
        atmo_state=atmo,
        band=band,
        output_variables=("xap", "xbp", "xcp", "tgasm"),
    )
    direct_case_count = call_sizes[-1]
    scene = scene_runner.compute_coefficients(
        geometry=geometry,
        atmo_state=atmo,
        band=band,
        output_variables=("xap", "xbp", "xcp", "tgasm"),
    )
    scene_case_count = call_sizes[-1]

    assert scene_case_count < direct_case_count
    for name in ("xap", "xbp", "xcp", "tgasm"):
        np.testing.assert_allclose(
            scene[name].values, direct[name].values, rtol=1.0e-10, atol=1.0e-10
        )


def test_joint_grid_search_lut_plan_preserves_aot_tcwv_axes() -> None:
    """The joint LUT plan must always preserve the caller-supplied aot/tcwv axes.

    Even when the case-count budget is tight, only the geometric axes should
    be coarsened — the (aot, tcwv) axes are the whole point of the joint LUT
    and must equal the grid-search points so the lookup is exact at each
    candidate.
    """
    case_arrays = {
        "sza_deg": np.linspace(20.0, 40.0, 64, dtype=np.float64),
        "saa_deg": np.linspace(100.0, 130.0, 64, dtype=np.float64),
        "vza_deg": np.linspace(0.0, 8.0, 64, dtype=np.float64),
        "vaa_deg": np.linspace(80.0, 110.0, 64, dtype=np.float64),
        "aot550": np.full(64, 0.15, dtype=np.float64),  # ignored — axis comes from caller
        "tcwv_cm": np.full(64, 2.0, dtype=np.float64),  # ignored — axis comes from caller
        "tco3_atmcm": np.linspace(0.28, 0.34, 64, dtype=np.float64),
        "elevation_km": np.linspace(0.0, 0.5, 64, dtype=np.float64),
    }
    aot_axis = np.linspace(0.05, 0.8, 11, dtype=np.float64)
    tcwv_axis = np.linspace(0.0, 6.0, 11, dtype=np.float64)

    plan = _build_joint_grid_search_lut_plan(
        case_arrays,
        aot_axis=aot_axis,
        tcwv_axis=tcwv_axis,
        max_nodes_per_axis=4,
        # Tight budget: 11*11 = 121 for aot×tcwv alone, leaves ~2 nodes per
        # geometric axis after trimming.
        max_cases=2048,
    )

    np.testing.assert_array_equal(plan.axes["aot550"], aot_axis)
    np.testing.assert_array_equal(plan.axes["tcwv_cm"], tcwv_axis)
    assert plan.lut_case_count <= 2048
    assert plan.direct_case_count == 64


def test_joint_grid_search_lut_evaluate_matches_direct_at_grid_points(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """At grid-search nodes, joint-LUT lookup must equal a fresh direct compute.

    The key correctness invariant: because the LUT's (aot, tcwv) axes are
    exactly the grid-search axes, evaluating at any (aot_node, tcwv_node)
    pair should yield identical coefficients to running the scene through
    the per-candidate (direct) path with that aot/tcwv. Only the geometric
    dimensions are linearly interpolated, but the synthetic kernel below is
    linear so even those are exact.
    """
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    geometry = GeometryAngles.from_degrees(
        xr.DataArray(np.full((4, 4), 25.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((4, 4), 110.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((4, 4), 5.0, dtype=np.float32), dims=("y", "x")),
        xr.DataArray(np.full((4, 4), 95.0, dtype=np.float32), dims=("y", "x")),
    )
    # Per-pixel aot/tcwv prior values are irrelevant to the joint LUT —
    # the LUT axes come from the caller, not from the prior.
    atmo = AtmosphericState(
        aot=xr.DataArray(0.15 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        tcwv=xr.DataArray(2.0 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        tco3=xr.DataArray(0.30 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        aot_unc=xr.DataArray(0.01 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        tcwv_unc=xr.DataArray(0.05 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        tco3_unc=xr.DataArray(0.01 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
        elevation=xr.DataArray(0.2 * np.ones((4, 4), dtype=np.float32), dims=("y", "x")),
    )
    band = SensorBand(
        name="B04",
        center_wavelength=665.0,
        bandwidth=30.0,
        resolution=10.0,
        band_index=3,
    )

    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            output_variables=("xap", "xbp", "xcp"),
            scene_lut_max_nodes_per_axis=2,
            scene_lut_max_cases=256,
            joint_grid_search_lut_max_nodes_per_axis=2,
            joint_grid_search_lut_max_cases=4096,
            # Force the openmp path so the _run_native_batch mock below
            # is exercised. The default worker_libraries path bypasses it
            # and loads real .so sessions.
            parallel_backend="openmp",
        )
    )
    runner = backend._runner

    def _fake_run_native_batch(**kwargs):
        n_cases = int(np.asarray(kwargs["sza_deg"]).size)
        aot = np.asarray(kwargs["aot550"], dtype=np.float64)
        tcwv = np.asarray(kwargs["tcwv_cm"], dtype=np.float64)
        sza = np.asarray(kwargs["sza_deg"], dtype=np.float64)
        # Linear in everything so interpolation is exact at any sample point.
        outputs = {
            "xap": 0.1 * aot + 0.05 * tcwv + 0.001 * sza,
            "xbp": 0.02 * aot + 0.001 * tcwv,
            "xcp": 0.01 * aot + 0.002 * tcwv,
        }
        outputs.update(
            {
                # All other native outputs unused by the test but required by
                # the bundle; supply zeros.
                name: np.zeros(n_cases, dtype=np.float64)
                for name in (
                    "tgasm",
                    "totg",
                    "rho_atm",
                    "rho_atm_pol",
                    "trans_solar",
                    "trans_view",
                    "spher_albedo",
                    "raylscatd",
                    "aerscatd",
                    "raylscatu",
                    "aerscatu",
                    "transm_h2o",
                    "transm_o3",
                    "transm_other",
                    "dwn_irr_dir",
                    "dwn_irr_dif",
                    "dwn_irr_env",
                    "rad_path",
                )
            }
        )
        # Make every payload-relevant output backed by float64.
        return _NativeBatchResult(
            outputs={
                name: np.ascontiguousarray(values, dtype=np.float64)
                for name, values in outputs.items()
            },
            status=np.zeros(n_cases, dtype=np.int32),
        )

    monkeypatch.setattr(runner, "_run_native_batch", _fake_run_native_batch)

    aot_axis = np.array([0.10, 0.20, 0.30], dtype=np.float64)
    tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float64)

    joint = backend.build_joint_grid_search_lut(
        geometry=geometry,
        atmo_state=atmo,
        aot_axis=aot_axis,
        tcwv_axis=tcwv_axis,
        bands=[band],
    )
    assert joint is not None

    # Compare joint-LUT outputs at every grid node against a fresh
    # compute_coefficients call with the same (aot, tcwv) value broadcast.
    for aot_val in aot_axis:
        for tcwv_val in tcwv_axis:
            atmo_at_node = AtmosphericState(
                aot=xr.DataArray(
                    float(aot_val) * np.ones((4, 4), dtype=np.float32),
                    dims=("y", "x"),
                ),
                tcwv=xr.DataArray(
                    float(tcwv_val) * np.ones((4, 4), dtype=np.float32),
                    dims=("y", "x"),
                ),
                tco3=atmo.tco3,
                aot_unc=atmo.aot_unc,
                tcwv_unc=atmo.tcwv_unc,
                tco3_unc=atmo.tco3_unc,
                elevation=atmo.elevation,
            )
            direct = backend.compute_coefficients(
                geometry, atmo_at_node, band, compute_jacobian=False
            )
            band_outputs = joint.evaluate(float(aot_val), float(tcwv_val))
            joint_xap = band_outputs[0]["xap"].values
            joint_xbp = band_outputs[0]["xbp"].values
            joint_xcp = band_outputs[0]["xcp"].values
            # Tolerance: per-pixel values pass through a float32 DataArray
            # template at one step in the pipeline, so absolute differences
            # at the float32 epsilon scale (~1e-7 for unit-magnitude values,
            # less for these ~0.01-scale coefficients) are expected.
            np.testing.assert_allclose(
                joint_xap, direct.xap.values, rtol=1.0e-6, atol=1.0e-9
            )
            np.testing.assert_allclose(
                joint_xbp, direct.xbp.values, rtol=1.0e-6, atol=1.0e-9
            )
            np.testing.assert_allclose(
                joint_xcp, direct.xcp.values, rtol=1.0e-6, atol=1.0e-9
            )


def test_joint_grid_search_lut_disabled_returns_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the config opt-out flag is False, build_joint_grid_search_lut returns None."""
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            joint_grid_search_lut_enabled=False,
        )
    )

    result = backend.build_joint_grid_search_lut(
        geometry=_sample_geometry((4, 4)),
        atmo_state=_sample_atmo((4, 4)),
        aot_axis=np.linspace(0.05, 0.8, 3),
        tcwv_axis=np.linspace(0.0, 6.0, 3),
        bands=[
            SensorBand(
                name="B04",
                center_wavelength=665.0,
                bandwidth=30.0,
                resolution=10.0,
                band_index=3,
            )
        ],
    )
    assert result is None


def test_joint_grid_search_lut_direct_mode_returns_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When sixs.mode == 'direct', joint LUT is disabled (consistent with scene LUT)."""
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    backend = SixSBackend(sixs_config=SixSAlgorithmConfig(mode="direct"))

    result = backend.build_joint_grid_search_lut(
        geometry=_sample_geometry((4, 4)),
        atmo_state=_sample_atmo((4, 4)),
        aot_axis=np.linspace(0.05, 0.8, 3),
        tcwv_axis=np.linspace(0.0, 6.0, 3),
        bands=[
            SensorBand(
                name="B04",
                center_wavelength=665.0,
                bandwidth=30.0,
                resolution=10.0,
                band_index=3,
            )
        ],
    )
    assert result is None


def test_joint_grid_search_lut_amortises_native_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One LUT build + N evaluate() calls must invoke 6S fewer times than N
    direct compute_coefficients invocations.

    This is the whole performance rationale of the optimization: the inner
    block-grid-search runs the coefficient provider once per (aot, tcwv)
    candidate. With the joint LUT, all candidates share a single 6S batch.
    """
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            output_variables=("xap", "xbp", "xcp"),
            scene_lut_max_nodes_per_axis=2,
            scene_lut_max_cases=256,
            joint_grid_search_lut_max_nodes_per_axis=2,
            joint_grid_search_lut_max_cases=4096,
            # The mock below patches _run_native_batch — force the openmp
            # path through it instead of the band-parallel worker_libraries
            # path that goes via session.run_batch directly.
            parallel_backend="openmp",
        )
    )
    runner = backend._runner

    call_counter = {"count": 0}

    def _counting_run_native_batch(**kwargs):
        call_counter["count"] += 1
        n_cases = int(np.asarray(kwargs["sza_deg"]).size)
        outputs = {
            name: np.zeros(n_cases, dtype=np.float64)
            for name in (
                "xap", "xbp", "xcp", "tgasm", "totg", "rho_atm", "rho_atm_pol",
                "trans_solar", "trans_view", "spher_albedo", "raylscatd",
                "aerscatd", "raylscatu", "aerscatu", "transm_h2o", "transm_o3",
                "transm_other", "dwn_irr_dir", "dwn_irr_dif", "dwn_irr_env",
                "rad_path",
            )
        }
        return _NativeBatchResult(
            outputs={
                name: np.ascontiguousarray(values, dtype=np.float64)
                for name, values in outputs.items()
            },
            status=np.zeros(n_cases, dtype=np.int32),
        )

    monkeypatch.setattr(runner, "_run_native_batch", _counting_run_native_batch)

    geometry = _sample_geometry((4, 4))
    atmo = _sample_atmo((4, 4))
    bands = [
        SensorBand(name=f"B{i:02d}", center_wavelength=400.0 + 50 * i,
                   bandwidth=30.0, resolution=10.0, band_index=i)
        for i in range(3)
    ]

    aot_axis = np.linspace(0.05, 0.8, 5, dtype=np.float64)
    tcwv_axis = np.linspace(0.0, 6.0, 5, dtype=np.float64)
    n_candidates = aot_axis.size * tcwv_axis.size  # 25

    # Joint path: 1 build × N_bands native calls, then N_candidates × 0 native calls.
    call_counter["count"] = 0
    joint = backend.build_joint_grid_search_lut(
        geometry=geometry,
        atmo_state=atmo,
        aot_axis=aot_axis,
        tcwv_axis=tcwv_axis,
        bands=bands,
    )
    assert joint is not None
    joint_build_calls = call_counter["count"]
    for aot_val in aot_axis:
        for tcwv_val in tcwv_axis:
            _ = joint.evaluate(float(aot_val), float(tcwv_val))
    joint_total_calls = call_counter["count"]
    assert joint_build_calls == len(bands), (
        f"Expected one 6S batch per band at build time, got {joint_build_calls}"
    )
    assert joint_total_calls == joint_build_calls, (
        f"evaluate() must not run native 6S; got {joint_total_calls - joint_build_calls} stray calls"
    )

    # Direct path: N_candidates × N_bands native calls.
    call_counter["count"] = 0
    for aot_val in aot_axis:
        for tcwv_val in tcwv_axis:
            atmo_at = AtmosphericState(
                aot=xr.DataArray(
                    float(aot_val) * np.ones((4, 4), dtype=np.float32), dims=("y", "x")
                ),
                tcwv=xr.DataArray(
                    float(tcwv_val) * np.ones((4, 4), dtype=np.float32), dims=("y", "x")
                ),
                tco3=atmo.tco3,
                aot_unc=atmo.aot_unc,
                tcwv_unc=atmo.tcwv_unc,
                tco3_unc=atmo.tco3_unc,
                elevation=atmo.elevation,
            )
            for band in bands:
                backend.compute_coefficients(
                    geometry, atmo_at, band, compute_jacobian=False
                )
    direct_calls = call_counter["count"]
    assert direct_calls == n_candidates * len(bands)

    # The headline assertion — the joint path must use far fewer 6S calls.
    assert joint_total_calls < direct_calls, (
        f"Joint LUT ({joint_total_calls}) should use fewer 6S calls than "
        f"direct ({direct_calls})."
    )


def test_joint_grid_search_lut_band_parallel_runs_each_band_on_a_session(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Wave 18: with parallel_backend='worker_libraries' the joint-LUT band
    loop dispatches each band's 6S batch to a different isolated library
    session in parallel, rather than re-using one OpenMP-shared session.

    The test verifies (a) each band's call goes through ``session.run_batch``
    rather than the OpenMP-shared ``_run_native_batch``, and (b) the
    final per-band outputs preserve band order.
    """
    monkeypatch.setattr(
        "siac.algorithms.rt.direct.sixs_native.ensure_native_sixs_module",
        lambda _config: Path("/tmp/fake_sixs_native.so"),
    )

    backend = SixSBackend(
        sixs_config=SixSAlgorithmConfig(
            mode="scene_lut",
            output_variables=("xap", "xbp", "xcp"),
            scene_lut_max_nodes_per_axis=2,
            scene_lut_max_cases=256,
            joint_grid_search_lut_max_nodes_per_axis=2,
            joint_grid_search_lut_max_cases=4096,
            parallel_backend="worker_libraries",
            worker_libraries=3,
        )
    )
    runner = backend._runner

    session_call_counts: list[int] = []

    class _StubSession:
        def __init__(self, marker: int) -> None:
            self._marker = marker
            self.calls = 0

        def run_batch(self, *, n_threads: int, **kwargs):
            _ = n_threads
            self.calls += 1
            n_cases = int(np.asarray(kwargs["sza_deg"]).size)
            # Encode the band's spectral_wlinf into xap so the test can
            # verify each result came from the band it was issued for.
            wlinf = float(kwargs["spectral_wlinf"])
            outputs = {
                name: np.zeros(n_cases, dtype=np.float64)
                for name in (
                    "xap", "xbp", "xcp", "tgasm", "totg", "rho_atm",
                    "rho_atm_pol", "trans_solar", "trans_view", "spher_albedo",
                    "raylscatd", "aerscatd", "raylscatu", "aerscatu",
                    "transm_h2o", "transm_o3", "transm_other", "dwn_irr_dir",
                    "dwn_irr_dif", "dwn_irr_env", "rad_path",
                )
            }
            outputs["xap"] = np.full(n_cases, wlinf, dtype=np.float64)
            return _NativeBatchResult(
                outputs={
                    name: np.ascontiguousarray(values, dtype=np.float64)
                    for name, values in outputs.items()
                },
                status=np.zeros(n_cases, dtype=np.int32),
            )

    stub_sessions = [_StubSession(i) for i in range(3)]

    def _fake_ensure_workers(worker_count: int):
        return stub_sessions[:worker_count]

    monkeypatch.setattr(runner, "_ensure_worker_sessions", _fake_ensure_workers)

    # If anything went wrong and the code falls back to _run_native_batch,
    # this raises and the test fails clearly rather than silently passing.
    def _explode_run_native_batch(**_kwargs):
        raise AssertionError(
            "_run_native_batch must not be called when band-parallel "
            "worker_libraries path is in use."
        )

    monkeypatch.setattr(runner, "_run_native_batch", _explode_run_native_batch)

    bands = [
        SensorBand(
            name=f"B{i:02d}",
            # Use distinct, increasing wavelengths so xap encodes a
            # band-unique value once spectral_wlinf is extracted.
            center_wavelength=400.0 + 50.0 * i,
            bandwidth=30.0,
            resolution=10.0,
            band_index=i,
        )
        for i in range(5)
    ]
    aot_axis = np.array([0.10, 0.20, 0.30], dtype=np.float64)
    tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float64)

    joint = backend.build_joint_grid_search_lut(
        geometry=_sample_geometry((4, 4)),
        atmo_state=_sample_atmo((4, 4)),
        aot_axis=aot_axis,
        tcwv_axis=tcwv_axis,
        bands=bands,
    )
    assert joint is not None
    assert joint.band_count == len(bands)

    # 5 bands across 3 sessions → at least each session should have run
    # at least once, and the total calls should equal the band count.
    total_session_calls = sum(s.calls for s in stub_sessions)
    assert total_session_calls == len(bands), (
        f"Expected one session call per band ({len(bands)}); "
        f"got {total_session_calls} (per-session: "
        f"{[s.calls for s in stub_sessions]})"
    )

    # Verify band-order is preserved in the returned outputs. xap was
    # encoded with each band's spectral_wlinf — extract back and check
    # bands appear in their input order.
    for band_idx, band in enumerate(bands):
        # All non-NaN xap values for this band should equal that band's wlinf.
        outputs_for_band = joint.evaluate(0.20, 2.0)[band_idx]
        xap_arr = np.asarray(outputs_for_band["xap"].values).ravel()
        finite = xap_arr[np.isfinite(xap_arr)]
        # Every finite value should be near the band's spectral lower bound.
        # Our SensorBand has bandwidth=30.0 (nm) so wlinf = center/1000 - 0.015 µm.
        expected_wlinf = (band.center_wavelength - 15.0) / 1000.0
        assert np.allclose(finite, expected_wlinf, rtol=1e-3, atol=1e-3), (
            f"Band {band.name}: expected xap≈{expected_wlinf:.4f}, "
            f"got {finite[:3]}"
        )


def test_worker_library_backend_slices_batches_and_merges_outputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    monkeypatch.setattr(
        runner, "_ensure_worker_sessions", lambda worker_count: sessions[:worker_count]
    )

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

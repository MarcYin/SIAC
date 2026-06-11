"""Fetch, patch, and build the native 6SV2.1 Python extension."""

from __future__ import annotations

import contextlib
import importlib.machinery
import importlib.metadata
import io
import logging
import os
import re
import shlex
import shutil
import site
import subprocess
import sys
import sysconfig
import tarfile
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from siac.algorithms.rt.direct._build_common import archive_sha256, fetch_archive
from siac.sixs_outputs import SIXS_CORE_OUTPUT_SPECS

if TYPE_CHECKING:
    from siac.config.algorithms import SixSAlgorithmConfig

logger = logging.getLogger(__name__)

_ARCHIVE_NAME = "6sV2.1.tar"
_UPSTREAM_DIRNAME = "upstream"
_PATCHED_DIRNAME = "patched"
_F2PY_BUILD_DIRNAME = "f2py_build"
_DIAGNOSTICS_DIRNAME = "diagnostics"
_MODULE_BASENAME = "_siac_rt6s_native"
_SIGNATURE_FILENAME = f"{_MODULE_BASENAME}.pyf"
_PATCHED_MAIN_SENTINEL = "subroutine sixs_case_core"
_COMMON_BLOCK_RE = re.compile(r"^\s*common\s*/([^/]+)/", re.IGNORECASE)

_FORTRAN_DIR = Path(__file__).with_name("fortran")
_BRIDGE_SOURCE_NAME = "siac_rt6s_bridge.f90"

_DISCOM_OLD = "\n".join(
    [
        "        if ((l.lt.20).and.(wldis(l).lt.wlinf).and.",
        "     a     (wldis(l+1).lt.wlinf)) goto 50",
        "        if ((l.gt.1).and.(wldis(l).gt.wlsup).and.",
        "     a      (wldis(l-1).gt.wlsup)) goto 50",
    ]
)
_DISCOM_NEW = "\n".join(
    [
        "        if (l.lt.20) then",
        "          if ((wldis(l).lt.wlinf).and.(wldis(l+1).lt.wlinf)) goto 50",
        "        endif",
        "        if (l.gt.1) then",
        "          if ((wldis(l).gt.wlsup).and.(wldis(l-1).gt.wlsup)) goto 50",
        "        endif",
    ]
)
_KERNELPOL_OLD = "\n".join(
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
_KERNELPOL_NEW = "\n".join(
    [
        "        psl(2,j)=xdb",
        "        psl(2,-j)=xdb",
        "        rsl(0,j)=0.0D+00",
        "        rsl(0,-j)=0.0D+00",
        "        rsl(1,j)=0.0D+00",
        "        rsl(1,-j)=0.0D+00",
        "        xdb=3.D+00*(1.D+00-c*c)/2.D+00/sqrt(6.D+00)",
        "\tif (abs(xdb).lt.1.e-30) xdb =0.0D+00",
        "        rsl(2,j)=xdb",
        "        rsl(2,-j)=xdb",
        "        tsl(0,j)=0.0D+00",
        "        tsl(0,-j)=0.0D+00",
        "        tsl(1,j)=0.0D+00",
        "        tsl(1,-j)=0.0D+00",
        "        tsl(2,j)=0.0D+00",
        "        tsl(2,-j)=0.0D+00",
    ]
)


@dataclass(frozen=True)
class SixSBuildPaths:
    """Resolved paths used for the native 6S build."""

    root_dir: Path
    archive_path: Path
    upstream_dir: Path
    patched_dir: Path
    f2py_build_dir: Path
    module_name: str
    module_hint_path: Path | None


class _EncodedStringIO(io.StringIO):
    """String buffer that reports a real encoding for numpy.distutils."""

    @property
    def encoding(self) -> str:
        return "utf-8"


def _core_output_assignment_lines() -> list[str]:
    lines = ["      ier_out=ier"]
    for index, (_name, expression, condition) in enumerate(SIXS_CORE_OUTPUT_SPECS, start=1):
        target = f"output_values_out({index})"
        if condition is None:
            lines.append(f"      {target}={expression}")
        else:
            lines.append(f"      if ({condition}) {target}={expression}")
    return lines


def _main_subroutine_signature() -> str:
    return "\n".join(
        [
            "      subroutine sixs_case_core(iwr_unit,asol_in,phi0_in,avis_in,",
            "     s phiv_in,month_in,jday_in,idatm_in,",
            "     s atmospheric_columns_mode_in,uw_in,uo3_in,",
            "     s radiosonde_altitude_in,radiosonde_pressure_in,",
            "     s radiosonde_temperature_in,radiosonde_water_in,",
            "     s radiosonde_ozone_in,iaer_in,aerosol_mixture_in,",
            "     s aerosol_dist_rmin_in,aerosol_dist_rmax_in,",
            "     s aerosol_dist_icp_in,aerosol_dist_x1_in,",
            "     s aerosol_dist_x2_in,aerosol_dist_x3_in,",
            "     s aerosol_dist_cij_in,aerosol_dist_rn_in,",
            "     s aerosol_dist_ri_in,aerosol_sun_count_in,",
            "     s aerosol_sun_radius_in,aerosol_sun_dvlogr_in,",
            "     s aerosol_layer_count_in,aerosol_layer_height_in,",
            "     s aerosol_layer_aot_in,aerosol_layer_type_in,",
            "     s aerosol_filename_in,taer55_in,elevation_km_in,",
            "     s wlinf_in,wlsup_in,spectral_response_in,",
            "     s surface_inhomo_in,surface_idirec_in,",
            "     s surface_target_mode_in,surface_target_constant_in,",
            "     s surface_target_spectrum_in,surface_env_mode_in,",
            "     s surface_env_constant_in,surface_env_spectrum_in,",
            "     s surface_radius_km_in,surface_brdf_model_in,",
            "     s surface_brdf_params_in,surface_brdf_options_in,",
            "     s surface_brdf_struct_in,surface_brdf_optics_in,",
            "     s surface_brdf_table_solar_in,surface_brdf_table_view_in,",
            "     s surface_brdf_spherical_albedo_in,",
            "     s surface_brdf_directional_reflectance_in,",
            "     s reference_reflectance_in,irapp_in,rapp_in,ier_out,",
            "     s output_values_out)",
        ]
    )


def _main_runtime_initialization() -> str:
    return "\n".join(
        [
            "      iwr=iwr_unit",
            "      ier=.FALSE.",
            "      ier_out=.FALSE.",
            "      iinf=1",
            "      rpfet_report=0.",
            "      plumet_report=0.",
            "      xpol_report=0.",
            "      rpfet_over_refet_report=0.",
            "      rocave_report=0.",
            "      robar1_over_xnorm1_report=0.",
            "      robar2_over_xnorm2_report=0.",
            "      rbard_report=0.",
            "      albbrdf_report=0.",
            "      rfoamave_report=0.",
            "      rwatave_report=0.",
            "      rglitave_report=0.",
        ]
    )


def shared_library_suffix() -> str:
    if sys.platform == "darwin":
        return ".dylib"
    return ".so"


def native_module_suffixes() -> tuple[str, ...]:
    return tuple(importlib.machinery.EXTENSION_SUFFIXES)


def default_build_root() -> Path:
    return Path.home().expanduser() / ".cache" / "siac" / "rt6s"


def resolve_build_paths(config: SixSAlgorithmConfig) -> SixSBuildPaths:
    default_root = default_build_root() / str(getattr(config, "build_profile", "release"))
    root_dir = Path(config.build_dir or default_root).expanduser()
    module_hint_path = None
    configured_module_path = getattr(config, "module_path", None)
    if configured_module_path is not None:
        module_hint_path = Path(configured_module_path).expanduser()
    return SixSBuildPaths(
        root_dir=root_dir,
        archive_path=root_dir / _ARCHIVE_NAME,
        upstream_dir=root_dir / _UPSTREAM_DIRNAME,
        patched_dir=root_dir / _PATCHED_DIRNAME,
        f2py_build_dir=root_dir / _F2PY_BUILD_DIRNAME,
        module_name=_MODULE_BASENAME,
        module_hint_path=module_hint_path,
    )


def _compiler_flags(build_profile: str) -> list[str]:
    optimization_flag = "-O0" if build_profile == "parity" else "-O3"
    return [
        optimization_flag,
        "-fopenmp",
        "-frecursive",
        "-ffixed-line-length-132",
        "-fallow-argument-mismatch",
        "-std=legacy",
    ]


def _resolve_f2py_backends() -> tuple[str, ...]:
    distutils_supported, distutils_reason = _distutils_backend_supported()
    override = os.environ.get("SIAC_SIXS_F2PY_BACKEND")
    if override:
        requested = tuple(part.strip().lower() for part in override.split(",") if part.strip())
        valid = {"meson", "distutils"}
        invalid = [backend for backend in requested if backend not in valid]
        if invalid:
            raise RuntimeError(
                "Invalid SIAC_SIXS_F2PY_BACKEND override. "
                f"Unsupported backends: {', '.join(invalid)}."
            )
        if "distutils" in requested and not distutils_supported:
            raise RuntimeError(
                "SIAC_SIXS_F2PY_BACKEND requested the distutils backend, "
                f"but it is unavailable: {distutils_reason}"
            )
        return requested or ("distutils",)

    backends: list[str] = []
    if shutil.which("meson") is not None and shutil.which("ninja") is not None:
        backends.append("meson")
    if distutils_supported:
        backends.append("distutils")
    elif distutils_reason is not None:
        logger.info("Skipping unavailable F2PY distutils backend: %s", distutils_reason)
    if not backends:
        raise RuntimeError(
            "No supported F2PY backends are available. "
            "Install meson+ninja for the Meson backend, or use Python 3.11 "
            "with setuptools < 60 to enable the distutils backend."
        )
    return tuple(backends)


def _distutils_backend_supported() -> tuple[bool, str | None]:
    try:
        from numpy.f2py.f2py2e import MESON_ONLY_VER
    except Exception as exc:  # pragma: no cover - defensive import guard
        return False, f"unable to inspect NumPy F2PY backend support ({exc})"

    if MESON_ONLY_VER:
        return False, "NumPy F2PY is Meson-only on Python 3.12+"

    try:
        setuptools_version = importlib.metadata.version("setuptools")
    except importlib.metadata.PackageNotFoundError:
        return False, "setuptools is not installed"

    match = re.match(r"(\d+)", setuptools_version)
    if match is not None and int(match.group(1)) >= 60:
        return False, (
            f"setuptools {setuptools_version} is too new for numpy.distutils; use setuptools < 60"
        )
    return True, None


def _is_fixed_form_continuation(line: str) -> bool:
    if not line:
        return False
    if line[0] in {"c", "C", "*", "!"}:
        return False
    return len(line) > 5 and line[5] not in {" ", "0"}


def patch_threadprivate_directives(text: str) -> str:
    """Append OpenMP THREADPRIVATE directives for named COMMON blocks."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        match = _COMMON_BLOCK_RE.match(line)
        if match:
            name = match.group(1).strip()
            if name:
                j = i + 1
                while j < len(lines) and _is_fixed_form_continuation(lines[j]):
                    out.append(lines[j])
                    j += 1
                directive = f"c$omp threadprivate(/{name}/)"
                already_present = j < len(lines) and lines[j].strip().lower() == directive.lower()
                if not already_present:
                    out.append(directive)
                i = j
                continue
        i += 1
    return "\n".join(out) + ("\n" if text.endswith("\n") or not text else "")


def patch_discom_source(text: str) -> str:
    """Apply the out-of-bounds fix to DISCOM.f."""
    if _DISCOM_NEW in text:
        return text
    if _DISCOM_OLD not in text:
        raise RuntimeError("Failed to patch DISCOM.f boundary checks.")
    return text.replace(_DISCOM_OLD, _DISCOM_NEW, 1)


def patch_kernelpol_source(text: str) -> str:
    """Apply the zero-order kernel initialization fix to KERNELPOL.f."""
    if _KERNELPOL_NEW in text:
        return text
    if _KERNELPOL_OLD not in text:
        raise RuntimeError("Failed to patch KERNELPOL.f polarization kernel initialization.")
    return text.replace(_KERNELPOL_OLD, _KERNELPOL_NEW, 1)


def patch_aeroso_source(text: str) -> str:
    """Widen the Mie-file path storage in AEROSO.f."""
    old = "      character FILE*80"
    new = "      character FILE*1024"
    if new in text:
        patched = text
    else:
        patched = _replace_once(text, old, new, "AEROSO.f Mie-file path length")

    save_old = "\n".join(
        [
            "      if (iaer.ge.8.and.iaer.le.11) then",
            "        open(10,file=FILE)",
        ]
    )
    save_new = "\n".join(
        [
            "      if (iaer.ge.8.and.iaer.le.11.and.len_trim(FILE).gt.0) then",
            "        open(10,file=FILE)",
        ]
    )
    if save_new not in patched:
        patched = _replace_once(
            patched, save_old, save_new, "AEROSO.f optional Mie-file output guard"
        )
    return patched


def patch_mie_source(text: str) -> str:
    """Make the Mie kernels safe for OpenMP-backed native execution."""
    if "double precision, allocatable ::" in text:
        return text

    text = _replace_once(
        text,
        "\n".join(
            [
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
            ]
        ),
        "\n".join(
            [
                "      double precision Ren,Imn,X,Y,Up,XnumRDnY,XnumIDnY",
                "      double precision XdenDnY,coxj,Qsca,Qext,xJonH,XdenGNX",
                "      double precision Xnum1An,Xnum2An,XdenAn,Xden1An,Xden2An,RAnb,IAnb",
                "      double precision Xnum1Bn,Xnum2Bn,XdenBn,Xden1Bn,Xden2Bn,RBnb,IBnb",
                "      double precision xmud,xpond,RS1,RS2,IS1,IS2,co_n,test",
                "      integer scratch_n",
                "      double precision, allocatable ::",
                "     & xj(:),xy(:),Rn(:),IDnY(:),RDnX(:),RDnY(:)",
                "      double precision, allocatable ::",
                "     & IGnX(:),RGnX(:),RAn(:),IAn(:),RBn(:),IBn(:)",
                "      double precision, allocatable ::",
                "     & TAUn(:),PIn(:)",
            ]
        ),
        "MIE.f EXSCPHASE scratch storage declarations",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      do i=1,icp",
                "        np(i)=0.D+00",
            ]
        ),
        "\n".join(
            [
                "      do i=1,icp",
                "        np(i)=0.D+00",
                "        vi(i)=0.0",
            ]
        ),
        "MIE.f particle-volume initialization",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
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
            ]
        ),
        "\n".join(
            [
                "      if (mu.le.0) then",
                '         write(6,*) " Error, mu is too small, mu is equal to : ",mu',
                "         Stop",
                "         endif",
                "      scratch_n=max0(mu+1,2)",
                "      allocate(xj(0:scratch_n),xy(-1:scratch_n),Rn(0:scratch_n),",
                "     & IDnY(0:scratch_n),RDnX(0:scratch_n),RDnY(0:scratch_n),",
                "     & IGnX(0:scratch_n),RGnX(0:scratch_n),RAn(0:scratch_n),",
                "     & IAn(0:scratch_n),RBn(0:scratch_n),IBn(0:scratch_n),",
                "     & TAUn(0:scratch_n),PIn(0:scratch_n))",
                "      xj=0.D+00",
                "      xy=0.D+00",
                "      Rn=0.D+00",
                "      IDnY=0.D+00",
                "      RDnX=0.D+00",
                "      RDnY=0.D+00",
                "      IGnX=0.D+00",
                "      RGnX=0.D+00",
                "      RAn=0.D+00",
                "      IAn=0.D+00",
                "      RBn=0.D+00",
                "      IBn=0.D+00",
                "      TAUn=0.D+00",
                "      PIn=0.D+00",
            ]
        ),
        "MIE.f EXSCPHASE dynamic scratch allocation",
    )
    tail_old = "\n".join(
        [
            "      enddo",
            "               ",
            "      return",
            "      end",
        ]
    )
    tail_new = "\n".join(
        [
            "      enddo",
            "      if (allocated(xj)) then",
            "        deallocate(xj,xy,Rn,IDnY,RDnX,RDnY,IGnX,RGnX,",
            "     &   RAn,IAn,RBn,IBn,TAUn,PIn)",
            "      endif",
            "      return",
            "      end",
        ]
    )
    if tail_old not in text:
        raise RuntimeError("Failed to patch MIE.f EXSCPHASE cleanup.")
    head, _, tail = text.rpartition(tail_old)
    text = f"{head}{tail_new}{tail}"
    return text


def patch_main_source(text: str) -> str:
    """Turn upstream main.f into a typed SIAC case subroutine."""
    if _PATCHED_MAIN_SENTINEL in text:
        return text
    text = _replace_once(
        text,
        "      program ssssss",
        _main_subroutine_signature(),
        "main.f program entrypoint",
    )
    text = _replace_once(
        text,
        "      character FILE*80,FILE2*80",
        "      character FILE*1024,FILE2*1024,aerosol_filename_in*1024",
        "main.f aerosol filename storage",
    )
    text = _replace_once(
        text,
        "        integer igrou1,igrou2,isort,irapp,ilut",
        "\n".join(
            [
                "        integer igrou1,igrou2,isort,irapp,ilut",
                "        integer iwr_unit,month_in,jday_in,idatm_in,iaer_in",
                "        integer atmospheric_columns_mode_in,idatm_abstra",
                "        integer aerosol_dist_icp_in,aerosol_sun_count_in",
                "        integer aerosol_layer_count_in,surface_inhomo_in",
                "        integer surface_idirec_in,surface_target_mode_in",
                "        integer surface_env_mode_in,surface_brdf_model_in",
                "        integer irapp_in",
                "        integer aerosol_layer_type_in(50)",
                "        integer surface_brdf_options_in(5)",
                "        integer sixs_output_count",
                f"        parameter (sixs_output_count={len(SIXS_CORE_OUTPUT_SPECS)})",
                "        logical ier_out",
                "        real asol_in,phi0_in,avis_in,phiv_in,uw_in,uo3_in",
                "        real radiosonde_altitude_in(34)",
                "        real radiosonde_pressure_in(34)",
                "        real radiosonde_temperature_in(34)",
                "        real radiosonde_water_in(34)",
                "        real radiosonde_ozone_in(34)",
                "        real taer55_in,elevation_km_in,wlinf_in,wlsup_in",
                "        real spectral_response_in(1501),aerosol_mixture_in(4)",
                "        real aerosol_dist_rmin_in,aerosol_dist_rmax_in",
                "        real aerosol_dist_x1_in(4),aerosol_dist_x2_in(4)",
                "        real aerosol_dist_x3_in(4),aerosol_dist_cij_in(4)",
                "        real aerosol_dist_rn_in(20,4),aerosol_dist_ri_in(20,4)",
                "        real aerosol_sun_radius_in(50)",
                "        real aerosol_sun_dvlogr_in(50)",
                "        real aerosol_layer_height_in(50)",
                "        real aerosol_layer_aot_in(50)",
                "        real surface_target_constant_in,surface_target_spectrum_in(1501)",
                "        real surface_env_constant_in,surface_env_spectrum_in(1501)",
                "        real surface_radius_km_in",
                "        real surface_brdf_params_in(12)",
                "        real surface_brdf_struct_in(4)",
                "        real surface_brdf_optics_in(3)",
                "        real surface_brdf_table_solar_in(10,13)",
                "        real surface_brdf_table_view_in(10,13)",
                "        real surface_brdf_spherical_albedo_in",
                "        real surface_brdf_directional_reflectance_in",
                "        real reference_reflectance_in,rapp_in",
                "        real output_values_out(sixs_output_count)",
                "        real rpfet_report,plumet_report,xpol_report",
                "        real rpfet_over_refet_report",
                "        real rocave_report,robar1_over_xnorm1_report",
                "        real robar2_over_xnorm2_report,rbard_report",
                "        real albbrdf_report,rfoamave_report",
                "        real rwatave_report,rglitave_report",
            ]
        ),
        "main.f output declarations",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      iwr=6",
                "      ier=.FALSE.",
                "      iinf=1",
            ]
        ),
        _main_runtime_initialization(),
        "main.f runtime initialization",
    )
    text = _replace_once(text, "      read(iread,*) igeom", "      igeom=0", "geometry mode input")

    # Wave 16 (performance): the upstream 6S hardcodes ``ipol=1`` to compute
    # the full polarized radiative transfer (Stokes I, Q, U via ospol_ /
    # kernelpol_). The polarization kernels dominate runtime when the
    # view geometry is significantly off-nadir. SIAC's downstream pipeline
    # uses only the scalar (I) component for atmospheric correction —
    # path reflectance, transmittance, spherical albedo — so the Q/U
    # work is wasted compute. Patch to ``ipol=0`` by default so the
    # unpolarized scalar RT runs ~3-5x faster.
    #
    # Set the environment variable ``SIAC_SIXS_COMPUTE_POLARIZATION=1``
    # at build time to opt back into the polarized path (e.g. for
    # science work that needs Stokes parameters; note the bridge
    # currently doesn't expose them downstream so this only affects
    # the scalar-I result via second-order polarization corrections,
    # which are typically <1% at typical S2 geometry).
    _polarization_env = os.environ.get("SIAC_SIXS_COMPUTE_POLARIZATION", "0").strip().lower()
    if _polarization_env in {"0", "false", "no", "off", ""}:
        text = _replace_once(
            text,
            "       ipol=1",
            "\n".join(
                [
                    "c      SIAC build: polarized RT disabled for performance",
                    "c      (set SIAC_SIXS_COMPUTE_POLARIZATION=1 at build time to revert)",
                    "       ipol=0",
                ]
            ),
            "ipol polarization toggle (SIAC build: default unpolarized)",
        )
    else:
        logger.info(
            "SIAC_SIXS_COMPUTE_POLARIZATION=%s — keeping upstream ipol=1; "
            "the build will run the full polarized RT (slower).",
            _polarization_env,
        )
    text = _replace_once(
        text,
        "      read(iread,*) asol,phi0,avis,phiv,month,jday",
        "\n".join(
            [
                "      asol=asol_in",
                "      phi0=phi0_in",
                "      avis=avis_in",
                "      phiv=phiv_in",
                "      month=month_in",
                "      jday=jday_in",
            ]
        ),
        "geometry scalar input",
    )
    text = _replace_once(
        text, "      read(iread,*) idatm", "      idatm=idatm_in", "atmosphere mode input"
    )
    text = _replace_once(
        text,
        "      if(idatm.eq.8) read(iread,*) uw,uo3",
        "\n".join(
            [
                "      if(idatm.eq.8) then",
                "        uw=uw_in",
                "        uo3=uo3_in",
                "      endif",
            ]
        ),
        "atmosphere water and ozone input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      do 7 k=1,34",
                "       read(iread,*) z(k),p(k),t(k),wh(k),wo(k)",
                "    7 continue",
            ]
        ),
        "\n".join(
            [
                "      do 7 k=1,34",
                "       z(k)=radiosonde_altitude_in(k)",
                "       p(k)=radiosonde_pressure_in(k)",
                "       t(k)=radiosonde_temperature_in(k)",
                "       wh(k)=radiosonde_water_in(k)",
                "       wo(k)=radiosonde_ozone_in(k)",
                "    7 continue",
            ]
        ),
        "user_profile atmospheric input",
    )
    text = _replace_once(
        text, "      read(iread,*) iaer", "      iaer=iaer_in", "aerosol mode input"
    )
    text = _replace_once(
        text,
        "      read(5,*) num_z",
        "      num_z=min(max(aerosol_layer_count_in,0),50)",
        "layered aerosol profile count input",
    )
    text = _replace_once(
        text,
        "       read(5,*) height_z(num_z-i),taer55_z(num_z-i),iaer",
        "\n".join(
            [
                "       height_z(num_z-i)=aerosol_layer_height_in(num_z-i)",
                "       taer55_z(num_z-i)=aerosol_layer_aot_in(num_z-i)",
                "       iaer=aerosol_layer_type_in(num_z-i)",
            ]
        ),
        "layered aerosol profile layer input",
    )
    text = _replace_once(
        text,
        "      if (iaer.ge.0.and.iaer.le.7) nquad=nqdef_p",
        "      if ((iaer.ge.0.and.iaer.le.7).or.(iaer.eq.12)) nquad=nqdef_p",
        "aerosol quadrature selection",
    )
    text = _replace_once(
        text,
        "      if(iaer.eq.4) read(iread,*) (c(n),n=1,4)",
        "\n".join(
            [
                "      if(iaer.eq.4) then",
                "        do n=1,4",
                "          c(n)=aerosol_mixture_in(n)",
                "        enddo",
                "      endif",
            ]
        ),
        "user aerosol mixture input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "   43 read(iread,*) rmin,rmax,icp",
                "      do i=1,icp",
                "       read(5,*)x1(i),x2(i),cij(i)",
                "       read(5,*)(rn(l,i),l=1,20)",
                "       read(5,*)(ri(l,i),l=1,20)",
                "      enddo",
                "        do i=1,icp",
                "         cij_out(i)=cij(i)",
                "        enddo",
                "      go to 49",
            ]
        ),
        "\n".join(
            [
                "   43 rmin=aerosol_dist_rmin_in",
                "      rmax=aerosol_dist_rmax_in",
                "      icp=max(1,min(aerosol_dist_icp_in,4))",
                "      do i=1,icp",
                "       x1(i)=aerosol_dist_x1_in(i)",
                "       x2(i)=aerosol_dist_x2_in(i)",
                "       cij(i)=aerosol_dist_cij_in(i)",
                "       cij_out(i)=cij(i)",
                "       do l=1,20",
                "        rn(l,i)=aerosol_dist_rn_in(l,i)",
                "        ri(l,i)=aerosol_dist_ri_in(l,i)",
                "       enddo",
                "      enddo",
                "      go to 49",
            ]
        ),
        "multimodal aerosol distribution input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "   44 read(iread,*) rmin,rmax",
                "      read(iread,*) x1(1),x2(1),x3(1)",
                "      read(5,*)(rn(l,1),l=1,20)",
                "      read(5,*)(ri(l,1),l=1,20)",
                "      go to 49",
            ]
        ),
        "\n".join(
            [
                "   44 rmin=aerosol_dist_rmin_in",
                "      rmax=aerosol_dist_rmax_in",
                "      icp=max(1,min(aerosol_dist_icp_in,4))",
                "      do i=1,icp",
                "       x1(i)=aerosol_dist_x1_in(i)",
                "       x2(i)=aerosol_dist_x2_in(i)",
                "       x3(i)=aerosol_dist_x3_in(i)",
                "       cij(i)=aerosol_dist_cij_in(i)",
                "       cij_out(i)=cij(i)",
                "       do l=1,20",
                "        rn(l,i)=aerosol_dist_rn_in(l,i)",
                "        ri(l,i)=aerosol_dist_ri_in(l,i)",
                "       enddo",
                "      enddo",
                "      go to 49",
            ]
        ),
        "modified gamma aerosol input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "   45 read(iread,*) rmin,rmax",
                "      read(iread,*) x1(1)",
                "      read(5,*)(rn(l,1),l=1,20)",
                "      read(5,*)(ri(l,1),l=1,20)",
                "      go to 49",
            ]
        ),
        "\n".join(
            [
                "   45 rmin=aerosol_dist_rmin_in",
                "      rmax=aerosol_dist_rmax_in",
                "      icp=max(1,min(aerosol_dist_icp_in,4))",
                "      do i=1,icp",
                "       x1(i)=aerosol_dist_x1_in(i)",
                "       cij(i)=aerosol_dist_cij_in(i)",
                "       cij_out(i)=cij(i)",
                "       do l=1,20",
                "        rn(l,i)=aerosol_dist_rn_in(l,i)",
                "        ri(l,i)=aerosol_dist_ri_in(l,i)",
                "       enddo",
                "      enddo",
                "      go to 49",
            ]
        ),
        "junge aerosol input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "   46 read(5,*)irsunph",
                "      do i=1,irsunph",
                "       read(5,*)rsunph(i),nrsunph(i)",
                "C       nrsunph(i)=nrsunph(i)/(rsunph(i)**4.)/(4*3.1415/3)",
                "      enddo",
                "      rmin=rsunph(1)",
                "      rmax=rsunph(irsunph)+1e-07",
                "      read(5,*)(rn(l,1),l=1,20)",
                "      read(5,*)(ri(l,1),l=1,20)",
                "      go to 49",
            ]
        ),
        "\n".join(
            [
                "   46 irsunph=max(1,min(aerosol_sun_count_in,50))",
                "      do i=1,irsunph",
                "       rsunph(i)=aerosol_sun_radius_in(i)",
                "       nrsunph(i)=aerosol_sun_dvlogr_in(i)",
                "      enddo",
                "      rmin=rsunph(1)",
                "      rmax=rsunph(irsunph)+1e-07",
                "      do l=1,20",
                "       rn(l,1)=aerosol_dist_rn_in(l,1)",
                "       ri(l,1)=aerosol_dist_ri_in(l,1)",
                "      enddo",
                "      go to 49",
            ]
        ),
        "sun photometer aerosol input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "   47 read(5,'(A80)')FILE2",
                "      i2=index(FILE2,' ')-1",
            ]
        ),
        "\n".join(
            [
                "   47 FILE2=aerosol_filename_in",
                "      i2=len_trim(FILE2)",
            ]
        ),
        "user aerosol model filename input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      if (iaer.ge.8.and.iaer.le.11)then",
                "       read(5,*)iaerp",
                "       if (iaerp.eq.1)read(5,'(A80)')FILE",
                "       i1=index(FILE,' ')-1",
                "       FILE2=FILE(1:I1)//'.mie'",
                "       i2=index(FILE2,' ')-1",
                "      endif",
            ]
        ),
        "\n".join(
            [
                "      if (iaer.ge.8.and.iaer.le.11)then",
                "       iaerp=0",
                "      endif",
            ]
        ),
        "skip interactive aerosol mie-save input",
    )
    text = _replace_once(
        text,
        "      read(iread,*) v",
        "\n".join(
            [
                "      if (iaer.eq.0) then",
                "        v=-1.",
                "      else",
                "        v=0.",
                "      endif",
            ]
        ),
        "aerosol optical thickness mode input",
    )
    text = _replace_once(
        text,
        "   10 read(iread,*) taer55",
        "   10 taer55=max(taer55_in,0.)",
        "aerosol optical thickness input",
    )
    text = _replace_once(
        text,
        " 771   read(iread,*) xps",
        " 771   xps=-max(elevation_km_in,0.)",
        "target elevation input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "         if (idatm.ne.8) then",
                "         call pressure(uw,uo3,xps)",
                "        else",
                "         call pressure(uwus,uo3us,xps)",
                "        endif",
            ]
        ),
        "\n".join(
            [
                "         if (idatm.ne.8) then",
                "         call pressure(uw,uo3,xps)",
                "         idatm_abstra=idatm",
                "         if (atmospheric_columns_mode_in.ne.0) then",
                "          uwus=max(uw,1.0e-06)",
                "          uo3us=max(uo3,1.0e-06)",
                "          uw=max(uw_in,0.)",
                "          uo3=max(uo3_in,0.)",
                "          idatm_abstra=8",
                "         endif",
                "        else",
                "         call pressure(uwus,uo3us,xps)",
                "         idatm_abstra=8",
                "        endif",
            ]
        ),
        "atmospheric column scaling mode",
    )
    text = _replace_once(
        text,
        "if (idatm.eq.8) then",
        "if (idatm_abstra.eq.8) then",
        "plane atmospheric column scaling mode",
    )
    text = _replace_once(
        text, "        read(iread,*) xpp", "        xpp=-1000.", "sensor altitude input"
    )
    text = _replace_once(
        text, "       s(l)=1.", "       s(l)=0.", "spectral response initialization"
    )
    text = _replace_once(text, "      read(iread,*) iwave", "      iwave=1", "spectral mode input")
    text = _replace_once(
        text,
        "  110 read(iread,*) wlinf,wlsup",
        "\n".join(["  110 wlinf=wlinf_in", "      wlsup=wlsup_in"]),
        "spectral limits input",
    )
    text = _replace_once(
        text,
        "      read(iread,*) (s(i),i=iinf,isup)",
        "\n".join(
            [
                "      do 1101 i=iinf,isup",
                "        s(i)=spectral_response_in(i)",
                " 1101 continue",
            ]
        ),
        "spectral response input",
    )
    text = _replace_once(
        text,
        "      read(iread,*) inhomo",
        "      inhomo=surface_inhomo_in",
        "surface homogeneity input",
    )
    text = _replace_once(
        text,
        "  30  read(iread,*) idirec",
        "  30  idirec=surface_idirec_in",
        "surface directional mode input",
    )
    text = _replace_once(
        text,
        " 25   read(iread,*) ibrdf",
        " 25   ibrdf=surface_brdf_model_in",
        "surface brdf model input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "        call abstra(idatm,wl,xmus,xmuv,uw/2.,uo3,uwus,uo3us,",
                "     a             idatmp,puw/2.,puo3,puwus,puo3us,",
                "     a      dtwava,dtozon,dtdica,dtoxyg,dtniox,dtmeth,dtmoca,",
                "     a      utwava,utozon,utdica,utoxyg,utniox,utmeth,utmoca,",
                "     a      attwava,ttozon,ttdica,ttoxyg,ttniox,ttmeth,ttmoca )",
                "        call abstra(idatm,wl,xmus,xmuv,uw,uo3,uwus,uo3us,",
                "     a             idatmp,puw,puo3,puwus,puo3us,",
                "     a      dtwava,dtozon,dtdica,dtoxyg,dtniox,dtmeth,dtmoca,",
                "     a      utwava,utozon,utdica,utoxyg,utniox,utmeth,utmoca,",
                "     a      ttwava,ttozon,ttdica,ttoxyg,ttniox,ttmeth,ttmoca )",
            ]
        ),
        "\n".join(
            [
                "        call abstra(idatm_abstra,wl,xmus,xmuv,uw/2.,uo3,uwus,uo3us,",
                "     a             idatmp,puw/2.,puo3,puwus,puo3us,",
                "     a      dtwava,dtozon,dtdica,dtoxyg,dtniox,dtmeth,dtmoca,",
                "     a      utwava,utozon,utdica,utoxyg,utniox,utmeth,utmoca,",
                "     a      attwava,ttozon,ttdica,ttoxyg,ttniox,ttmeth,ttmoca )",
                "        call abstra(idatm_abstra,wl,xmus,xmuv,uw,uo3,uwus,uo3us,",
                "     a             idatmp,puw,puo3,puwus,puo3us,",
                "     a      dtwava,dtozon,dtdica,dtoxyg,dtniox,dtmeth,dtmoca,",
                "     a      utwava,utozon,utdica,utoxyg,utniox,utmeth,utmoca,",
                "     a      ttwava,ttozon,ttdica,ttoxyg,ttniox,ttmeth,ttmoca )",
            ]
        ),
        "abstra atmospheric mode override",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "  23  do 900 k=1,13",
                "        read(iread,*) (brdfdats(10-j+1,k),j=1,10)",
                "  900 continue",
                "      do 901 k=1,13",
                "        read(iread,*) (brdfdatv(10-j+1,k),j=1,10)",
                "  901 continue",
                "      read(iread,*) albbrdf",
                "      read(iread,*) rodir",
            ]
        ),
        "\n".join(
            [
                "  23  do 900 k=1,13",
                "        do 901 j=1,10",
                "          brdfdats(j,k)=surface_brdf_table_solar_in(j,k)",
                "  901   continue",
                "  900 continue",
                "      do 902 k=1,13",
                "        do 903 j=1,10",
                "          brdfdatv(j,k)=surface_brdf_table_view_in(j,k)",
                "  903   continue",
                "  902 continue",
                "      albbrdf=surface_brdf_spherical_albedo_in",
                "      rodir=surface_brdf_directional_reflectance_in",
            ]
        ),
        "user-defined brdf tables input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) par1,par2,par3,par4",
        "\n".join(
            [
                "        par1=surface_brdf_params_in(1)",
                "        par2=surface_brdf_params_in(2)",
                "        par3=surface_brdf_params_in(3)",
                "        par4=surface_brdf_params_in(4)",
            ]
        ),
        "hapke and walthall parameters input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) par1,par2,par3,par4",
        "\n".join(
            [
                "        par1=surface_brdf_params_in(1)",
                "        par2=surface_brdf_params_in(2)",
                "        par3=surface_brdf_params_in(3)",
                "        par4=surface_brdf_params_in(4)",
            ]
        ),
        "walthall parameters input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "        read(iread,*) (options(i),i=3,5)",
                "        options(1)=1",
                "        options(2)=1",
                "        read(iread,*) (struct(i),i=1,4)",
                "        read(iread,*) (optics(i),i=1,3)",
            ]
        ),
        "\n".join(
            [
                "        options(1)=1",
                "        options(2)=1",
                "        do i=3,5",
                "          options(i)=surface_brdf_options_in(i)",
                "        enddo",
                "        do i=1,4",
                "          struct(i)=surface_brdf_struct_in(i)",
                "        enddo",
                "        do i=1,3",
                "          optics(i)=surface_brdf_optics_in(i)",
                "        enddo",
            ]
        ),
        "verstraete parameters input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "        read(iread,*) pild,pihs",
                "        read(iread,*) pxLt,pc",
                "        read(iread,*) pRl,pTl,pRs",
            ]
        ),
        "\n".join(
            [
                "        pild=surface_brdf_options_in(1)",
                "        pihs=surface_brdf_options_in(2)",
                "        pxLt=surface_brdf_params_in(1)",
                "        pc=surface_brdf_params_in(2)",
                "        pRl=surface_brdf_params_in(3)",
                "        pTl=surface_brdf_params_in(4)",
                "        pRs=surface_brdf_params_in(5)",
            ]
        ),
        "iaquinta-pinty input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) par1,par2,par3",
        "\n".join(
            [
                "        par1=surface_brdf_params_in(1)",
                "        par2=surface_brdf_params_in(2)",
                "        par3=surface_brdf_params_in(3)",
            ]
        ),
        "roujean parameters input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) par1,par2,par3",
        "\n".join(
            [
                "        par1=surface_brdf_params_in(1)",
                "        par2=surface_brdf_params_in(2)",
                "        par3=surface_brdf_params_in(3)",
            ]
        ),
        "rahman parameters input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) pws,phi_wind,xsal,pcl",
        "\n".join(
            [
                "        pws=surface_brdf_params_in(1)",
                "        phi_wind=surface_brdf_params_in(2)",
                "        xsal=surface_brdf_params_in(3)",
                "        pcl=surface_brdf_params_in(4)",
            ]
        ),
        "ocean brdf input",
    )
    text = _replace_once(
        text,
        "        read(iread,*) par1,par2",
        "\n".join(
            [
                "        par1=surface_brdf_params_in(1)",
                "        par2=surface_brdf_params_in(2)",
            ]
        ),
        "minnaert parameters input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "         read(iread,*) uli,eei,thmi,sli",
                "         read(iread,*) cabi,cwi,vaii,rnci,rsl1i",
            ]
        ),
        "\n".join(
            [
                "         uli=surface_brdf_params_in(1)",
                "         eei=surface_brdf_params_in(2)",
                "         thmi=surface_brdf_params_in(3)",
                "         sli=surface_brdf_params_in(4)",
                "         cabi=surface_brdf_params_in(5)",
                "         cwi=surface_brdf_params_in(6)",
                "         vaii=surface_brdf_params_in(7)",
                "         rnci=surface_brdf_params_in(8)",
                "         rsl1i=surface_brdf_params_in(9)",
            ]
        ),
        "kuusk input",
    )
    text = _replace_once(
        text,
        "         read(iread,*)p1,p2,p3",
        "\n".join(
            [
                "         p1=surface_brdf_params_in(1)",
                "         p2=surface_brdf_params_in(2)",
                "         p3=surface_brdf_params_in(3)",
            ]
        ),
        "modis brdf input",
    )
    text = _replace_once(
        text,
        "         read(iread,*)p1,p2,p3",
        "\n".join(
            [
                "         p1=surface_brdf_params_in(1)",
                "         p2=surface_brdf_params_in(2)",
                "         p3=surface_brdf_params_in(3)",
            ]
        ),
        "ross-li-maignan brdf input",
    )
    text = _replace_once(
        text,
        "  21  read(iread,*) igroun",
        "  21  igroun=surface_target_mode_in",
        "surface group input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "  29  read(iread,*) nwlinf,nwlsup",
                "      niinf=(nwlinf-.25)/0.0025+1.5",
                "      nisup=(nwlsup-.25)/0.0025+1.5",
                "      read(iread,*) (rocl(i),i=niinf,nisup)",
            ]
        ),
        "\n".join(
            [
                "  29  do i=1,1501",
                "        rocl(i)=surface_target_spectrum_in(i)",
                "      enddo",
            ]
        ),
        "target spectral reflectance input",
    )
    text = _replace_once(
        text,
        "  32  read(iread,*) ro",
        "  32  ro=surface_target_constant_in",
        "target constant reflectance input",
    )
    text = _replace_once(
        text,
        " 31   read(iread,*) igrou1,igrou2,rad",
        "\n".join(
            [
                " 31   igrou1=surface_target_mode_in",
                "      igrou2=surface_env_mode_in",
                "      rad=surface_radius_km_in",
            ]
        ),
        "heterogeneous surface input",
    )
    text = _replace_once(
        text,
        "  59  read(iread,*) (rocl(i),i=iinf,isup)",
        "\n".join(
            [
                "  59  do i=1,1501",
                "        rocl(i)=surface_target_spectrum_in(i)",
                "      enddo",
            ]
        ),
        "heterogeneous target spectral reflectance input",
    )
    text = _replace_once(
        text,
        "  60  read(iread,*) roc",
        "  60  roc=surface_target_constant_in",
        "heterogeneous target constant reflectance input",
    )
    text = _replace_once(
        text,
        "  66  read(iread,*) (roel(i),i=iinf,isup)",
        "\n".join(
            [
                "  66  do i=1,1501",
                "        roel(i)=surface_env_spectrum_in(i)",
                "      enddo",
            ]
        ),
        "heterogeneous environment spectral reflectance input",
    )
    text = _replace_once(
        text,
        "  62  read(iread,*) roe",
        "  62  roe=surface_env_constant_in",
        "heterogeneous environment constant reflectance input",
    )
    text = _replace_once(
        text,
        "      read(iread,*) irapp",
        "      irapp=irapp_in",
        "atmospheric correction mode input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      if (irapp.ge.0) then",
                "         irapp=1",
                "         read(iread,*) rapp",
                "         endif",
            ]
        ),
        "\n".join(
            [
                "      if (irapp.ge.0) then",
                "         rapp=rapp_in",
                "         if (abs(rapp).lt.1.e-12) rapp=-reference_reflectance_in",
                "         irapp=1",
                "      endif",
            ]
        ),
        "atmospheric correction value input",
    )
    text = _replace_once(
        text,
        "       read(iread,*,end=37) irop",
        "\n".join(["       irop=0", "       goto 37"]),
        "surface polarization input",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      if (ipol.eq.1)then",
                "        rpfet=sqrt(rqfet*rqfet+rufet*rufet)",
                "\tplumet=sqrt(qlumet*qlumet+ulumet*ulumet)",
                "\txpol=atan2(rufet,rqfet)*180.0/3.14159/2.",
                "        write(iwr, 429 )rpfet,plumet,xpol,rpfet/refet",
                "C       write(iwr, 428 )rpfet1,rpfet2,rpfet3",
                "      endif",
            ]
        ),
        "\n".join(
            [
                "      if (ipol.eq.1)then",
                "        rpfet=sqrt(rqfet*rqfet+rufet*rufet)",
                "\tplumet=sqrt(qlumet*qlumet+ulumet*ulumet)",
                "\txpol=atan2(rufet,rqfet)*180.0/3.14159/2.",
                "        rpfet_report=rpfet",
                "        plumet_report=plumet",
                "        xpol_report=xpol",
                "        rpfet_over_refet_report=rpfet/refet",
                "        write(iwr, 429 )rpfet,plumet,xpol,rpfet/refet",
                "C       write(iwr, 428 )rpfet1,rpfet2,rpfet3",
                "      endif",
            ]
        ),
        "polarized report capture",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "      rpfet=0.",
                "      rpfet1=0.",
                "      rpfet2=0.",
                "      rpfet3=0.",
                "      alumet=0.",
                "      plumet=0.",
                "      tgasm=0.",
            ]
        ),
        "\n".join(
            [
                "      rpfet=0.",
                "      rpfet1=0.",
                "      rpfet2=0.",
                "      rpfet3=0.",
                "      alumet=0.",
                "      plumet=0.",
                "      rqfet=0.",
                "      rufet=0.",
                "      qlumet=0.",
                "      ulumet=0.",
                "      tgasm=0.",
            ]
        ),
        "polarization accumulator initialization",
    )
    text = _replace_once(
        text,
        "\n".join(
            [
                "        rocave=rocave/seb",
                "\trfoamave=rfoamave/seb",
                "\trwatave=rwatave/seb",
                "\trglitave=rglitave/seb",
                "\t",
                "         goto(2000,2001,2002,2003,2004,2005,2006,2007,2008,2010,2011,2012)",
            ]
        ),
        "\n".join(
            [
                "        rocave=rocave/seb",
                "\trfoamave=rfoamave/seb",
                "\trwatave=rwatave/seb",
                "\trglitave=rglitave/seb",
                "        rocave_report=rocave",
                "        robar1_over_xnorm1_report=robar1/xnorm1",
                "        robar2_over_xnorm2_report=robar2/xnorm2",
                "        rbard_report=rbard",
                "        albbrdf_report=albbrdf",
                "        rfoamave_report=rfoamave",
                "        rwatave_report=rwatave",
                "        rglitave_report=rglitave",
                "\t",
                "         goto(2000,2001,2002,2003,2004,2005,2006,2007,2008,2010,2011,2012)",
            ]
        ),
        "brdf surface report capture",
    )
    text = text.replace(
        "      if(ier) stop",
        "\n".join(
            [
                "      if(ier) then",
                "       ier_out=.TRUE.",
                "       return",
                "      endif",
            ]
        ),
    )
    text = re.sub(r"\bwrite\s*\(\s*6\s*,", "write(iwr,", text, flags=re.IGNORECASE)
    text = re.sub(r"\bprint\s*\*", "write(iwr,*)", text, flags=re.IGNORECASE)

    tail_old = "\n      stop\n"
    tail_new = "\n".join(["", *_core_output_assignment_lines(), "      return", ""])
    if tail_old not in text:
        raise RuntimeError("Failed to patch main.f final return block.")
    head, _, tail = text.rpartition(tail_old)
    return f"{head}{tail_new}{tail}"


def parse_makefile_sources(makefile_path: Path) -> list[Path]:
    """Return the fixed-form source files used by the upstream 6S build."""
    lines = makefile_path.read_text(encoding="utf-8").splitlines()
    collecting = False
    object_names: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not collecting and stripped.startswith("OBJECTS0"):
            collecting = True
            fragment = stripped.split("=", 1)[1].strip()
        elif collecting:
            fragment = stripped
        else:
            continue

        has_continuation = fragment.endswith("\\")
        fragment = fragment[:-1].strip() if has_continuation else fragment
        object_names.extend(token for token in fragment.split() if token.endswith(".o"))
        if not has_continuation:
            break

    if not object_names:
        raise RuntimeError(f"Could not parse OBJECTS0 from {makefile_path}.")

    sources = [makefile_path.parent / name.replace(".o", ".f") for name in object_names]
    for path in sources:
        if not path.exists():
            raise RuntimeError(f"Expected upstream 6S source {path} was not found.")
    return sources


def ensure_native_sixs_module(config: SixSAlgorithmConfig) -> Path:
    """Ensure the native 6SV2.1 Python extension exists and return its path."""
    paths = resolve_build_paths(config)
    existing = find_built_extension(paths)
    if existing is not None:
        return existing
    if not config.auto_build:
        raise RuntimeError(
            "6S native Python extension is not available. "
            f"Expected a built module under {paths.root_dir} and auto_build is disabled."
        )
    return build_native_sixs_module(config)


def build_native_sixs_module(config: SixSAlgorithmConfig) -> Path:
    """Fetch, patch, and compile the 6SV2.1 Python extension."""
    paths = resolve_build_paths(config)
    paths.root_dir.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which(config.compiler) or os.fspath(Path(config.compiler).expanduser())
    if shutil.which(config.compiler) is None and not Path(compiler).exists():
        raise RuntimeError(
            "6S native build requires a Fortran compiler. "
            f"Could not resolve {config.compiler!r}. "
            "Use the Pixi `rt6s` environment or set `algorithms.rt.sixs.compiler`."
        )

    source_dir = _prepare_source_tree(config, paths)
    _compile_f2py_extension(
        source_dir=source_dir,
        build_paths=paths,
        compiler=compiler,
        build_profile=str(getattr(config, "build_profile", "release")),
    )
    built_module = find_built_extension(paths)
    if built_module is None:
        raise RuntimeError("F2PY completed but the compiled 6S Python extension was not found.")
    return built_module


def find_built_extension(paths: SixSBuildPaths) -> Path | None:
    """Return the compiled extension path if it exists."""
    return _find_built_extension(paths)


def _find_built_extension(
    paths: SixSBuildPaths,
    *,
    extra_roots: tuple[Path, ...] = (),
    min_mtime: float | None = None,
) -> Path | None:
    """Return the newest compiled extension path if it exists."""
    if (
        paths.module_hint_path is not None
        and paths.module_hint_path.exists()
        and (min_mtime is None or paths.module_hint_path.stat().st_mtime >= min_mtime)
    ):
        return paths.module_hint_path
    roots = (paths.root_dir, paths.f2py_build_dir, *extra_roots)
    candidates: dict[Path, float] = {}
    for root in roots:
        if not root.exists():
            continue
        resolved_root = root.resolve()
        for suffix in native_module_suffixes():
            matches = {
                *resolved_root.glob(f"{paths.module_name}*{suffix}"),
                *resolved_root.rglob(f"{paths.module_name}*{suffix}"),
            }
            for candidate in matches:
                stat = candidate.stat()
                if min_mtime is not None and stat.st_mtime < min_mtime:
                    continue
                candidates[candidate] = stat.st_mtime
    if candidates:
        return max(candidates, key=candidates.get)
    return None


def _prepare_source_tree(config: SixSAlgorithmConfig, paths: SixSBuildPaths) -> Path:
    source_dir = Path(config.source_dir).expanduser() if config.source_dir else None
    if source_dir is None:
        source_dir = _fetch_and_unpack_source(
            config.source_url, paths.archive_path, paths.upstream_dir
        )
    else:
        if not source_dir.exists():
            raise RuntimeError(f"Configured 6S source_dir does not exist: {source_dir}")
        if not (source_dir / "main.f").exists():
            raise RuntimeError(
                f"Configured 6S source_dir does not look like a 6S source tree: {source_dir}"
            )

    if paths.patched_dir.exists():
        shutil.rmtree(paths.patched_dir)
    shutil.copytree(
        source_dir,
        paths.patched_dir,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("*.o", "*.a", "*.so", "*.dylib", "__pycache__"),
    )

    for source_path in paths.patched_dir.glob("*.f"):
        text = source_path.read_text(encoding="utf-8")
        if source_path.name == "DISCOM.f":
            text = patch_discom_source(text)
        if source_path.name == "AEROSO.f":
            text = patch_aeroso_source(text)
        if source_path.name == "KERNELPOL.f":
            text = patch_kernelpol_source(text)
        if source_path.name == "MIE.f":
            text = patch_mie_source(text)
        if source_path.name == "main.f":
            text = patch_main_source(text)
        text = patch_threadprivate_directives(text)
        _write_text(source_path, text)

    bridge_source = _FORTRAN_DIR / _BRIDGE_SOURCE_NAME
    if not bridge_source.exists():
        raise RuntimeError(f"Missing checked-in 6S bridge source: {bridge_source}")
    shutil.copy2(bridge_source, paths.patched_dir / _BRIDGE_SOURCE_NAME)
    return paths.patched_dir


#: SHA-256 of the upstream 6SV2.1 source tarball that SIAC's bridge expects.
#: Compared against the locally-cached archive after download. If the
#: upstream source is ever re-published (or this constant drifts from
#: the actual archive), update it here in one place rather than at
#: every callsite. ``None`` disables the check (set
#: ``SIAC_SIXS_SOURCE_SHA256`` env var to override per-run).
_SIXS_SOURCE_SHA256: str | None = None

#: Hard cap for the streamed source download (the real tarball is ~3 MB).
_MAX_SOURCE_DOWNLOAD_BYTES = 1 * 1024**3


def _archive_sha256(path: Path) -> str:
    """SHA-256 of *path* in fixed-size chunks (avoid loading the whole tarball)."""
    return archive_sha256(path)


def _fetch_and_unpack_source(source_url: str, archive_path: Path, upstream_dir: Path) -> Path:
    """Fetch the 6SV2.1 source tarball and unpack it under ``upstream_dir``.

    REVIEW.md §3.4 sixs_build.py:1576-1597 flagged two robustness issues
    in the original implementation: the download was not atomic (a
    partially-written tarball could survive a crash and be re-used on
    the next run), and there was no integrity check against a known
    hash. This version fixes both:

    - Download to ``<archive>.tmp`` and ``Path.replace`` on success so
      a crash mid-transfer leaves no usable file.
    - Verify the SHA-256 of the cached archive against
      ``_SIXS_SOURCE_SHA256`` (or ``SIAC_SIXS_SOURCE_SHA256`` env
      override). When the constant is ``None``, the check is skipped
      with a one-line warning.
    """
    fetch_archive(
        source_url,
        archive_path,
        expected_sha256=_SIXS_SOURCE_SHA256,
        sha_env_var="SIAC_SIXS_SOURCE_SHA256",
        max_bytes=_MAX_SOURCE_DOWNLOAD_BYTES,
        timeout_s=120.0,
        what="6SV2.1",
    )

    if upstream_dir.exists() and (upstream_dir / "main.f").exists():
        return upstream_dir

    if upstream_dir.exists():
        shutil.rmtree(upstream_dir)
    upstream_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path) as archive:
        _safe_extract_tar(archive, upstream_dir)

    if (upstream_dir / "main.f").exists():
        return upstream_dir

    candidates = [
        path for path in upstream_dir.iterdir() if path.is_dir() and (path / "main.f").exists()
    ]
    if len(candidates) == 1:
        return candidates[0]
    raise RuntimeError(f"Unable to locate unpacked 6S source tree under {upstream_dir}.")


def _safe_extract_tar(archive: tarfile.TarFile, dest: Path) -> None:
    """Extract a tar archive only if every member stays under ``dest``."""
    root = dest.resolve()
    for member in archive.getmembers():
        target = (root / member.name).resolve()
        if target != root and root not in target.parents:
            raise RuntimeError(
                f"Refusing to extract unsafe 6S source archive member: {member.name!r}"
            )
        if member.islnk() or member.issym():
            raise RuntimeError(
                f"Refusing to extract linked 6S source archive member: {member.name!r}"
            )
    archive.extractall(path=root)


def _compile_f2py_extension(
    *,
    source_dir: Path,
    build_paths: SixSBuildPaths,
    compiler: str,
    build_profile: str,
) -> None:
    source_dir = source_dir.resolve()
    build_root = build_paths.root_dir.resolve()
    f2py_build_dir = build_paths.f2py_build_dir.resolve()
    if f2py_build_dir.exists():
        shutil.rmtree(f2py_build_dir)
    f2py_build_dir.mkdir(parents=True, exist_ok=True)

    existing = find_built_extension(build_paths)
    if existing is not None:
        existing.unlink()

    source_files = parse_makefile_sources(source_dir / "Makefile")
    compile_sources = list(
        dict.fromkeys([*source_files, source_dir / "main.f", source_dir / _BRIDGE_SOURCE_NAME])
    )
    flags = _compiler_flags(build_profile)
    signature_path = _generate_f2py_signature(
        source_dir=source_dir,
        module_name=build_paths.module_name,
        f2py_build_dir=f2py_build_dir,
    )

    env = os.environ.copy()
    env["FC"] = compiler
    env["F77"] = compiler
    build_started_at = time.time()
    caller_cwd = Path.cwd().resolve()
    extension_search_roots = (caller_cwd, *_environment_extension_roots())
    failures: list[tuple[str, list[str], subprocess.CompletedProcess[str], str]] = []
    for backend in _resolve_f2py_backends():
        if backend == "distutils":
            logger.info("Compiling native 6S Python extension with distutils backend.")
            completed = _run_distutils_backend(
                module_name=build_paths.module_name,
                source_dir=source_dir,
                signature_path=signature_path,
                compile_sources=compile_sources,
                flags=flags,
                f2py_build_dir=f2py_build_dir,
                build_root=build_root,
                compiler=compiler,
            )
            cmd = list(completed.args)
        else:
            cmd = _build_f2py_command(
                backend=backend,
                module_name=build_paths.module_name,
                source_dir=source_dir,
                signature_path=signature_path,
                compile_sources=compile_sources,
                flags=flags,
                f2py_build_dir=f2py_build_dir,
            )
            logger.info(
                "Compiling native 6S Python extension with %s backend: %s", backend, " ".join(cmd)
            )
            completed = subprocess.run(
                cmd,
                cwd=build_root,
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

        built_extension = _find_built_extension(
            build_paths,
            extra_roots=extension_search_roots,
            min_mtime=build_started_at,
        )
        attempt_status = "success-with-extension"
        if built_extension is None and completed.returncode != 0:
            attempt_status = "failed"
        elif built_extension is None:
            attempt_status = "missing-extension"
        _write_backend_diagnostics(
            build_paths=build_paths,
            backend=backend,
            cmd=cmd,
            completed=completed,
            status=attempt_status,
            built_extension=built_extension,
        )
        if built_extension is not None:
            validation = _validate_extension_import(build_paths.module_name, built_extension)
            if validation.returncode != 0:
                _write_backend_diagnostics(
                    build_paths=build_paths,
                    backend=backend,
                    cmd=cmd,
                    completed=completed,
                    status="import-failed",
                    built_extension=built_extension,
                    validation=validation,
                )
                with contextlib.suppress(OSError):
                    built_extension.unlink()
                log_summary = _summarize_build_output(validation.stdout, validation.stderr)
                failures.append(
                    (
                        backend,
                        cmd,
                        validation,
                        f"built extension failed import validation: {log_summary}",
                    )
                )
                logger.warning(
                    "F2PY %s backend built extension %s, but import validation failed: %s",
                    backend,
                    built_extension,
                    log_summary,
                )
                continue
            if completed.returncode != 0:
                logger.warning(
                    "F2PY %s backend reported exit code %s but built extension %s was found; "
                    "continuing with the produced module. Summary: %s",
                    backend,
                    completed.returncode,
                    built_extension,
                    _summarize_build_output(completed.stdout, completed.stderr),
                )
            if (
                build_paths.module_hint_path is not None
                and built_extension != build_paths.module_hint_path
            ):
                build_paths.module_hint_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(os.fspath(built_extension), os.fspath(build_paths.module_hint_path))
            return

        if completed.returncode != 0:
            log_summary = _summarize_build_output(completed.stdout, completed.stderr)
            failures.append((backend, cmd, completed, log_summary))
            logger.warning(
                "F2PY %s backend failed: %s",
                backend,
                log_summary,
            )
            continue

        log_summary = _summarize_build_output(completed.stdout, completed.stderr)
        failures.append((backend, cmd, completed, log_summary))
        logger.warning(
            "F2PY %s backend reported success but no extension was found: %s",
            backend,
            log_summary,
        )
        continue

    details = []
    for backend, cmd, completed, log_summary in failures:
        details.append(
            f"Backend {backend}: {log_summary}\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    _write_build_failure_summary(build_paths, details)
    raise RuntimeError("6S native Python extension build failed.\n" + "\n\n".join(details))


def _build_f2py_command(
    *,
    backend: str,
    module_name: str,
    source_dir: Path,
    signature_path: Path,
    compile_sources: list[Path],
    flags: list[str],
    f2py_build_dir: Path,
) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-c",
        "--backend",
        backend,
    ]
    if backend == "meson":
        backend_build_dir = f2py_build_dir / backend
        if backend_build_dir.exists():
            shutil.rmtree(backend_build_dir)
        backend_build_dir.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--build-dir", os.fspath(backend_build_dir)])
    cmd.extend(
        [
            "-m",
            module_name,
            f"--f77flags={' '.join(flags)}",
            f"--f90flags={' '.join(flags)}",
            f"-I{os.fspath(source_dir)}",
        ]
    )
    cmd.append(os.fspath(signature_path.resolve()))
    cmd.extend(os.fspath(path.resolve()) for path in compile_sources)
    return cmd


def _validate_extension_import(
    module_name: str, module_path: Path
) -> subprocess.CompletedProcess[str]:
    validation_code = """
import ctypes
import importlib.machinery
import importlib.util
import os
import pathlib
import sys

module_name = sys.argv[1]
module_path = pathlib.Path(sys.argv[2]).resolve()
runtime_candidates = []
env_runtime = os.getenv("SIAC_SIXS_OPENMP_RUNTIME")
if env_runtime:
    runtime_candidates.append(pathlib.Path(env_runtime).expanduser())
for name in (
    "libgomp.dylib",
    "libgomp.1.dylib",
    "libgomp.so.1",
    "libgomp.so",
    "libomp.dylib",
    "libomp.so",
):
    runtime_candidates.append(pathlib.Path(sys.prefix) / "lib" / name)
for runtime in runtime_candidates:
    if runtime.exists():
        try:
            ctypes.CDLL(str(runtime), mode=ctypes.RTLD_GLOBAL)
            break
        except OSError:
            pass
loader = importlib.machinery.ExtensionFileLoader(module_name, str(module_path))
spec = importlib.util.spec_from_file_location(module_name, module_path, loader=loader)
if spec is None or spec.loader is None:
    raise ImportError(f"Unable to create import spec for {module_path}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
print(module_path)
"""
    return subprocess.run(
        [sys.executable, "-c", validation_code, module_name, os.fspath(module_path.resolve())],
        capture_output=True,
        text=True,
        check=False,
    )


def _distutils_extra_link_args(flags: list[str]) -> list[str]:
    if "-fopenmp" not in flags:
        return []
    return ["-lgomp", "-lquadmath", "-lm"]


def _run_distutils_backend(
    *,
    module_name: str,
    source_dir: Path,
    signature_path: Path,
    compile_sources: list[Path],
    flags: list[str],
    f2py_build_dir: Path,
    build_root: Path,
    compiler: str,
) -> subprocess.CompletedProcess[str]:
    from numpy.distutils.core import Extension, setup

    backend_build_dir = f2py_build_dir / "distutils"
    if backend_build_dir.exists():
        shutil.rmtree(backend_build_dir)
    backend_build_dir.mkdir(parents=True, exist_ok=True)

    ext = Extension(
        name=module_name,
        sources=[
            os.fspath(signature_path.resolve()),
            *(os.fspath(path.resolve()) for path in compile_sources),
        ],
        include_dirs=[os.fspath(source_dir)],
        extra_link_args=_distutils_extra_link_args(flags),
    )

    script_name = os.fspath(build_root / "setup.py")
    script_args = [
        "build",
        "--build-temp",
        os.fspath(backend_build_dir),
        "--build-base",
        os.fspath(backend_build_dir),
        "--build-platlib",
        os.fspath(build_root),
        "--disable-optimization",
        "config_fc",
        f"--f77flags={' '.join(flags)}",
        f"--f90flags={' '.join(flags)}",
    ]
    argv = [script_name, *script_args]

    stdout_buffer = _EncodedStringIO()
    stderr_buffer = _EncodedStringIO()
    old_cwd = Path.cwd()
    previous_fc = os.environ.get("FC")
    previous_f77 = os.environ.get("F77")
    os.environ["FC"] = compiler
    os.environ["F77"] = compiler

    try:
        os.chdir(build_root)
        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            setup(script_name=script_name, script_args=script_args, ext_modules=[ext])
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
        return subprocess.CompletedProcess(
            args=argv,
            returncode=code,
            stdout=stdout_buffer.getvalue(),
            stderr=stderr_buffer.getvalue(),
        )
    except Exception:
        traceback.print_exc(file=stderr_buffer)
        return subprocess.CompletedProcess(
            args=argv,
            returncode=1,
            stdout=stdout_buffer.getvalue(),
            stderr=stderr_buffer.getvalue(),
        )
    finally:
        os.chdir(old_cwd)
        if previous_fc is None:
            os.environ.pop("FC", None)
        else:
            os.environ["FC"] = previous_fc
        if previous_f77 is None:
            os.environ.pop("F77", None)
        else:
            os.environ["F77"] = previous_f77

    return subprocess.CompletedProcess(
        args=argv,
        returncode=0,
        stdout=stdout_buffer.getvalue(),
        stderr=stderr_buffer.getvalue(),
    )


def _generate_f2py_signature(
    *,
    source_dir: Path,
    module_name: str,
    f2py_build_dir: Path,
) -> Path:
    signature_path = (f2py_build_dir / _SIGNATURE_FILENAME).resolve()
    bridge_path = (source_dir / _BRIDGE_SOURCE_NAME).resolve()
    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-h",
        os.fspath(signature_path),
        "-m",
        module_name,
        os.fspath(bridge_path),
        "only:",
        "sixs_f2py_run_batch",
        ":",
        "--overwrite-signature",
    ]
    logger.info("Generating native 6S F2PY signature: %s", " ".join(cmd))
    completed = subprocess.run(
        cmd,
        cwd=f2py_build_dir.resolve(),
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0 or not signature_path.exists():
        details = _summarize_build_output(completed.stdout, completed.stderr)
        raise RuntimeError(
            "Failed to generate the 6S F2PY signature. "
            f"Command: {' '.join(cmd)}. Details: {details}"
        )
    return signature_path


def _diagnostics_dir(build_paths: SixSBuildPaths) -> Path:
    path = build_paths.root_dir / _DIAGNOSTICS_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def _environment_extension_roots() -> tuple[Path, ...]:
    roots: list[Path] = []
    seen: set[Path] = set()

    def _add(candidate: str | os.PathLike[str] | None) -> None:
        if not candidate:
            return
        path = Path(candidate).expanduser().resolve()
        if path in seen or not path.exists():
            return
        seen.add(path)
        roots.append(path)

    paths = sysconfig.get_paths()
    _add(paths.get("platlib"))
    _add(paths.get("purelib"))
    with contextlib.suppress(AttributeError):
        for candidate in site.getsitepackages():
            _add(candidate)
    with contextlib.suppress(AttributeError):
        _add(site.getusersitepackages())
    return tuple(roots)


def _write_backend_diagnostics(
    *,
    build_paths: SixSBuildPaths,
    backend: str,
    cmd: list[str],
    completed: subprocess.CompletedProcess[str],
    status: str,
    built_extension: Path | None,
    validation: subprocess.CompletedProcess[str] | None = None,
) -> None:
    diagnostics_dir = _diagnostics_dir(build_paths)
    prefix = f"f2py-{backend}"
    _write_text(diagnostics_dir / f"{prefix}.command.txt", shlex.join(cmd))
    _write_text(diagnostics_dir / f"{prefix}.stdout.txt", completed.stdout or "")
    _write_text(diagnostics_dir / f"{prefix}.stderr.txt", completed.stderr or "")
    summary_lines = [
        f"backend={backend}",
        f"status={status}",
        f"returncode={completed.returncode}",
        f"built_extension={built_extension if built_extension is not None else ''}",
        f"command={shlex.join(cmd)}",
        f"summary={_summarize_build_output(completed.stdout, completed.stderr)}",
    ]
    if validation is not None:
        _write_text(diagnostics_dir / f"{prefix}.import_check.stdout.txt", validation.stdout or "")
        _write_text(diagnostics_dir / f"{prefix}.import_check.stderr.txt", validation.stderr or "")
        summary_lines.extend(
            [
                f"import_check_returncode={validation.returncode}",
                f"import_check_summary={_summarize_build_output(validation.stdout, validation.stderr)}",
            ]
        )
    _write_text(diagnostics_dir / f"{prefix}.summary.txt", "\n".join(summary_lines) + "\n")


def _write_build_failure_summary(build_paths: SixSBuildPaths, details: list[str]) -> None:
    diagnostics_dir = _diagnostics_dir(build_paths)
    _write_text(
        diagnostics_dir / "build_failure_summary.txt",
        "6S native Python extension build failed.\n\n" + "\n\n".join(details),
    )


def _replace_once(text: str, old: str, new: str, description: str) -> str:
    if old not in text:
        raise RuntimeError(f"Failed to patch {description}.")
    return text.replace(old, new, 1)


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _summarize_build_output(stdout: str, stderr: str, *, max_chars: int = 1200) -> str:
    combined = (stderr.strip() or stdout.strip() or "no build output captured").replace("\r", "\n")
    lines = [line.strip() for line in combined.splitlines() if line.strip()]
    if not lines:
        return "no build output captured"
    summary = " | ".join(lines[-12:])
    if len(summary) > max_chars:
        return f"...{summary[-max_chars:]}"
    return summary

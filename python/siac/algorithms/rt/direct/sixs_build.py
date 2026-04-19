"""Fetch, patch, and build the native 6SV2.1 Python extension."""

from __future__ import annotations

import importlib.machinery
import importlib.metadata
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import requests

from siac.sixs_outputs import SIXS_CORE_OUTPUT_SPECS

if TYPE_CHECKING:
    from siac.config.schema import SixSAlgorithmConfig

logger = logging.getLogger(__name__)

_ARCHIVE_NAME = "6sV2.1.tar"
_UPSTREAM_DIRNAME = "upstream"
_PATCHED_DIRNAME = "patched"
_F2PY_BUILD_DIRNAME = "f2py_build"
_MODULE_BASENAME = "_siac_rt6s_native"
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


def _core_output_assignment_lines() -> list[str]:
    lines = ["      ier_out=ier"]
    for index, (_name, expression, condition) in enumerate(SIXS_CORE_OUTPUT_SPECS, start=1):
        target = f"output_values_out({index})"
        if condition is None:
            lines.append(f"      {target}={expression}")
        else:
            lines.append(f"      if ({condition}) {target}={expression}")
    return lines


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
    if configured_module_path is None and getattr(config, "library_path", None) is not None:
        candidate = Path(config.library_path).expanduser()
        if candidate.suffix in native_module_suffixes():
            configured_module_path = candidate
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
            f"setuptools {setuptools_version} is too new for numpy.distutils; "
            "use setuptools < 60"
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
        patched = _replace_once(patched, save_old, save_new, "AEROSO.f optional Mie-file output guard")
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
        "\n".join(
            [
                "      subroutine sixs_case_core(iwr_unit,asol_in,phi0_in,avis_in,",
                "     s phiv_in,month_in,jday_in,idatm_in,uw_in,uo3_in,",
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
        ),
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
        "\n".join(
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
        ),
        "main.f runtime initialization",
    )
    text = _replace_once(text, "      read(iread,*) igeom", "      igeom=0", "geometry mode input")
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
    text = _replace_once(text, "      read(iread,*) idatm", "      idatm=idatm_in", "atmosphere mode input")
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
    text = _replace_once(text, "      read(iread,*) iaer", "      iaer=iaer_in", "aerosol mode input")
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
    text = _replace_once(text, "   10 read(iread,*) taer55", "   10 taer55=max(taer55_in,0.)", "aerosol optical thickness input")
    text = _replace_once(text, " 771   read(iread,*) xps", " 771   xps=-max(elevation_km_in,0.)", "target elevation input")
    text = _replace_once(text, "        read(iread,*) xpp", "        xpp=-1000.", "sensor altitude input")
    text = _replace_once(text, "       s(l)=1.", "       s(l)=0.", "spectral response initialization")
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
    text = _replace_once(text, "      read(iread,*) inhomo", "      inhomo=surface_inhomo_in", "surface homogeneity input")
    text = _replace_once(text, "  30  read(iread,*) idirec", "  30  idirec=surface_idirec_in", "surface directional mode input")
    text = _replace_once(text, " 25   read(iread,*) ibrdf", " 25   ibrdf=surface_brdf_model_in", "surface brdf model input")
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
    text = _replace_once(text, "  21  read(iread,*) igroun", "  21  igroun=surface_target_mode_in", "surface group input")
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
    text = _replace_once(text, "  32  read(iread,*) ro", "  32  ro=surface_target_constant_in", "target constant reflectance input")
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
    text = _replace_once(text, "  60  read(iread,*) roc", "  60  roc=surface_target_constant_in", "heterogeneous target constant reflectance input")
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
    text = _replace_once(text, "  62  read(iread,*) roe", "  62  roe=surface_env_constant_in", "heterogeneous environment constant reflectance input")
    text = _replace_once(text, "      read(iread,*) irapp", "      irapp=irapp_in", "atmospheric correction mode input")
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


def ensure_native_sixs_library(config: SixSAlgorithmConfig) -> Path:
    """Backward-compatible alias for the compiled Python extension path."""
    return ensure_native_sixs_module(config)


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


def build_native_sixs_library(config: SixSAlgorithmConfig) -> Path:
    """Backward-compatible alias for the compiled Python extension path."""
    return build_native_sixs_module(config)


def find_built_extension(paths: SixSBuildPaths) -> Path | None:
    """Return the compiled extension path if it exists."""
    if paths.module_hint_path is not None and paths.module_hint_path.exists():
        return paths.module_hint_path
    for suffix in native_module_suffixes():
        candidates = sorted(
            {
                *paths.root_dir.glob(f"{paths.module_name}*{suffix}"),
                *paths.root_dir.rglob(f"{paths.module_name}*{suffix}"),
                *paths.f2py_build_dir.glob(f"{paths.module_name}*{suffix}"),
                *paths.f2py_build_dir.rglob(f"{paths.module_name}*{suffix}"),
            }
        )
        if candidates:
            return max(candidates, key=lambda path: path.stat().st_mtime)
    return None


def _prepare_source_tree(config: SixSAlgorithmConfig, paths: SixSBuildPaths) -> Path:
    source_dir = Path(config.source_dir).expanduser() if config.source_dir else None
    if source_dir is None:
        source_dir = _fetch_and_unpack_source(config.source_url, paths.archive_path, paths.upstream_dir)
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


def _fetch_and_unpack_source(source_url: str, archive_path: Path, upstream_dir: Path) -> Path:
    if not archive_path.exists():
        logger.info("Downloading 6SV2.1 source from %s", source_url)
        response = requests.get(source_url, timeout=120)
        response.raise_for_status()
        archive_path.write_bytes(response.content)

    if upstream_dir.exists() and (upstream_dir / "main.f").exists():
        return upstream_dir

    if upstream_dir.exists():
        shutil.rmtree(upstream_dir)
    upstream_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path) as archive:
        archive.extractall(path=upstream_dir)

    if (upstream_dir / "main.f").exists():
        return upstream_dir

    candidates = [path for path in upstream_dir.iterdir() if path.is_dir() and (path / "main.f").exists()]
    if len(candidates) == 1:
        return candidates[0]
    raise RuntimeError(f"Unable to locate unpacked 6S source tree under {upstream_dir}.")


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
    compile_sources = list(dict.fromkeys([*source_files, source_dir / "main.f", source_dir / _BRIDGE_SOURCE_NAME]))
    flags = _compiler_flags(build_profile)

    env = os.environ.copy()
    env["FC"] = compiler
    env["F77"] = compiler
    failures: list[tuple[str, list[str], subprocess.CompletedProcess[str], str]] = []
    for backend in _resolve_f2py_backends():
        cmd = _build_f2py_command(
            backend=backend,
            module_name=build_paths.module_name,
            source_dir=source_dir,
            compile_sources=compile_sources,
            flags=flags,
            f2py_build_dir=f2py_build_dir,
        )
        logger.info("Compiling native 6S Python extension with %s backend: %s", backend, " ".join(cmd))
        completed = subprocess.run(
            cmd,
            cwd=build_root,
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        if completed.returncode != 0:
            log_summary = _summarize_build_output(completed.stdout, completed.stderr)
            failures.append((backend, cmd, completed, log_summary))
            logger.warning("F2PY %s backend failed: %s", backend, log_summary)
            continue

        built_extension = find_built_extension(build_paths)
        if built_extension is None:
            log_summary = _summarize_build_output(completed.stdout, completed.stderr)
            failures.append((backend, cmd, completed, log_summary))
            logger.warning(
                "F2PY %s backend reported success but no extension was found: %s",
                backend,
                log_summary,
            )
            continue
        if build_paths.module_hint_path is not None and built_extension != build_paths.module_hint_path:
            build_paths.module_hint_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(os.fspath(built_extension), os.fspath(build_paths.module_hint_path))
        return

    details = []
    for backend, cmd, completed, log_summary in failures:
        details.append(
            f"Backend {backend}: {log_summary}\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    raise RuntimeError("6S native Python extension build failed.\n" + "\n\n".join(details))


def _build_f2py_command(
    *,
    backend: str,
    module_name: str,
    source_dir: Path,
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
    if backend == "distutils" and "-fopenmp" in flags:
        cmd.append("-lgomp")
    cmd.extend(os.fspath(path.resolve()) for path in compile_sources)
    cmd.extend(["only:", "sixs_f2py_run_batch", ":"])
    return cmd


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


__all__ = [
    "SixSBuildPaths",
    "build_native_sixs_library",
    "build_native_sixs_module",
    "default_build_root",
    "ensure_native_sixs_library",
    "ensure_native_sixs_module",
    "find_built_extension",
    "native_module_suffixes",
    "parse_makefile_sources",
    "patch_aeroso_source",
    "patch_discom_source",
    "patch_kernelpol_source",
    "patch_main_source",
    "patch_mie_source",
    "patch_threadprivate_directives",
    "resolve_build_paths",
    "shared_library_suffix",
]

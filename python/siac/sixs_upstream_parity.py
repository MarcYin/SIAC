"""Parity validation between the original 6SV2.1 executable and the native SIAC bridge."""

from __future__ import annotations

import base64
import json
import math
import os
import pickle
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs_build import build_native_sixs_module
from siac.algorithms.rt.direct.sixs_native import (
    SixSNativeRunner,
    _build_spectral_response,
    _resolve_aerosol_inputs,
    _resolve_atmospheric_correction_inputs,
    _resolve_atmospheric_inputs,
    _resolve_surface_inputs,
)
from siac.config import RTSetupConfig, SixSAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles
from siac.sixs_outputs import SIXS_DEFAULT_OUTPUT_VARIABLES

_FLOAT_TOKEN = r"(?:NaN|[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
_FLOAT_RE = re.compile(_FLOAT_TOKEN)


@dataclass(frozen=True)
class SixSParityCase:
    """Single original-vs-native comparison case."""

    name: str
    description: str
    sixs_config: SixSAlgorithmConfig
    rt_setup: RTSetupConfig
    band: SensorBand
    sza_deg: float
    saa_deg: float
    vza_deg: float
    vaa_deg: float
    aot550: float
    tcwv_cm: float
    tco3_atmcm: float
    elevation_km: float


@dataclass(frozen=True)
class ParsedMetric:
    """Single parsed numeric value from original 6S stdout."""

    value: float
    tolerance: float
    digits: int
    source: str


@dataclass(frozen=True)
class ComparedMetric:
    """Comparison record between original stdout and native output."""

    name: str
    original: float | None
    native: float | None
    tolerance: float
    abs_diff: float | None
    matched: bool
    source: str


@dataclass(frozen=True)
class CaseReport:
    """Full parity report for one case."""

    name: str
    description: str
    input_path: str
    stdout_path: str
    stderr_path: str
    compared_variable_count: int
    matched_variable_count: int
    mismatches: tuple[ComparedMetric, ...]
    compared: tuple[ComparedMetric, ...]

    @property
    def matched(self) -> bool:
        return len(self.mismatches) == 0


_TRIPLET_LINE_SPECS: tuple[tuple[str, tuple[str, str, str], int], ...] = (
    ("global gas. trans.", ("dgasm", "ugasm", "tgasm"), 5),
    ('water   "     "', ("sdwava", "suwava", "stwava"), 5),
    ('ozone   "     "', ("sdozon", "suozon", "stozon"), 5),
    ('co2     "     "', ("sddica", "sudica", "stdica"), 5),
    ('oxyg    "     "', ("sdoxyg", "suoxyg", "stoxyg"), 5),
    ('no2     "     "', ("sdniox", "suniox", "stniox"), 5),
    ('ch4     "     "', ("sdmeth", "sumeth", "stmeth"), 5),
    ('co      "     "', ("sdmoca", "sumoca", "stmoca"), 5),
    ("rayl.  sca. trans.", ("sdtotr", "sutotr", "sttotr"), 5),
    ('aeros. sca.   "   ', ("sdtota", "sutota", "sttota"), 5),
    ('total  sca.   "   ', ("sdtott", "sutott", "sttott"), 5),
    ("spherical albedo", ("sasr", "sasa", "sast"), 5),
    ("optical depth total", ("sodray", "sodaer", "sodtot"), 5),
    ("optical depth plane", ("sodrayp", "sodaerp", "sodtotp"), 5),
    ("reflectance I", ("sroray", "sroaer", "srotot"), 5),
    ("reflectance Q", ("srqray", "srqaer", "srqtot"), 5),
    ("reflectance U", ("sruray", "sruaer", "srutot"), 5),
    ("polarized reflect.", ("srpray", "srpaer", "srptot"), 5),
    ("degree of polar.", ("sdpray", "sdpaer", "sdptot"), 2),
    ("dir. plane polar.", ("sdppray", "sdppaer", "sdpptot"), 2),
    ("phase function I", ("fophsr", "fophsa", "fophst"), 5),
    ("phase function Q", ("foqhsr", "foqhsa", "foqhst"), 5),
    ("phase function U", ("fouhsr", "fouhsa", "fouhst"), 5),
    ("primary deg. of pol", ("spdpray", "spdpaer", "spdptot"), 5),
    ("sing. scat. albedo", ("pizerr", "pizera", "pizert"), 5),
)


def _print_step_for_digits(digits: int) -> float:
    return 10.0 ** (-digits)


def _parse_float(token: str) -> float:
    if token.strip().lower() == "nan":
        return math.nan
    return float(token)


def _extract_floats(line: str) -> list[float]:
    return [_parse_float(token) for token in _FLOAT_RE.findall(line)]


def _label_pattern(label: str) -> str:
    return re.escape(label).replace(r"\ ", r"\s+")


def _match_line(line: str, label: str, count: int) -> list[float] | None:
    pattern = re.compile(
        rf"^{_label_pattern(label)}\s*:\s*"
        + r"\s+".join([rf"({_FLOAT_TOKEN})" for _ in range(count)])
        + r"$",
        flags=re.IGNORECASE,
    )
    match = pattern.match(line)
    if match is None:
        return None
    return [_parse_float(item) for item in match.groups()]


def _store_metric(metrics: dict[str, ParsedMetric], name: str, value: float, digits: int, source: str) -> None:
    metric = ParsedMetric(
        value=value,
        tolerance=_print_step_for_digits(digits),
        digits=digits,
        source=source,
    )
    previous = metrics.get(name)
    if previous is None or metric.digits > previous.digits:
        metrics[name] = metric


def _clean_stdout_line(line: str) -> str:
    return line.replace("*", " ").strip()


def parse_original_sixs_stdout(stdout: str) -> dict[str, ParsedMetric]:
    """Parse the numeric report emitted by the unmodified 6SV2.1 executable."""
    metrics: dict[str, ParsedMetric] = {}
    expect_triplet: tuple[tuple[str, str, str], int] | None = None
    expect_pair: tuple[tuple[str, int], tuple[str, int]] | None = None
    expect_brdf = False

    for raw_line in stdout.splitlines():
        line = _clean_stdout_line(raw_line)
        if not line:
            continue

        if expect_triplet is not None:
            values = _extract_floats(line)
            if len(values) >= 3:
                for name, value in zip(expect_triplet[0], values[:3], strict=True):
                    _store_metric(metrics, name, value, expect_triplet[1], line)
                expect_triplet = None
                continue

        if expect_pair is not None:
            values = _extract_floats(line)
            if len(values) >= 2:
                for (name, digits), value in zip(expect_pair, values[:2], strict=True):
                    _store_metric(metrics, name, value, digits, line)
                expect_pair = None
                continue

        if expect_brdf:
            values = _extract_floats(line)
            if len(values) >= 5:
                names = ("rocave", "robar1_over_xnorm1", "robar2_over_xnorm2", "rbard", "albbrdf")
                for name, value in zip(names, values[:5], strict=True):
                    _store_metric(metrics, name, value, 5, line)
                expect_brdf = False
                continue

        if "rodir    robar    ropbar    robarbar  albedo" in line:
            expect_brdf = True
            continue
        if "% of direct  irr." in line:
            expect_triplet = (("aini_1_1", "aini_1_2", "aini_1_3"), 3)
            continue
        if "atm. intrin. ref." in line:
            expect_triplet = (("ainr_1_1", "ainr_1_2", "ainr_1_3"), 3)
            continue
        if "direct solar irr." in line:
            expect_triplet = (("aini_2_1", "aini_2_2", "aini_2_3"), 3)
            continue
        if "atm. intrin. rad." in line:
            expect_triplet = (("ainr_2_1", "ainr_2_2", "ainr_2_3"), 3)
            continue
        if "int. funct filter (in mic)" in line:
            expect_pair = (("sb", 7), ("seb", 3))
            continue

        match = re.search(
            rf"apparent reflectance\s+({_FLOAT_TOKEN})\s+appar\.\s+rad\.\(w/m2/sr/mic\)\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "refet", _parse_float(match.group(1)), 7, line)
            _store_metric(metrics, "alumet", _parse_float(match.group(2)), 3, line)
            continue

        match = re.search(
            rf"total gaseous transmittance\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "tgasm", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(
            rf"wv above aerosol\s*:\s*({_FLOAT_TOKEN}).*wv mixed with aerosol\s*:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "refet1", _parse_float(match.group(1)), 3, line)
            _store_metric(metrics, "refet2", _parse_float(match.group(2)), 3, line)
            continue

        match = re.search(
            rf"wv under aerosol\s*:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "refet3", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(
            rf"app\.\s+polarized refl\.\s+({_FLOAT_TOKEN}).*app\.\s+pol\.\s+rad\.\s+\(w/m2/sr/mic\)\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "rpfet", _parse_float(match.group(1)), 4, line)
            _store_metric(metrics, "plumet", _parse_float(match.group(2)), 3, line)
            continue

        match = re.search(
            rf"direction of the plane of polarization\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "xpol", _parse_float(match.group(1)), 2, line)
            continue

        match = re.search(
            rf"total polarization ratio\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "rpfet_over_refet", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(
            rf"Foam:\s*({_FLOAT_TOKEN})\s+Water:\s*({_FLOAT_TOKEN})\s+Glint:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "rfoamave", _parse_float(match.group(1)), 5, line)
            _store_metric(metrics, "rwatave", _parse_float(match.group(2)), 5, line)
            _store_metric(metrics, "rglitave", _parse_float(match.group(3)), 5, line)
            continue

        match = re.search(rf"roocean water\s+({_FLOAT_TOKEN})", line, flags=re.IGNORECASE)
        if match is not None:
            _store_metric(metrics, "rooceaw", _parse_float(match.group(1)), 7, line)
            continue

        match = re.search(
            rf"input apparent reflectance\s*:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "rapp", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(
            rf"measured radiance\s+\[w/m2/sr/mic\]\s*:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "xrad", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(rf"Lambertian case\s*:\s*({_FLOAT_TOKEN})", line, flags=re.IGNORECASE)
        if match is not None:
            _store_metric(metrics, "rog", _parse_float(match.group(1)), 5, line)
            continue

        match = re.search(
            rf"atmospherically corrected reflectance\s*:\s*({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "rog", _parse_float(match.group(1)), 3, line)
            continue

        match = re.search(rf"BRDF\s+case\s*:\s*({_FLOAT_TOKEN})", line, flags=re.IGNORECASE)
        if match is not None:
            _store_metric(metrics, "rogbrdf", _parse_float(match.group(1)), 5, line)
            continue

        match = re.search(
            rf"coefficients xa xb xc\s*:\s*({_FLOAT_TOKEN})\s+({_FLOAT_TOKEN})\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "xa", _parse_float(match.group(1)), 5, line)
            _store_metric(metrics, "xb", _parse_float(match.group(2)), 5, line)
            _store_metric(metrics, "xc", _parse_float(match.group(3)), 5, line)
            continue

        match = re.search(
            rf"coefficients xap xb xc\s*:\s*({_FLOAT_TOKEN})\s+({_FLOAT_TOKEN})\s+({_FLOAT_TOKEN})",
            line,
            flags=re.IGNORECASE,
        )
        if match is not None:
            _store_metric(metrics, "xap", _parse_float(match.group(1)), 6, line)
            _store_metric(metrics, "xb", _parse_float(match.group(2)), 6, line)
            _store_metric(metrics, "xc", _parse_float(match.group(3)), 6, line)
            continue

        matched_triplet = False
        for label, names, digits in _TRIPLET_LINE_SPECS:
            values = _match_line(line, label, 3)
            if values is None:
                continue
            for name, value in zip(names, values, strict=True):
                _store_metric(metrics, name, value, digits, line)
            matched_triplet = True
            break
        if matched_triplet:
            continue

    return metrics


def _scalar_array(value: float) -> xr.DataArray:
    return xr.DataArray(np.array([[value]], dtype=np.float64), dims=("y", "x"))


def _native_case_scalars(case: SixSParityCase, module_path: Path, source_dir: Path) -> dict[str, float]:
    config = case.sixs_config.model_copy(
        update={
            "source_dir": source_dir,
            "module_path": module_path,
            "auto_build": False,
            "native_threads": 1,
            "output_variables": SIXS_DEFAULT_OUTPUT_VARIABLES,
        }
    )
    runner = SixSNativeRunner(sixs_config=config, rt_setup=case.rt_setup)
    geometry = GeometryAngles.from_degrees(
        _scalar_array(case.sza_deg),
        _scalar_array(case.saa_deg),
        _scalar_array(case.vza_deg),
        _scalar_array(case.vaa_deg),
    )
    atmo = AtmosphericState(
        aot=_scalar_array(case.aot550),
        tcwv=_scalar_array(case.tcwv_cm),
        tco3=_scalar_array(case.tco3_atmcm),
        aot_unc=_scalar_array(0.0),
        tcwv_unc=_scalar_array(0.0),
        tco3_unc=_scalar_array(0.0),
        elevation=_scalar_array(case.elevation_km),
    )
    outputs = runner.compute_coefficients(
        geometry=geometry,
        atmo_state=atmo,
        band=case.band,
        output_variables=SIXS_DEFAULT_OUTPUT_VARIABLES,
    )
    return {
        name: float(np.asarray(data.values, dtype=np.float64).reshape(-1)[0])
        for name, data in outputs.items()
    }


def _native_case_scalars_subprocess(case: SixSParityCase, module_path: Path, source_dir: Path) -> dict[str, float]:
    payload = base64.b64encode(
        pickle.dumps(
            {
                "case": case,
                "module_path": os.fspath(module_path),
                "source_dir": os.fspath(source_dir),
            }
        )
    ).decode("ascii")
    script = """
import base64
import json
import pickle
import sys
from pathlib import Path

from siac.sixs_upstream_parity import _native_case_scalars

payload = pickle.loads(base64.b64decode(sys.argv[1].encode("ascii")))
metrics = _native_case_scalars(
    payload["case"],
    Path(payload["module_path"]),
    Path(payload["source_dir"]),
)
print(json.dumps(metrics, sort_keys=True, allow_nan=True))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script, payload],
        capture_output=True,
        text=True,
        check=False,
        cwd=Path(__file__).resolve().parents[2],
        env=os.environ.copy(),
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Native parity subprocess failed for case {case.name}.\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    json_line = next((line for line in reversed(completed.stdout.splitlines()) if line.strip().startswith("{")), None)
    if json_line is None:
        raise RuntimeError(
            f"Native parity subprocess for case {case.name} did not emit JSON.\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return json.loads(json_line)


def _format_number(value: float) -> str:
    if math.isnan(value):
        return "NaN"
    return f"{value:.12g}"


def _render_surface_reflectance_lines(mode: int, constant: float, spectrum: np.ndarray, wlinf: float, wlsup: float) -> list[str]:
    if mode == 0:
        return [_format_number(constant)]
    if mode > 0:
        return []
    if mode == -1:
        iinf = int((wlinf - 0.25) / 0.0025 + 1.5)
        isup = int((wlsup - 0.25) / 0.0025 + 1.5)
        values = spectrum[iinf - 1 : isup]
        return [
            f"{_format_number(wlinf)} {_format_number(wlsup)}",
            " ".join(_format_number(float(item)) for item in values),
        ]
    raise ValueError(f"Unsupported surface reflectance mode: {mode}")


def render_original_sixs_input(case: SixSParityCase) -> str:
    """Render an original 6S stdin deck using the same resolved inputs as the native bridge."""
    config = case.sixs_config
    response, wlinf, wlsup = _build_spectral_response(case.band)
    atmospheric_mode, radiosonde = _resolve_atmospheric_inputs(case.rt_setup, month=int(config.month))
    aerosol_mode, aerosol_inputs = _resolve_aerosol_inputs(case.rt_setup)
    surface_inputs = _resolve_surface_inputs(case.rt_setup)
    atmospheric_correction_mode, atmospheric_correction_value = _resolve_atmospheric_correction_inputs(case.rt_setup)

    iinf = int((wlinf - 0.25) / 0.0025 + 1.5)
    isup = int((wlsup - 0.25) / 0.0025 + 1.5)
    response_values = response[iinf - 1 : isup]

    lines: list[str] = [
        "0",
        " ".join(
            [
                _format_number(case.sza_deg),
                _format_number(case.saa_deg),
                _format_number(case.vza_deg),
                _format_number(case.vaa_deg),
                str(config.month),
                str(config.day),
            ]
        ),
        str(atmospheric_mode),
    ]

    if atmospheric_mode == 8:
        lines.append(f"{_format_number(case.tcwv_cm)} {_format_number(case.tco3_atmcm)}")
    elif atmospheric_mode == 7:
        for index in range(34):
            lines.append(
                " ".join(
                    [
                        _format_number(float(radiosonde["altitude_km"][index])),
                        _format_number(float(radiosonde["pressure_mb"][index])),
                        _format_number(float(radiosonde["temperature_k"][index])),
                        _format_number(float(radiosonde["water_g_m3"][index])),
                        _format_number(float(radiosonde["ozone_g_m3"][index])),
                    ]
                )
            )

    lines.append(str(aerosol_mode))
    if aerosol_mode == 4:
        lines.append(" ".join(_format_number(float(item)) for item in aerosol_inputs["mixture"]))
    elif aerosol_mode in {8, 9, 10}:
        component_count = int(aerosol_inputs["dist_component_count"])
        lines.append(
            " ".join(
                [
                    _format_number(float(aerosol_inputs["dist_rmin"])),
                    _format_number(float(aerosol_inputs["dist_rmax"])),
                    str(component_count),
                ]
            )
        )
        for component_index in range(component_count):
            x1 = float(aerosol_inputs["dist_x1"][component_index])
            x2 = float(aerosol_inputs["dist_x2"][component_index])
            x3 = float(aerosol_inputs["dist_x3"][component_index])
            cij = float(aerosol_inputs["dist_cij"][component_index])
            if aerosol_mode == 8:
                lines.append(" ".join([_format_number(x1), _format_number(x2), _format_number(cij)]))
            elif aerosol_mode == 9:
                lines.append(" ".join([_format_number(x1), _format_number(x2), _format_number(x3)]))
            else:
                lines.append(_format_number(x1))
            lines.append(
                " ".join(
                    _format_number(float(item))
                    for item in aerosol_inputs["dist_rn"][:, component_index]
                )
            )
            lines.append(
                " ".join(
                    _format_number(float(item))
                    for item in aerosol_inputs["dist_ri"][:, component_index]
                )
            )
    elif aerosol_mode == 11:
        count = int(aerosol_inputs["sun_count"])
        lines.append(str(count))
        for index in range(count):
            lines.append(
                " ".join(
                    [
                        _format_number(float(aerosol_inputs["sun_radius"][index])),
                        _format_number(float(aerosol_inputs["sun_dvlogr"][index])),
                    ]
                )
            )
        lines.append(" ".join(_format_number(float(item)) for item in aerosol_inputs["dist_rn"][:, 0]))
        lines.append(" ".join(_format_number(float(item)) for item in aerosol_inputs["dist_ri"][:, 0]))
    elif aerosol_mode == -1:
        count = int(aerosol_inputs["layer_count"])
        lines.append(str(count))
        for offset in range(count):
            index = count - offset - 1
            lines.append(
                " ".join(
                    [
                        _format_number(float(aerosol_inputs["layer_height"][index])),
                        _format_number(float(aerosol_inputs["layer_aot"][index])),
                        str(int(aerosol_inputs["layer_type"][index])),
                    ]
                )
            )
    elif aerosol_mode == 12:
        aerosol_path = Path(case.rt_setup.aerosol.model_path or "")
        if not aerosol_path:
            raise ValueError("user_model parity cases require `aerosol_model_path`.")
        lines.append(str(aerosol_path))

    visibility_mode = "-1" if aerosol_mode == 0 else "0"
    lines.append(visibility_mode)
    if aerosol_mode != 0:
        lines.append(_format_number(case.aot550))

    lines.extend(
        [
            _format_number(-case.elevation_km),
            "-1000",
            "1",
            f"{_format_number(wlinf)} {_format_number(wlsup)}",
            " ".join(_format_number(float(item)) for item in response_values),
            str(int(surface_inputs["inhomo"])),
        ]
    )

    if int(surface_inputs["inhomo"]) == 0:
        lines.append(str(int(surface_inputs["idirec"])))
        if int(surface_inputs["idirec"]) == 0:
            lines.append(str(int(surface_inputs["target_mode"])))
            lines.extend(
                _render_surface_reflectance_lines(
                    int(surface_inputs["target_mode"]),
                    float(surface_inputs["target_constant"]),
                    np.asarray(surface_inputs["target_spectrum"], dtype=np.float64),
                    wlinf,
                    wlsup,
                )
            )
        else:
            lines.append(str(int(surface_inputs["brdf_model"])))
            params = np.asarray(surface_inputs["brdf_params"], dtype=np.float64)
            options = np.asarray(surface_inputs["brdf_options"], dtype=np.int32)
            struct = np.asarray(surface_inputs["brdf_struct"], dtype=np.float64)
            optics = np.asarray(surface_inputs["brdf_optics"], dtype=np.float64)
            model = int(surface_inputs["brdf_model"])
            if model in {1, 4}:
                lines.append(" ".join(_format_number(float(item)) for item in params[:4]))
            elif model == 2:
                lines.append(" ".join(str(int(item)) for item in options[2:5]))
                lines.append(" ".join(_format_number(float(item)) for item in struct[:4]))
                lines.append(" ".join(_format_number(float(item)) for item in optics[:3]))
            elif model == 3:
                lines.append(" ".join(_format_number(float(item)) for item in params[:3]))
            elif model == 5:
                lines.append(" ".join(_format_number(float(item)) for item in params[:2]))
            elif model == 6:
                lines.append(" ".join(_format_number(float(item)) for item in params[:4]))
            elif model == 7:
                lines.append(
                    " ".join(
                        [
                            str(int(options[0])),
                            str(int(options[1])),
                        ]
                    )
                )
                lines.append(" ".join(_format_number(float(item)) for item in params[:2]))
                lines.append(" ".join(_format_number(float(item)) for item in params[2:5]))
            elif model == 8:
                lines.append(" ".join(_format_number(float(item)) for item in params[:3]))
            elif model == 9:
                lines.append(" ".join(_format_number(float(item)) for item in params[:9]))
            elif model in {10, 11}:
                lines.append(" ".join(_format_number(float(item)) for item in params[:3]))
            elif model == 0:
                solar_table = np.asarray(surface_inputs["brdf_table_solar"], dtype=np.float64)
                view_table = np.asarray(surface_inputs["brdf_table_view"], dtype=np.float64)
                for azimuth_index in range(13):
                    lines.append(
                        " ".join(
                            _format_number(float(item))
                            for item in solar_table[:, azimuth_index][::-1]
                        )
                    )
                for azimuth_index in range(13):
                    lines.append(
                        " ".join(
                            _format_number(float(item))
                            for item in view_table[:, azimuth_index][::-1]
                        )
                    )
                lines.append(_format_number(float(surface_inputs["brdf_spherical_albedo"])))
                lines.append(_format_number(float(surface_inputs["brdf_directional_reflectance"])))
            else:
                raise ValueError(f"Unsupported BRDF model for upstream parity deck: {model}")
    else:
        lines.append(
            " ".join(
                [
                    str(int(surface_inputs["target_mode"])),
                    str(int(surface_inputs["env_mode"])),
                    _format_number(float(surface_inputs["radius_km"])),
                ]
            )
        )
        lines.extend(
            _render_surface_reflectance_lines(
                int(surface_inputs["target_mode"]),
                float(surface_inputs["target_constant"]),
                np.asarray(surface_inputs["target_spectrum"], dtype=np.float64),
                wlinf,
                wlsup,
            )
        )
        lines.extend(
            _render_surface_reflectance_lines(
                int(surface_inputs["env_mode"]),
                float(surface_inputs["env_constant"]),
                np.asarray(surface_inputs["env_spectrum"], dtype=np.float64),
                wlinf,
                wlsup,
            )
        )

    lines.append(str(int(atmospheric_correction_mode)))
    if atmospheric_correction_mode >= 0:
        lines.append(_format_number(float(atmospheric_correction_value)))

    return "\n".join(lines) + "\n"


def build_original_sixs_executable(source_dir: Path, build_dir: Path, compiler: str = "gfortran") -> Path:
    """Build the untouched upstream 6SV2.1 executable in an isolated copy."""
    source_dir = source_dir.expanduser().resolve()
    build_dir = build_dir.expanduser().resolve()
    if not (source_dir / "main.f").exists():
        raise RuntimeError(f"{source_dir} does not look like a 6S source tree.")

    if build_dir.exists():
        shutil.rmtree(build_dir)
    shutil.copytree(
        source_dir,
        build_dir,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("*.o", "*.a", "*.so", "*.dylib", "__pycache__"),
    )

    clean = subprocess.run(
        ["make", "clean"],
        cwd=build_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if clean.returncode != 0:
        raise RuntimeError(f"Upstream 6S clean failed:\nstdout:\n{clean.stdout}\nstderr:\n{clean.stderr}")

    build = subprocess.run(
        [
            "make",
            "sixs",
            f"FC={compiler}",
            "EXTRA=-O -ffixed-line-length-132 -std=legacy",
        ],
        cwd=build_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if build.returncode != 0:
        raise RuntimeError(f"Upstream 6S build failed:\nstdout:\n{build.stdout}\nstderr:\n{build.stderr}")

    executable = build_dir / "sixsV2.1"
    if not executable.exists():
        raise RuntimeError(f"Expected upstream executable was not produced: {executable}")
    return executable


def default_parity_cases() -> tuple[SixSParityCase, ...]:
    """Representative parity cases covering multiple 6S branches."""
    uniform_band = SensorBand(
        name="PARITY_BAND",
        center_wavelength=665.0,
        bandwidth=30.0,
        resolution=10.0,
        band_index=1,
        rsrf_wavelengths_nm=np.arange(650.0, 680.0 + 2.5, 2.5, dtype=np.float64),
        rsrf_response=np.ones(13, dtype=np.float64),
    )
    return (
        SixSParityCase(
            name="lambertian_user_water_ozone_continental",
            description="Homogeneous Lambertian surface with user water/ozone and continental aerosol.",
            sixs_config=SixSAlgorithmConfig(
                month=7,
                day=12,
            ),
            rt_setup=RTSetupConfig(
                atmosphere={"profile": "user_water_ozone", "columns_mode": "input_columns"},
                aerosol={"profile": "continental"},
                surface={"mode": "homogeneous_lambertian", "target": {"kind": "constant", "constant": 0.0}},
                atmospheric_correction={"mode": "lambertian_reflectance", "value": 0.1},
                reference_reflectance=0.1,
            ),
            band=uniform_band,
            sza_deg=30.0,
            saa_deg=0.0,
            vza_deg=5.0,
            vaa_deg=40.0,
            aot550=0.15,
            tcwv_cm=2.0,
            tco3_atmcm=0.3,
            elevation_km=0.1,
        ),
        SixSParityCase(
            name="lambertian_builtin_surface_user_mixture",
            description="Homogeneous Lambertian built-in vegetation surface with user aerosol mixture.",
            sixs_config=SixSAlgorithmConfig(
                month=5,
                day=20,
            ),
            rt_setup=RTSetupConfig(
                atmosphere={"profile": "midlatitude_summer", "columns_mode": "input_columns"},
                aerosol={
                    "profile": "user_mixture",
                    "mixture": {"dust": 0.3, "water": 0.4, "oceanic": 0.2, "soot": 0.1},
                },
                surface={
                    "mode": "homogeneous_lambertian",
                    "target": {"kind": "built_in", "built_in": "green_vegetation"},
                },
                atmospheric_correction={"mode": "lambertian_reflectance", "value": 0.12},
                reference_reflectance=0.12,
            ),
            band=uniform_band,
            sza_deg=24.0,
            saa_deg=20.0,
            vza_deg=10.0,
            vaa_deg=70.0,
            aot550=0.22,
            tcwv_cm=1.9,
            tco3_atmcm=0.34,
            elevation_km=0.3,
        ),
        SixSParityCase(
            name="rahman_brdf_biomass_burning",
            description="Rahman BRDF surface with biomass-burning aerosol.",
            sixs_config=SixSAlgorithmConfig(
                month=7,
                day=12,
            ),
            rt_setup=RTSetupConfig(
                atmosphere={"profile": "user_water_ozone", "columns_mode": "input_columns"},
                aerosol={"profile": "biomass_burning"},
                surface={
                    "mode": "homogeneous_brdf",
                    "target": {"kind": "constant", "constant": 0.0},
                    "brdf": {
                        "model": "rahman",
                        "parameters": {
                            "intensity": 0.12,
                            "asymmetry_factor": -0.1,
                            "structural_parameter": 0.3,
                        },
                    },
                },
                atmospheric_correction={"mode": "brdf_reflectance", "value": 0.1},
                reference_reflectance=0.1,
            ),
            band=uniform_band,
            sza_deg=30.0,
            saa_deg=0.0,
            vza_deg=5.0,
            vaa_deg=40.0,
            aot550=0.15,
            tcwv_cm=2.0,
            tco3_atmcm=0.3,
            elevation_km=0.1,
        ),
        SixSParityCase(
            name="ocean_brdf_maritime",
            description="Ocean BRDF surface with maritime aerosol, including ocean-only outputs.",
            sixs_config=SixSAlgorithmConfig(
                month=7,
                day=12,
            ),
            rt_setup=RTSetupConfig(
                atmosphere={"profile": "user_water_ozone", "columns_mode": "input_columns"},
                aerosol={"profile": "maritime"},
                surface={
                    "mode": "homogeneous_brdf",
                    "target": {"kind": "constant", "constant": 0.0},
                    "brdf": {
                        "model": "ocean",
                        "parameters": {
                            "wind_speed": 5.0,
                            "wind_azimuth": 60.0,
                            "salinity": 35.0,
                            "pigment_concentration": 0.2,
                        },
                    },
                },
                atmospheric_correction={"mode": "brdf_reflectance", "value": 0.1},
                reference_reflectance=0.1,
            ),
            band=uniform_band,
            sza_deg=30.0,
            saa_deg=0.0,
            vza_deg=5.0,
            vaa_deg=40.0,
            aot550=0.15,
            tcwv_cm=2.0,
            tco3_atmcm=0.3,
            elevation_km=0.1,
        ),
    )


def _compare_metrics(original: dict[str, ParsedMetric], native: dict[str, float]) -> tuple[ComparedMetric, ...]:
    compared: list[ComparedMetric] = []
    for name in sorted(original):
        native_value = native.get(name)
        if native_value is None:
            continue
        original_value = original[name].value
        if math.isnan(original_value) and math.isnan(native_value):
            matched = True
            diff = math.nan
        elif math.isnan(original_value) or math.isnan(native_value):
            matched = False
            diff = math.nan
        else:
            diff = abs(native_value - original_value)
            matched = diff <= (original[name].tolerance + 1e-12)
        compared.append(
            ComparedMetric(
                name=name,
                original=None if math.isnan(original_value) else original_value,
                native=None if math.isnan(native_value) else native_value,
                tolerance=original[name].tolerance,
                abs_diff=None if math.isnan(diff) else diff,
                matched=matched,
                source=original[name].source,
            )
        )
    return tuple(compared)


def _run_original_case(case: SixSParityCase, executable: Path, case_dir: Path) -> tuple[dict[str, ParsedMetric], Path, Path, Path]:
    case_dir.mkdir(parents=True, exist_ok=True)
    input_path = case_dir / f"{case.name}.in"
    stdout_path = case_dir / f"{case.name}.stdout.txt"
    stderr_path = case_dir / f"{case.name}.stderr.txt"

    deck = render_original_sixs_input(case)
    input_path.write_text(deck, encoding="utf-8")

    completed = subprocess.run(
        [os.fspath(executable)],
        input=deck,
        capture_output=True,
        text=True,
        check=False,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Original 6S case {case.name} failed with exit code {completed.returncode}.\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return parse_original_sixs_stdout(completed.stdout), input_path, stdout_path, stderr_path


def run_upstream_parity_suite(
    *,
    source_dir: Path,
    upstream_build_dir: Path,
    native_build_dir: Path,
    output_dir: Path,
    compiler: str = "gfortran",
    cases: tuple[SixSParityCase, ...] | None = None,
) -> dict[str, Any]:
    """Run the full original-vs-native parity suite and return a JSON-serializable report."""
    source_dir = source_dir.expanduser().resolve()
    upstream_build_dir = upstream_build_dir.expanduser().resolve()
    native_build_dir = native_build_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    executable = build_original_sixs_executable(source_dir, upstream_build_dir, compiler=compiler)
    build_config = SixSAlgorithmConfig(
        source_dir=source_dir,
        build_dir=native_build_dir,
        compiler=compiler,
        build_profile="parity",
    )
    native_module_path = build_native_sixs_module(build_config)
    active_cases = cases or default_parity_cases()

    case_reports: list[CaseReport] = []
    for case in active_cases:
        case_dir = output_dir / case.name
        original_metrics, input_path, stdout_path, stderr_path = _run_original_case(case, executable, case_dir)
        native_metrics = _native_case_scalars_subprocess(case, native_module_path, source_dir)
        compared = _compare_metrics(original_metrics, native_metrics)
        mismatches = tuple(item for item in compared if not item.matched)
        case_reports.append(
            CaseReport(
                name=case.name,
                description=case.description,
                input_path=os.fspath(input_path),
                stdout_path=os.fspath(stdout_path),
                stderr_path=os.fspath(stderr_path),
                compared_variable_count=len(compared),
                matched_variable_count=len(compared) - len(mismatches),
                mismatches=mismatches,
                compared=compared,
            )
        )

    report = {
        "source_dir": os.fspath(source_dir),
        "upstream_executable": os.fspath(executable),
        "native_module": os.fspath(native_module_path),
        "all_cases_matched": all(case.matched for case in case_reports),
        "cases": [
            {
                **asdict(case_report),
                "matched": case_report.matched,
            }
            for case_report in case_reports
        ],
    }
    return report


def write_upstream_parity_report(report: dict[str, Any], path: Path) -> None:
    """Write the parity report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")


__all__ = [
    "CaseReport",
    "ComparedMetric",
    "ParsedMetric",
    "SixSParityCase",
    "build_original_sixs_executable",
    "default_parity_cases",
    "parse_original_sixs_stdout",
    "render_original_sixs_input",
    "run_upstream_parity_suite",
    "write_upstream_parity_report",
]

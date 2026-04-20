"""Build and smoke-test the native 6SV2.1 backend."""

from __future__ import annotations

import platform
import sys
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import _DIAGNOSTICS_DIRNAME, build_native_sixs_module
from siac.config import SixSAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles


def _field(value: float) -> xr.DataArray:
    return xr.DataArray(np.full((2, 2), value, dtype=np.float32), dims=("y", "x"))


def _geometry() -> GeometryAngles:
    return GeometryAngles.from_degrees(
        _field(30.0),
        _field(150.0),
        _field(5.0),
        _field(110.0),
    )


def _atmosphere() -> AtmosphericState:
    return AtmosphericState(
        aot=_field(0.15),
        tcwv=_field(2.0),
        tco3=_field(0.3),
        aot_unc=_field(0.01),
        tcwv_unc=_field(0.05),
        tco3_unc=_field(0.01),
        elevation=_field(0.1),
    )


def _band() -> SensorBand:
    return SensorBand(
        name="B04",
        center_wavelength=665.0,
        bandwidth=30.0,
        resolution=10.0,
        band_index=3,
        rsrf_wavelengths_nm=np.array([650.0, 657.5, 665.0, 672.5, 680.0], dtype=np.float64),
        rsrf_response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float64),
    )


def _local_source_dir() -> Path | None:
    root = Path("tmp/6s_upstream").resolve()
    if (root / "main.f").exists():
        return root
    matches = sorted(path.parent for path in root.rglob("main.f"))
    return matches[0] if matches else None


def _render_tree(root: Path, *, limit: int = 2000) -> str:
    if not root.exists():
        return f"{root} does not exist\n"
    lines: list[str] = []
    for entries, path in enumerate(sorted(root.rglob("*")), start=1):
        relative = path.relative_to(root)
        suffix = "/" if path.is_dir() else ""
        lines.append(f"{relative}{suffix}")
        if entries >= limit:
            lines.append(f"... truncated after {limit} entries")
            break
    if not lines:
        return ".\n"
    return "\n".join(lines) + "\n"


def _tail_text(path: Path, *, lines: int = 120) -> str:
    try:
        content = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return f"[failed to read {path}: {exc}]"
    return "\n".join(content[-lines:])


def _annotation_text(text: str, *, limit: int = 1800) -> str:
    summary = " | ".join(line.strip() for line in text.splitlines() if line.strip())
    if len(summary) > limit:
        return f"{summary[:limit]}..."
    return summary


def _extract_key_diagnostics(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return ""

    keywords = (
        "error:",
        "undefined reference",
        "collect2:",
        "ld:",
        "cannot find",
        "no such file",
        "traceback",
        "fatal",
        "undefined symbol",
    )
    matches = [line for line in lines if any(keyword in line.lower() for keyword in keywords)]
    if matches:
        return "\n".join(matches[:80])

    if len(lines) <= 80:
        return "\n".join(lines)
    return "\n".join(lines[-80:])


def _write_failure_diagnostics(build_dir: Path, exc: BaseException) -> list[Path]:
    diagnostics_dir = build_dir / _DIAGNOSTICS_DIRNAME
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    summary = " ".join(str(exc).split())
    if len(summary) > 4000:
        summary = f"{summary[:4000]}..."

    written_paths: list[Path] = []
    summary_path = diagnostics_dir / "smoke_failure_summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                f"exception_type={type(exc).__name__}",
                f"summary={summary}",
                f"platform={platform.platform()}",
                f"python={sys.version}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    written_paths.append(summary_path)

    traceback_path = diagnostics_dir / "smoke_traceback.txt"
    traceback_path.write_text(traceback.format_exc(), encoding="utf-8")
    written_paths.append(traceback_path)

    tree_path = diagnostics_dir / "build_tree.txt"
    tree_path.write_text(_render_tree(build_dir), encoding="utf-8")
    written_paths.append(tree_path)

    candidate_paths = sorted(
        {
            *build_dir.rglob("*.so"),
            *build_dir.rglob("*.dylib"),
            *build_dir.rglob("*.pyd"),
        }
    )
    candidates_path = diagnostics_dir / "module_candidates.txt"
    candidates_path.write_text(
        ("\n".join(str(path) for path in candidate_paths) + "\n") if candidate_paths else "",
        encoding="utf-8",
    )
    written_paths.append(candidates_path)
    return written_paths


def main() -> int:
    build_dir = Path("tmp/rt6s_ci_smoke").resolve()
    build_dir.mkdir(parents=True, exist_ok=True)
    local_source_dir = _local_source_dir()
    requested_outputs = ("xap", "xbp", "xcp", "tgasm", "sutott", "sast")

    build_config = SixSAlgorithmConfig(
        source_dir=local_source_dir,
        build_dir=build_dir,
        build_profile="release",
        native_threads=2,
        output_variables=requested_outputs,
    )
    module_path = build_native_sixs_module(build_config)
    runtime_config = build_config.model_copy(update={"module_path": module_path, "auto_build": False})

    backend = SixSBackend(sixs_config=runtime_config)
    backend.set_observation_time(datetime(2025, 7, 12, 10, 30))
    coeffs = backend.compute_coefficients(_geometry(), _atmosphere(), _band())

    assert coeffs.output_names == requested_outputs
    for name in requested_outputs:
        values = coeffs.get_output(name).values
        assert values.shape == (2, 2)
        assert np.all(np.isfinite(values)), f"Output {name} contains non-finite values"

    print(module_path)
    print(build_config.source_dir or build_config.source_url)
    print("native 6S smoke OK")
    return 0


if __name__ == "__main__":
    build_dir = Path("tmp/rt6s_ci_smoke").resolve()
    try:
        raise SystemExit(main())
    except Exception as exc:
        written_paths = _write_failure_diagnostics(build_dir, exc)
        summary = " ".join(str(exc).split())
        if len(summary) > 1800:
            summary = f"{summary[:1800]}..."
        print(
            "::error title=Native 6S smoke failed::"
            f"{type(exc).__name__}: {summary} [platform={platform.platform()}]"
        )
        print("Native 6S smoke diagnostics:")
        for path in written_paths:
            print(f"  - {path}")
        for log_name in (
            "build_failure_summary.txt",
            "f2py-meson.summary.txt",
            "f2py-meson.stderr.txt",
            "f2py-meson.stdout.txt",
            "f2py-distutils.summary.txt",
            "f2py-distutils.stderr.txt",
            "f2py-distutils.stdout.txt",
        ):
            log_path = build_dir / _DIAGNOSTICS_DIRNAME / log_name
            if log_path.exists():
                print(f"===== tail: {log_path} =====")
                log_text = log_path.read_text(encoding="utf-8", errors="replace")
                print(_tail_text(log_path))
                if log_name in {
                    "build_failure_summary.txt",
                    "f2py-meson.stderr.txt",
                    "f2py-meson.stdout.txt",
                    "f2py-distutils.stderr.txt",
                    "f2py-distutils.stdout.txt",
                }:
                    print(
                        f"::error title=Native 6S {log_name}::"
                        f"{_annotation_text(_extract_key_diagnostics(log_text))}"
                    )
        candidates_path = build_dir / _DIAGNOSTICS_DIRNAME / "module_candidates.txt"
        if candidates_path.exists():
            candidates_text = candidates_path.read_text(encoding="utf-8", errors="replace")
            if candidates_text.strip():
                print(
                    "::notice title=Native 6S module candidates::"
                    f"{_annotation_text(candidates_text)}"
                )
            else:
                build_tree_path = build_dir / _DIAGNOSTICS_DIRNAME / "build_tree.txt"
                tree_tail = _tail_text(build_tree_path, lines=80)
                top_level = "\n".join(
                    str(path.relative_to(build_dir)) + ("/" if path.is_dir() else "")
                    for path in sorted(build_dir.iterdir())
                )
                tree_excerpt = top_level + ("\n" if top_level else "") + tree_tail
                print(
                    "::notice title=Native 6S build tree tail::"
                    f"{_annotation_text(tree_excerpt)}"
                )
        traceback.print_exc()
        raise

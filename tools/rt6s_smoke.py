"""Build and smoke-test the native 6SV2.1 backend."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.algorithms.rt.direct.sixs_build import build_native_sixs_module
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


def main() -> int:
    build_dir = Path("tmp/rt6s_ci_smoke").resolve()
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
    raise SystemExit(main())

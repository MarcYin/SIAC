"""Unit tests for the configured output writer adapter."""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

from siac.adapters.output import ConfiguredOutputWriter
from siac.config.system import OutputDefaultsConfig
from siac.runtime import CorrectionDiagnostics, CorrectionResult


def _spatial_da(value: float) -> xr.DataArray:
    data = xr.DataArray(
        np.full((4, 4), value, dtype=np.float32),
        dims=["y", "x"],
        coords={
            "y": np.array([40.0, 30.0, 20.0, 10.0], dtype=np.float32),
            "x": np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float32),
        },
    )
    return data.rio.write_crs("EPSG:32632")


def test_configured_output_writer_writes_netcdf_and_quicklook(tmp_path):
    result = CorrectionResult(
        boa=xr.Dataset(
            {
                "B02": _spatial_da(0.10),
                "B03": _spatial_da(0.15),
                "B04": _spatial_da(0.20),
            }
        ),
        boa_unc=None,
        aot=_spatial_da(0.12),
        tcwv=_spatial_da(1.8),
        cloud_mask=_spatial_da(0.0).astype(bool),
        diagnostics=CorrectionDiagnostics(processing_time_s=1.5),
    )
    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_auxiliary=True,
            include_rgb=True,
            include_uncertainty=False,
        )
    )

    artifacts = writer.write(result, tmp_path / "products")

    assert artifacts["boa"].exists()
    assert artifacts["auxiliary"].exists()
    assert artifacts["quicklook.rgb"].exists()

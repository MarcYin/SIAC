"""
Layer 9 — Performance benchmarks.

These are not strict pass/fail tests but ensure operations complete in
reasonable time and allocate reasonable memory.  They are marked with
``@pytest.mark.slow`` so they can be skipped in quick CI runs.
"""

import time

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.grid.assembler import assemble_grids
from siac.domain import SensorBand, SensorConfig
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SurfacePrior,
)

# ── Helpers ──────────────────────────────────────────────────────────


def _build_large_obs(n: int = 256) -> ObservationBundle:
    """Build a large ObservationBundle (n×n)."""
    from datetime import datetime

    rng = np.random.RandomState(7)
    shape = (n, n)
    config = SensorConfig(
        sensor_id="PERF",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
        ),
    )
    toa = xr.Dataset(
        {
            b.name: xr.DataArray(rng.uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"])
            for b in config.bands
        }
    )
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )
    cloud = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
    return ObservationBundle(
        toa=toa,
        geometry=geometry,
        cloud_mask=cloud,
        sensor_config=config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 300000.0 + n * 10.0, 5500000.0 + n * 10.0),
    )


def _build_large_atmo(n: int = 256) -> AtmosphericState:
    shape = (n, n)
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )


def _build_large_surface(n: int = 256) -> SurfacePrior:
    shape = (n, n)
    brdf = BRDFKernelWeights(
        f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
        f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
    )
    return SurfacePrior(
        boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        kernels=brdf,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )


# ── Benchmark tests ──────────────────────────────────────────────────


@pytest.mark.slow
class TestGridAssemblerPerf:
    def test_256x256_completes_in_time(self, mock_rt_model):
        """Grid assembler on 256×256 should complete in < 5 seconds."""
        obs = _build_large_obs(256)
        atmo = _build_large_atmo(256)
        surface = _build_large_surface(256)

        t0 = time.perf_counter()
        sib = assemble_grids(obs, atmo, surface, mock_rt_model)
        elapsed = time.perf_counter() - t0

        assert sib is not None
        assert elapsed < 5.0, f"Grid assembler took {elapsed:.2f}s (limit: 5s)"

    def test_512x512_completes_in_time(self, mock_rt_model):
        """Grid assembler on 512×512 should complete in < 15 seconds."""
        obs = _build_large_obs(512)
        atmo = _build_large_atmo(512)
        surface = _build_large_surface(512)

        t0 = time.perf_counter()
        sib = assemble_grids(obs, atmo, surface, mock_rt_model)
        elapsed = time.perf_counter() - t0

        assert sib is not None
        assert elapsed < 15.0, f"Grid assembler took {elapsed:.2f}s (limit: 15s)"

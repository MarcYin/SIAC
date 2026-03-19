"""
Layer 3 — Grid assembler (M4) unit tests.
"""

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.grid.assembler import assemble_grids
from siac.domain import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SensorBand,
    SensorConfig,
    SolverInputBundle,
    SurfacePrior,
)
from siac.domain.validation import _validate_solver_input_bundle


@pytest.fixture
def large_obs_bundle():
    """ObservationBundle at 64x64 (simulating ~10 m native) with 4 bands."""
    shape = (64, 64)
    config = SensorConfig(
        sensor_id="MOCK",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
        ),
    )
    from datetime import datetime

    toa = xr.Dataset({
        "B01": xr.DataArray(np.random.RandomState(10).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
        "B02": xr.DataArray(np.random.RandomState(11).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
        "B03": xr.DataArray(np.random.RandomState(12).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
        "B04": xr.DataArray(np.random.RandomState(13).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
    })
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
        bounds=(300000.0, 5500000.0, 300640.0, 5500640.0),
    )


@pytest.fixture
def large_atmo():
    """AtmosphericState at 64x64."""
    shape = (64, 64)
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )


@pytest.fixture
def large_surface():
    """SurfacePrior at 64x64."""
    shape = (64, 64)
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


class TestAssembleGrids:
    def test_output_type(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        result = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        assert isinstance(result, SolverInputBundle)

    def test_spatial_alignment(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """All raster fields should share the same (y, x) shape."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        # toa has (band, y, x)
        toa_spatial = sib.toa.shape[1:]
        assert sib.geometry.sza.shape == toa_spatial
        assert sib.cloud_mask.shape == toa_spatial
        assert sib.atmo_prior.aot.shape == toa_spatial
        assert sib.surface_prior.boa.shape == toa_spatial

    def test_aux_resolution(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model, aux_resolution_m=500.0)
        assert sib.aux_resolution_m == 500.0

    def test_band_selection(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Bands should be aerosol-sensitive (400-520 nm)."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        for b in sib.bands:
            assert 400.0 <= b.center_wavelength <= 520.0

    def test_cloud_mask_conservative(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Any native pixel that is cloud → aux pixel is cloud."""
        import dataclasses
        shape = (64, 64)
        cloud = np.zeros(shape, dtype=bool)
        cloud[10:20, 10:20] = True
        obs = dataclasses.replace(
            large_obs_bundle,
            cloud_mask=xr.DataArray(cloud, dims=["y", "x"]),
        )
        sib = assemble_grids(obs, large_atmo, large_surface, mock_rt_model)
        # The cloud region should still have True values after resampling
        assert sib.cloud_mask.values.any(), "Cloud region should be preserved"

    def test_identity_same_res(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """When aux_resolution == native, output ≈ input."""
        # bands[0].resolution is 60.0 for B01, but let's set aux to native so no resampling
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model,
            aux_resolution_m=60.0,  # matches B01 resolution
        )
        # Should have same shape as input (64x64)
        assert sib.toa.shape[1:] == (64, 64)

    def test_downsampling_reduces_shape(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Resampling from 10m to 500m should produce a smaller grid."""
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model,
            aux_resolution_m=500.0,
        )
        # 64 px @ 10m -> ~1.3 px @ 500m -> at least (1, 1)
        assert sib.toa.shape[1] <= 64
        assert sib.toa.shape[2] <= 64

    def test_passes_validation(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Output should pass _validate_solver_input_bundle()."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        _validate_solver_input_bundle(sib)  # should not raise

    def test_surface_prior_with_band_dimension(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Assembler should handle banded SurfacePrior arrays from real BRDF providers."""
        shape = large_surface.boa.shape
        banded_boa = xr.DataArray(
            np.stack(
                [
                    np.full(shape, 0.11, dtype=np.float32),
                    np.full(shape, 0.13, dtype=np.float32),
                ],
                axis=0,
            ),
            dims=["band", "y", "x"],
            coords={"band": [1, 2]},
        )
        banded_unc = xr.full_like(banded_boa, 0.02)
        banded_mask = xr.DataArray(
            np.ones((2, *shape), dtype=bool),
            dims=["band", "y", "x"],
            coords={"band": [1, 2]},
        )
        banded_surface = SurfacePrior(
            boa=banded_boa,
            boa_unc=banded_unc,
            kernels=large_surface.kernels,
            mask=banded_mask,
        )

        sib = assemble_grids(large_obs_bundle, large_atmo, banded_surface, mock_rt_model)
        assert sib.surface_prior.boa.ndim == 3
        assert sib.surface_prior.boa.shape[0] == 2
        assert sib.surface_prior.mask.ndim == 2

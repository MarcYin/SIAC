"""
Layer 3 — Grid assembler (M4) unit tests.
"""

import dataclasses

import numpy as np
import pytest
import xarray as xr
from rasterio.enums import Resampling as RasterioResampling

from siac.algorithms.grid.assembler import assemble_grids
from siac.domain import SensorBand, SensorConfig
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.validation import validate_solver_input_bundle


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

    def test_aux_resolution_falls_back_to_solver_grid_for_legacy_callers(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model, aux_resolution_m=320.0)

        assert sib.toa.shape[1:] == (2, 2)
        assert sib.aerosol_resolution_m == 320.0

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
        """When aerosol_resolution matches the scene pixel size, output ≈ input."""
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model,
            aerosol_resolution_m=10.0,
        )
        # Should have same shape as input (64x64)
        assert sib.toa.shape[1:] == (64, 64)

    def test_downsampling_reduces_shape(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Resampling to a coarser aerosol grid should produce a smaller grid."""
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model,
            aerosol_resolution_m=500.0,
        )
        # 64 px @ 10m -> ~1.3 px @ 500m -> at least (1, 1)
        assert sib.toa.shape[1] <= 64
        assert sib.toa.shape[2] <= 64

    def test_toa_downsampling_uses_gdal_average_when_georeferenced(
        self,
        monkeypatch: pytest.MonkeyPatch,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        import rioxarray  # noqa: F401
        from affine import Affine
        from rioxarray.raster_array import RasterArray

        x = 300000.0 + (np.arange(64, dtype=np.float64) + 0.5) * 10.0
        y = 5500640.0 - (np.arange(64, dtype=np.float64) + 0.5) * 10.0
        transform = Affine(10.0, 0.0, 300000.0, 0.0, -10.0, 5500640.0)

        def _georef(da: xr.DataArray) -> xr.DataArray:
            out = da.assign_coords({"x": x, "y": y})
            out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
            out = out.rio.write_crs("EPSG:32632")
            return out.rio.write_transform(transform)

        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset({name: _georef(data) for name, data in large_obs_bundle.toa.data_vars.items()}),
        )

        calls: list[RasterioResampling] = []
        original_reproject_match = RasterArray.reproject_match

        def _fake_reproject_match(self, target, *, resampling, **kwargs):  # type: ignore[no-untyped-def]
            calls.append(resampling)
            return original_reproject_match(self, target, resampling=resampling, **kwargs)

        monkeypatch.setattr(RasterArray, "reproject_match", _fake_reproject_match)

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=320.0,
        )

        assert sib.toa.shape[1:] == (2, 2)
        assert calls == [RasterioResampling.average, RasterioResampling.average]

    def test_aerosol_resolution_controls_grid_and_georeference(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=320.0,
        )

        assert sib.toa.shape[1:] == (2, 2)
        assert sib.atmo_prior.aot.shape == (2, 2)
        assert sib.atmo_prior.aot.rio.crs is not None
        assert sib.atmo_prior.aot.rio.crs.to_string() == "EPSG:32632"
        assert sib.atmo_prior.aot.coords["x"].values.tolist() == pytest.approx([300160.0, 300480.0])
        assert sib.atmo_prior.aot.coords["y"].values.tolist() == pytest.approx([5500480.0, 5500160.0])

    def test_passes_validation(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Output should pass validate_solver_input_bundle()."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        validate_solver_input_bundle(sib)  # should not raise

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

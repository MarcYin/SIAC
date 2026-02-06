"""
Integration smoke test: solver -> corrector pipeline.
"""

import numpy as np
import pytest
import xarray as xr

from siac.core.types import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    RTCoefficients,
    SensorBand,
    SurfacePrior,
    SENTINEL2A_CONFIG,
)
from siac.correction.atmospheric import AtmosphericCorrector
from siac.solver.multigrid import MultiGridSolver, MultiGridConfig


@pytest.mark.integration
class TestPipelineSmoke:
    """End-to-end smoke test with synthetic data."""

    @pytest.fixture
    def synthetic_scene(self):
        """Create a synthetic 32x32 scene."""
        shape = (32, 32)
        np.random.seed(99)

        # TOA reflectance (3 bands as Dataset for corrector, DataArray for solver)
        bands_list = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
            SensorBand("B04", 665.0, 30.0, 10.0, 2),
        ]

        toa_vals = np.random.RandomState(42).uniform(0.05, 0.35, (3, *shape)).astype(np.float32)
        toa_da = xr.DataArray(toa_vals, dims=["band", "y", "x"])
        toa_ds = xr.Dataset({
            b.name: xr.DataArray(toa_vals[i], dims=["y", "x"])
            for i, b in enumerate(bands_list)
        })

        # Geometry
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        # Atmospheric state prior
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        # Surface prior
        brdf_weights = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )

        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            kernels=brdf_weights,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])

        return {
            "toa_da": toa_da,
            "toa_ds": toa_ds,
            "geometry": geometry,
            "atmo_prior": atmo_prior,
            "surface_prior": surface_prior,
            "cloud_mask": cloud_mask,
            "bands": bands_list,
        }

    def test_solver_then_corrector(self, synthetic_scene, mock_rt_model):
        """Solver output feeds into corrector and produces valid BOA."""
        scene = synthetic_scene

        # Run solver (minimal config for speed)
        config = MultiGridConfig(n_levels=2, min_grid_size=8, max_iter_per_level=5)
        solver = MultiGridSolver(config)

        solver_result = solver.solve(
            toa=scene["toa_da"],
            surface_prior=scene["surface_prior"],
            geometry=scene["geometry"],
            atmo_prior=scene["atmo_prior"],
            rt_model=mock_rt_model,
            cloud_mask=scene["cloud_mask"],
            bands=scene["bands"],
        )

        # Basic solver result checks
        assert solver_result.aot.shape == (32, 32)
        assert solver_result.tcwv.shape == (32, 32)
        assert solver_result.success or solver_result.n_iterations > 0

        # Build updated atmospheric state from solver output
        solved_atmo = scene["atmo_prior"].with_updated_aot_tcwv(
            aot=solver_result.aot,
            tcwv=solver_result.tcwv,
        )

        # Run corrector
        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        correction_result = corrector.correct(
            scene["toa_ds"], scene["geometry"], solved_atmo
        )

        # BOA should be in physically valid range
        for band_name in ["B02", "B03", "B04"]:
            boa = correction_result.boa[band_name].values
            valid = np.isfinite(boa) & (boa > 0)
            assert valid.sum() > 0, f"No valid BOA pixels for {band_name}"
            assert np.nanmax(boa[valid]) < 1.5, f"BOA too high for {band_name}"
            assert np.nanmin(boa[valid]) > -0.1, f"BOA too low for {band_name}"

    def test_corrector_standalone(self, synthetic_scene, mock_rt_model):
        """Corrector alone should produce valid BOA from synthetic inputs."""
        scene = synthetic_scene

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(
            scene["toa_ds"], scene["geometry"], scene["atmo_prior"]
        )

        assert len(result.boa.data_vars) == 3
        for band_name in result.boa.data_vars:
            boa = result.boa[band_name].values
            valid = np.isfinite(boa)
            assert valid.any()

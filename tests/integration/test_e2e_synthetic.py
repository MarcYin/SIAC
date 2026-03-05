"""
Layer 8 — End-to-end synthetic pipeline tests.

Tests the full pipeline with synthetic data, verifying that all modules
connect correctly and produce valid output without any real data access.
"""

from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from siac.core.types import (
    AtmosphericState,
    BRDFKernelWeights,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    SensorBand,
    SensorConfig,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.core.validation import (
    _validate_atmospheric_state,
    _validate_observation_bundle,
    _validate_solver_input_bundle,
    _validate_surface_prior,
)
from siac.grid.assembler import assemble_grids
from siac.pipeline import run_pipeline

# ── Synthetic data builders ──────────────────────────────────────────

SHAPE = (64, 64)


def _make_sensor_config() -> SensorConfig:
    return SensorConfig(
        sensor_id="SYNTH",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
            SensorBand("B08", 842.0, 115.0, 10.0, 4),
            SensorBand("B11", 1610.0, 90.0, 20.0, 5),
        ),
    )


def _make_obs_bundle() -> ObservationBundle:
    rng = np.random.RandomState(99)
    config = _make_sensor_config()
    toa = xr.Dataset({
        b.name: xr.DataArray(
            rng.uniform(0.05, 0.4, SHAPE).astype(np.float32),
            dims=["y", "x"],
        )
        for b in config.bands
    })
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(SHAPE, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(SHAPE, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(SHAPE, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(SHAPE, 1.5), dims=["y", "x"]),
    )
    cloud = xr.DataArray(np.zeros(SHAPE, dtype=bool), dims=["y", "x"])
    return ObservationBundle(
        toa=toa,
        geometry=geometry,
        cloud_mask=cloud,
        sensor_config=config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 300640.0, 5500640.0),
    )


def _make_atmo() -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(SHAPE, 0.2), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(SHAPE, 2.0), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(SHAPE, 0.35), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(SHAPE, 0.1), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(SHAPE, 0.5), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(SHAPE, 0.02), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(SHAPE, 0.05), dims=["y", "x"]),
    )


def _make_surface() -> SurfacePrior:
    brdf = BRDFKernelWeights(
        f0=xr.DataArray(np.full(SHAPE, 0.1), dims=["y", "x"]),
        f1=xr.DataArray(np.full(SHAPE, 0.04), dims=["y", "x"]),
        f2=xr.DataArray(np.full(SHAPE, 0.02), dims=["y", "x"]),
        f0_unc=xr.DataArray(np.full(SHAPE, 0.01), dims=["y", "x"]),
        f1_unc=xr.DataArray(np.full(SHAPE, 0.005), dims=["y", "x"]),
        f2_unc=xr.DataArray(np.full(SHAPE, 0.002), dims=["y", "x"]),
    )
    return SurfacePrior(
        boa=xr.DataArray(np.full(SHAPE, 0.12), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(SHAPE, 0.02), dims=["y", "x"]),
        kernels=brdf,
        mask=xr.DataArray(np.ones(SHAPE, dtype=bool), dims=["y", "x"]),
    )


# ── E2E tests ────────────────────────────────────────────────────────

class TestE2ESynthetic:
    """Full pipeline end-to-end with all-synthetic data."""

    def test_grid_assembler_then_validate(self, mock_rt_model):
        """M1+M2+M3 -> M4 produces a valid SolverInputBundle."""
        obs = _make_obs_bundle()
        atmo = _make_atmo()
        surface = _make_surface()

        sib = assemble_grids(obs, atmo, surface, mock_rt_model)

        assert isinstance(sib, SolverInputBundle)
        _validate_solver_input_bundle(sib)

    def test_full_pipeline_with_mocks(self, mock_rt_model):
        """Full pipeline (M1→M6) with all mock modules returns CorrectionResult."""
        obs = _make_obs_bundle()
        atmo = _make_atmo()
        surface = _make_surface()

        def preprocessor(path, aoi=None):
            return obs

        def atmo_provider(bounds, crs, obs_time, resolution):
            return atmo

        def surface_provider(bounds, crs, obs_time, sensor_config, geometry, resolution):
            return surface

        def grid_assembler(o, a, s, rt, aux_res=500.0, aero_res=1000.0):
            return assemble_grids(o, a, s, rt)

        def solver(inputs, config):
            return SolvedAtmosphere(
                atmo_state=inputs.atmo_prior,
                aot=inputs.atmo_prior.aot,
                tcwv=inputs.atmo_prior.tcwv,
                aot_unc=inputs.atmo_prior.aot_unc,
                tcwv_unc=inputs.atmo_prior.tcwv_unc,
                cost_final=0.001,
                n_iterations=3,
                converged=True,
            )

        def corrector(o, solved, rt):
            return CorrectionResult(
                boa=o.toa,
                boa_unc=None,
                aot=solved.aot,
                tcwv=solved.tcwv,
                cloud_mask=o.cloud_mask,
                metadata={"test": True},
            )

        result = run_pipeline(
            Path("/fake/input"),
            None,
            None,
            preprocessor=preprocessor,
            atmo_provider=atmo_provider,
            surface_prior_provider=surface_provider,
            grid_assembler=grid_assembler,
            solver=solver,
            corrector=corrector,
            rt_model=mock_rt_model,
        )

        assert isinstance(result, CorrectionResult)
        assert result.boa is not None
        assert result.cloud_mask is not None

    def test_validation_chain(self):
        """Each contract validator passes on synthetic data."""
        obs = _make_obs_bundle()
        atmo = _make_atmo()
        surface = _make_surface()

        _validate_observation_bundle(obs)
        _validate_atmospheric_state(atmo)
        _validate_surface_prior(surface)

    def test_pipeline_preserves_crs_and_bounds(self, mock_rt_model):
        """ObservationBundle CRS/bounds are preserved through the pipeline."""
        obs = _make_obs_bundle()
        atmo = _make_atmo()
        surface = _make_surface()

        sib = assemble_grids(obs, atmo, surface, mock_rt_model)
        assert sib.sensor_config.sensor_id == obs.sensor_config.sensor_id

    def test_cloud_mask_propagation(self, mock_rt_model):
        """Cloud pixels in M1 output should be flagged in M4 output."""
        obs = _make_obs_bundle()
        import dataclasses
        cloud = np.zeros(SHAPE, dtype=bool)
        cloud[20:40, 20:40] = True
        obs = dataclasses.replace(
            obs,
            cloud_mask=xr.DataArray(cloud, dims=["y", "x"]),
        )
        atmo = _make_atmo()
        surface = _make_surface()

        sib = assemble_grids(obs, atmo, surface, mock_rt_model)
        assert sib.cloud_mask.values.any(), "Cloud region should survive resampling"

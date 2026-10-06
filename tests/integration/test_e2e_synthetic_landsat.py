"""Synthetic OLI M1-M6 integration and Collection-2 preprocessing seam tests.

Like the Sentinel-style synthetic E2E test, M2/M3/M5/M6 and RT are injected;
M4 is real. No emulator forward pass, BRDF kernels, or network are required.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pytest
import xarray as xr

from siac.adapters.auth import CredentialManager
from siac.adapters.satellite.landsat import LandsatPreprocessor
from siac.adapters.satellite.observation import preprocess_to_observation_bundle
from siac.algorithms.grid.assembler import assemble_grids
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.catalog.sensors import get_sensor_config
from siac.config import SIACConfig
from siac.runtime import (
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.validation import (
    validate_atmospheric_state,
    validate_observation_bundle,
    validate_solver_input_bundle,
    validate_surface_prior,
)
from siac.workflows.pipeline import run_pipeline
from siac.workflows.scene import execute_plan
from tests.integration.test_e2e_synthetic import _make_atmo, _make_surface
from tests.landsat_scene import CRS as SCENE_CRS
from tests.landsat_scene import HEIGHT, ORIGIN_X, ORIGIN_Y, WIDTH, _make_scene

if TYPE_CHECKING:
    from siac.adapters.satellite.observation import SatellitePreprocessorLike
    from siac.domain import SensorConfig

pytestmark = pytest.mark.integration
SHAPE = (64, 64)
BANDS = {f"B{n}" for n in range(1, 8)}
CRS = "EPSG:32632"
BOUNDS = (300000.0, 5500000.0, 301920.0, 5501920.0)


@pytest.fixture(autouse=True)
def _forbid_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_request(*args: Any, **kwargs: Any) -> Any:
        pytest.fail("Synthetic Landsat tests must not make network requests")

    monkeypatch.setattr("requests.sessions.Session.request", fail_request)


def _make_observation(config: SensorConfig) -> ObservationBundle:
    rng = np.random.RandomState(99)
    coords = {
        "x": BOUNDS[0] + (np.arange(SHAPE[1]) + 0.5) * 30.0,
        "y": BOUNDS[3] - (np.arange(SHAPE[0]) + 0.5) * 30.0,
    }

    def raster(values: Any) -> xr.DataArray:
        return cast(
            "xr.DataArray", xr.DataArray(values, dims=("y", "x"), coords=coords).rio.write_crs(CRS)
        )

    toa = xr.Dataset(
        {
            band.name: raster(rng.uniform(0.05, 0.4, SHAPE).astype(np.float32))
            for band in config.bands
        }
    )
    cloud = np.zeros(SHAPE, dtype=bool)
    cloud[20:40, 20:40] = True
    return ObservationBundle(
        toa=toa,
        geometry=GeometryAngles(
            sza=raster(np.full(SHAPE, 0.5)),
            saa=raster(np.full(SHAPE, 2.5)),
            vza=raster(np.full(SHAPE, 0.1)),
            vaa=raster(np.full(SHAPE, 1.5)),
        ),
        cloud_mask=raster(cloud),
        sensor_config=config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs=CRS,
        bounds=BOUNDS,
    )


def _injected_stages(rt_model: Any, grids: list[SolverInputBundle]) -> dict[str, Any]:
    atmo = _make_atmo()
    surface = _make_surface()
    validate_atmospheric_state(atmo)
    validate_surface_prior(surface)

    def atmo_provider(
        bounds: tuple[float, float, float, float], crs: str, obs_time: datetime, resolution: float
    ) -> AtmosphericState:
        return atmo

    def surface_provider(
        observation: ObservationBundle,
        atmo_state: AtmosphericState,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
        return surface

    def grid_assembler(
        observation: ObservationBundle,
        atmo_state: AtmosphericState,
        surface_prior: SurfacePrior,
        rt: Any,
        **kwargs: Any,
    ) -> SolverInputBundle:
        # The default pipeline water mask is remote; synthetic tests omit that optional input.
        kwargs["water_mask_path"] = None
        # PSF setup uses native BRDF kernels; omit it as in the original synthetic E2E.
        kwargs["toa_psf_config"] = None
        bundle = assemble_grids(observation, atmo_state, surface_prior, rt, **kwargs)
        validate_solver_input_bundle(bundle)
        grids.append(bundle)
        return bundle

    def solver(inputs: SolverInputBundle, config: Any) -> SolvedAtmosphere:
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

    def corrector(
        observation: ObservationBundle, solved: SolvedAtmosphere, rt: Any
    ) -> CorrectionResult:
        return CorrectionResult(
            boa=observation.toa,
            boa_unc=None,
            aot=solved.aot,
            tcwv=solved.tcwv,
            cloud_mask=observation.cloud_mask,
            metadata={"corrector": "synthetic_identity"},
        )

    return {
        "atmo_provider": atmo_provider,
        "surface_prior_provider": surface_provider,
        "grid_assembler": grid_assembler,
        "solver": solver,
        "corrector": corrector,
        "rt_model": rt_model,
    }


def _assert_grid(bundle: SolverInputBundle, config: SensorConfig) -> None:
    assert bundle.sensor_config == config
    assert set(bundle.toa.band.values) == {
        band.name for band in config.default_aerosol_solver_bands()
    }
    assert bundle.cloud_mask.values.any()
    assert not bundle.cloud_mask.values.all()


@pytest.mark.parametrize("satellite", ["L8", "L9"])
def test_full_oli_pipeline_with_injected_stages(satellite: str, mock_rt_model: Any) -> None:
    config = get_sensor_config("OLI", satellite)
    assert {band.name for band in config.bands} == BANDS
    assert all(band.resolution == 30.0 for band in config.bands)
    observation = _make_observation(config)
    validate_observation_bundle(observation)
    grids: list[SolverInputBundle] = []

    def preprocessor(path: Path, aoi: Any = None) -> ObservationBundle:
        return observation

    result = run_pipeline(
        Path("synthetic-oli"),
        None,
        None,
        preprocessor=preprocessor,
        **_injected_stages(mock_rt_model, grids),
    )

    assert isinstance(result, CorrectionResult)
    assert set(result.boa.data_vars) == BANDS
    assert len(grids) == 1
    _assert_grid(grids[0], config)
    assert grids[0].toa.shape[1:] == (16, 16)  # 64 native 30 m pixels -> 120 m grid
    assert str(grids[0].toa.rio.crs) == CRS
    assert grids[0].toa.rio.bounds() == pytest.approx(BOUNDS)
    assert np.all(grids[0].cloud_mask.values[5:10, 5:10])
    for name in BANDS:
        assert result.boa[name].shape == SHAPE
        assert str(result.boa[name].rio.crs) == observation.crs
        assert result.boa[name].rio.bounds() == pytest.approx(observation.bounds)
        assert result.boa[name].rio.resolution() == pytest.approx((30.0, -30.0))
    xr.testing.assert_identical(result.cloud_mask, observation.cloud_mask)


@pytest.mark.parametrize("spacecraft", ["LANDSAT_8", "LANDSAT_9"])
def test_collection2_preprocessing_and_execution_plan(
    tmp_path: Path, spacecraft: str, mock_rt_model: Any
) -> None:
    scene = _make_scene(tmp_path, spacecraft=spacecraft)
    observation = preprocess_to_observation_bundle(
        # The adapter only reads sensor_config, which Landsat exposes as a property.
        cast(
            "SatellitePreprocessorLike",
            LandsatPreprocessor(config={"cloud_mask": {"provider": "product_qa"}}),
        ),
        scene,
    )
    validate_observation_bundle(observation)
    assert set(observation.toa.data_vars) == BANDS
    assert observation.crs == SCENE_CRS
    assert observation.bounds == pytest.approx(
        (ORIGIN_X, ORIGIN_Y - HEIGHT * 30.0, ORIGIN_X + WIDTH * 30.0, ORIGIN_Y)
    )
    assert bool(observation.cloud_mask.values[12, 12])
    assert not bool(observation.cloud_mask.values[7, 0])

    grids: list[SolverInputBundle] = []
    config = SIACConfig.model_validate(
        {"sensor": "auto", "algorithms": {"cloud_mask": {"provider": "product_qa"}}}
    )
    plan = build_execution_plan(
        SceneProcessRequest(config=config, input_path=scene, auth=CredentialManager()),
        **_injected_stages(mock_rt_model, grids),
    )
    planned_observation = plan.preprocessor(scene, None)
    validate_observation_bundle(planned_observation)
    assert planned_observation.sensor_config.sensor_id == "OLI"
    assert planned_observation.sensor_config.satellite_id == spacecraft.replace("LANDSAT_", "L")
    assert planned_observation.bounds == observation.bounds
    xr.testing.assert_identical(planned_observation.cloud_mask, observation.cloud_mask)

    # Planning eagerly loads cloud/reference bands; real M4 loads remaining solver bands lazily.
    assert "B4" in planned_observation.toa
    result = execute_plan(plan)
    assert len(grids) == 1
    _assert_grid(grids[0], planned_observation.sensor_config)
    assert str(result.boa.rio.crs) == SCENE_CRS
    assert result.boa.rio.bounds() == pytest.approx(observation.bounds)
    xr.testing.assert_identical(result.cloud_mask, observation.cloud_mask)

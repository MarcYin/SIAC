"""Workflow tests for the real scene planning and execution path."""

from __future__ import annotations

import pytest

from siac.adapters.auth import CredentialManager
from siac.adapters.output import ConfiguredOutputWriter
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.domain.aoi import AOI
from siac.workflows.scene import process_scene


def _empty_auth() -> CredentialManager:
    return CredentialManager()


def test_build_execution_plan_resolves_runtime_aoi_and_output_writer(
    tmp_path,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
):
    cfg = SIACConfig(sensor="s2", output={"defaults": {"format": "netcdf"}})
    request = SceneProcessRequest(
        config=cfg,
        input_path=tmp_path / "in.SAFE",
        output_path=tmp_path / "products",
        aoi=[1.0, 2.0, 3.0, 4.0],
        auth=_empty_auth(),
    )

    plan = build_execution_plan(
        request,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert isinstance(plan.runtime_aoi, AOI)
    assert plan.runtime_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert isinstance(plan.output_writer, ConfiguredOutputWriter)
    assert plan.output_path == tmp_path / "products"


def test_build_execution_plan_accepts_geojson_dict_aoi(
    tmp_path,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
):
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(
        config=cfg,
        input_path=tmp_path / "in.SAFE",
        aoi={
            "type": "Polygon",
            "coordinates": [
                [
                    [1.0, 2.0],
                    [3.0, 2.0],
                    [3.0, 4.0],
                    [1.0, 4.0],
                    [1.0, 2.0],
                ]
            ],
        },
        auth=_empty_auth(),
    )

    plan = build_execution_plan(
        request,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert isinstance(plan.runtime_aoi, AOI)
    assert plan.runtime_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))


def test_process_scene_executes_real_plan_building(
    tmp_path,
    mock_observation_bundle,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
    mock_solved_atmosphere,
):
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(
        config=cfg,
        input_path=tmp_path / "in.SAFE",
        aoi=[1.0, 2.0, 3.0, 4.0],
        auth=_empty_auth(),
    )

    result = process_scene(
        request,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert result.boa.identical(mock_observation_bundle.toa)
    assert result.cloud_mask.identical(mock_observation_bundle.cloud_mask)
    assert result.aot.identical(mock_solved_atmosphere.aot)
    assert result.tcwv.identical(mock_solved_atmosphere.tcwv)
    assert result.diagnostics.processing_time_s == pytest.approx(0.01)

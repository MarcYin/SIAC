"""Workflow tests for the real scene planning and execution path."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

import siac.workflows.scene as scene_workflow
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
    mock_observation_bundle,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
):
    input_path = tmp_path / "in.SAFE"
    input_path.mkdir()
    cfg = SIACConfig(sensor="s2", output={"defaults": {"format": "netcdf"}})
    request = SceneProcessRequest(
        config=cfg,
        input_path=input_path,
        output_path=tmp_path / "products",
        aoi=[1.0, 2.0, 3.0, 4.0],
        auth=_empty_auth(),
    )

    def _build_runtime(*_args, **_kwargs):
        return SimpleNamespace(
            preprocessor=mock_preprocessor,
            sensor_config=mock_observation_bundle.sensor_config,
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
        build_preprocessor_runtime_fn=_build_runtime,
    )

    assert isinstance(plan.runtime_aoi, AOI)
    assert plan.runtime_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert isinstance(plan.output_writer, ConfiguredOutputWriter)
    assert plan.output_path == tmp_path / "products"


def test_build_execution_plan_accepts_geojson_dict_aoi(
    tmp_path,
    mock_observation_bundle,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
):
    input_path = tmp_path / "in.SAFE"
    input_path.mkdir()
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(
        config=cfg,
        input_path=input_path,
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

    def _build_runtime(*_args, **_kwargs):
        return SimpleNamespace(
            preprocessor=mock_preprocessor,
            sensor_config=mock_observation_bundle.sensor_config,
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
        build_preprocessor_runtime_fn=_build_runtime,
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
    input_path = tmp_path / "in.SAFE"
    input_path.mkdir()
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(
        config=cfg,
        input_path=input_path,
        aoi=[1.0, 2.0, 3.0, 4.0],
        auth=_empty_auth(),
    )

    def _build_runtime(*_args, **_kwargs):
        return SimpleNamespace(
            preprocessor=mock_preprocessor,
            sensor_config=mock_observation_bundle.sensor_config,
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
        build_preprocessor_runtime_fn=_build_runtime,
    )

    assert result.boa.identical(mock_observation_bundle.toa)
    assert result.cloud_mask.identical(mock_observation_bundle.cloud_mask)
    assert result.aot.identical(mock_solved_atmosphere.aot)
    assert result.tcwv.identical(mock_solved_atmosphere.tcwv)
    assert result.diagnostics.processing_time_s == pytest.approx(0.01)


def test_process_scene_show_progress_derives_report_path(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
):
    input_path = tmp_path / "in.SAFE"
    input_path.mkdir()
    cfg = SIACConfig(
        sensor="s2",
        runtime={
            "execution": {
                "show_progress": True,
                "profiling_sample_interval_s": 0.01,
                "progress_heartbeat_s": 0.01,
            }
        },
    )
    request = SceneProcessRequest(
        config=cfg,
        input_path=input_path,
        output_path=tmp_path / "products",
        auth=_empty_auth(),
    )

    def _skip_write_output(*_args, **_kwargs) -> None:
        return None

    monkeypatch.setattr(scene_workflow, "write_output", _skip_write_output)

    def _build_runtime(*_args, **_kwargs):
        return SimpleNamespace(
            preprocessor=mock_preprocessor,
            sensor_config=mock_observation_bundle.sensor_config,
        )

    process_scene(
        request,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
        build_preprocessor_runtime_fn=_build_runtime,
    )

    report_path = tmp_path / "products" / "reports" / "execution_profile.json"
    assert report_path.exists()
    summary = json.loads(report_path.read_text(encoding="utf-8"))
    assert summary["run"]["status"] == "success"

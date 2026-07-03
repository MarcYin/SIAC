"""Workflow tests for the real scene planning and execution path."""

from __future__ import annotations

import json
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import siac.workflows.scene as scene_workflow
from siac.adapters.auth import CredentialManager
from siac.adapters.output import ConfiguredOutputWriter
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.domain.aoi import AOI
from siac.runtime import CorrectionResult
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


def _correction_result_with_aod(value: float) -> CorrectionResult:
    arr = xr.DataArray(np.full((2, 2), value, dtype=np.float32), dims=("y", "x"))
    return CorrectionResult(
        boa=xr.Dataset(),
        boa_unc=None,
        aot=arr,
        tcwv=xr.zeros_like(arr),
        cloud_mask=xr.zeros_like(arr, dtype=bool),
        metadata={"base": "metadata"},
    )


def test_aod_quality_score_uses_solver_fit_uncertainty_and_spatial_terms():
    candidate = _correction_result_with_aod(0.4)
    stronger_candidate = _correction_result_with_aod(0.4)
    candidate = replace(
        candidate,
        metadata={
            "solver": {
                "cost_final": 3.0,
                "aot_unc_median": 0.08,
                "aot_std": 0.3,
                "aot_finite_fraction": 1.0,
                "converged": True,
            }
        },
    )
    stronger_candidate = replace(
        stronger_candidate,
        metadata={
            "solver": {
                "cost_final": 1.0,
                "aot_unc_median": 0.02,
                "aot_std": 0.05,
                "aot_finite_fraction": 1.0,
                "converged": True,
            }
        },
    )

    assert scene_workflow.aod_quality_score(stronger_candidate) < scene_workflow.aod_quality_score(
        candidate
    )


def test_aod_scene_mean_returns_finite_aod_mean():
    candidate = _correction_result_with_aod(0.4)
    values = candidate.aot.values.copy()
    values[0, 0] = np.nan
    candidate = replace(candidate, aot=candidate.aot.copy(data=values))

    assert scene_workflow.aod_scene_mean(candidate) == pytest.approx(0.4)


def test_aod_quality_score_prefers_per_band_cost_when_available():
    candidate = replace(
        _correction_result_with_aod(0.4),
        metadata={
            "solver": {
                "cost_final": 1.0,
                "cost_final_per_band": 0.5,
                "aot_unc_median": 0.02,
                "aot_std": 0.05,
                "aot_finite_fraction": 1.0,
                "converged": True,
            }
        },
    )
    raw_cost_candidate = replace(
        candidate,
        metadata={
            "solver": {
                "cost_final": 0.75,
                "aot_unc_median": 0.02,
                "aot_std": 0.05,
                "aot_finite_fraction": 1.0,
                "converged": True,
            }
        },
    )

    assert scene_workflow.aod_quality_score(candidate) < scene_workflow.aod_quality_score(
        raw_cost_candidate
    )


def test_process_scene_with_aod_selector_keeps_base_when_selector_rejects_fallback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[object] = []
    base_result = _correction_result_with_aod(0.4)
    fallback_result = _correction_result_with_aod(0.9)

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append(surface_prior_provider)
        return fallback_result if surface_prior_provider == "fallback-provider" else base_result

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_selector(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        surface_prior_provider="base-provider",
        fallback_surface_prior_provider="fallback-provider",
        selector=lambda _base, _fallback: False,
    )

    assert calls == ["base-provider", "fallback-provider"]
    assert result.aot.identical(base_result.aot)
    assert result.metadata["aod_selector"] == {
        "enabled": True,
        "selected": "base",
        "selector": "<lambda>",
    }


def test_process_scene_with_aod_selector_uses_fallback_request_and_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[tuple[SceneProcessRequest, object]] = []
    base_cfg = SIACConfig(sensor="s2")
    fallback_cfg = base_cfg.with_overrides(runtime={"log_level": "DEBUG"})
    base_request = SceneProcessRequest(config=base_cfg, input_path=tmp_path / "in.SAFE")
    fallback_request = SceneProcessRequest(config=fallback_cfg, input_path=tmp_path / "in.SAFE")
    base_result = _correction_result_with_aod(1.3)
    fallback_result = _correction_result_with_aod(0.9)

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append((_request, surface_prior_provider))
        return fallback_result if surface_prior_provider == "fallback-provider" else base_result

    def production_rule(_base: CorrectionResult, _fallback: CorrectionResult):
        return True, {"base_aod": 1.3, "fallback_aod": 0.9, "reason": "test"}

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_selector(
        base_request,
        fallback_request=fallback_request,
        surface_prior_provider="base-provider",
        fallback_surface_prior_provider="fallback-provider",
        selector=production_rule,
    )

    assert [provider for _, provider in calls] == ["base-provider", "fallback-provider"]
    assert calls[0][0].config is base_cfg
    assert calls[0][0].output_path is None
    assert calls[1][0].config is fallback_cfg
    assert calls[1][0].output_path is None
    assert result.aot.identical(fallback_result.aot)
    assert result.metadata["aod_selector"] == {
        "enabled": True,
        "selected": "fallback",
        "selector": "production_rule",
        "base_aod": 1.3,
        "fallback_aod": 0.9,
        "reason": "test",
    }


def test_process_scene_with_aod_ensemble_returns_fixed_mean(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[object] = []
    results = {
        "base-provider": _correction_result_with_aod(0.3),
        "pred-provider": _correction_result_with_aod(0.6),
        "floor03-provider": _correction_result_with_aod(0.9),
    }

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append(surface_prior_provider)
        return results[surface_prior_provider]

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_ensemble(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        candidate_surface_prior_providers={
            "base": "base-provider",
            "predictor": "pred-provider",
            "predictor_unc03": "floor03-provider",
        },
    )

    assert calls == ["base-provider", "pred-provider", "floor03-provider"]
    assert np.allclose(result.aot.values, 0.6)
    assert result.metadata["aod_ensemble"]["method"] == "mean"
    assert result.metadata["aod_ensemble"]["candidates"] == [
        "base",
        "predictor",
        "predictor_unc03",
    ]
    assert result.metadata["aod_ensemble"]["candidate_statistics"]["base"]["median"] == (
        pytest.approx(0.3)
    )


def test_process_scene_with_aod_ensemble_supports_power_mean(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    results = {
        "low": _correction_result_with_aod(0.2),
        "high": _correction_result_with_aod(0.8),
    }

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        return results[surface_prior_provider]

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_ensemble(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        candidate_surface_prior_providers={"low": "low", "high": "high"},
        ensemble_power=4.0,
    )

    expected = ((0.2**4 + 0.8**4) / 2.0) ** 0.25
    assert np.allclose(result.aot.values, expected)
    assert result.metadata["aod_ensemble"]["method"] == "power_mean"
    assert result.metadata["aod_ensemble"]["power"] == pytest.approx(4.0)


def test_process_scene_with_aod_ensemble_supports_fixed_weights(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    results = {
        "base": _correction_result_with_aod(0.2),
        "predictor": _correction_result_with_aod(0.8),
    }

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        return results[surface_prior_provider]

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_ensemble(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        candidate_surface_prior_providers={"base": "base", "predictor": "predictor"},
        candidate_weights={"base": 3.0, "predictor": 7.0},
    )

    assert np.allclose(result.aot.values, 0.62)
    assert result.metadata["aod_ensemble"]["weights"] == {
        "base": pytest.approx(0.3),
        "predictor": pytest.approx(0.7),
    }


def test_process_scene_with_aod_ensemble_supports_median(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    results = {
        "low": _correction_result_with_aod(0.2),
        "middle": _correction_result_with_aod(0.5),
        "high": _correction_result_with_aod(1.4),
    }

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        return results[surface_prior_provider]

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_ensemble(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        candidate_surface_prior_providers={"low": "low", "middle": "middle", "high": "high"},
        ensemble_method="median",
    )

    assert np.allclose(result.aot.values, 0.5)
    assert result.metadata["aod_ensemble"]["method"] == "median"
    assert result.metadata["aod_ensemble"]["weights"] is None


def test_process_scene_with_aod_ensemble_rejects_bad_weights(tmp_path):
    with pytest.raises(ValueError, match="missing entries"):
        scene_workflow.process_scene_with_aod_ensemble(
            SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
            candidate_surface_prior_providers={"base": "base", "predictor": "predictor"},
            candidate_weights={"base": 1.0},
        )
    with pytest.raises(ValueError, match="not supported"):
        scene_workflow.process_scene_with_aod_ensemble(
            SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
            candidate_surface_prior_providers={"base": "base", "predictor": "predictor"},
            candidate_weights={"base": 1.0, "predictor": 1.0},
            ensemble_method="median",
        )


def test_process_scene_with_aod_rail_fallback_keeps_base_below_threshold(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[object] = []
    base_result = _correction_result_with_aod(0.5)

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append(surface_prior_provider)
        return base_result

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_rail_fallback(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        surface_prior_provider="base-provider",
        fallback_surface_prior_provider="fallback-provider",
        aod_rail_threshold=1.7,
    )

    assert calls == ["base-provider"]
    assert result.aot.identical(base_result.aot)
    assert result.metadata["base"] == "metadata"
    assert result.metadata["aod_rail_fallback"] == {
        "enabled": True,
        "triggered": False,
        "threshold": 1.7,
        "statistic": "max",
        "statistic_value": pytest.approx(0.5),
        "source": "base",
    }


def test_process_scene_with_aod_rail_fallback_reruns_with_fallback_provider(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[object] = []
    base_result = _correction_result_with_aod(1.8)
    fallback_result = _correction_result_with_aod(1.4)

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append(surface_prior_provider)
        return fallback_result if surface_prior_provider == "fallback-provider" else base_result

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_rail_fallback(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        surface_prior_provider="base-provider",
        fallback_surface_prior_provider="fallback-provider",
        aod_rail_threshold=1.7,
    )

    assert calls == ["base-provider", "fallback-provider"]
    assert result.aot.identical(fallback_result.aot)
    assert result.metadata["aod_rail_fallback"]["triggered"] is True
    assert result.metadata["aod_rail_fallback"]["source"] == "fallback"
    assert result.metadata["aod_rail_fallback"]["statistic_value"] == pytest.approx(1.8)


def test_process_scene_with_aod_rail_fallback_accepts_callable_statistic(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
):
    calls: list[object] = []
    base_result = _correction_result_with_aod(0.4)
    fallback_result = _correction_result_with_aod(0.9)

    def _fake_process_scene(_request, *, surface_prior_provider=None, **_kwargs):  # noqa: ANN001
        calls.append(surface_prior_provider)
        return fallback_result if surface_prior_provider == "fallback-provider" else base_result

    def station_aod(_result: CorrectionResult) -> float:
        return 1.75

    monkeypatch.setattr(scene_workflow, "process_scene", _fake_process_scene)

    result = scene_workflow.process_scene_with_aod_rail_fallback(
        SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE"),
        surface_prior_provider="base-provider",
        fallback_surface_prior_provider="fallback-provider",
        aod_rail_threshold=1.7,
        aod_rail_statistic=station_aod,
    )

    assert calls == ["base-provider", "fallback-provider"]
    assert result.aot.identical(fallback_result.aot)
    assert result.metadata["aod_rail_fallback"]["statistic"] == "station_aod"
    assert result.metadata["aod_rail_fallback"]["statistic_value"] == pytest.approx(1.75)

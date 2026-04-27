from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import xarray as xr

import siac.workflows.scene as scene_mod
from siac.adapters.auth import CredentialManager
from siac.app.planning import (
    ExecutionPlan,
    build_execution_plan,
    coerce_aoi_spec,
    resolve_run_config,
)
from siac.app.requests import SceneProcessRequest
from siac.domain.aoi import AOI
from siac.workflows.scene import execute_plan, process_scene, write_output


def _purge_modules(*module_names: str) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)


def test_resolve_run_config_converts_string_output_path_to_path() -> None:
    captured: dict[str, object] = {}

    class _Config:
        def resolve(self, request):  # noqa: ANN001
            captured["request"] = request
            return "resolved"

    result = resolve_run_config(
        _Config(),
        input_path=Path("/tmp/input.SAFE"),
        output_path="/tmp/out",
        sensor="s2",
        aoi=[1.0, 2.0, 3.0, 4.0],
        s2_query="T31UDQ_20240101",
    )

    assert result == "resolved"
    assert captured["request"].input_path == Path("/tmp/input.SAFE")
    assert captured["request"].output_path == Path("/tmp/out")
    assert captured["request"].sensor == "s2"
    assert captured["request"].aoi == [1.0, 2.0, 3.0, 4.0]
    assert captured["request"].s2_query == Path("T31UDQ_20240101")


def test_coerce_aoi_spec_covers_supported_inputs_and_invalid_type(tmp_path: Path) -> None:
    existing_aoi = AOI.from_bounds((1.0, 2.0, 3.0, 4.0))
    assert coerce_aoi_spec(None) is None
    assert coerce_aoi_spec(existing_aoi) is existing_aoi

    bounds_aoi = coerce_aoi_spec([1.0, 2.0, 3.0, 4.0])
    assert isinstance(bounds_aoi, AOI)
    assert bounds_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))

    geojson_payload = {
        "type": "Polygon",
        "coordinates": [[[1.0, 2.0], [3.0, 2.0], [3.0, 4.0], [1.0, 4.0], [1.0, 2.0]]],
    }
    dict_aoi = coerce_aoi_spec(geojson_payload)
    assert isinstance(dict_aoi, AOI)
    assert dict_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))

    geojson_path = tmp_path / "aoi.geojson"
    geojson_path.write_text(
        '{"type":"Polygon","coordinates":[[[1,2],[3,2],[3,4],[1,4],[1,2]]]}',
        encoding="utf-8",
    )
    path_aoi = coerce_aoi_spec(geojson_path)
    assert isinstance(path_aoi, AOI)
    assert path_aoi.get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))

    with pytest.raises(ValueError, match="longitude bounds must look like degrees"):
        coerce_aoi_spec([500000.0, 4100000.0, 501000.0, 4101000.0])

    with pytest.raises(ValueError, match="Could not parse AOI"):
        coerce_aoi_spec(123)


def test_build_execution_plan_uses_auth_fallback_and_skips_output_writer_when_unset(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    auth_obj = CredentialManager()
    captured: dict[str, object] = {}
    resolved_config = SimpleNamespace(sensor="s2", aoi=None, run=SimpleNamespace(output_path=None))

    class _Config:
        sensor = "s2"
        aoi = None

        def resolve(self, _request):  # noqa: ANN001
            return resolved_config

    request = SceneProcessRequest(config=_Config(), input_path=tmp_path / "in.SAFE")

    def _fake_build_preprocessor_runtime(config, **kwargs):  # noqa: ANN001
        captured["build_preprocessor_runtime"] = {"config": config, **kwargs}
        default_aoi = kwargs["default_aoi_resolver"]
        dataset = xr.Dataset({"B02": xr.DataArray([[1.0]])})
        captured["default_aoi"] = default_aoi(dataset)
        return SimpleNamespace(preprocessor="preprocessor", sensor_config="sensor-config")

    monkeypatch.setattr(
        "siac.app.planning.CredentialManager.from_config",
        lambda _config: auth_obj,
    )

    plan = build_execution_plan(
        request,
        build_preprocessor_runtime_fn=_fake_build_preprocessor_runtime,
        resolve_atmo_provider_fn=lambda config, auth=None: ("atmo", config, auth),
        resolve_surface_prior_provider_fn=lambda config, auth=None: ("surface", config, auth),
        resolve_grid_assembler_fn=lambda: "grid",
        resolve_solver_fn=lambda config: ("solver", config),
        resolve_corrector_fn=lambda config: ("corrector", config),
        resolve_rt_model_fn=lambda config, auth=None, sensor_config=None: (
            "rt",
            config,
            auth,
            sensor_config,
        ),
        resolve_output_writer_fn=lambda config: ("writer", config),
        aoi_resolver=lambda toa: ("aoi", tuple(toa.data_vars)),
    )

    assert isinstance(plan, ExecutionPlan)
    assert plan.auth is auth_obj
    assert plan.runtime_aoi is None
    assert plan.output_writer is None
    assert captured["build_preprocessor_runtime"]["sensor"] == "s2"
    assert captured["default_aoi"] == ("aoi", ("B02",))
    assert plan.atmo_provider == ("atmo", resolved_config, auth_obj)
    assert plan.surface_prior_provider == ("surface", resolved_config, auth_obj)
    assert plan.grid_assembler == "grid"
    assert plan.solver == ("solver", resolved_config)
    assert plan.corrector == ("corrector", resolved_config)
    assert plan.rt_model == ("rt", resolved_config, auth_obj, "sensor-config")


def test_app_and_workflows_lazy_submodule_access_imports_modules() -> None:
    module_names = (
        "siac.app",
        "siac.app.planning",
        "siac.workflows",
        "siac.workflows.scene",
    )
    original_modules = {module_name: sys.modules.get(module_name) for module_name in module_names}

    try:
        _purge_modules(*module_names)

        app_pkg = importlib.import_module("siac.app")
        workflows_pkg = importlib.import_module("siac.workflows")

        assert "siac.app.planning" not in sys.modules
        assert "siac.workflows.scene" not in sys.modules

        planning_mod = app_pkg.planning
        scene_submodule = workflows_pkg.scene

        assert planning_mod.__name__ == "siac.app.planning"
        assert scene_submodule.__name__ == "siac.workflows.scene"
        assert "siac.app.planning" in sys.modules
        assert "siac.workflows.scene" in sys.modules
    finally:
        _purge_modules(*module_names)
        for module_name, module in original_modules.items():
            if module is not None:
                sys.modules[module_name] = module


def test_write_output_and_execute_plan_delegate_to_writer_and_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    writer_calls: dict[str, object] = {}

    class _Writer:
        def write(self, result, path):  # noqa: ANN001
            writer_calls["result"] = result
            writer_calls["path"] = path

    write_output("result", tmp_path / "out", output_writer=_Writer())
    assert writer_calls == {"result": "result", "path": tmp_path / "out"}

    pipeline_calls: dict[str, object] = {}
    monkeypatch.setitem(
        scene_mod.execute_plan.__globals__,
        "run_pipeline",
        lambda *args, **kwargs: pipeline_calls.setdefault("call", (args, kwargs)) or "ignored",
    )
    plan = SimpleNamespace(
        input_path=tmp_path / "in.SAFE",
        runtime_aoi="aoi",
        config="config",
        preprocessor="preprocessor",
        atmo_provider="atmo",
        surface_prior_provider="surface",
        grid_assembler="grid",
        solver="solver",
        corrector="corrector",
        rt_model="rt",
    )

    execute_plan(plan)
    args, kwargs = pipeline_calls["call"]
    assert args == (tmp_path / "in.SAFE", "aoi", "config")
    assert kwargs["preprocessor"] == "preprocessor"
    assert kwargs["atmo_provider"] == "atmo"
    assert kwargs["surface_prior_provider"] == "surface"
    assert kwargs["grid_assembler"] == "grid"
    assert kwargs["solver"] == "solver"
    assert kwargs["corrector"] == "corrector"
    assert kwargs["rt_model"] == "rt"


def test_process_scene_only_writes_when_output_path_and_writer_are_present(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    write_calls: list[tuple[object, object, object]] = []
    request = SceneProcessRequest(config=SimpleNamespace(sensor="s2", aoi=None), input_path="input")
    plan = SimpleNamespace(config=request.config, output_path="out", output_writer="writer")
    plan_without_writer = SimpleNamespace(output_path="out", output_writer=None)
    results = iter(["result-1", "result-2"])

    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "build_execution_plan",
        lambda _request, **kwargs: plan if kwargs["preprocessor"] == "pre" else plan_without_writer,
    )

    def _fake_execute_plan(_built_plan, execution=None):  # noqa: ANN001
        _ = execution
        return next(results)

    def _fake_resolve_execution_settings(_config, execution=None, max_workers=None):  # noqa: ANN001
        _ = execution, max_workers
        return {
            "show_progress": False,
            "performance_report_path": None,
        }

    monkeypatch.setitem(scene_mod.process_scene.__globals__, "execute_plan", _fake_execute_plan)
    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "_resolve_execution_settings",
        _fake_resolve_execution_settings,
    )
    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "write_output",
        lambda result, output_path, *, output_writer: write_calls.append(
            (result, output_path, output_writer)
        ),
    )

    result_written = process_scene(request, preprocessor="pre")
    result_skipped = process_scene(request, preprocessor="skip")

    assert result_written == "result-1"
    assert result_skipped == "result-2"
    assert write_calls == [("result-1", "out", "writer")]


def test_process_scene_skips_post_write_when_pipeline_streamed_outputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    write_calls: list[tuple[object, object, object]] = []
    request = SceneProcessRequest(config=SimpleNamespace(sensor="s2", aoi=None), input_path="input")
    streamed_result = SimpleNamespace(metadata={"_siac_outputs_written": True})
    plan = SimpleNamespace(config=request.config, output_path="out", output_writer="writer")

    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "build_execution_plan",
        lambda *_a, **_kw: plan,
    )
    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "execute_plan",
        lambda *_a, **_kw: streamed_result,
    )
    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "_resolve_execution_settings",
        lambda *_a, **_kw: {"show_progress": False, "performance_report_path": None},
    )
    monkeypatch.setitem(
        scene_mod.process_scene.__globals__,
        "write_output",
        lambda result, output_path, *, output_writer: write_calls.append(
            (result, output_path, output_writer)
        ),
    )

    result = process_scene(request, preprocessor="pre")

    assert result is streamed_result
    assert write_calls == []

from __future__ import annotations

import importlib
import sys
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.app.planning import build_execution_plan, coerce_aoi_spec
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.domain.aoi import AOI
from siac.workflows.scene import write_output

if TYPE_CHECKING:
    from pathlib import Path


def _purge_modules(*module_names: str) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)


@pytest.mark.parametrize(
    ("package_name", "submodule_name"),
    [
        ("siac.app", "planning"),
        ("siac.workflows", "scene"),
    ],
)
def test_package_submodules_are_lazy_and_listed_in_dir(
    package_name: str,
    submodule_name: str,
) -> None:
    module_name = f"{package_name}.{submodule_name}"
    original_package = sys.modules.get(package_name)
    original_submodule = sys.modules.get(module_name)

    try:
        _purge_modules(package_name, module_name)

        package = importlib.import_module(package_name)

        assert module_name not in sys.modules
        assert "__getattr__" in dir(package)
        assert submodule_name in dir(package)
        assert submodule_name not in package.__all__

        submodule = getattr(package, submodule_name)

        assert submodule is sys.modules[module_name]
    finally:
        _purge_modules(package_name, module_name)
        if original_package is not None:
            sys.modules[package_name] = original_package
        if original_submodule is not None:
            sys.modules[module_name] = original_submodule


def test_coerce_aoi_spec_returns_existing_aoi_and_rejects_unknown() -> None:
    aoi = AOI.from_bounds((1.0, 2.0, 3.0, 4.0))

    assert coerce_aoi_spec(aoi) is aoi

    with pytest.raises(ValueError, match="Could not parse AOI from int"):
        coerce_aoi_spec(5)


def test_build_execution_plan_default_aoi_resolver_uses_callable_override(tmp_path: Path) -> None:
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(config=cfg, input_path=tmp_path / "in.SAFE")
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"])}
    )
    expected_aoi = AOI.from_bounds((5.0, 6.0, 7.0, 8.0))

    def _fake_build_preprocessor_runtime(config, *, input_path, sensor, default_aoi_resolver):  # noqa: ANN001
        assert config.sensor == "s2"
        assert input_path == tmp_path / "in.SAFE"
        assert sensor == "s2"

        def _preprocessor(_path: Path) -> SimpleNamespace:
            resolved = default_aoi_resolver(toa)
            return SimpleNamespace(aoi=resolved)

        return SimpleNamespace(preprocessor=_preprocessor, sensor_config="sensor-config")

    plan = build_execution_plan(
        request,
        build_preprocessor_runtime_fn=_fake_build_preprocessor_runtime,
        resolve_atmo_provider_fn=lambda _config, auth=None: ("atmo", auth),
        resolve_surface_prior_provider_fn=lambda _config, auth=None: ("surface", auth),
        resolve_grid_assembler_fn=lambda: "grid",
        resolve_solver_fn=lambda _config: "solver",
        resolve_corrector_fn=lambda _config: "corrector",
        resolve_rt_model_fn=lambda _config, auth=None, sensor_config=None: ("rt", auth, sensor_config),
        resolve_output_writer_fn=lambda _config: "writer",
        aoi_resolver=lambda _toa: expected_aoi,
    )

    result = plan.preprocessor(tmp_path / "in.SAFE")

    assert result.aoi is expected_aoi


def test_build_execution_plan_default_aoi_resolver_falls_back_to_raster(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cfg = SIACConfig(sensor="s2")
    request = SceneProcessRequest(config=cfg, input_path=tmp_path / "in.SAFE")
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"])}
    )
    expected_aoi = AOI.from_bounds((9.0, 10.0, 11.0, 12.0))

    def _fake_build_preprocessor_runtime(config, *, input_path, sensor, default_aoi_resolver):  # noqa: ANN001
        assert config.sensor == "s2"
        assert input_path == tmp_path / "in.SAFE"
        assert sensor == "s2"

        def _preprocessor(_path: Path) -> SimpleNamespace:
            resolved = default_aoi_resolver(toa)
            return SimpleNamespace(aoi=resolved)

        return SimpleNamespace(preprocessor=_preprocessor, sensor_config="sensor-config")

    monkeypatch.setattr("siac.app.planning.AOI.from_raster", lambda _raster: expected_aoi)

    plan = build_execution_plan(
        request,
        build_preprocessor_runtime_fn=_fake_build_preprocessor_runtime,
        resolve_atmo_provider_fn=lambda _config, auth=None: ("atmo", auth),
        resolve_surface_prior_provider_fn=lambda _config, auth=None: ("surface", auth),
        resolve_grid_assembler_fn=lambda: "grid",
        resolve_solver_fn=lambda _config: "solver",
        resolve_corrector_fn=lambda _config: "corrector",
        resolve_rt_model_fn=lambda _config, auth=None, sensor_config=None: ("rt", auth, sensor_config),
        resolve_output_writer_fn=lambda _config: "writer",
    )

    result = plan.preprocessor(tmp_path / "in.SAFE")

    assert result.aoi is expected_aoi


def test_write_output_passes_a_path_to_output_writer(tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    class _FakeWriter:
        def write(self, result, output_path):  # noqa: ANN001
            captured["result"] = result
            captured["output_path"] = output_path

    result = object()
    write_output(result, str(tmp_path / "out"), output_writer=_FakeWriter())

    assert captured["result"] is result
    assert captured["output_path"] == tmp_path / "out"


def test_process_scene_writes_outputs_when_plan_has_writer(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import siac.workflows.scene as scene_mod

    request = SceneProcessRequest(config=SIACConfig(sensor="s2"), input_path=tmp_path / "in.SAFE")
    result = object()
    plan = SimpleNamespace(
        config=request.config,
        output_path=tmp_path / "products",
        output_writer="writer",
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(scene_mod, "build_execution_plan", lambda *_args, **_kwargs: plan)
    monkeypatch.setattr(
        scene_mod,
        "_resolve_execution_settings",
        lambda _config, **_kwargs: {
            "show_progress": False,
            "performance_report_path": None,
        },
    )
    monkeypatch.setattr(scene_mod, "execute_plan", lambda _plan, **_kwargs: result)

    def _fake_write_output(result_arg, output_path, *, output_writer):  # noqa: ANN001
        captured["result"] = result_arg
        captured["output_path"] = output_path
        captured["output_writer"] = output_writer

    monkeypatch.setattr(scene_mod, "write_output", _fake_write_output)

    out = scene_mod.process_scene(request)

    assert out is result
    assert captured == {
        "result": result,
        "output_path": tmp_path / "products",
        "output_writer": "writer",
    }

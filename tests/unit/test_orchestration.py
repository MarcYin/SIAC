"""Public-orchestration tests for the greenfield API layer."""

from __future__ import annotations

import numpy as np
import xarray as xr

from siac.api.public import SIAC, process_landsat8, process_sentinel2, siac_process
from siac.api.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.runtime import CorrectionResult


def _da(value: float, shape: tuple[int, int] = (2, 2)) -> xr.DataArray:
    return xr.DataArray(np.full(shape, value, dtype=np.float32), dims=["y", "x"])


def _result() -> CorrectionResult:
    return CorrectionResult(
        boa=xr.Dataset({"B02": _da(0.2), "B03": _da(0.25)}),
        boa_unc=None,
        aot=_da(0.2),
        tcwv=_da(2.3),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        metadata={"ok": True},
    )


def _empty_auth():
    from siac.adapters.auth import CredentialManager

    return CredentialManager()


def test_siac_from_helpers(monkeypatch):
    cfg = SIACConfig(sensor="s2")
    monkeypatch.setattr(
        "siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth()
    )
    monkeypatch.setattr(
        "siac.api.public.SIACConfig.from_file", classmethod(lambda _cls, _path: cfg)
    )

    siac_from_file = SIAC.from_config("dummy.toml")
    siac_default = SIAC.from_defaults(sensor="l8")

    assert isinstance(siac_from_file, SIAC)
    assert siac_from_file.config is cfg
    assert siac_default.config.sensor == "l8"


def test_siac_process_delegates_with_scene_request(monkeypatch, tmp_path):
    cfg = SIACConfig(sensor="s2")
    monkeypatch.setattr(
        "siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth()
    )

    captured: dict[str, object] = {}

    def _fake_process_scene(request):  # noqa: ANN001
        captured["request"] = request
        return _result()

    monkeypatch.setattr("siac.api.public.workflow_process_scene", _fake_process_scene)

    siac_obj = SIAC(cfg)
    result = siac_obj.process(
        tmp_path / "in.SAFE", output_path=tmp_path / "out", aoi=[1.0, 2.0, 3.0, 4.0]
    )

    assert result.metadata["ok"] is True
    request = captured["request"]
    assert isinstance(request, SceneProcessRequest)
    assert request.config is cfg
    assert request.input_path == tmp_path / "in.SAFE"
    assert request.output_path == tmp_path / "out"
    assert request.aoi == [1.0, 2.0, 3.0, 4.0]


def test_siac_process_helper_delegates_with_scene_request(monkeypatch, tmp_path):
    cfg = SIACConfig(sensor="s2")
    captured: dict[str, object] = {}

    def _fake_process_scene(request):  # noqa: ANN001
        captured["request"] = request
        return "pipeline-result"

    monkeypatch.setattr("siac.api.public.workflow_process_scene", _fake_process_scene)

    out = siac_process(cfg, tmp_path / "in.SAFE", aoi="/tmp/aoi.geojson")
    assert out == "pipeline-result"
    request = captured["request"]
    assert isinstance(request, SceneProcessRequest)
    assert request.input_path == tmp_path / "in.SAFE"
    assert request.aoi == "/tmp/aoi.geojson"


def test_sensor_wrappers_delegate(monkeypatch):
    calls: list[object] = []

    class _Runner:
        def process(self, input_path, output_path=None, **kwargs):  # noqa: ANN001
            calls.append((input_path, output_path, kwargs))
            return "ok"

    monkeypatch.setattr(
        "siac.api.public.SIAC.from_defaults", lambda sensor: calls.append(sensor) or _Runner()
    )

    assert process_sentinel2("in1", "out1", aoi="a.geojson") == "ok"
    assert process_landsat8("in2", "out2") == "ok"
    assert calls[0] == "s2"
    assert calls[2] == "l8"
    assert calls[1][2]["aoi"] == "a.geojson"

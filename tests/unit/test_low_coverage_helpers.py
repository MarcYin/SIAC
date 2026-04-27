"""Targeted coverage lifts for low-coverage helper modules."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from siac.adapters.rt import build_rt_model
from siac.algorithms.rt.lut import srf_kernel
from siac.api.requests import SceneProcessRequest, Sentinel2ProcessRequest, Sentinel2ResolveRequest
from siac.config import RunRequest, SIACConfig
from siac.config.load import (
    CONFIG_PATH_ENV,
    load_system_config_from_default,
    write_default_system_config,
)
from siac.config.public import get_default_config
from siac.config.schema import RTSetupConfig
from siac.config.snapshot import (
    _normalize,
    snapshot_resolved_config,
    snapshot_system_config,
    write_runtime_snapshot,
)
from siac.workflows.sentinel2 import process_s2


def test_process_s2_resolves_config_auth_and_delegates(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    config = SIACConfig(sensor="auto", aoi={"type": "Point", "coordinates": [0.0, 0.0]})
    request = Sentinel2ProcessRequest(
        config=config,
        query="T31UDQ_20240101",
        output_path=tmp_path / "out",
    )
    resolved_config = SimpleNamespace(token="resolved")
    auth_obj = object()
    resolved_input = tmp_path / "scene.SAFE"
    resolved_input.mkdir()
    captured: dict[str, object] = {}

    def _fake_resolve_run_config(config_arg, **kwargs):  # noqa: ANN001
        captured["resolve_run_config"] = {"config": config_arg, **kwargs}
        return resolved_config

    def _fake_from_config(config_arg):  # noqa: ANN001
        captured["credential_manager_config"] = config_arg
        return auth_obj

    def _fake_resolve_s2_input(request_arg):  # noqa: ANN001
        captured["resolve_s2_input_request"] = request_arg
        return resolved_input

    def _fake_process_scene(*, request):  # noqa: ANN001
        captured["process_scene_request"] = request
        return "done"

    monkeypatch.setitem(process_s2.__globals__, "resolve_run_config", _fake_resolve_run_config)
    monkeypatch.setitem(
        process_s2.__globals__, "CredentialManager", SimpleNamespace(from_config=_fake_from_config)
    )
    monkeypatch.setitem(process_s2.__globals__, "process_scene", _fake_process_scene)

    out = process_s2(request, resolve_s2_input_fn=_fake_resolve_s2_input)

    assert out == "done"
    assert captured["resolve_run_config"]["sensor"] == "s2"
    assert captured["resolve_run_config"]["aoi"] == config.aoi
    assert captured["resolve_run_config"]["s2_query"] == "T31UDQ_20240101"
    assert captured["credential_manager_config"] is resolved_config
    resolve_request = captured["resolve_s2_input_request"]
    assert isinstance(resolve_request, Sentinel2ResolveRequest)
    assert resolve_request.config is config
    assert resolve_request.query == "T31UDQ_20240101"
    assert resolve_request.auth is auth_obj
    scene_request = captured["process_scene_request"]
    assert isinstance(scene_request, SceneProcessRequest)
    assert scene_request.config is config
    assert scene_request.input_path == resolved_input
    assert scene_request.output_path == tmp_path / "out"
    assert scene_request.auth is auth_obj


def test_process_s2_uses_default_scene_resolver(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    config = SIACConfig(sensor="s2")
    auth_obj = object()
    request = Sentinel2ProcessRequest(
        config=config,
        query="T31UDQ_20240101",
        output_path=tmp_path / "out",
        auth=auth_obj,
    )
    resolved_config = SimpleNamespace(token="resolved")
    resolved_input = tmp_path / "scene.SAFE"
    resolved_input.mkdir()
    captured: dict[str, object] = {}

    def _fake_resolve_run_config(config_arg, **kwargs):  # noqa: ANN001
        captured["resolve_run_config"] = {"config": config_arg, **kwargs}
        return resolved_config

    def _fake_resolve_s2_scene_input(request_arg):  # noqa: ANN001
        captured["resolve_s2_scene_input_request"] = request_arg
        return resolved_input

    def _fake_process_scene(*, request):  # noqa: ANN001
        captured["process_scene_request"] = request
        return "done"

    monkeypatch.setitem(process_s2.__globals__, "resolve_run_config", _fake_resolve_run_config)
    monkeypatch.setitem(
        process_s2.__globals__, "resolve_s2_scene_input", _fake_resolve_s2_scene_input
    )
    monkeypatch.setitem(process_s2.__globals__, "process_scene", _fake_process_scene)

    out = process_s2(request)

    assert out == "done"
    assert captured["resolve_run_config"]["sensor"] == "s2"
    assert "credential_manager_config" not in captured
    resolve_request = captured["resolve_s2_scene_input_request"]
    assert isinstance(resolve_request, Sentinel2ResolveRequest)
    assert resolve_request.auth is auth_obj
    assert captured["process_scene_request"].input_path == resolved_input


def test_snapshot_helpers_cover_redaction_and_path_states(tmp_path: Path) -> None:
    dem = tmp_path / "dem.tif"
    dem.write_text("dem")
    emulator_dir = tmp_path / "emulator"
    emulator_dir.mkdir()
    lut_path = "s3://bucket/lut.zarr"
    output_dir = tmp_path / "products"
    output_dir.mkdir()

    config = SIACConfig(
        auth={
            "cdse": {"username": "alice", "password": "secret"},
            "aws": {"access_key_id": "AKIA", "secret_access_key": "shh"},
            "gcs": {"credentials_file": tmp_path / "gcs.json"},
        },
        paths={
            "dem": dem,
            "emulator_dir": emulator_dir,
            "lut_path": lut_path,
            "cache_root": tmp_path / "cache",
        },
        output={"defaults": {"output_dir": output_dir}},
    )
    resolved = config.resolve(RunRequest(input_path=tmp_path / "scene.SAFE"))

    system_snapshot = snapshot_system_config(config, source_path=Path("~/siac.toml"))
    resolved_snapshot = snapshot_resolved_config(resolved, source_path=tmp_path / "runtime.toml")

    assert system_snapshot["source_path"].endswith("siac.toml")
    assert system_snapshot["config"]["auth"]["cdse"]["username"] == "<redacted>"
    assert system_snapshot["config"]["auth"]["aws"]["secret_access_key"] == "<redacted>"
    assert system_snapshot["config"]["auth"]["gcs"]["credentials_file"] == "<redacted>"
    assert resolved_snapshot["path_states"]["paths.dem"]["kind"] == "file"
    assert resolved_snapshot["path_states"]["paths.emulator_dir"]["kind"] == "directory"
    assert resolved_snapshot["path_states"]["paths.lut_path"]["kind"] == "remote"
    assert resolved_snapshot["path_states"]["output.defaults.output_dir"]["kind"] == "directory"


def test_snapshot_helpers_without_redaction_and_normalization(tmp_path: Path) -> None:
    config = SIACConfig(
        auth={"cdse": {"username": "alice", "password": "secret"}},
        paths={"lut_path": "s3://bucket/lut.zarr"},
    )
    resolved = config.resolve(RunRequest(input_path=tmp_path / "scene.SAFE"))

    system_snapshot = snapshot_system_config(config, redact_secrets=False)
    resolved_snapshot = snapshot_resolved_config(resolved, redact_secrets=False)

    assert system_snapshot["config"]["auth"]["cdse"]["username"] == "alice"
    assert resolved_snapshot["config"]["auth"]["cdse"]["username"] == "alice"
    assert _normalize((Path("/tmp/demo"), ["a", {"b": Path("/tmp/path")}])) == [
        "/tmp/demo",
        ["a", {"b": "/tmp/path"}],
    ]


def test_write_runtime_snapshot_serializes_state(tmp_path: Path) -> None:
    config = SIACConfig(output={"defaults": {"output_dir": tmp_path / "products"}})
    resolved = config.resolve(RunRequest(input_path=tmp_path / "scene.SAFE"))
    snapshot_path = tmp_path / "snapshot" / "runtime.yaml"

    write_runtime_snapshot(resolved, snapshot_path, source_path=tmp_path / "config.toml")

    assert snapshot_path.exists()
    import yaml

    content = yaml.safe_load(snapshot_path.read_text(encoding="utf-8"))
    assert content["state"]["source_path"].endswith("config.toml")
    assert "generated_at" in content
    assert content["state"]["config"]["run"]["input_path"].endswith("scene.SAFE")


def test_config_wrapper_helpers_and_loaders(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    config = SIACConfig(sensor="s2")
    toml_path = tmp_path / "config.toml"
    config.to_toml(toml_path)

    loaded = SIACConfig.load(path=toml_path, sensor="l8")
    assert loaded.sensor == "l8"

    monkeypatch.setattr("siac.config.public.load_system_config_from_default", lambda: config)
    default_loaded = SIACConfig.load()
    assert default_loaded.sensor == "s2"

    assert SIACConfig.from_toml(toml_path).sensor == "auto"
    assert get_default_config().sensor == "auto"

    with pytest.raises(ValueError, match="TOML"):
        SIACConfig.from_yaml(toml_path)
    with pytest.raises(ValueError, match="TOML"):
        config.to_yaml(toml_path)


def test_load_helpers_cover_default_path_and_suffixless_output(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "config.toml"
    SIACConfig().to_toml(config_path)
    monkeypatch.setenv(CONFIG_PATH_ENV, str(config_path))

    loaded = load_system_config_from_default()
    assert loaded.paths.lut_path is not None

    output_path = tmp_path / "default-config"
    written = write_default_system_config(output_path)
    assert written.suffix == ".toml"
    assert written.exists()


def test_write_state_snapshot_delegates_to_runtime_snapshot(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    config = SIACConfig(sensor="l8", aoi={"kind": "bbox"})
    snapshot_path = tmp_path / "runtime.yaml"
    captured: dict[str, object] = {}

    resolved_config = SimpleNamespace(token="resolved")

    def _fake_resolve_config(config_arg, request_arg):  # noqa: ANN001
        captured["resolve_config"] = {"config": config_arg, "request": request_arg}
        return resolved_config

    def _fake_write_runtime_snapshot(config_arg, path_arg, *, redact_secrets=True):  # noqa: ANN001
        captured["write_runtime_snapshot"] = {
            "config": config_arg,
            "path": path_arg,
            "redact_secrets": redact_secrets,
        }

    monkeypatch.setattr("siac.config.public.resolve_config", _fake_resolve_config)
    monkeypatch.setattr("siac.config.public.write_runtime_snapshot", _fake_write_runtime_snapshot)

    config.write_state_snapshot(snapshot_path, redact_secrets=False)

    resolved_input = captured["resolve_config"]["config"]
    assert isinstance(resolved_input, SIACConfig)
    assert resolved_input.sensor == config.sensor
    assert resolved_input.aoi == config.aoi

    request = captured["resolve_config"]["request"]
    assert request.sensor == "l8"
    assert request.aoi == {"kind": "bbox"}
    assert captured["write_runtime_snapshot"]["config"] is resolved_config
    assert captured["write_runtime_snapshot"]["path"] == snapshot_path
    assert captured["write_runtime_snapshot"]["redact_secrets"] is False


def test_build_rt_model_covers_emulator_lut_and_error_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeEmulator:
        instances: list[dict[str, object | None]] = []

        def __init__(self, emulator_dir=None, sensor_id=None, satellite_id=None):  # noqa: ANN001
            self.kwargs = {
                "emulator_dir": emulator_dir,
                "sensor_id": sensor_id,
                "satellite_id": satellite_id,
            }
            _FakeEmulator.instances.append(self.kwargs)

    class _FakeLUT:
        instances: list[dict[str, object]] = []

        def __init__(self, lut_path, interpolation_method, storage_options, rt_setup=None):  # noqa: ANN001
            self.kwargs = {
                "lut_path": lut_path,
                "interpolation_method": interpolation_method,
                "storage_options": storage_options,
                "rt_setup": rt_setup,
            }
            _FakeLUT.instances.append(self.kwargs)

    class _FakeSixS:
        instances: list[dict[str, object | None]] = []

        def __init__(self, *, sixs_config=None, sensor_config=None, rt_setup=None, runner=None):  # noqa: ANN001
            self.kwargs = {
                "sixs_config": sixs_config,
                "sensor_config": sensor_config,
                "rt_setup": rt_setup,
                "runner": runner,
            }
            _FakeSixS.instances.append(self.kwargs)

    class _FakeAws:
        def has_credentials(self) -> bool:
            return True

        def storage_options(self) -> dict[str, str]:
            return {"key": "abc", "secret": "def"}

    class _FakeAuth:
        def aws(self) -> _FakeAws:
            return _FakeAws()

    monkeypatch.setitem(build_rt_model.__globals__, "TwoLayerNNEmulator", _FakeEmulator)
    monkeypatch.setitem(build_rt_model.__globals__, "ZarrLUTBackend", _FakeLUT)
    monkeypatch.setitem(build_rt_model.__globals__, "SixSBackend", _FakeSixS)

    sensor_config = SimpleNamespace(sensor_id="MSI", satellite_id="S2C")
    emulator_config = SimpleNamespace(
        sensor="auto",
        rt_model=SimpleNamespace(
            backend="emulator",
            fallback_to_lut=True,
            lut_interpolation="linear",
            lut_storage_options={},
        ),
        paths=SimpleNamespace(emulator_dir=Path("/tmp/emu"), lut_path="s3://bucket/lut.zarr"),
    )

    emulator = build_rt_model(emulator_config, sensor_config=sensor_config)
    assert emulator.kwargs["sensor_id"] == "MSI"
    assert emulator.kwargs["satellite_id"] == "S2C"
    assert emulator.kwargs["emulator_dir"] == Path("/tmp/emu")

    fallback_config = SimpleNamespace(
        sensor="l8",
        rt_model=SimpleNamespace(
            backend="emulator",
            fallback_to_lut=True,
            lut_interpolation="nearest",
            lut_storage_options={"anon": True},
        ),
        paths=SimpleNamespace(lut_path="s3://bucket/lut.zarr"),
    )

    def _boom(*args, **kwargs):  # noqa: ANN001
        raise RuntimeError("missing emulator")

    monkeypatch.setitem(build_rt_model.__globals__, "TwoLayerNNEmulator", _boom)

    lut = build_rt_model(fallback_config, auth=_FakeAuth())
    assert lut.kwargs["lut_path"] == "s3://bucket/lut.zarr"
    assert lut.kwargs["interpolation_method"] == "nearest"
    assert lut.kwargs["storage_options"] == {"anon": True, "key": "abc", "secret": "def"}
    assert lut.kwargs["rt_setup"].atmosphere.profile == "us_standard_62"
    assert lut.kwargs["rt_setup"].aerosol.profile == "continental_average"
    assert lut.kwargs["rt_setup"].surface.mode == "homogeneous_lambertian"

    sixs_config = SimpleNamespace(
        sensor="auto",
        rt_model=SimpleNamespace(
            backend="sixs",
            setup=RTSetupConfig(
                atmosphere={"profile": "tropical", "columns_mode": "input_columns"},
                aerosol={"profile": "continental"},
            ),
            sixs=SimpleNamespace(
                mode="direct",
                output_variables=("xap", "xbp", "xcp", "tgasm"),
            ),
            fallback_to_lut=False,
            lut_interpolation="linear",
            lut_storage_options={},
        ),
        paths=SimpleNamespace(lut_path=None, emulator_dir=None),
    )

    sixs = build_rt_model(sixs_config, sensor_config=sensor_config)
    assert sixs.kwargs["sensor_config"] is sensor_config
    assert sixs.kwargs["sixs_config"].output_variables == ("xap", "xbp", "xcp", "tgasm")
    assert sixs.kwargs["rt_setup"].atmosphere.profile == "tropical"
    assert sixs.kwargs["rt_setup"].aerosol.profile == "continental"

    no_fallback_config = SimpleNamespace(
        sensor="l8",
        rt_model=SimpleNamespace(
            backend="emulator",
            fallback_to_lut=False,
            lut_interpolation="nearest",
            lut_storage_options={},
        ),
        paths=SimpleNamespace(lut_path=None),
    )

    monkeypatch.setattr("siac.adapters.rt.TwoLayerNNEmulator", _boom)
    with pytest.raises(ValueError, match="fallback is unavailable"):
        build_rt_model(no_fallback_config)

    unknown_config = SimpleNamespace(
        sensor="auto",
        rt_model=SimpleNamespace(
            backend="bogus",
            fallback_to_lut=False,
            lut_interpolation="linear",
            lut_storage_options={},
        ),
        paths=SimpleNamespace(lut_path=None, emulator_dir=None),
    )

    with pytest.raises(ValueError, match="Cannot resolve RT model"):
        build_rt_model(unknown_config)


def test_srf_kernel_shim_reexports() -> None:
    from siac.algorithms.rt.lut.rsrf_kernel import build_aligned_rsrf_kernel

    assert srf_kernel.build_aligned_rsrf_kernel is build_aligned_rsrf_kernel
    assert srf_kernel.build_aligned_srf_kernel is build_aligned_rsrf_kernel
    assert srf_kernel.AlignedSRFKernel is srf_kernel.AlignedRSRFKernel

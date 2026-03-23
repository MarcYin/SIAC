"""Unit tests for the SIAC CLI surface."""

from __future__ import annotations

import logging
from types import SimpleNamespace
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    import pytest

import siac.cli as cli
from siac.config import SIACConfig


def test_version_flag_prints_version(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = cli.main(["--version"])

    captured = capsys.readouterr()
    assert exit_code == 0
    assert captured.out.strip() == f"siac {cli.__version__}"


def test_process_s2_delegates_with_parsed_arguments(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    config = SIACConfig(sensor="s2")
    captured: dict[str, object] = {}

    def _fake_from_file(cls, _path: Path) -> SIACConfig:
        _ = cls
        captured["config_path"] = _path
        return config

    def _fake_process(config_obj, query, **kwargs):  # noqa: ANN001
        captured["config_obj"] = config_obj
        captured["query"] = query
        captured["kwargs"] = kwargs
        return SimpleNamespace(aot=SimpleNamespace(mean=lambda: 0.1234))

    monkeypatch.setattr(cli.SIACConfig, "from_file", classmethod(_fake_from_file))
    monkeypatch.setattr("siac.api.public.siac_process_s2", _fake_process)

    config_path = tmp_path / "config.toml"
    exit_code = cli.main(
        [
            "process-s2",
            "T31UDQ_20240101",
            "--config",
            str(config_path),
            "--output-path",
            str(tmp_path / "out"),
            "--aoi-bbox",
            "1",
            "2",
            "3",
            "4",
        ]
    )

    captured_out = capsys.readouterr()
    assert exit_code == 0
    assert captured["config_path"] == config_path
    assert captured["config_obj"] is config
    assert captured["query"] == "T31UDQ_20240101"
    assert captured["kwargs"]["output_path"] == tmp_path / "out"
    assert captured["kwargs"]["aoi"] == (1.0, 2.0, 3.0, 4.0)
    assert "Sentinel-2 processing complete" in captured_out.out


def test_process_s2_reports_errors_non_zero(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def _boom(config_obj, query, **kwargs):  # noqa: ANN001
        _ = (config_obj, query, kwargs)
        raise ValueError("boom")

    monkeypatch.setattr("siac.api.public.siac_process_s2", _boom)

    exit_code = cli.main(["process-s2", "T31UDQ_20240101"])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "siac: error: boom" in captured.err


def test_process_s2_honors_runtime_log_level_per_invocation(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    configs = {
        "quiet.toml": SIACConfig(sensor="s2", runtime={"log_level": "ERROR"}),
        "info.toml": SIACConfig(sensor="s2", runtime={"log_level": "INFO"}),
    }

    def _fake_from_file(cls, path: Path) -> SIACConfig:
        _ = cls
        return configs[path.name]

    def _fake_process(config_obj, query, **kwargs):  # noqa: ANN001
        _ = (config_obj, query, kwargs)
        logging.getLogger("siac.test.workflow").info("workflow info")
        return SimpleNamespace(aot=SimpleNamespace(mean=lambda: 0.1234))

    monkeypatch.setattr(cli.SIACConfig, "from_file", classmethod(_fake_from_file))
    monkeypatch.setattr("siac.api.public.siac_process_s2", _fake_process)

    quiet_code = cli.main(["process-s2", "T31UDQ_20240101", "--config", str(tmp_path / "quiet.toml")])
    quiet_output = capsys.readouterr()

    info_code = cli.main(["process-s2", "T31UDQ_20240101", "--config", str(tmp_path / "info.toml")])
    info_output = capsys.readouterr()

    assert quiet_code == 0
    assert "workflow info" not in quiet_output.err
    assert "Sentinel-2 processing complete" in quiet_output.out
    assert info_code == 0
    assert "INFO:siac.test.workflow:workflow info" in info_output.err
    assert "Sentinel-2 processing complete" in info_output.out


def test_process_s2_uses_default_info_logging_without_config(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def _fake_load(cls) -> SIACConfig:
        _ = cls
        return SIACConfig(sensor="s2")

    def _fake_process(config_obj, query, **kwargs):  # noqa: ANN001
        _ = (config_obj, query, kwargs)
        logging.getLogger("siac.test.default").info("default info")
        return SimpleNamespace(aot=SimpleNamespace(mean=lambda: 0.1234))

    monkeypatch.setattr(cli.SIACConfig, "load", classmethod(_fake_load))
    monkeypatch.setattr("siac.api.public.siac_process_s2", _fake_process)

    exit_code = cli.main(["process-s2", "T31UDQ_20240101"])

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "INFO:siac.test.default:default info" in captured.err
    assert "Sentinel-2 processing complete" in captured.out


def test_process_s2_requires_query(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = cli.main(["process-s2"])

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "the following arguments are required: query" in captured.err

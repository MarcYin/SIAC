"""Unit tests for the SIAC CLI surface."""

from __future__ import annotations

import logging
import re
from types import SimpleNamespace
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

import siac.cli as cli
from siac.config import SIACConfig
from siac.domain.aoi import AOI


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
    assert isinstance(captured["kwargs"]["aoi"], AOI)
    assert captured["kwargs"]["aoi"].get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert "Sentinel-2 processing complete" in captured_out.out


def test_process_s2_accepts_aoi_file_argument(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _fake_process(config_obj, query, **kwargs):  # noqa: ANN001
        captured["config_obj"] = config_obj
        captured["query"] = query
        captured["kwargs"] = kwargs
        return SimpleNamespace(aot=SimpleNamespace(mean=lambda: 0.1234))

    monkeypatch.setattr("siac.api.public.siac_process_s2", _fake_process)

    aoi_path = tmp_path / "aoi.geojson"
    aoi_path.write_text(
        '{"type":"Polygon","coordinates":[[[1,2],[3,2],[3,4],[1,4],[1,2]]]}',
        encoding="utf-8",
    )
    exit_code = cli.main(["process-s2", "T31UDQ_20240101", "--aoi-file", str(aoi_path)])

    captured_out = capsys.readouterr()
    assert exit_code == 0
    assert captured["query"] == "T31UDQ_20240101"
    assert isinstance(captured["kwargs"]["aoi"], AOI)
    assert captured["kwargs"]["aoi"].get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert "Sentinel-2 processing complete" in captured_out.out


def test_prepare_monthly_composites_delegates_with_period_selection(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _fake_prepare(config_obj, **kwargs):  # noqa: ANN001
        captured["config_obj"] = config_obj
        captured["kwargs"] = kwargs
        return SimpleNamespace(
            store_path=tmp_path / "prepared_store",
            period_count=2,
            periods=("2022-07", "2023-08"),
            representation="kernel_weights",
        )

    monkeypatch.setattr("siac.api.public.prepare_monthly_composites", _fake_prepare)

    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            str(tmp_path / "prepared_store"),
            "--aoi-bbox",
            "1",
            "2",
            "3",
            "4",
            "--year-month",
            "2023-08",
            "--year-month",
            "2022-07",
        ]
    )

    captured_out = capsys.readouterr()
    assert exit_code == 0
    assert isinstance(captured["kwargs"]["aoi"], AOI)
    assert captured["kwargs"]["aoi"].get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert captured["kwargs"]["year_months"] == ((2022, 7), (2023, 8))
    assert captured["kwargs"]["resolution"] is None
    assert captured["kwargs"]["output_path"] == tmp_path / "prepared_store"
    assert "Prepared monthly composites written to" in captured_out.out


def test_prepare_monthly_composites_rejects_duplicate_period_selection(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            "/tmp/prepared-store",
            "--aoi-bbox",
            "1",
            "2",
            "3",
            "4",
            "--year-month",
            "2023-07",
            "--year-month",
            "2023-07",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "Duplicate year/month selections" in captured.err


def test_prepare_monthly_composites_accepts_year_month_cross_product_and_wkt(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _fake_prepare(config_obj, **kwargs):  # noqa: ANN001
        captured["config_obj"] = config_obj
        captured["kwargs"] = kwargs
        return SimpleNamespace(
            store_path=tmp_path / "prepared_store",
            period_count=4,
            periods=("2022-07", "2022-08", "2023-07", "2023-08"),
            representation="kernel_weights",
        )

    monkeypatch.setattr("siac.api.public.prepare_monthly_composites", _fake_prepare)

    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            str(tmp_path / "prepared_store"),
            "--aoi-wkt",
            "POLYGON ((1 2, 3 2, 3 4, 1 4, 1 2))",
            "--year",
            "2023",
            "--year",
            "2022",
            "--month",
            "8",
            "--month",
            "7",
        ]
    )

    captured_out = capsys.readouterr()
    assert exit_code == 0
    assert isinstance(captured["kwargs"]["aoi"], AOI)
    assert captured["kwargs"]["aoi"].get_bounds() == pytest.approx((1.0, 2.0, 3.0, 4.0))
    assert str(captured["kwargs"]["aoi"].crs) == "EPSG:4326"
    assert captured["kwargs"]["year_months"] == ((2022, 7), (2022, 8), (2023, 7), (2023, 8))
    assert "Prepared monthly composites written to" in captured_out.out


def test_cli_applies_aoi_crs_to_wkt_inputs(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _fake_prepare(config_obj, **kwargs):  # noqa: ANN001
        captured["config_obj"] = config_obj
        captured["kwargs"] = kwargs
        return SimpleNamespace(
            store_path=tmp_path / "prepared_store",
            period_count=1,
            periods=("2023-08",),
            representation="kernel_weights",
        )

    monkeypatch.setattr("siac.api.public.prepare_monthly_composites", _fake_prepare)

    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            str(tmp_path / "prepared_store"),
            "--aoi-wkt",
            "POLYGON ((500000 4100000, 501000 4100000, 501000 4101000, 500000 4101000, 500000 4100000))",
            "--aoi-crs",
            "EPSG:32632",
            "--year-month",
            "2023-08",
        ]
    )

    captured_out = capsys.readouterr()
    assert exit_code == 0
    assert isinstance(captured["kwargs"]["aoi"], AOI)
    assert str(captured["kwargs"]["aoi"].crs) == "EPSG:32632"
    assert "Prepared monthly composites written to" in captured_out.out


def test_cli_rejects_degree_like_bbox_when_crs_is_omitted_and_values_are_projected(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            "/tmp/prepared-store",
            "--aoi-bbox",
            "500000",
            "4100000",
            "501000",
            "4101000",
            "--year-month",
            "2023-08",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "longitude bounds must look like degrees" in captured.err


def test_prepare_monthly_composites_rejects_mixed_period_selection(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = cli.main(
        [
            "prepare-monthly-composites",
            "--output-path",
            "/tmp/prepared-store",
            "--aoi-bbox",
            "1",
            "2",
            "3",
            "4",
            "--year-month",
            "2023-07",
            "--year",
            "2023",
            "--month",
            "8",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "Use either --year-month or --year/--month" in captured.err


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

    quiet_code = cli.main(
        ["process-s2", "T31UDQ_20240101", "--config", str(tmp_path / "quiet.toml")]
    )
    quiet_output = capsys.readouterr()

    info_code = cli.main(["process-s2", "T31UDQ_20240101", "--config", str(tmp_path / "info.toml")])
    info_output = capsys.readouterr()

    assert quiet_code == 0
    assert "workflow info" not in quiet_output.err
    assert "Sentinel-2 processing complete" in quiet_output.out
    assert info_code == 0
    assert re.search(
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} INFO:siac\.test\.workflow:workflow info",
        info_output.err,
    )
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
    assert re.search(
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} INFO:siac\.test\.default:default info",
        captured.err,
    )
    assert "Sentinel-2 processing complete" in captured.out


def test_process_s2_requires_query(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = cli.main(["process-s2"])

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "the following arguments are required: query" in captured.err

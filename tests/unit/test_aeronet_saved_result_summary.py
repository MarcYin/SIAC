from __future__ import annotations

import json
from typing import TYPE_CHECKING

from tools.aeronet_validation.summarize_saved_results import (
    expand_result_patterns,
    main,
    summarize_result_files,
)

if TYPE_CHECKING:
    from pathlib import Path


def _write_json(path: Path, payload: dict[str, object]) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_summarize_result_files_supports_seasonal_flag_shape(tmp_path: Path) -> None:
    paths = [
        _write_json(tmp_path / "a.json", {"matchup_id": "a", "flag": "OK"}),
        _write_json(tmp_path / "b.json", {"matchup_id": "b", "flag": "OUT"}),
        _write_json(tmp_path / "c.json", {"matchup_id": "c", "flag": "OK"}),
    ]

    summary = summarize_result_files(paths)

    assert summary.total == 3
    assert summary.within_ee == 2
    assert round(summary.pct, 1) == 66.7
    assert summary.out_sites == ("b",)


def test_summarize_result_files_prefers_phase_d_within_ee(tmp_path: Path) -> None:
    path = _write_json(
        tmp_path / "site.json",
        {"matchup_id": "site", "status": "OK", "within_ee": False, "flag": "OK"},
    )

    summary = summarize_result_files([path])

    assert summary.within_ee == 0
    assert summary.out_sites == ("site",)


def test_summarize_result_files_counts_phase_d_failure_as_out(tmp_path: Path) -> None:
    path = _write_json(
        tmp_path / "failed.json",
        {"matchup_id": "failed-site", "status": "FAILED", "reason": "missing prior"},
    )

    summary = summarize_result_files([path])

    assert summary.total == 1
    assert summary.within_ee == 0
    assert summary.out_sites == ("failed-site",)


def test_expand_result_patterns_accepts_globs(tmp_path: Path) -> None:
    first = _write_json(tmp_path / "a.json", {"flag": "OK"})
    second = _write_json(tmp_path / "b.json", {"flag": "OK"})

    assert expand_result_patterns([str(tmp_path / "*.json")]) == [first, second]


def test_main_passes_when_saved_results_meet_target(tmp_path: Path, capsys) -> None:
    _write_json(tmp_path / "a.json", {"matchup_id": "a", "flag": "OK"})
    _write_json(tmp_path / "b.json", {"matchup_id": "b", "flag": "OK"})
    _write_json(tmp_path / "c.json", {"matchup_id": "c", "flag": "OUT"})

    code = main(
        [
            str(tmp_path / "*.json"),
            "--target-pct",
            "60",
            "--expect-count",
            "3",
            "--label",
            "toy",
        ]
    )

    captured = capsys.readouterr()
    assert code == 0
    assert "toy: 2/3 = 66.7% within EE" in captured.out
    assert "PASSED" in captured.out


def test_main_fails_when_saved_results_miss_target(tmp_path: Path) -> None:
    _write_json(tmp_path / "a.json", {"matchup_id": "a", "flag": "OK"})
    _write_json(tmp_path / "b.json", {"matchup_id": "b", "flag": "OUT"})

    code = main([str(tmp_path / "*.json"), "--target-pct", "80"])

    assert code == 1

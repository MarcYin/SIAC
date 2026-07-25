from __future__ import annotations

import json
from typing import TYPE_CHECKING

from tools.aeronet_validation.summarize_phase_d_options import evaluate

if TYPE_CHECKING:
    from pathlib import Path


def _write(directory: Path, matchup_id: str, record: dict[str, object]) -> None:
    directory.mkdir(exist_ok=True)
    (directory / f"{matchup_id}.json").write_text(json.dumps(record), encoding="utf-8")


def test_evaluate_counts_no_valid_and_missing_as_strict_misses(tmp_path: Path) -> None:
    mids = ["a", "b", "c", "d"]
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    _write(baseline, "a", {"status": "OK", "retrieved": 0.1, "truth": 0.1, "within_ee": True})
    _write(baseline, "b", {"status": "OK", "retrieved": 0.5, "truth": 0.1, "within_ee": False})
    _write(baseline, "c", {"status": "NO_VALID_OBSERVATION", "retrieved": None, "truth": 0.2})
    _write(candidate, "a", {"status": "OK", "retrieved": 0.5, "truth": 0.1, "within_ee": False})
    _write(candidate, "b", {"status": "OK", "retrieved": 0.1, "truth": 0.1, "within_ee": True})
    _write(candidate, "c", {"status": "NO_VALID_OBSERVATION", "retrieved": None, "truth": 0.2})
    _write(candidate, "d", {"status": "FAILED", "truth": 0.3})

    result = evaluate(
        mids=mids,
        baseline_label="baseline",
        baseline_dir=baseline,
        candidates=[("candidate", candidate)],
    )

    base_summary, candidate_summary = result["summaries"]
    assert base_summary["hits"] == 1
    assert base_summary["valid"] == 2
    assert base_summary["no_valid_observation"] == 1
    assert base_summary["missing"] == 1
    assert candidate_summary["hits"] == 1
    assert candidate_summary["valid"] == 2
    assert candidate_summary["failed"] == 1
    assert result["paired"][0]["gains"] == 1
    assert result["paired"][0]["losses"] == 1
    assert result["paired"][0]["net"] == 0


def test_evaluate_recomputes_ee_for_alternative_retrieval_field(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    _write(
        baseline,
        "a",
        {
            "status": "OK",
            "retrieved": 0.5,
            "retrieved_winmed": 0.1,
            "truth": 0.1,
            "within_ee": False,
        },
    )

    result = evaluate(
        mids=["a"],
        baseline_label="window",
        baseline_dir=baseline,
        candidates=[],
        retrieved_key="retrieved_winmed",
    )

    assert result["retrieved_key"] == "retrieved_winmed"
    assert result["summaries"][0]["hits"] == 1

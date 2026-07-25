from __future__ import annotations

from tools.aeronet_validation.build_low_cloud_aod_report import (
    REQUIRED_CANDIDATES,
    _failure_analysis_completion_issues,
    _hit,
    _required_completion_issues,
    _single_prior_completion_issues,
)


def _summary(expected: int, *, present: int | None = None, failed: int = 0) -> dict:
    statuses = {"OK": expected - failed}
    if failed:
        statuses["FAILED"] = failed
    return {
        "expected": expected,
        "present": expected if present is None else present,
        "statuses": statuses,
    }


def test_required_completion_accepts_terminal_required_experiments() -> None:
    candidates = {}
    for name, expected in REQUIRED_CANDIDATES.items():
        summary = _summary(expected)
        candidates[name] = {"cohort": summary, "screen": summary}

    assert _required_completion_issues({"candidates": candidates}) == []


def test_required_completion_reports_missing_and_failed_records() -> None:
    candidates = {}
    for name, expected in REQUIRED_CANDIDATES.items():
        summary = _summary(expected)
        candidates[name] = {"cohort": summary, "screen": summary}
    candidates["modis_monthly_anchor_calibrated"]["screen"] = _summary(24, present=23, failed=1)

    issues = _required_completion_issues({"candidates": candidates})

    assert "modis_monthly_anchor_calibrated: 23/24 records" in issues
    assert "modis_monthly_anchor_calibrated: 1 failed records" in issues


def test_expected_error_boundary_is_inclusive() -> None:
    truth = 0.2
    assert _hit(truth + 0.05 + 0.15 * truth, truth)


def test_single_prior_completion_requires_full_ok_screens_and_promoted_cohort() -> None:
    analysis = {
        "screen_count": 24,
        "variants": {
            "b03_chi2": {
                "screen": _summary(24),
                "cohort": _summary(152),
            },
            "profile": {
                "screen": _summary(24),
                "cohort": _summary(24),
            },
        },
    }

    assert _single_prior_completion_issues(analysis) == []

    analysis["variants"]["profile"]["screen"] = _summary(24, present=23, failed=1)
    assert _single_prior_completion_issues(analysis) == [
        "single-prior profile screen: 23/24 records",
        "single-prior profile screen: non-OK statuses {'FAILED': 1}",
    ]


def test_failure_analysis_completion_requires_reconciled_rows_and_diagnostics() -> None:
    analysis = {
        "summary": {
            "cohort_count": 152,
            "current_hits": 110,
            "current_misses": 42,
            "tested_one_prior_unresolved": 12,
        },
        "failure_cases": [{} for _ in range(42)],
        "unresolved_cases": [{} for _ in range(12)],
        "data_quality": {
            "current_non_ok": 0,
            "required_diagnostic_complete": 152,
        },
    }

    assert _failure_analysis_completion_issues(analysis) == []

    analysis["failure_cases"].pop()
    analysis["data_quality"]["required_diagnostic_complete"] = 151
    assert _failure_analysis_completion_issues(analysis) == [
        "failure analysis rows: 41/42",
        "failure analysis diagnostics: 151/152",
    ]

from __future__ import annotations

from tools.aeronet_validation.build_final_aod_dashboard import (
    binned_summaries,
    summarize_rows,
)


def _row(truth: float, baseline: float, unadjusted: float, adjusted: float) -> dict[str, float]:
    return {
        "truth": truth,
        "baseline": baseline,
        "unadjusted": unadjusted,
        "adjusted": adjusted,
        "cloud": truth * 10,
    }


def test_summarize_rows_uses_requested_prediction_field() -> None:
    rows = [_row(0.1, 0.5, 0.1, 0.1), _row(0.2, 0.6, 0.2, 0.2)]
    assert summarize_rows(rows, "baseline")["hits"] == 0
    assert summarize_rows(rows, "adjusted")["hits"] == 2


def test_binned_summaries_reconcile_counts_and_metrics() -> None:
    rows = [_row(0.05, 0.05, 0.05, 0.05), _row(0.25, 0.8, 0.25, 0.25)]
    result = binned_summaries(rows, "truth", (0, 0.1, 0.5), ("low", "high"))
    assert [row["count"] for row in result] == [1, 1]
    assert sum(row["final"]["hits"] for row in result) == 2

from __future__ import annotations

import numpy as np
from tools.aeronet_validation.evaluate_global_aod_offset import (
    apply_shifted_log_offset,
    evaluate_prediction,
    rank_offsets,
)


def _rows() -> list[dict[str, object]]:
    return [
        {
            "matchup_id": f"m{index}",
            "site": f"s{index}",
            "site_fold": index,
            "truth": truth,
            "baseline": baseline,
        }
        for index, (truth, baseline) in enumerate(
            ((0.1, 0.1), (0.4, 0.1), (0.8, 0.8))
        )
    ]


def test_apply_shifted_log_offset_is_identity_at_zero() -> None:
    prediction = np.asarray([0.0, 0.1, 0.8])
    assert np.allclose(apply_shifted_log_offset(prediction, 0.0), prediction)
    assert np.all(apply_shifted_log_offset(prediction, -0.1) <= prediction)


def test_rank_offsets_uses_only_requested_development_rows() -> None:
    rows = _rows()
    prediction = np.asarray([0.1, 0.4, 0.2])
    ranking = rank_offsets(rows, prediction, np.asarray([0, 1]), (-0.2, 0.0, 0.2))
    assert ranking[0]["offset"] == 0.0
    assert ranking[0]["metrics"]["hits"] == 2


def test_evaluate_prediction_counts_gain_and_loss() -> None:
    rows = _rows()
    prediction = np.asarray([0.5, 0.4, 0.0])
    result = evaluate_prediction(rows, prediction, np.asarray([0, 1, 2]))
    assert result["gains"] == 1
    assert result["losses"] == 2

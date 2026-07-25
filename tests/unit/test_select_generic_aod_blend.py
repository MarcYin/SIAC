from __future__ import annotations

import numpy as np
from tools.aeronet_validation.select_generic_aod_blend import select_convex_blend


def test_select_convex_blend_can_choose_complementary_models() -> None:
    truth = np.asarray([0.2, 0.2, 0.8, 0.8])
    predictions = {
        "low": np.asarray([0.0, 0.0, 0.6, 0.6]),
        "high": np.asarray([0.4, 0.4, 1.0, 1.0]),
    }
    selected, ranking = select_convex_blend(truth, predictions)
    assert selected["left"] == "low"
    assert selected["right"] == "high"
    assert selected["left_weight"] == 0.5
    assert selected["metrics"]["hits"] == 4
    assert ranking[0] == selected


def test_select_convex_blend_retains_best_individual_when_blends_do_not_help() -> None:
    truth = np.asarray([0.1, 0.3, 0.5])
    predictions = {
        "exact": truth.copy(),
        "biased": truth + 0.4,
    }
    selected, _ = select_convex_blend(truth, predictions)
    assert selected["left"] == "exact"
    assert selected["right"] is None

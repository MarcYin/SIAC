from __future__ import annotations

import numpy as np
from tools.aeronet_validation.evaluate_target_domain_recipe_selection import (
    rank_prediction_matrix,
)


def test_rank_prediction_matrix_uses_only_selected_indices() -> None:
    truth = np.asarray([0.1, 0.3, 0.8])
    predictions = {
        "development_winner": np.asarray([0.1, 0.3, 0.0]),
        "confirmation_winner": np.asarray([0.0, 0.0, 0.8]),
    }
    ranking = rank_prediction_matrix(truth, predictions, np.asarray([0, 1]))
    assert ranking[0]["recipe"] == "development_winner"
    assert ranking[0]["metrics"]["hits"] == 2

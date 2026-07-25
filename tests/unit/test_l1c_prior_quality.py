from __future__ import annotations

import numpy as np
from tools.aeronet_validation.evaluate_l1c_prior_quality import (
    _paired_summary,
    _window_mask,
)


def test_window_mask_uses_pixel_centers_and_square_radius():
    mask = _window_mask(
        shape=(3, 3),
        transform=np.asarray([60.0, 0.0, 499880.0, 0.0, -60.0, 120.0]),
        epsg=32631,
        lon=3.0,
        lat=0.0,
        radius_m=61.0,
    )
    assert mask.shape == (3, 3)
    assert int(mask.sum()) == 4


def test_paired_summary_uses_only_common_matchups():
    def row(matchup_id: str, b02: float, b04: float):
        return {
            "matchup_id": matchup_id,
            "B02": {"residual": b02},
            "B04": {"residual": b04},
        }

    result = _paired_summary(
        "old",
        [row("a", 0.02, -0.01), row("b", 0.04, 0.03), row("old-only", 9.0, 9.0)],
        "new",
        [row("a", 0.01, -0.02), row("b", 0.04, 0.01), row("new-only", 0.0, 0.0)],
    )

    assert result["common"] == 2
    assert result["bands"]["B02"]["improved"] == 1
    assert result["bands"]["B02"]["tied"] == 1
    assert result["bands"]["B04"]["improved"] == 1
    assert result["bands"]["B04"]["worse"] == 1
    assert np.isclose(result["bands"]["B02"]["baseline_median_abs_error"], 0.03)
    assert np.isclose(result["bands"]["B02"]["candidate_median_abs_error"], 0.025)

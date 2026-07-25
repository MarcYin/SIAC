from __future__ import annotations

from tools.aeronet_validation.evaluate_band_consensus_aod import (
    aggregation_consensus_aod,
    band_consensus_aod,
)


def _row(spread: float, delta: float) -> dict:
    return {
        "aod": 0.4,
        "species_candidates": [
            {
                "surface_cost": 1.0,
                "band_argmin_spread": spread,
                "curve_relative_second_delta": delta,
                "band_argmin_aod": {"B02": 0.15, "B04": 0.475},
            },
            {
                "surface_cost": 2.0,
                "band_argmin_spread": 0.0,
                "curve_relative_second_delta": 0.1,
                "band_argmin_aod": {"B02": 0.1, "B04": 0.1},
            },
        ],
    }


def test_band_consensus_uses_winning_family_median_on_disagreement():
    aod, used = band_consensus_aod(_row(0.325, 0.002), spread_threshold=0.2)
    assert used
    assert aod == (0.15 + 0.475) / 2.0


def test_band_consensus_can_require_flat_curve():
    aod, used = band_consensus_aod(
        _row(0.325, 0.02),
        spread_threshold=0.2,
        require_flat_delta=0.005,
    )
    assert not used
    assert aod == 0.4


def test_aggregation_consensus_uses_agreeing_minima_far_from_final_node():
    row = _row(0.03, 0.002)
    row["species_candidates"][0]["band_argmin_aod"] = {"B02": 0.02, "B04": 0.05}
    row["aod"] = 0.25

    aod, used = aggregation_consensus_aod(row, spread_max=0.1, deviation_min=0.1)

    assert used
    assert aod == 0.035


def test_aggregation_consensus_rejects_disagreeing_minima():
    aod, used = aggregation_consensus_aod(
        _row(0.325, 0.002),
        spread_max=0.1,
        deviation_min=0.1,
    )

    assert not used
    assert aod == 0.4

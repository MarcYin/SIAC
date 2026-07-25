from __future__ import annotations

import math

from tools.aeronet_validation.analyze_local_cost_diagnostics import (
    aggregation_consensus_aod,
)


def test_aggregation_consensus_uses_agreeing_band_minima() -> None:
    aod, used = aggregation_consensus_aod(0.25, [0.02, 0.05])

    assert used is True
    assert aod == 0.035


def test_aggregation_consensus_rejects_band_disagreement() -> None:
    aod, used = aggregation_consensus_aod(0.40, [0.15, 0.475])

    assert used is False
    assert aod == 0.40


def test_aggregation_consensus_requires_material_aggregation_deviation() -> None:
    aod, used = aggregation_consensus_aod(0.10, [0.04, 0.08])

    assert used is False
    assert aod == 0.10


def test_aggregation_consensus_requires_two_finite_bands() -> None:
    aod, used = aggregation_consensus_aod(0.10, [0.04, math.nan])

    assert used is False
    assert aod == 0.10

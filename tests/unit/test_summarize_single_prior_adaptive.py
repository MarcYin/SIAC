from __future__ import annotations

from tools.aeronet_validation.summarize_single_prior_adaptive import _hit_set


def test_hit_set_returns_unique_matchup_ids() -> None:
    assert _hit_set({"hit_matchup_ids": ["a", "b", "a"]}) == {"a", "b"}

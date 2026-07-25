from __future__ import annotations

import pytest
from tools.aeronet_validation.build_low_cloud_cohort import build_cohort


def test_build_cohort_preserves_campaign_order_and_uses_strict_threshold() -> None:
    campaign_ids = ["a", "b", "c"]
    records = {
        "a": {"cloud_frac": 0.19},
        "b": {"cloud_frac": 0.20},
        "c": {"cloud_frac": 0.0},
    }

    assert build_cohort(campaign_ids, records, threshold=0.20) == ["a", "c"]


def test_build_cohort_rejects_missing_or_invalid_cloud_values() -> None:
    with pytest.raises(ValueError, match="missing cloud records"):
        build_cohort(["a", "b"], {"a": {"cloud_frac": 0.1}}, threshold=0.20)

    with pytest.raises(ValueError, match="invalid cloud_frac"):
        build_cohort(["a"], {"a": {"cloud_frac": None}}, threshold=0.20)

    assert (
        build_cohort(
            ["a"],
            {"a": {"cloud_frac": None}},
            threshold=0.20,
            max_unclassified=1,
        )
        == []
    )

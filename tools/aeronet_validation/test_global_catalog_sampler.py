"""Tests for stratified global AOI sampling."""

from __future__ import annotations

import pytest
from tools.aeronet_validation.global_catalog_sampler import (
    Candidate,
    composition,
    latitude_band,
    sample_catalog,
)


def _pool(**counts: int) -> list[Candidate]:
    out, index = [], 0
    for land_cover, count in counts.items():
        for i in range(count):
            index += 1
            out.append(
                Candidate(
                    sample_id=f"s{index:05d}",
                    longitude=float(i % 180),
                    latitude=float((i % 60) - 30),
                    land_cover=land_cover,
                    continent=("africa", "asia", "europe")[i % 3],
                )
            )
    return out


def test_latitude_bands_partition_by_magnitude() -> None:
    assert latitude_band(10.0) == "tropical"
    assert latitude_band(-30.0) == "subtropical"
    assert latitude_band(50.0) == "temperate"
    assert latitude_band(-75.0) == "polar"
    assert latitude_band(90.0) == "polar"


def test_bare_sparse_is_lifted_to_its_target_share() -> None:
    # Naturally rare (5% of the pool) but the failure mode, so it must reach 25%.
    pool = _pool(bare_sparse=50, tree_cover=500, cropland=450)
    selected = sample_catalog(pool, total=200)
    share = composition(selected)["land_cover"]["bare_sparse"]
    assert share == pytest.approx(0.25, abs=0.02)


def test_shortfall_is_filled_rather_than_producing_a_short_catalog() -> None:
    # Only 5 bare_sparse exist but 25% of 200 is 50; the rest must come from
    # elsewhere rather than the catalogue silently coming up short.
    pool = _pool(bare_sparse=5, tree_cover=500, cropland=500)
    selected = sample_catalog(pool, total=200)
    assert composition(selected)["land_cover"]["bare_sparse"] * len(selected) == 5
    assert len(selected) >= 190


def test_selection_is_deterministic_for_a_seed() -> None:
    pool = _pool(bare_sparse=100, tree_cover=300)
    a = sample_catalog(pool, total=80, seed=7)
    b = sample_catalog(pool, total=80, seed=7)
    c = sample_catalog(pool, total=80, seed=8)
    assert [x.sample_id for x in a] == [x.sample_id for x in b]
    assert [x.sample_id for x in a] != [x.sample_id for x in c]


def test_continents_are_balanced_within_a_class() -> None:
    pool = _pool(bare_sparse=300)
    selected = sample_catalog(pool, total=60, targets={"bare_sparse": 1.0})
    shares = composition(selected)["continent"]
    assert max(shares.values()) - min(shares.values()) < 0.15


def test_targets_are_validated() -> None:
    pool = _pool(bare_sparse=10)
    with pytest.raises(ValueError, match="exceed 1.0"):
        sample_catalog(pool, total=5, targets={"bare_sparse": 0.8, "grassland": 0.5})
    with pytest.raises(ValueError, match="non-negative"):
        sample_catalog(pool, total=5, targets={"bare_sparse": -0.1})


def test_empty_inputs_are_rejected() -> None:
    with pytest.raises(ValueError, match="no candidates"):
        sample_catalog([], total=5)
    with pytest.raises(ValueError, match="total must be positive"):
        sample_catalog(_pool(bare_sparse=5), total=0)

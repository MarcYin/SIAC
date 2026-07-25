from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.select_canonical_aerosol_family import (
    band_edge_fraction,
    family_priors_from_percentages,
    scale_surface_costs,
    selection_score,
)


def test_family_priors_follow_dominant_cci_component():
    dust = family_priors_from_percentages(np.asarray([80.0, 5.0, 5.0, 10.0]))
    smoke = family_priors_from_percentages(np.asarray([5.0, 5.0, 80.0, 10.0]))
    maritime = family_priors_from_percentages(np.asarray([5.0, 80.0, 5.0, 10.0]))
    assert dust["desert"] == max(dust.values())
    assert smoke["biomass_burning"] == max(smoke.values())
    assert maritime["maritime"] == max(maritime.values())
    assert sum(dust.values()) == pytest.approx(1.0)


def test_selection_score_penalizes_unlikely_flat_disagreeing_family():
    clean = selection_score(
        surface_cost=1.0,
        prior_probability=0.7,
        band_spread=0.05,
        curve_delta=0.05,
        prior_weight=1.0,
        spread_weight=1.0,
        flat_weight=1.0,
    )
    poor = selection_score(
        surface_cost=1.0,
        prior_probability=0.05,
        band_spread=0.5,
        curve_delta=0.001,
        prior_weight=1.0,
        spread_weight=1.0,
        flat_weight=1.0,
    )
    assert poor > clean


def test_selection_score_penalizes_axis_clipped_band_minima():
    interior = selection_score(
        surface_cost=1.0,
        prior_probability=0.25,
        band_spread=0.1,
        curve_delta=0.03,
        prior_weight=0.0,
        spread_weight=0.0,
        flat_weight=0.0,
        edge_fraction=0.0,
        edge_weight=2.0,
    )
    clipped = selection_score(
        surface_cost=1.0,
        prior_probability=0.25,
        band_spread=0.1,
        curve_delta=0.03,
        prior_weight=0.0,
        spread_weight=0.0,
        flat_weight=0.0,
        edge_fraction=0.5,
        edge_weight=2.0,
    )
    assert clipped - interior == pytest.approx(1.0)
    assert band_edge_fraction({"B02": 0.01, "B04": 0.8}) == 0.5
    assert band_edge_fraction({"B02": 0.2, "B04": 4.0}) == 0.5


def test_relative_surface_costs_use_within_scene_evidence_scale():
    assert scale_surface_costs([1000.0, 1100.0, 2000.0], "relative") == pytest.approx(
        [0.0, 0.1, 1.0]
    )
    assert scale_surface_costs([1.0, 2.0], "raw") == [1.0, 2.0]

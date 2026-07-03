from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.rt.aerosol_species import (
    candidate_cci_aerosol_setups,
    candidate_fraction_sets,
    cci_distribution,
    climatology_fraction_percentages,
)


def test_climatology_fraction_percentages_sum_to_100() -> None:
    fractions = climatology_fraction_percentages(lon=13.4, lat=41.9, month=6)

    assert fractions.shape == (4,)
    assert np.all(np.isfinite(fractions))
    assert float(np.sum(fractions)) == pytest.approx(100.0)


def test_candidate_fraction_sets_return_nearest_requested_count() -> None:
    candidates = candidate_fraction_sets(lon=13.4, lat=41.9, month=6, n=3)

    assert len(candidates) == 3
    for candidate in candidates:
        assert set(candidate) == {"dust", "sea_salt", "fine_strong", "fine_weak"}
        assert sum(candidate.values()) == pytest.approx(1.0)


def test_cci_distribution_builds_native_sixs_multimodal_config() -> None:
    distribution = cci_distribution(
        {"dust": 0.25, "sea_salt": 0.0, "fine_strong": 0.25, "fine_weak": 0.50}
    )

    assert distribution.rmin == pytest.approx(0.005)
    assert distribution.rmax == pytest.approx(20.0)
    assert len(distribution.components) == 3
    assert sum(c.percentage_density for c in distribution.components) == pytest.approx(100.0)
    assert all(len(c.refr_real) == 20 and len(c.refr_imag) == 20 for c in distribution.components)


def test_candidate_cci_aerosol_setups_are_rt_setup_payloads() -> None:
    setups = candidate_cci_aerosol_setups(lon=13.4, lat=41.9, month=6, n=2)

    assert len(setups) == 2
    assert all(setup.profile == "multimodal_log_normal" for setup in setups)
    assert all(setup.distribution is not None for setup in setups)

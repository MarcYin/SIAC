from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.rt.aerosol_species import (
    candidate_cci_aerosol_setups,
    candidate_fraction_sets,
    cci_distribution,
    climatology_cci_aerosol_setup,
    climatology_fraction_percentages,
    climatology_fraction_set,
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


def test_exact_climatology_fraction_set_and_setup() -> None:
    fractions = climatology_fraction_set(lon=13.4, lat=41.9, month=6)
    setup = climatology_cci_aerosol_setup(lon=13.4, lat=41.9, month=6)

    assert set(fractions) == {"dust", "sea_salt", "fine_strong", "fine_weak"}
    assert sum(fractions.values()) == pytest.approx(1.0)
    assert setup.profile == "multimodal_log_normal"
    assert setup.distribution is not None


def test_climatology_setup_caps_fine_strong(monkeypatch: pytest.MonkeyPatch) -> None:
    # ICIMOD in April: uncapped fine_strong is ~0.32, above the 0.25 default cap.
    uncapped = climatology_fraction_set(lon=85.32, lat=27.65, month=4)
    assert uncapped["fine_strong"] > 0.25

    monkeypatch.setenv("SIAC_SPECIES_FINE_STRONG_CAP", "0.15")
    capped = climatology_cci_aerosol_setup(lon=85.32, lat=27.65, month=4)
    assert capped.distribution is not None

    monkeypatch.setenv("SIAC_SPECIES_FINE_STRONG_CAP", "none")
    disabled = climatology_cci_aerosol_setup(lon=85.32, lat=27.65, month=4)
    assert disabled.distribution is not None
    # Disabling the cap must reproduce the raw climatology mixture, which
    # differs from the capped one.
    assert [c.percentage_density for c in disabled.distribution.components] != [
        c.percentage_density for c in capped.distribution.components
    ]


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

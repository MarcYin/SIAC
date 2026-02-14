"""Unit tests for Earthaccess-backed provider scaffolds (Phase M0)."""

from __future__ import annotations

from datetime import datetime

import pytest

from siac.priors.atmospheric.mcd19_earthaccess import MCD19AODProvider
from siac.priors.atmospheric.merra2 import MERRA2Provider
from siac.priors.brdf.mcd43_earthaccess import (
    MCD43EarthAccessProvider,
    VNP43EarthAccessProvider,
)


def test_merra2_provider_returns_default_prior_without_probe():
    provider = MERRA2Provider(probe_earthdata=False)
    state = provider.get_prior(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=500.0,
    )

    assert provider.source_name == "MERRA-2"
    assert state.aot.shape == (2, 2)
    assert float(state.aot.mean()) == pytest.approx(0.15)
    assert state.aot_unc.shape == (2, 2)


def test_mcd19_provider_returns_default_prior_without_probe():
    provider = MCD19AODProvider(probe_earthdata=False)
    state = provider.get_prior(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=500.0,
    )

    assert provider.source_name == "MCD19"
    assert state.aot.shape == (2, 2)
    assert float(state.aot.mean()) == pytest.approx(0.12)
    assert state.tco3_unc.shape == (2, 2)


def test_mcd43_provider_returns_default_weights_without_probe():
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    weights = provider.get_brdf_parameters(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        target_resolution=500.0,
        bands=[1, 2],
        temporal_window=16,
    )

    assert provider.source_name == "MCD43"
    assert weights.f0.shape == (2, 2, 2)
    assert list(weights.f0.coords["band"].values) == [1, 2]
    assert float(weights.f0.mean()) == pytest.approx(0.20)


def test_vnp43_provider_returns_default_weights_without_probe():
    provider = VNP43EarthAccessProvider(probe_earthdata=False)
    weights = provider.get_brdf_parameters(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        target_resolution=500.0,
        bands=[3],
        temporal_window=16,
    )

    assert provider.source_name == "VNP43"
    assert weights.f0.shape == (1, 2, 2)
    assert float(weights.f2.mean()) == pytest.approx(0.02)

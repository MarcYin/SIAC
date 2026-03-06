"""Unit tests for canonical spectral response function handling."""

from __future__ import annotations

import numpy as np
import pytest

from siac.core.srf import SpectralResponseFunction
from siac.core.srf_repository import SRFRepository


def test_srf_from_tabulated_trims_normalizes_and_derives_metadata():
    wavelengths = np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32)
    response = np.array([0.0, 0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32)

    srf = SpectralResponseFunction.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=wavelengths,
        response=response,
        source_id="s2-srf",
        source_version="1.0",
    )

    np.testing.assert_allclose(
        srf.wavelengths_nm,
        np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
    )
    np.testing.assert_allclose(
        srf.response_raw,
        np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )
    assert np.all(np.diff(srf.wavelengths_nm) > 0.0)
    assert np.all(srf.response >= 0.0)
    assert np.isfinite(srf.response).all()
    assert np.trapezoid(srf.response, srf.wavelengths_nm) == pytest.approx(1.0)
    assert srf.effective_wavelength_nm == pytest.approx(460.0)
    assert srf.centre_wavelength_nm == pytest.approx(460.0)
    assert srf.fwhm_nm == pytest.approx(20.0)


def test_srf_from_tabulated_rejects_all_zero_response():
    with pytest.raises(ValueError, match="non-zero"):
        SpectralResponseFunction.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([450.0, 460.0, 470.0], dtype=np.float32),
            response=np.zeros(3, dtype=np.float32),
        )


def test_repository_returns_platform_specific_srf():
    s2a = SpectralResponseFunction.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=np.array([450.0, 460.0, 470.0], dtype=np.float32),
        response=np.array([0.0, 1.0, 0.0], dtype=np.float32),
    )
    s2b = SpectralResponseFunction.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2B",
        band_name="B02",
        wavelengths_nm=np.array([452.0, 462.0, 472.0], dtype=np.float32),
        response=np.array([0.0, 1.0, 0.0], dtype=np.float32),
    )
    repo = SRFRepository([s2a, s2b])

    loaded = repo.get_band_srf("MSI", "S2B", "B02")
    sensor_srfs = repo.get_sensor_srfs("MSI", "S2A")

    assert loaded.satellite_id == "S2B"
    assert loaded.band_name == "B02"
    assert set(sensor_srfs) == {"B02"}
    assert sensor_srfs["B02"].satellite_id == "S2A"


def test_repository_missing_key_raises_keyerror():
    repo = SRFRepository([])

    with pytest.raises(KeyError, match="MSI"):
        repo.get_band_srf("MSI", "S2C", "B08")

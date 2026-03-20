"""Unit tests for canonical spectral response function handling."""

from __future__ import annotations

import numpy as np
import pytest

from siac.domain.spectral import SpectralResponseFunction


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


def test_srf_support_bounds_follow_trimmed_nonzero_extent():
    srf = SpectralResponseFunction.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=np.array([430.0, 440.0, 450.0, 460.0, 470.0], dtype=np.float32),
        response=np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32),
    )
    assert srf.support_min_nm == pytest.approx(440.0)
    assert srf.support_max_nm == pytest.approx(460.0)

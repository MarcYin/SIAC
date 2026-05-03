"""Unit tests for canonical relative spectral response handling."""

from __future__ import annotations

import numpy as np
import pytest

import siac.domain.spectral as spectral_mod
from siac.domain.spectral import RelativeSpectralResponse


def test_rsrf_from_tabulated_trims_normalizes_and_derives_metadata():
    wavelengths = np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32)
    response = np.array([0.0, 0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32)

    rsrf = RelativeSpectralResponse.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=wavelengths,
        response=response,
        source_id="s2-srf",
        source_version="1.0",
    )

    np.testing.assert_allclose(
        rsrf.wavelengths_nm,
        np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
    )
    np.testing.assert_allclose(
        rsrf.response_raw,
        np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )
    assert np.all(np.diff(rsrf.wavelengths_nm) > 0.0)
    assert np.all(rsrf.response >= 0.0)
    assert np.isfinite(rsrf.response).all()
    assert np.trapezoid(rsrf.response, rsrf.wavelengths_nm) == pytest.approx(1.0)
    assert rsrf.effective_wavelength_nm == pytest.approx(460.0)
    assert rsrf.center_wavelength_nm == pytest.approx(460.0)
    assert rsrf.fwhm_nm == pytest.approx(20.0)


def test_rsrf_from_tabulated_rejects_all_zero_response():
    with pytest.raises(ValueError, match="non-zero"):
        RelativeSpectralResponse.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([450.0, 460.0, 470.0], dtype=np.float32),
            response=np.zeros(3, dtype=np.float32),
        )


def test_rsrf_support_bounds_follow_trimmed_nonzero_extent():
    rsrf = RelativeSpectralResponse.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=np.array([430.0, 440.0, 450.0, 460.0, 470.0], dtype=np.float32),
        response=np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32),
    )
    assert rsrf.support_min_nm == pytest.approx(440.0)
    assert rsrf.support_max_nm == pytest.approx(460.0)


def test_trapezoid_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delattr(spectral_mod.np, "trapezoid", raising=False)

    area = spectral_mod._trapezoid(  # noqa: SLF001
        np.array([0.0, 1.0, 0.0], dtype=np.float32),
        np.array([450.0, 460.0, 470.0], dtype=np.float32),
    )

    assert area == pytest.approx(10.0)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {
                "wavelengths_nm": np.array([[450.0, 460.0]], dtype=np.float32),
                "response": np.array([0.0, 1.0], dtype=np.float32),
            },
            "must be 1-D",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, 460.0], dtype=np.float32),
                "response": np.array([0.0], dtype=np.float32),
            },
            "equal length",
        ),
        (
            {
                "wavelengths_nm": np.array([460.0], dtype=np.float32),
                "response": np.array([1.0], dtype=np.float32),
            },
            "at least two samples",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, np.nan], dtype=np.float32),
                "response": np.array([0.0, 1.0], dtype=np.float32),
            },
            "must be finite",
        ),
        (
            {
                "wavelengths_nm": np.array([460.0, 450.0], dtype=np.float32),
                "response": np.array([0.0, 1.0], dtype=np.float32),
            },
            "strictly increasing",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, 460.0], dtype=np.float32),
                "response": np.array([-1.0, 1.0], dtype=np.float32),
            },
            "non-negative",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, 460.0], dtype=np.float32),
                "response": np.array([0.0, 1.0], dtype=np.float32),
                "response_raw": np.array([0.0], dtype=np.float32),
            },
            "raw response must match wavelengths",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, 460.0], dtype=np.float32),
                "response": np.array([0.0, 0.0], dtype=np.float32),
            },
            "integrate to a positive value",
        ),
        (
            {
                "wavelengths_nm": np.array([450.0, 460.0, 470.0], dtype=np.float32),
                "response": np.array([0.0, 1.0, 0.0], dtype=np.float32),
            },
            "area-normalized",
        ),
    ],
)
def test_rsrf_post_init_validation_errors(kwargs: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        RelativeSpectralResponse(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            **kwargs,
        )


@pytest.mark.parametrize(
    ("wavelengths", "response", "message"),
    [
        (
            np.array([450.0, 460.0], dtype=np.float32),
            np.array([0.0], dtype=np.float32),
            "equal length",
        ),
        (
            np.array([460.0], dtype=np.float32),
            np.array([1.0], dtype=np.float32),
            "at least two samples",
        ),
        (
            np.array([450.0, np.nan], dtype=np.float32),
            np.array([0.0, 1.0], dtype=np.float32),
            "must be finite",
        ),
    ],
)
def test_rsrf_from_tabulated_validates_basic_inputs(
    wavelengths: np.ndarray,
    response: np.ndarray,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        RelativeSpectralResponse.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=wavelengths,
            response=response,
        )


def test_fwhm_returns_none_for_small_or_invalid_curves() -> None:
    assert (
        spectral_mod._fwhm(  # noqa: SLF001
            np.array([460.0], dtype=np.float32),
            np.array([1.0], dtype=np.float32),
        )
        is None
    )
    assert (
        spectral_mod._fwhm(  # noqa: SLF001
            np.array([450.0, 460.0], dtype=np.float32),
            np.array([-1.0, -0.5], dtype=np.float32),
        )
        is None
    )

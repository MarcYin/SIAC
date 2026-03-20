"""Unit tests for LUT-aligned SRF kernel generation."""

from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.rt.lut.srf_kernel import AlignedSRFKernel, build_aligned_srf_kernel
from siac.domain.spectral import SpectralResponseFunction


def _srf() -> SpectralResponseFunction:
    return SpectralResponseFunction.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
        response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )


def test_build_aligned_srf_kernel_uses_support_slice_and_padding():
    lut_wavelengths = np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32)

    kernel = build_aligned_srf_kernel(
        _srf(),
        lut_wavelengths_nm=lut_wavelengths,
        lut_id="lut-v1",
        support_padding=1,
    )

    assert isinstance(kernel, AlignedSRFKernel)
    assert kernel.start_index == 0
    assert kernel.end_index == 7
    np.testing.assert_allclose(kernel.wavelengths_nm, lut_wavelengths)
    assert kernel.response_on_lut.shape == lut_wavelengths.shape
    assert np.trapezoid(kernel.response_on_lut, kernel.wavelengths_nm) == pytest.approx(1.0)
    assert kernel.response_on_lut[0] == pytest.approx(0.0)
    assert kernel.response_on_lut[-1] == pytest.approx(0.0)


def test_build_aligned_srf_kernel_solar_weighting_normalizes_weighted_response():
    lut_wavelengths = np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32)
    solar = np.array([1.0, 2.0, 4.0, 2.0, 1.0], dtype=np.float32)

    kernel = build_aligned_srf_kernel(
        _srf(),
        lut_wavelengths_nm=lut_wavelengths,
        lut_id="lut-v1",
        solar_irradiance=solar,
        support_padding=0,
    )

    assert kernel.start_index == 0
    assert kernel.end_index == 5
    assert kernel.solar_weighted_response_on_lut is not None
    assert np.trapezoid(kernel.solar_weighted_response_on_lut, kernel.wavelengths_nm) == pytest.approx(1.0)
    peak_idx = int(np.argmax(kernel.solar_weighted_response_on_lut))
    assert kernel.wavelengths_nm[peak_idx] == pytest.approx(460.0)


def test_build_aligned_srf_kernel_rejects_non_monotonic_lut_axis():
    with pytest.raises(ValueError, match="strictly increasing"):
        build_aligned_srf_kernel(
            _srf(),
            lut_wavelengths_nm=np.array([440.0, 460.0, 450.0], dtype=np.float32),
            lut_id="lut-v1",
        )

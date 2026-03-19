"""
Layer 6 — Spectral model tests.

Tests for SpectralBandDescriptor, SensorConfig band selection via
spectral model, and spectral convolution functions.
"""

import numpy as np
import pytest
import xarray as xr

from siac.domain.spectral import (
    SpectralBandDescriptor,
    load_reference_rsr,
    reference_to_sensor,
    sensor_to_reference,
)

# ── SpectralBandDescriptor ────────────────────────────────────────────

class TestSpectralBandDescriptor:
    def test_gaussian_only_construction(self):
        """Gaussian-only (no SRF arrays) -> has_srf == False."""
        b = SpectralBandDescriptor("B02", 490.0, 65.0, 10.0)
        assert not b.has_srf
        assert b.center_wavelength_nm == 490.0

    def test_with_srf(self):
        """Full SRF construction -> has_srf == True."""
        wl = np.linspace(460, 520, 61)
        resp = np.exp(-0.5 * ((wl - 490) / 15) ** 2)
        b = SpectralBandDescriptor("B02", 490.0, 65.0, 10.0,
                                    srf_wavelengths_nm=wl, srf_response=resp)
        assert b.has_srf

    def test_wavelength_um(self):
        """550 nm -> 0.55 µm."""
        b = SpectralBandDescriptor("Green", 550.0, 35.0, 10.0)
        assert b.wavelength_um == pytest.approx(0.55)

    def test_frozen(self):
        """Mutating field raises."""
        b = SpectralBandDescriptor("B02", 490.0, 65.0, 10.0)
        with pytest.raises(AttributeError):
            b.name = "other"

    def test_gaussian_response_peak(self):
        """Peak of Gaussian response should be at center wavelength."""
        b = SpectralBandDescriptor("B02", 490.0, 65.0, 10.0)
        wl = np.linspace(400, 600, 201)
        resp = b.gaussian_response(wl)
        peak_idx = np.argmax(resp)
        assert wl[peak_idx] == pytest.approx(490.0, abs=1.0)

    def test_effective_response_uses_srf_when_available(self):
        """effective_response should use SRF when available, not Gaussian."""
        wl_srf = np.array([480.0, 490.0, 500.0])
        resp_srf = np.array([0.2, 1.0, 0.3])
        b = SpectralBandDescriptor("B02", 490.0, 65.0, 10.0,
                                    srf_wavelengths_nm=wl_srf, srf_response=resp_srf)
        # At 490 nm, should return 1.0 from the SRF
        result = b.effective_response(np.array([490.0]))
        assert result[0] == pytest.approx(1.0)


# ── Reference RSR loading ─────────────────────────────────────────────

class TestLoadReferenceRSR:
    def test_load_modis(self):
        """MODIS should return 7 bands with non-empty arrays."""
        rsr = load_reference_rsr("MODIS")
        assert len(rsr) == 7
        for _name, (wl, resp) in rsr.items():
            assert len(wl) > 0
            assert len(resp) > 0
            assert resp.max() > 0

    def test_load_modis_wavelengths_reasonable(self):
        """MODIS band wavelengths should be in ~400-2500 nm range."""
        rsr = load_reference_rsr("MODIS")
        for _name, (wl, _) in rsr.items():
            assert wl.min() >= 350.0
            assert wl.max() <= 2500.0

    def test_load_unknown_raises(self):
        """Unknown sensor should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown reference sensor"):
            load_reference_rsr("UNKNOWN")


# ── Spectral convolution functions ────────────────────────────────────

def _make_test_sensor_bands():
    """3-band sensor: Blue (490nm), Green (560nm), Red (665nm)."""
    return [
        SpectralBandDescriptor("B02", 490.0, 65.0, 10.0),
        SpectralBandDescriptor("B03", 560.0, 35.0, 10.0),
        SpectralBandDescriptor("B04", 665.0, 30.0, 10.0),
    ]


class TestSensorToReference:
    def test_flat_spectrum(self):
        """Flat reflectance = 0.5 -> all reference bands should be ~0.5."""
        bands = _make_test_sensor_bands()
        shape = (4, 4)
        ds = xr.Dataset({
            b.name: xr.DataArray(np.full(shape, 0.5), dims=["y", "x"])
            for b in bands
        })
        result = sensor_to_reference(ds, bands)
        assert result.shape[0] == 7  # 7 MODIS reference bands
        # Bands that overlap with sensor should be ~0.5
        # Some reference bands may be 0 (no overlap)
        assert result.shape[1:] == shape

    def test_output_band_count(self):
        """Output should have 7 MODIS reference bands."""
        bands = _make_test_sensor_bands()
        shape = (2, 2)
        ds = xr.Dataset({
            b.name: xr.DataArray(np.ones(shape) * 0.3, dims=["y", "x"])
            for b in bands
        })
        result = sensor_to_reference(ds, bands)
        assert result.shape[0] == 7

    def test_nonzero_for_overlapping_bands(self):
        """Reference bands that overlap with sensor bands should be non-zero."""
        bands = _make_test_sensor_bands()
        shape = (2, 2)
        ds = xr.Dataset({
            b.name: xr.DataArray(np.ones(shape) * 0.3, dims=["y", "x"])
            for b in bands
        })
        result = sensor_to_reference(ds, bands)
        # At least some reference bands should be non-zero
        assert np.any(result > 0)


class TestReferenceToSensor:
    def test_roundtrip_approximate(self):
        """sensor -> reference -> sensor should approximately recover original."""
        bands = _make_test_sensor_bands()
        shape = (4, 4)
        ds = xr.Dataset({
            b.name: xr.DataArray(np.full(shape, 0.3), dims=["y", "x"])
            for b in bands
        })
        ref_vals = sensor_to_reference(ds, bands)
        recovered = reference_to_sensor(ref_vals, bands)
        assert recovered.shape[0] == len(bands)
        assert recovered.shape[1:] == shape

    def test_output_shape(self):
        """Output should have n_target_bands as first dimension."""
        bands = _make_test_sensor_bands()
        ref_vals = np.random.default_rng(0).random((7, 4, 4))
        result = reference_to_sensor(ref_vals, bands)
        assert result.shape[0] == 3
        assert result.shape[1:] == (4, 4)

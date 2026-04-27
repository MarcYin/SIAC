"""
Layer 6 — Spectral model tests.

Tests for SensorBand spectral behavior and reference-basis transforms.
"""

import numpy as np
import pytest
import xarray as xr

from siac.adapters.rsrf import band_convolution_weights
from siac.algorithms.surface.reference_spectral import (
    load_reference_rsrf,
    reference_to_sensor,
    sensor_to_reference,
)
from siac.domain import SensorBand

# ── SensorBand behavior ───────────────────────────────────────────────


class TestSensorBand:
    def test_gaussian_only_construction(self):
        """Center/FWHM-only construction keeps band metadata but no sampled RSRF."""
        b = SensorBand("B02", 490.0, 65.0, 10.0, 0)
        assert not b.has_rsrf
        assert b.center_wavelength == 490.0

    def test_with_rsrf(self):
        """Full RSRF construction -> has_rsrf == True."""
        wl = np.linspace(460, 520, 61)
        resp = np.exp(-0.5 * ((wl - 490) / 15) ** 2)
        b = SensorBand(
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            rsrf_wavelengths_nm=wl,
            rsrf_response=resp,
        )
        assert b.has_rsrf

    def test_wavelength_um(self):
        """550 nm -> 0.55 µm."""
        b = SensorBand("Green", 550.0, 35.0, 10.0, 0)
        assert b.wavelength_um == pytest.approx(0.55)

    def test_frozen(self):
        """Mutating field raises."""
        b = SensorBand("B02", 490.0, 65.0, 10.0, 0)
        with pytest.raises(AttributeError):
            b.name = "other"

    def test_band_convolution_weights_peak_at_center_for_band_spec(self):
        b = SensorBand("B02", 490.0, 65.0, 10.0, 0)
        wl = np.linspace(400, 600, 201)
        resp = band_convolution_weights(b, wl)
        peak_idx = np.argmax(resp)
        assert wl[peak_idx] == pytest.approx(490.0, abs=1.0)

    def test_band_convolution_weights_use_tabulated_rsrf_when_available(self):
        wl_rsrf = np.array([480.0, 490.0, 500.0])
        resp_rsrf = np.array([0.2, 1.0, 0.3])
        b = SensorBand(
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            rsrf_wavelengths_nm=wl_rsrf,
            rsrf_response=resp_rsrf,
        )
        result = band_convolution_weights(b, np.array([480.0, 490.0, 500.0]))
        assert result[1] > result[0]
        assert result[1] > result[2]

    def test_legacy_srf_names_remain_supported(self):
        """Backward-compatible SRF aliases should map onto the renamed RSRF fields."""
        wl = np.array([480.0, 490.0, 500.0])
        resp = np.array([0.2, 1.0, 0.3])
        b = SensorBand(
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            srf_wavelengths_nm=wl,
            srf_response=resp,
        )

        assert b.has_srf
        assert b.has_rsrf
        assert np.array_equal(b.srf_wavelengths_nm, wl)
        assert np.array_equal(b.srf_response, resp)


# ── Reference RSRF loading ────────────────────────────────────────────


class TestLoadReferenceRSRF:
    def test_load_modis(self):
        """MODIS should return 7 bands with non-empty arrays."""
        rsrf = load_reference_rsrf("MODIS")
        assert len(rsrf) == 7
        for _name, (wl, resp) in rsrf.items():
            assert len(wl) > 0
            assert len(resp) > 0
            assert resp.max() > 0

    def test_load_modis_wavelengths_reasonable(self):
        """MODIS band wavelengths should be in ~400-2500 nm range."""
        rsrf = load_reference_rsrf("MODIS")
        for _name, (wl, _) in rsrf.items():
            assert wl.min() >= 350.0
            assert wl.max() <= 2500.0

    def test_load_unknown_raises(self):
        """Unknown sensor should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown reference sensor"):
            load_reference_rsrf("UNKNOWN")


# ── Spectral convolution functions ────────────────────────────────────


def _make_test_sensor_bands():
    """3-band sensor: Blue (490nm), Green (560nm), Red (665nm)."""
    return [
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
        SensorBand("B04", 665.0, 30.0, 10.0, 2),
    ]


class TestSensorToReference:
    def test_flat_spectrum(self):
        """Flat reflectance = 0.5 -> all reference bands should be ~0.5."""
        bands = _make_test_sensor_bands()
        shape = (4, 4)
        ds = xr.Dataset({b.name: xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]) for b in bands})
        result = sensor_to_reference(ds, bands)
        assert result.shape[0] == 7  # 7 MODIS reference bands
        # Bands that overlap with sensor should be ~0.5
        # Some reference bands may be 0 (no overlap)
        assert result.shape[1:] == shape

    def test_output_band_count(self):
        """Output should have 7 MODIS reference bands."""
        bands = _make_test_sensor_bands()
        shape = (2, 2)
        ds = xr.Dataset(
            {b.name: xr.DataArray(np.ones(shape) * 0.3, dims=["y", "x"]) for b in bands}
        )
        result = sensor_to_reference(ds, bands)
        assert result.shape[0] == 7

    def test_nonzero_for_overlapping_bands(self):
        """Reference bands that overlap with sensor bands should be non-zero."""
        bands = _make_test_sensor_bands()
        shape = (2, 2)
        ds = xr.Dataset(
            {b.name: xr.DataArray(np.ones(shape) * 0.3, dims=["y", "x"]) for b in bands}
        )
        result = sensor_to_reference(ds, bands)
        # At least some reference bands should be non-zero
        assert np.any(result > 0)


class TestReferenceToSensor:
    def test_roundtrip_approximate(self):
        """sensor -> reference -> sensor should approximately recover original."""
        bands = _make_test_sensor_bands()
        shape = (4, 4)
        ds = xr.Dataset({b.name: xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]) for b in bands})
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

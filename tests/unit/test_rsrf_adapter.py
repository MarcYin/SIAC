"""Tests for the thin RSRF adapter."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rsrf.models import BandSpec, SampledCurve

from siac.adapters.rsrf import load_band_srf_from_rsrf, load_sensor_config_from_rsrf
from siac.catalog import SENTINEL2C_CONFIG
from siac.domain import SensorBand, SensorConfig


def test_load_band_srf_from_rsrf_uses_sampled_curve_directly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, object] = {}

    def _load(sensor_unit_id: str, band_id: str, representation_variant: str, *, root=None):
        seen["sensor_unit_id"] = sensor_unit_id
        seen["band_id"] = band_id
        seen["representation_variant"] = representation_variant
        seen["root"] = root
        return SampledCurve(
            band_id=band_id,
            wavelength_nm=[440.0, 450.0, 460.0, 470.0, 480.0],
            response=[0.0, 0.5, 1.0, 0.5, 0.0],
            source_variant=representation_variant,
        )

    monkeypatch.setattr("siac.adapters.rsrf.rsrf.load_response_definition", _load)

    band = SENTINEL2C_CONFIG.get_band("B02")
    srf = load_band_srf_from_rsrf(
        band,
        sensor_id="MSI",
        satellite_id="S2C",
        rsrf_root="/tmp/rsrf-root",
    )

    assert seen == {
        "sensor_unit_id": "sentinel-2c_msi",
        "band_id": "B02",
        "representation_variant": "band_average",
        "root": Path("/tmp/rsrf-root").resolve(),
    }
    assert srf.band_name == "B02"
    assert np.trapezoid(srf.response, srf.wavelengths_nm) == pytest.approx(1.0)


def test_load_sensor_config_from_rsrf_realizes_band_spec(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base_band = SENTINEL2C_CONFIG.get_band("B02")
    base_config = SensorConfig(
        sensor_id="MSI",
        satellite_id="S2C",
        bands=(base_band,),
    )

    monkeypatch.setattr(
        "siac.adapters.rsrf.rsrf.load_response_definition",
        lambda *_args, **_kwargs: BandSpec(
            band_id="B02",
            center_wavelength_nm=491.0,
            fwhm_nm=66.0,
        ),
    )
    monkeypatch.setattr(
        "siac.adapters.rsrf.rsrf.realize_curve",
        lambda band_spec, *, source_variant=None, **_kwargs: SampledCurve(
            band_id=band_spec.band_id,
            wavelength_nm=[460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0],
            response=[0.0, 0.2, 0.7, 1.0, 0.7, 0.2, 0.0],
            source_variant=source_variant,
        ),
    )

    config = load_sensor_config_from_rsrf(
        "MSI",
        "S2C",
        base_config=base_config,
    )
    band = config.get_band("B02")

    assert band.has_srf
    assert band.srf_wavelengths_nm is not None
    assert band.srf_response is not None
    assert band.center_wavelength == pytest.approx(490.0)
    assert band.bandwidth == pytest.approx(28.0)


def test_load_band_srf_from_rsrf_rejects_missing_identity() -> None:
    band = SensorBand(
        name="B02",
        center_wavelength=490.0,
        bandwidth=65.0,
        resolution=10.0,
        band_index=1,
    )

    with pytest.raises(ValueError, match="RSRF sensor unit id"):
        load_band_srf_from_rsrf(band, sensor_id="MSI", satellite_id="S2A")

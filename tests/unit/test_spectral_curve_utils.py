from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.surface import spectral_mapping as spectral_mapping_mod
from siac.algorithms.surface._spectral_curve_utils import (
    classify_band_region,
    normalized_band_response,
    primary_nir_band_index,
)
from siac.domain import SensorBand


def test_normalized_band_response_integrates_to_one() -> None:
    band = SensorBand("B03", 560.0, 35.0, 10.0, 0)
    wavelengths = np.arange(500.0, 621.0, 1.0, dtype=np.float32)

    response = normalized_band_response(band, wavelengths)
    area = np.trapezoid(response, wavelengths)

    assert area == pytest.approx(1.0, rel=1e-4)


def test_classify_band_region_marks_visible_and_infrared() -> None:
    assert classify_band_region(SensorBand("VIS", 560.0, 20.0, 10.0, 0)) == "visible"
    assert classify_band_region(SensorBand("NIR", 842.0, 115.0, 10.0, 1)) == "infrared"


def test_primary_nir_band_index_prefers_band_closest_to_865nm() -> None:
    bands = (
        SensorBand("B06", 740.0, 15.0, 20.0, 0),
        SensorBand("B08", 842.0, 115.0, 10.0, 1),
        SensorBand("B8A", 865.0, 20.0, 20.0, 2),
        SensorBand("B09", 945.0, 20.0, 60.0, 3),
    )

    assert primary_nir_band_index(bands) == 2


def test_spectral_mapping_module_keeps_private_curve_helpers_available() -> None:
    wavelengths, response = spectral_mapping_mod._canonicalize_curve(
        np.array([500.0, 490.0, 490.0, 510.0], dtype=np.float32),
        np.array([0.0, 0.2, 0.5, 0.0], dtype=np.float32),
    )

    np.testing.assert_allclose(wavelengths, np.array([490.0, 500.0], dtype=np.float32))
    np.testing.assert_allclose(response, np.array([0.2, 0.0], dtype=np.float32))

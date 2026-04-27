from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import siac.algorithms.rt.lut._spectral_math as spectral_math
from siac.algorithms.rt.lut._spectral_math import (
    build_spectral_integration_weights,
    derive_standard_rt_coefficients,
    finite_mean,
    finite_range,
    spectral_scene_cache_key,
    summarize_spectral_scene,
)
from siac.domain import SensorBand


def test_finite_helpers_use_fallbacks_for_nonfinite_inputs() -> None:
    values = np.array([np.nan, np.inf, -np.inf], dtype=np.float32)

    assert finite_mean(values, fallback=1.5) == 1.5
    assert finite_range(values, fallback=(2.0, 3.0)) == (2.0, 3.0)


def test_summarize_spectral_scene_rounds_and_ignores_nonfinite_values() -> None:
    summary = summarize_spectral_scene(
        sza=np.array([20.1234, np.nan], dtype=np.float32),
        vza=np.array([10.8764, 10.8766], dtype=np.float32),
        raa=np.array([90.4444, np.inf], dtype=np.float32),
        tco3=np.array([0.3219, np.nan, 0.3221], dtype=np.float32),
        elevation=np.array([0.1234, 0.9876, -np.inf], dtype=np.float32),
    )

    assert summary == (20.123, 10.877, 90.444, (0.322, 0.322), (0.123, 0.988))


def test_spectral_scene_cache_key_flattens_summary() -> None:
    key = spectral_scene_cache_key(
        sza=np.array([20.0], dtype=np.float32),
        vza=np.array([10.0], dtype=np.float32),
        raa=np.array([90.0], dtype=np.float32),
        tco3=np.array([0.2, 0.4], dtype=np.float32),
        elevation=np.array([0.1, 0.3], dtype=np.float32),
    )

    assert key == (20.0, 10.0, 90.0, 0.2, 0.4, 0.1, 0.3)


def test_derive_standard_rt_coefficients_handles_degenerate_denominators() -> None:
    xap, xbp, xcp = derive_standard_rt_coefficients(
        toa_rho1=np.array([0.25], dtype=np.float32),
        toa_rho2=np.array([0.25], dtype=np.float32),
        eg_rho1=np.array([0.5], dtype=np.float32),
        eg_rho2=np.array([0.15], dtype=np.float32),
        rho1=0.15,
        rho2=0.5,
    )

    assert np.isfinite(xap).all()
    assert np.isfinite(xbp).all()
    assert np.isfinite(xcp).all()
    assert xap.dtype == np.float32
    assert xbp.dtype == np.float32
    assert xcp.dtype == np.float32


def test_build_spectral_integration_weights_covers_rsrf_and_error_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = xr.Dataset(
        {
            "solar_irradiance": (
                ["sample"],
                np.array([10.0, 12.0], dtype=np.float32),
            )
        },
        coords={"wavelength": np.array([440.0, 450.0, 460.0, 470.0], dtype=np.float32)},
    )
    band = SensorBand(
        "B02",
        460.0,
        20.0,
        10.0,
        0,
        rsrf_wavelengths_nm=np.array([445.0, 455.0, 465.0], dtype=np.float32),
        rsrf_response=np.array([0.0, 1.0, 0.0], dtype=np.float32),
    )

    def _fake_build_aligned_rsrf_kernel(*args, **kwargs):  # noqa: ANN002, ANN003
        return SimpleNamespace(
            start_index=1,
            end_index=3,
            solar_weighted_response_on_lut=np.array([0.4, 0.6], dtype=np.float32),
            response_on_lut=np.array([0.25, 0.75], dtype=np.float32),
        )

    monkeypatch.setattr(spectral_math, "build_aligned_rsrf_kernel", _fake_build_aligned_rsrf_kernel)

    weights = build_spectral_integration_weights(
        dataset,
        band,
        lut_path="unused.zarr",
        solar_irradiance_names=("solar_irradiance",),
        band_rsrf=lambda _: object(),
    )
    np.testing.assert_allclose(weights.values, np.array([0.0, 0.4, 0.6, 0.0], dtype=np.float32))

    with pytest.raises(ValueError, match="wavelength coordinate"):
        build_spectral_integration_weights(
            xr.Dataset(),
            SensorBand("B03", 560.0, 35.0, 10.0, 0),
            lut_path="unused.zarr",
            solar_irradiance_names=("solar_irradiance",),
            band_rsrf=lambda _: object(),
        )

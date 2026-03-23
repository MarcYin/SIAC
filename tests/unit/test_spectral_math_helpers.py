from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.lut._spectral_math import (
    build_point_interpolation_coords,
    build_spectral_integration_weights,
    derive_standard_rt_coefficients,
    finite_mean,
    finite_range,
    spectral_scene_cache_key,
    summarize_spectral_scene,
    weighted_spectral_mean,
)
from siac.domain import SensorBand
from siac.domain.spectral import RelativeSpectralResponse


def test_build_point_interpolation_coords_requires_finite_and_clips_axes() -> None:
    lut = xr.Dataset(
        coords={
            "aot": np.array([0.01, 0.5, 1.0], dtype=np.float32),
            "tcwv": np.array([0.1, 2.0, 6.0], dtype=np.float32),
        }
    )

    coords = build_point_interpolation_coords(
        lut,
        aot=np.array([-1.0, 0.2, 5.0], dtype=np.float32),
        tcwv=np.array([-3.0, 3.0, 10.0], dtype=np.float32),
        require_finite_values=lambda values, *, name: np.asarray(values, dtype=np.float32),  # noqa: ARG005
        sanitize_point_values=lambda values, axis: np.clip(values, axis.min(), axis.max()),
    )

    np.testing.assert_allclose(coords["aot"].values, np.array([0.01, 0.2, 1.0], dtype=np.float32))
    np.testing.assert_allclose(coords["tcwv"].values, np.array([0.1, 3.0, 6.0], dtype=np.float32))


def test_build_spectral_integration_weights_supports_gaussian_and_rsrf_paths() -> None:
    wavelengths = np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32)
    lut = xr.Dataset(
        {
            "solar_irradiance": ("wavelength", np.linspace(1.0, 2.0, wavelengths.size, dtype=np.float32)),
        },
        coords={"wavelength": wavelengths},
    )

    gaussian_band = SensorBand("B03", 460.0, 80.0, 10.0, 0)
    gaussian_weights = build_spectral_integration_weights(
        lut,
        gaussian_band,
        lut_path="dummy.zarr",
        solar_irradiance_names=("solar_irradiance",),
        band_rsrf=lambda band: band,  # pragma: no cover - not used in gaussian branch
    )
    assert gaussian_weights.dims == ("wavelength",)
    assert float(gaussian_weights.sel(wavelength=460.0).values) > float(
        gaussian_weights.sel(wavelength=430.0).values
    )

    rsrf_band = SensorBand(
        "B02",
        460.0,
        200.0,
        10.0,
        1,
        rsrf_wavelengths_nm=np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
        rsrf_response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )
    rsrf_weights = build_spectral_integration_weights(
        lut,
        rsrf_band,
        lut_path="dummy.zarr",
        solar_irradiance_names=("solar_irradiance",),
        band_rsrf=lambda band: RelativeSpectralResponse.from_tabulated(
            sensor_id="S2A",
            satellite_id="S2A",
            band_name=band.name,
            wavelengths_nm=np.asarray(band.rsrf_wavelengths_nm, dtype=np.float32),
            response=np.asarray(band.rsrf_response, dtype=np.float32),
        ),
    )
    assert float(rsrf_weights.sel(wavelength=430.0).values) == pytest.approx(0.0)
    assert float(rsrf_weights.sel(wavelength=460.0).values) > float(rsrf_weights.sel(wavelength=450.0).values)


def test_weighted_spectral_mean_handles_single_sample_and_coordinate_alignment() -> None:
    single = xr.DataArray(
        np.array([4.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": np.array([460.0], dtype=np.float32)},
    )
    single_weights = xr.DataArray(
        np.array([0.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": np.array([460.0], dtype=np.float32)},
    )
    assert float(weighted_spectral_mean(single, single_weights).values) == pytest.approx(0.0)

    data = xr.DataArray(
        np.array([10.0, 1.0, 2.0, 4.0, 2.0, 1.0, 9.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32)},
    )
    weights = xr.DataArray(
        np.array([0.0, 0.0, 0.25, 0.5, 0.25, 0.0, 0.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32)},
    )
    assert float(weighted_spectral_mean(data, weights).values) == pytest.approx(3.0)


def test_coefficient_and_scene_summary_helpers_return_stable_values() -> None:
    xap, xbp, xcp = derive_standard_rt_coefficients(
        toa_rho1=np.array([0.1, 0.2], dtype=np.float32),
        toa_rho2=np.array([0.2, 0.3], dtype=np.float32),
        eg_rho1=np.array([0.4, 0.5], dtype=np.float32),
        eg_rho2=np.array([0.6, 0.7], dtype=np.float32),
        rho1=0.15,
        rho2=0.5,
    )
    assert xap.shape == (2,)
    assert xbp.shape == (2,)
    assert xcp.shape == (2,)
    assert np.all(np.isfinite(xap))
    assert np.all(np.isfinite(xbp))
    assert np.all(np.isfinite(xcp))

    assert finite_range(np.array([np.nan, 1.5, 3.5], dtype=np.float32), fallback=(0.0, 0.0)) == (1.5, 3.5)
    assert finite_mean(np.array([np.nan, 1.0, 3.0], dtype=np.float32), fallback=0.0) == pytest.approx(2.0)

    summary = summarize_spectral_scene(
        sza=np.array([10.1111, 10.2222], dtype=np.float32),
        vza=np.array([5.0, 7.0], dtype=np.float32),
        raa=np.array([90.0, 92.0], dtype=np.float32),
        tco3=np.array([280.0, 300.0], dtype=np.float32),
        elevation=np.array([0.5, 1.5], dtype=np.float32),
    )
    key = spectral_scene_cache_key(
        sza=np.array([10.1111, 10.2222], dtype=np.float32),
        vza=np.array([5.0, 7.0], dtype=np.float32),
        raa=np.array([90.0, 92.0], dtype=np.float32),
        tco3=np.array([280.0, 300.0], dtype=np.float32),
        elevation=np.array([0.5, 1.5], dtype=np.float32),
    )
    assert summary == (
        10.167,
        6.0,
        91.0,
        (280.0, 300.0),
        (0.5, 1.5),
    )
    assert key == (10.167, 6.0, 91.0, 280.0, 300.0, 0.5, 1.5)

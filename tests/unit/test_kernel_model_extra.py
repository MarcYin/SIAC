from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.surface import kernel_model as km
from siac.algorithms.surface.kernel_model import KernelModelDeriver, PSFConvolver
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights, GeometryAngles


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    coords = {
        "y": np.arange(shape[0], dtype=np.float32),
        "x": np.arange(shape[1], dtype=np.float32),
    }
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"], coords=coords),
        vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"], coords=coords),
    )


def _weights(shape: tuple[int, int] = (2, 2)) -> BRDFKernelWeights:
    coords = {
        "band": ["SRC_A", "SRC_B"],
        "y": np.arange(shape[0], dtype=np.float32),
        "x": np.arange(shape[1], dtype=np.float32),
    }
    f0 = xr.DataArray(
        np.full((2, *shape), 0.2, dtype=np.float32),
        dims=["band", "y", "x"],
        coords=coords,
    )
    return BRDFKernelWeights(
        f0=f0,
        f1=xr.full_like(f0, 0.03),
        f2=xr.full_like(f0, 0.01),
        f0_unc=xr.full_like(f0, 0.02),
        f1_unc=xr.full_like(f0, 0.01),
        f2_unc=xr.full_like(f0, 0.005),
    )


def test_surface_prior_experiment_toggles_are_gated_and_apply_expected_formulas(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    boa = xr.DataArray(
        np.array([[[0.03]], [[0.04]], [[0.10]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02", "B04", "B8A"], "y": [0], "x": [0]},
    )
    boa_unc = xr.full_like(boa, 0.02)

    monkeypatch.setattr(km, "_APPLY_PRIOR_DEBIAS", False)
    monkeypatch.setattr(km, "_APPLY_DARK_TARGET_UNC", False)
    assert km.debias_visible_prior(boa) is boa
    assert km.inflate_dark_target_uncertainty(boa, boa_unc) is boa_unc

    monkeypatch.setattr(km, "_APPLY_PRIOR_DEBIAS", True)
    debiased = km.debias_visible_prior(boa)
    a = np.array([0.0133, 0.0076, 0.0], dtype=np.float32)[:, None, None]
    b = np.array([0.928, 0.950, 1.0], dtype=np.float32)[:, None, None]
    full_correction = (boa.values + a) / b
    expected_debias = boa.values + np.float32(0.3) * (full_correction - boa.values)
    np.testing.assert_allclose(debiased.values, expected_debias, rtol=1e-6)

    monkeypatch.setattr(km, "_APPLY_DARK_TARGET_UNC", True)
    inflated = km.inflate_dark_target_uncertainty(boa, boa_unc)
    deficit = np.maximum(0.0, 0.06 - boa.values)
    expected_unc = np.sqrt(boa_unc.values**2 + (0.6 * deficit) ** 2)
    np.testing.assert_allclose(inflated.values, expected_unc, rtol=1e-6)
    assert inflated.values[2, 0, 0] == pytest.approx(boa_unc.values[2, 0, 0])


def test_compute_surface_prior_maps_multispectral_output_before_psf(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[tuple[str, ...], tuple[str, ...], int]] = []

    def fake_map(
        boa: xr.DataArray,
        *,
        source_bands,
        target_bands,
        source_uncertainty: xr.DataArray,
        spectral_library,
        k_neighbors: int,
    ) -> tuple[xr.DataArray, xr.DataArray]:
        calls.append(
            (
                tuple(band.name for band in source_bands),
                tuple(band.name for band in target_bands),
                k_neighbors,
            )
        )
        coords = {
            "band": [band.name for band in target_bands],
            "y": boa.coords["y"],
            "x": boa.coords["x"],
        }
        mapped = xr.DataArray(
            np.asarray(boa.values, dtype=np.float32) + 0.1,
            dims=boa.dims,
            coords=coords,
        )
        mapped_unc = xr.DataArray(
            np.asarray(source_uncertainty.values, dtype=np.float32) + 0.05,
            dims=source_uncertainty.dims,
            coords=coords,
        )
        return mapped, mapped_unc

    monkeypatch.setattr(km, "map_multispectral_reflectance", fake_map)

    deriver = KernelModelDeriver(psf_sigma_x=0.6, psf_sigma_y=0.4, apply_psf=True)
    source_bands = (
        SensorBand("SRC_A", 500.0, 20.0, 10.0, 0),
        SensorBand("SRC_B", 860.0, 20.0, 10.0, 1),
    )
    target_bands = (
        SensorBand("DST_A", 520.0, 20.0, 10.0, 0),
        SensorBand("DST_B", 870.0, 20.0, 10.0, 1),
    )

    prior = deriver.compute_surface_prior(
        _weights(),
        _geometry((2, 2)),
        source_bands=source_bands,
        target_bands=target_bands,
        spectral_k_neighbors=7,
    )

    assert calls == [(("SRC_A", "SRC_B"), ("DST_A", "DST_B"), 7)]
    assert tuple(prior.boa.coords["band"].values.tolist()) == ("DST_A", "DST_B")
    assert prior.boa_unc.shape == prior.boa.shape
    assert bool(np.isfinite(prior.boa.values).all())
    assert bool(prior.mask.values.all())


def test_align_kernels_to_brdf_grid_covers_passthrough_and_resize_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    deriver = KernelModelDeriver(apply_psf=False)
    ref = xr.DataArray(
        np.ones((2, 3), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [10.0, 20.0], "x": [30.0, 40.0, 50.0]},
    )

    raw_vol = np.zeros((2, 3), dtype=np.float32)
    raw_geo = np.ones((2, 3), dtype=np.float32)
    passthrough = deriver._align_kernels_to_brdf_grid(raw_vol, raw_geo, ref)
    assert passthrough == (raw_vol, raw_geo)

    no_xy_ref = xr.DataArray(np.array([1.0, 2.0], dtype=np.float32), dims=["band"])
    same_vol = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    same_geo = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    assert deriver._align_kernels_to_brdf_grid(same_vol, same_geo, no_xy_ref) == (
        same_vol,
        same_geo,
    )

    aligned_vol, aligned_geo = deriver._align_kernels_to_brdf_grid(
        xr.DataArray(np.ones((2, 3), dtype=np.float32), dims=["y", "x"], coords=ref.coords),
        xr.DataArray(np.ones((2, 3), dtype=np.float32), dims=["y", "x"], coords=ref.coords),
        ref,
    )
    assert aligned_vol.identical(
        xr.DataArray(np.ones((2, 3), dtype=np.float32), dims=["y", "x"], coords=ref.coords)
    )
    assert aligned_geo.identical(
        xr.DataArray(np.ones((2, 3), dtype=np.float32), dims=["y", "x"], coords=ref.coords)
    )

    def raising_interp(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        raise ValueError("force resize path")

    monkeypatch.setattr(xr.DataArray, "interp", raising_interp)
    fallback_vol, fallback_geo = deriver._align_kernels_to_brdf_grid(
        xr.DataArray(np.arange(4, dtype=np.float32).reshape(2, 2), dims=["y", "x"]),
        xr.DataArray(np.arange(4, dtype=np.float32).reshape(2, 2), dims=["y", "x"]),
        ref,
    )
    assert fallback_vol.shape == ref.shape
    assert fallback_geo.shape == ref.shape
    assert np.array_equal(fallback_vol.coords["x"].values, ref.coords["x"].values)


def test_resize_kernel_grid_handles_degenerate_inputs(monkeypatch: pytest.MonkeyPatch) -> None:
    target_y = xr.DataArray(np.array([0.0, 1.0], dtype=np.float32), dims=["y"])
    target_x = xr.DataArray(np.array([0.0, 1.0, 2.0], dtype=np.float32), dims=["x"])

    one_dim = xr.DataArray(np.array([1.0, 2.0], dtype=np.float32), dims=["z"])
    assert KernelModelDeriver._resize_kernel_grid(one_dim, (2, 3), target_y, target_x) is one_dim

    same_shape = xr.DataArray(np.ones((2, 3), dtype=np.float32), dims=["y", "x"])
    resized_same = KernelModelDeriver._resize_kernel_grid(same_shape, (2, 3), target_y, target_x)
    assert resized_same.shape == (2, 3)
    assert np.array_equal(resized_same.coords["y"].values, target_y.values)

    empty = xr.DataArray(np.empty((0, 2), dtype=np.float32), dims=["y", "x"])
    resized_empty = KernelModelDeriver._resize_kernel_grid(empty, (2, 3), target_y, target_x)
    assert resized_empty.shape == (2, 3)
    assert bool(np.isnan(resized_empty.values).all())

    def fake_zoom(_src, _zoom, order):  # type: ignore[no-untyped-def]
        del order
        return np.ones((1, 1), dtype=np.float32)

    monkeypatch.setattr(km.ndimage, "zoom", fake_zoom)
    padded = KernelModelDeriver._resize_kernel_grid(
        xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"]),
        (2, 3),
        target_y,
        target_x,
    )
    assert padded.shape == (2, 3)
    assert np.isnan(padded.values[1, 2])


def test_apply_psf_and_convolver_support_2d_inputs() -> None:
    deriver = KernelModelDeriver(apply_psf=True)
    data = xr.DataArray(
        np.array([[1.0, np.nan], [3.0, 5.0]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1, 2], "x": [3, 4]},
    )

    out = deriver._apply_psf(data, sigma_x=0.8, sigma_y=0.6)
    assert out.dims == ("y", "x")
    assert out.shape == (2, 2)

    convolved = PSFConvolver(sigma_x=0.8, sigma_y=0.6).convolve(
        np.array([[np.nan, 1.0], [2.0, 3.0]], dtype=np.float32)
    )
    assert isinstance(convolved, np.ndarray)
    assert convolved.shape == (2, 2)
    assert np.isfinite(convolved[1, 1])


def test_scale_psf_sigmas_preserves_physical_footprint() -> None:
    """PSF sigma (pixels @ target_resolution_m) rescales to the grid resolution."""
    deriver = KernelModelDeriver(psf_sigma_x=29.75, psf_sigma_y=39.0, target_resolution_m=10.0)

    # Coarser grid (120 m) -> fewer pixels for the same physical blur.
    sx, sy = deriver._scale_psf_sigmas(29.75, 39.0, 120.0)
    assert sx == pytest.approx(29.75 * 10.0 / 120.0)
    assert sy == pytest.approx(39.0 * 10.0 / 120.0)

    # Calibration resolution -> identity (scale == 1).
    assert deriver._scale_psf_sigmas(29.75, 39.0, 10.0) == pytest.approx((29.75, 39.0))


def test_scale_psf_sigmas_passthrough_on_missing_or_invalid_resolution() -> None:
    deriver = KernelModelDeriver(psf_sigma_x=2.0, psf_sigma_y=3.0, target_resolution_m=10.0)
    for bad in (None, 0.0, -5.0, float("nan"), float("inf")):
        assert deriver._scale_psf_sigmas(2.0, 3.0, bad) == (2.0, 3.0)


def test_compute_surface_prior_applies_grid_scaled_sigma(monkeypatch: pytest.MonkeyPatch) -> None:
    """``grid_resolution_m`` must reach ``_apply_psf`` as a rescaled sigma."""
    deriver = KernelModelDeriver(psf_sigma_x=29.75, psf_sigma_y=39.0, target_resolution_m=10.0)
    captured: list[tuple[float, float]] = []

    def fake_apply(data: xr.DataArray, sigma_x: float, sigma_y: float) -> xr.DataArray:
        captured.append((sigma_x, sigma_y))
        return data

    monkeypatch.setattr(deriver, "_apply_psf", fake_apply)

    deriver.compute_surface_prior(_weights(), _geometry((2, 2)), grid_resolution_m=120.0)

    assert captured, "PSF was not applied"
    sx, sy = captured[0]
    assert sx == pytest.approx(29.75 * 10.0 / 120.0)
    assert sy == pytest.approx(39.0 * 10.0 / 120.0)

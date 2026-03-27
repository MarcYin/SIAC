from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import rasterio
import xarray as xr

import siac.algorithms.surface.brdf_whittaker as brdf_whittaker_mod
import siac.algorithms.surface.kernel_model as kernel_model_mod
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


@pytest.fixture(autouse=True)
def _stub_kernel_model_kernels(monkeypatch: pytest.MonkeyPatch) -> None:
    class _FakeBRDFKernels:
        def __init__(self, hb: float = 2.0, br: float = 1.0) -> None:
            self.hb = hb
            self.br = br

        def compute(self, vza, sza, raa):  # noqa: ANN001
            del sza, raa
            if isinstance(vza, xr.DataArray):
                zeros = xr.zeros_like(vza, dtype=np.float32)
                return zeros, zeros.copy(deep=False)
            shape = np.asarray(vza).shape
            zeros = np.zeros(shape, dtype=np.float32)
            return zeros, zeros.copy()

    monkeypatch.setattr(kernel_model_mod, "BRDFKernels", _FakeBRDFKernels)
    monkeypatch.setattr(brdf_whittaker_mod, "BRDFKernels", _FakeBRDFKernels)

    def _fake_whittaker_smooth_cube(cube, weights, temporal_lambda):  # noqa: ANN001
        del weights, temporal_lambda
        values = np.asarray(cube, dtype=np.float32).copy()
        smoothed = values.copy()
        time_axis = np.arange(values.shape[0], dtype=np.float32)
        for band_index in range(values.shape[1]):
            for y_index in range(values.shape[2]):
                for x_index in range(values.shape[3]):
                    series = values[:, band_index, y_index, x_index]
                    finite = np.isfinite(series)
                    if not np.any(finite):
                        continue
                    if finite.sum() == 1:
                        smoothed[:, band_index, y_index, x_index] = series[finite][0]
                        continue
                    smoothed[:, band_index, y_index, x_index] = np.interp(
                        time_axis,
                        time_axis[finite],
                        series[finite],
                    ).astype(np.float32)
        return smoothed

    monkeypatch.setattr(brdf_whittaker_mod, "whittaker_smooth_cube", _fake_whittaker_smooth_cube)


def _weights(
    f0_values: np.ndarray,
    f0_unc_values: np.ndarray,
    *,
    times: np.ndarray | None = None,
    georeferenced: bool = False,
) -> BRDFKernelWeights:
    def _geo_metadata(height: int, width: int) -> tuple[np.ndarray, np.ndarray, rasterio.Affine]:
        xmin = 399960.0
        ymax = 4700040.0
        resolution = 500.0
        xmax = xmin + width * resolution
        ymin = ymax - height * resolution
        transform = rasterio.transform.from_bounds(xmin, ymin, xmax, ymax, width, height)
        x = np.linspace(xmin + resolution / 2.0, xmax - resolution / 2.0, width, dtype=np.float64)
        y = np.linspace(ymax - resolution / 2.0, ymin + resolution / 2.0, height, dtype=np.float64)
        return y, x, transform

    def _with_geo(data: xr.DataArray) -> xr.DataArray:
        y, x, transform = _geo_metadata(int(data.sizes["y"]), int(data.sizes["x"]))
        out = data.assign_coords({"y": y, "x": x}).rio.set_spatial_dims(x_dim="x", y_dim="y")
        return out.rio.write_crs("EPSG:32615").rio.write_transform(transform)

    if times is None:
        coords = {"band": ["B02"], "y": np.arange(f0_values.shape[1]), "x": np.arange(f0_values.shape[2])}
        arrays = {
            "f0": xr.DataArray(f0_values, dims=["band", "y", "x"], coords=coords),
            "f1": xr.DataArray(np.zeros_like(f0_values), dims=["band", "y", "x"], coords=coords),
            "f2": xr.DataArray(np.zeros_like(f0_values), dims=["band", "y", "x"], coords=coords),
            "f0_unc": xr.DataArray(f0_unc_values, dims=["band", "y", "x"], coords=coords),
            "f1_unc": xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["band", "y", "x"], coords=coords),
            "f2_unc": xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["band", "y", "x"], coords=coords),
        }
        if georeferenced:
            arrays = {name: _with_geo(data) for name, data in arrays.items()}
        return BRDFKernelWeights(**arrays)

    coords = {
        "time": times,
        "band": ["B02"],
        "y": np.arange(f0_values.shape[2]),
        "x": np.arange(f0_values.shape[3]),
    }
    zeros = np.zeros_like(f0_values)
    arrays = {
        "f0": xr.DataArray(f0_values, dims=["time", "band", "y", "x"], coords=coords),
        "f1": xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords=coords),
        "f2": xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords=coords),
        "f0_unc": xr.DataArray(f0_unc_values, dims=["time", "band", "y", "x"], coords=coords),
        "f1_unc": xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["time", "band", "y", "x"], coords=coords),
        "f2_unc": xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["time", "band", "y", "x"], coords=coords),
    }
    if georeferenced:
        arrays = {name: _with_geo(data) for name, data in arrays.items()}
    return BRDFKernelWeights(**arrays)


def test_whittaker_deriver_fills_center_gap_and_returns_surface_prior() -> None:
    times = np.array(
        [(datetime(2024, 7, 15) + timedelta(days=offset)).date().isoformat() for offset in (-2, -1, 0, 1, 2)],
        dtype="datetime64[D]",
    )
    f0 = np.array([0.10, 0.12, np.nan, 0.16, 0.18], dtype=np.float32).reshape(5, 1, 1, 1)
    unc = np.array([0.03, 0.03, np.nan, 0.03, 0.03], dtype=np.float32).reshape(5, 1, 1, 1)

    weights = _weights(f0, unc, times=times)

    deriver = BRDFWhittakerDeriver(
        temporal_lambda=10.0,
        apply_psf=False,
    )
    prior = deriver.compute_surface_prior(
        weights,
        _geometry((1, 1)),
        obs_time=datetime(2024, 7, 15),
    )

    assert isinstance(prior, SurfacePrior)
    assert prior.boa.dims == ("band", "y", "x")
    assert float(prior.boa.sel(band="B02", y=0, x=0)) == pytest.approx(0.14, abs=0.03)
    assert float(prior.boa_unc.sel(band="B02", y=0, x=0)) > 0.0
    assert bool(prior.mask.sel(y=0, x=0))


def test_whittaker_deriver_defaults_when_all_temporal_samples_missing() -> None:
    times = np.array(
        [(datetime(2024, 7, 15) + timedelta(days=offset)).date().isoformat() for offset in (-1, 0, 1)],
        dtype="datetime64[D]",
    )
    shape = (3, 1, 1, 1)
    nan_cube = np.full(shape, np.nan, dtype=np.float32)

    weights = _weights(nan_cube, nan_cube, times=times)

    prior = BRDFWhittakerDeriver(temporal_lambda=10.0, apply_psf=False).compute_surface_prior(
        weights,
        _geometry((1, 1)),
        obs_time=datetime(2024, 7, 15),
    )

    assert float(prior.boa.sel(band="B02", y=0, x=0)) == pytest.approx(0.20)
    assert float(prior.boa_unc.sel(band="B02", y=0, x=0)) == pytest.approx(0.08)


def test_whittaker_deriver_requires_time_dimension() -> None:
    weights = _weights(
        np.array([[[0.10]]], dtype=np.float32),
        np.array([[[0.03]]], dtype=np.float32),
    )

    with pytest.raises(ValueError, match="requires BRDF weights with a 'time' dimension"):
        BRDFWhittakerDeriver(apply_psf=False).compute_surface_prior(
            weights,
            _geometry((1, 1)),
        )


def test_whittaker_deriver_rejects_unexpected_kwargs() -> None:
    times = np.array(["2024-07-14", "2024-07-15", "2024-07-16"], dtype="datetime64[D]")
    values = np.array([0.10, 0.12, 0.14], dtype=np.float32).reshape(3, 1, 1, 1)
    unc = np.full_like(values, 0.03)

    with pytest.raises(TypeError, match="Unexpected keyword argument\\(s\\): unexpected_flag"):
        BRDFWhittakerDeriver(apply_psf=False).compute_surface_prior(
            _weights(values, unc, times=times),
            _geometry((1, 1)),
            unexpected_flag=True,
        )


def test_target_index_defaults_to_center_and_uses_nearest_date() -> None:
    times = np.array(["2024-07-10", "2024-07-14", "2024-07-21"], dtype="datetime64[D]")

    assert BRDFWhittakerDeriver._target_index(times, None) == 1
    assert BRDFWhittakerDeriver._target_index(times, datetime(2024, 7, 20)) == 2


def test_whittaker_deriver_uses_neighbor_uncertainty_when_target_is_missing() -> None:
    times = np.array(["2024-07-14", "2024-07-15", "2024-07-16"], dtype="datetime64[D]")
    values = np.array([0.10, 0.10, 0.10], dtype=np.float32).reshape(3, 1, 1, 1)
    unc = np.array([0.04, np.nan, 0.02], dtype=np.float32).reshape(3, 1, 1, 1)

    prior = BRDFWhittakerDeriver(apply_psf=False).compute_surface_prior(
        _weights(values, unc, times=times),
        _geometry((1, 1)),
        obs_time=datetime(2024, 7, 15),
    )
    target_unc = BRDFWhittakerDeriver._target_uncertainty(unc, 1)
    derived_unc = float(prior.boa_unc.sel(band="B02", y=0, x=0))

    assert float(prior.boa.sel(band="B02", y=0, x=0)) == pytest.approx(0.10)
    assert float(target_unc[0, 0, 0]) == pytest.approx(0.02, abs=1e-6)
    assert np.isfinite(derived_unc)
    assert 0.019999 <= derived_unc < 0.08


def test_whittaker_deriver_preserves_spatial_metadata() -> None:
    times = np.array(["2024-07-14", "2024-07-15", "2024-07-16"], dtype="datetime64[D]")
    values = np.full((3, 1, 2, 2), 0.10, dtype=np.float32)
    unc = np.full_like(values, 0.03)
    weights = _weights(values, unc, times=times, georeferenced=True)

    prior = BRDFWhittakerDeriver(apply_psf=False).compute_surface_prior(
        weights,
        _geometry((2, 2)),
        obs_time=datetime(2024, 7, 15),
    )

    expected = weights.f0.isel(time=0, drop=True)
    assert prior.boa.rio.crs is not None
    assert prior.boa.rio.crs.to_string() == "EPSG:32615"
    assert prior.boa_unc.rio.crs is not None
    assert prior.boa_unc.rio.crs.to_string() == "EPSG:32615"
    assert prior.mask.rio.crs is not None
    assert prior.mask.rio.crs.to_string() == "EPSG:32615"
    assert prior.boa.rio.transform(recalc=True) == expected.rio.transform(recalc=True)
    assert prior.mask.rio.transform(recalc=True) == expected.rio.transform(recalc=True)

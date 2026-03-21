from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


def _weights(
    f0_values: np.ndarray,
    f0_unc_values: np.ndarray,
    *,
    times: np.ndarray | None = None,
) -> BRDFKernelWeights:
    if times is None:
        return BRDFKernelWeights(
            f0=xr.DataArray(f0_values, dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
            f1=xr.DataArray(np.zeros_like(f0_values), dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
            f2=xr.DataArray(np.zeros_like(f0_values), dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
            f0_unc=xr.DataArray(f0_unc_values, dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
            f1_unc=xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
            f2_unc=xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["band", "y", "x"], coords={"band": ["B02"], "y": [0], "x": [0]}),
        )

    coords = {"time": times, "band": ["B02"], "y": [0], "x": [0]}
    zeros = np.zeros_like(f0_values)
    return BRDFKernelWeights(
        f0=xr.DataArray(f0_values, dims=["time", "band", "y", "x"], coords=coords),
        f1=xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords=coords),
        f2=xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords=coords),
        f0_unc=xr.DataArray(f0_unc_values, dims=["time", "band", "y", "x"], coords=coords),
        f1_unc=xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["time", "band", "y", "x"], coords=coords),
        f2_unc=xr.DataArray(np.full_like(f0_unc_values, 0.01), dims=["time", "band", "y", "x"], coords=coords),
    )


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
    assert 0.02 <= derived_unc < 0.08

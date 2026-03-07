from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from siac.core.types import BRDFKernelWeights, GeometryAngles, SurfacePrior
from siac.priors.surface.brdf_whittaker import BRDFWhittakerDeriver


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


def test_whittaker_deriver_fills_center_gap_and_returns_surface_prior() -> None:
    times = np.array(
        [(datetime(2024, 7, 15) + timedelta(days=offset)).date().isoformat() for offset in (-2, -1, 0, 1, 2)],
        dtype="datetime64[D]",
    )
    f0 = np.array([0.10, 0.12, np.nan, 0.16, 0.18], dtype=np.float32).reshape(5, 1, 1, 1)
    zeros = np.zeros_like(f0)
    unc = np.array([0.03, 0.03, np.nan, 0.03, 0.03], dtype=np.float32).reshape(5, 1, 1, 1)

    weights = BRDFKernelWeights(
        f0=xr.DataArray(f0, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f1=xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f2=xr.DataArray(zeros, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f1_unc=xr.DataArray(np.full_like(unc, 0.01), dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f2_unc=xr.DataArray(np.full_like(unc, 0.01), dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
    )

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

    weights = BRDFKernelWeights(
        f0=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f1=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f2=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f0_unc=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f1_unc=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
        f2_unc=xr.DataArray(nan_cube, dims=["time", "band", "y", "x"], coords={"time": times, "band": ["B02"], "y": [0], "x": [0]}),
    )

    prior = BRDFWhittakerDeriver(temporal_lambda=10.0, apply_psf=False).compute_surface_prior(
        weights,
        _geometry((1, 1)),
        obs_time=datetime(2024, 7, 15),
    )

    assert float(prior.boa.sel(band="B02", y=0, x=0)) == pytest.approx(0.20)
    assert float(prior.boa_unc.sel(band="B02", y=0, x=0)) == pytest.approx(0.08)

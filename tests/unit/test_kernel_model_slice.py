from __future__ import annotations

import numpy as np
import xarray as xr

from siac.algorithms.surface.kernel_model import KernelModelDeriver, PSFConvolver
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


def _weights(shape: tuple[int, int]) -> BRDFKernelWeights:
    coords = {
        "y": np.linspace(10.0, 10.0 + shape[0] - 1, shape[0], dtype=np.float32),
        "x": np.linspace(20.0, 20.0 + shape[1] - 1, shape[1], dtype=np.float32),
    }
    f0 = xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"], coords=coords)
    return BRDFKernelWeights(
        f0=f0,
        f1=xr.full_like(f0, 0.05),
        f2=xr.full_like(f0, 0.02),
        f0_unc=xr.full_like(f0, 0.01),
        f1_unc=xr.full_like(f0, 0.005),
        f2_unc=xr.full_like(f0, 0.002),
    )


def test_kernel_model_aligns_geometry_kernels_to_brdf_grid() -> None:
    deriver = KernelModelDeriver(apply_psf=False)
    weights = _weights((4, 4))

    prior = deriver.compute_surface_prior(weights, _geometry((8, 8)))

    assert prior.boa.shape == (4, 4)
    assert tuple(prior.boa.dims) == ("y", "x")
    assert np.array_equal(prior.boa.coords["y"].values, weights.f0.coords["y"].values)
    assert np.array_equal(prior.boa.coords["x"].values, weights.f0.coords["x"].values)


def test_psf_convolver_preserves_dataarray_template() -> None:
    data = xr.DataArray(
        np.array([[1.0, np.nan], [2.0, 3.0]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [3, 4], "x": [7, 8]},
    )

    out = PSFConvolver(sigma_x=0.8, sigma_y=1.2).convolve(data)

    assert isinstance(out, xr.DataArray)
    assert out.dims == data.dims
    assert np.array_equal(out.coords["y"].values, data.coords["y"].values)
    assert np.array_equal(out.coords["x"].values, data.coords["x"].values)
    assert out.shape == data.shape

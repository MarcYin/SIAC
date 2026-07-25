"""Unit tests for observation-side PSF convolution and integer shift fitting."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.grid.toa_psf import (
    ShiftFitResult,
    ToaPsfConfig,
    _pearson,
    _scaled_sigmas,
    apply_integer_shift,
    convolve_toa_band,
    convolve_toa_cube,
    fit_integer_shift,
    psf_convolve_and_align_toa,
)


def _da(arr: np.ndarray, bands: list[str] | None = None) -> xr.DataArray:
    if arr.ndim == 3:
        ny, nx = arr.shape[1], arr.shape[2]
        return xr.DataArray(
            arr,
            dims=("band", "y", "x"),
            coords={
                "band": bands or [f"B{i}" for i in range(arr.shape[0])],
                "y": np.arange(ny),
                "x": np.arange(nx),
            },
        )
    ny, nx = arr.shape
    return xr.DataArray(arr, dims=("y", "x"), coords={"y": np.arange(ny), "x": np.arange(nx)})


def _textured(ny: int = 80, nx: int = 90, seed: int = 0) -> np.ndarray:
    """Spatially-correlated random field (numpy-only 5-point box blur ×3)."""
    rng = np.random.RandomState(seed)
    a = rng.rand(ny, nx).astype("float64")
    for _ in range(3):
        a = (a + np.roll(a, 1, 0) + np.roll(a, -1, 0) + np.roll(a, 1, 1) + np.roll(a, -1, 1)) / 5.0
    return a


def test_scaled_sigmas_matches_resolution_ratio():
    sx, sy = _scaled_sigmas(ToaPsfConfig(), 120.0)
    assert sx == pytest.approx(29.75 * 10.0 / 120.0)
    assert sy == pytest.approx(39.0 * 10.0 / 120.0)


def test_scaled_sigmas_passthrough_when_no_grid():
    sx, sy = _scaled_sigmas(ToaPsfConfig(), 0.0)
    assert (sx, sy) == (29.75, 39.0)


def test_convolve_constant_field_preserved():
    arr = np.full((40, 50), 0.2, dtype="float64")
    out = convolve_toa_band(arr, 2.0, 3.0)
    assert np.allclose(out[5:-5, 5:-5], 0.2, atol=1e-5)


def test_convolve_reduces_variance_and_is_nan_aware():
    arr = _textured(40, 50)
    out = convolve_toa_band(arr, 2.0, 2.0)
    assert np.nanvar(out) < np.nanvar(arr)
    # a single NaN must not poison the whole neighbourhood
    arr2 = arr.copy()
    arr2[20, 25] = np.nan
    out2 = convolve_toa_band(arr2, 2.0, 2.0)
    assert np.isfinite(out2[10:30, 10:40]).all()


def test_convolve_zero_sigma_is_identity():
    arr = _textured(30, 30)
    out = convolve_toa_band(arr, 0.0, 5.0)
    assert np.allclose(out, arr.astype("float32"))


def test_convolve_cube_preserves_bands_and_coords():
    cube = _da(np.stack([_textured(40, 50, s) for s in (0, 1)]), bands=["B02", "B04"])
    out = convolve_toa_cube(cube, 1.5, 1.5)
    assert out.dims == ("band", "y", "x")
    assert list(out.coords["band"].values) == ["B02", "B04"]
    assert out.shape == cube.shape


def test_fit_integer_shift_recovers_known_shift():
    prior = _textured(80, 90)
    dy0, dx0 = 3, -2
    toa = np.roll(prior, (dy0, dx0), (0, 1)) + 0.005 * _textured(80, 90, seed=7)
    valid = np.ones((80, 90), bool)
    res = fit_integer_shift(
        toa.astype("float32"), prior, valid, search_radius_px=5, min_correlation=0.6
    )
    # the fit returns the shift that, applied to the TOA, undoes the roll
    assert res.accepted
    assert (res.dx, res.dy) == (-dx0, -dy0)
    assert res.correlation > 0.9


def test_fit_integer_shift_rejects_uncorrelated():
    rng = np.random.RandomState(1)
    toa = rng.rand(80, 90).astype("float32")
    prior = rng.rand(80, 90).astype("float64")
    valid = np.ones((80, 90), bool)
    res = fit_integer_shift(toa, prior, valid, search_radius_px=4, min_correlation=0.6)
    assert not res.accepted
    assert (res.dx, res.dy) == (0, 0)


def test_apply_integer_shift_blanks_wrapped_edges():
    cube = _da(np.ones((1, 20, 24)), bands=["B02"])
    out = apply_integer_shift(cube, dx=3, dy=-2)
    vals = out.isel(band=0).values
    # dy=-2 blanks the last 2 rows; dx=3 blanks the first 3 cols
    assert np.isnan(vals[-2:, :]).all()
    assert np.isnan(vals[:, :3]).all()
    # interior is preserved (no wrap-around bleed)
    assert np.isfinite(vals[2:-2, 3:]).all()


def test_apply_integer_shift_zero_is_noop():
    cube = _da(np.stack([_textured(20, 22)]), bands=["B02"])
    out = apply_integer_shift(cube, 0, 0)
    assert out is cube


def test_orchestrator_disabled_is_noop():
    cube = _da(np.stack([_textured(20, 22)]), bands=["B02"])
    out, res = psf_convolve_and_align_toa(
        cube, None, None, None, grid_resolution_m=120.0, cfg=ToaPsfConfig(enabled=False)
    )
    assert out is cube
    assert res == ShiftFitResult()


def test_orchestrator_convolves_without_refs():
    cube = _da(np.stack([_textured(40, 50)]), bands=["B02"])
    out, res = psf_convolve_and_align_toa(
        cube, None, None, None, grid_resolution_m=120.0, cfg=ToaPsfConfig()
    )
    # convolved (variance reduced) but no shift attempted
    assert np.nanvar(out.values) < np.nanvar(cube.values)
    assert not res.accepted


def test_convolve_cube_2d_input():
    arr = _da(_textured(40, 50))  # (y, x), no band dim
    out = convolve_toa_cube(arr, 1.5, 1.5)
    assert out.dims == ("y", "x")
    assert np.nanvar(out.values) < np.nanvar(arr.values)


def test_pearson_small_array_returns_nan():
    assert np.isnan(_pearson(np.ones(10), np.ones(10)))


def test_pearson_degenerate_constant_returns_nan():
    a = np.zeros(100)
    b = np.arange(100, dtype="float64")
    assert np.isnan(_pearson(a, b))  # zero-variance -> nan


def test_fit_constant_prior_rejects():
    # constant prior -> every candidate degenerate -> no finite correlation -> (0,0)
    prior = np.full((80, 90), 0.1)
    toa = _textured(80, 90).astype("float32")
    res = fit_integer_shift(
        toa, prior, np.ones((80, 90), bool), search_radius_px=3, min_correlation=0.6
    )
    assert not res.accepted
    assert (res.dx, res.dy) == (0, 0)


def test_fit_tiny_grid_insufficient_overlap():
    # grid smaller than the minimum-overlap guard -> every candidate skipped -> (0,0)
    prior = _textured(6, 7)
    toa = prior.astype("float32")
    res = fit_integer_shift(
        toa, prior, np.ones((6, 7), bool), search_radius_px=2, min_correlation=0.6
    )
    assert not res.accepted
    assert (res.dx, res.dy) == (0, 0)


def test_apply_integer_shift_single_axis():
    cube = _da(np.ones((1, 16, 18)), bands=["B02"])
    out_y = apply_integer_shift(cube, dx=0, dy=2).isel(band=0).values
    assert np.isnan(out_y[:2, :]).all() and np.isfinite(out_y[2:, :]).all()
    out_x = apply_integer_shift(cube, dx=-3, dy=0).isel(band=0).values
    assert np.isnan(out_x[:, -3:]).all() and np.isfinite(out_x[:, :-3]).all()


def test_orchestrator_shape_mismatch_skips_shift():
    cube = _da(np.stack([_textured(40, 50)]), bands=["B02"])
    out, res = psf_convolve_and_align_toa(
        cube,
        _da(_textured(40, 50)),
        _da(_textured(30, 30)),
        None,
        grid_resolution_m=120.0,
        cfg=ToaPsfConfig(),
    )
    assert not res.accepted
    assert out.shape == cube.shape


def test_orchestrator_unaccepted_shift_returns_convolved():
    rng = np.random.RandomState(3)
    cube = _da(np.stack([rng.rand(80, 90)]), bands=["B02"])
    out, res = psf_convolve_and_align_toa(
        cube,
        _da(rng.rand(80, 90)),
        _da(rng.rand(80, 90)),
        np.ones((80, 90), bool),
        grid_resolution_m=120.0,
        cfg=ToaPsfConfig(),
    )
    assert not res.accepted
    assert out.shape == cube.shape


def _msi_scene_with_b12(shape=(64, 64), banded_mask=False):
    """Minimal ObservationBundle/SurfacePrior/AtmosphericState with a B12 ref band."""
    from datetime import datetime

    from siac.domain import SensorBand, SensorConfig
    from siac.runtime import AtmosphericState, GeometryAngles, ObservationBundle, SurfacePrior

    cfg = SensorConfig(
        sensor_id="MSI",
        satellite_id="S2A",
        bands=(
            SensorBand("B02", 492.0, 65.0, 10.0, 1),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
            SensorBand("B12", 2186.0, 175.0, 20.0, 12),
        ),
    )
    base = _textured(*shape, seed=5).astype("float32")

    def ones(v):
        return xr.DataArray(np.full(shape, v), dims=["y", "x"])

    toa = xr.Dataset(
        {
            "B02": xr.DataArray(base * 0.5 + 0.05, dims=["y", "x"]),
            "B04": xr.DataArray(base * 0.7 + 0.05, dims=["y", "x"]),
            # B12 TOA shifted vs the prior so a co-registration shift is fittable
            "B12": xr.DataArray(np.roll(base, (1, -1), (0, 1)), dims=["y", "x"]),
        }
    )
    geom = GeometryAngles(sza=ones(0.5), saa=ones(2.5), vza=ones(0.1), vaa=ones(1.5))
    obs = ObservationBundle(
        toa=toa,
        geometry=geom,
        cloud_mask=xr.DataArray(np.zeros(shape, bool), dims=["y", "x"]),
        sensor_config=cfg,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 300640.0, 5500640.0),
    )
    atmo = AtmosphericState(
        aot=ones(0.15),
        tcwv=ones(2.5),
        tco3=ones(0.3),
        aot_unc=ones(0.05),
        tcwv_unc=ones(0.3),
        tco3_unc=ones(0.01),
        elevation=ones(0.1),
    )
    bands = ["B02", "B04", "B12"]
    boa = xr.DataArray(
        np.stack([base, base * 0.9, base]), dims=["band", "y", "x"], coords={"band": bands}
    )
    boa_unc = xr.DataArray(
        np.full((3, *shape), 0.02), dims=["band", "y", "x"], coords={"band": bands}
    )
    if banded_mask:
        mask = xr.DataArray(
            np.ones((3, *shape), bool), dims=["band", "y", "x"], coords={"band": bands}
        )
    else:
        mask = xr.DataArray(np.ones(shape, bool), dims=["y", "x"])
    surface = SurfacePrior(boa=boa, boa_unc=boa_unc, kernels=None, mask=mask)
    return obs, atmo, surface


def test_assemble_grids_observation_psf_drops_b12_and_keeps_solver_bands():
    from siac.algorithms.grid.assembler import assemble_grids

    obs, atmo, surface = _msi_scene_with_b12(banded_mask=True)
    bundle = assemble_grids(
        obs,
        atmo,
        surface,
        object(),  # rt_model passed through, unused by grid assembly
        aerosol_resolution_m=10.0,
        solver_band_names=("B02", "B04"),
        toa_psf_config=ToaPsfConfig(enabled=True, reference_bands=("B12",)),
    )
    # the B12 reference band is used for the shift only — dropped from the solver bundle
    assert list(bundle.toa.coords["band"].values) == ["B02", "B04"]
    assert list(bundle.surface_prior.boa.coords["band"].values) == ["B02", "B04"]
    assert bundle.toa.sizes["y"] == 64 and bundle.toa.sizes["x"] == 64


def test_resolve_psf_shift_reference_branches():
    from dataclasses import replace

    from siac.algorithms.grid.assembler import _resolve_psf_shift_reference

    obs, _, surface = _msi_scene_with_b12((16, 16))
    template = xr.DataArray(
        np.zeros((16, 16)), dims=["y", "x"], coords={"y": np.arange(16), "x": np.arange(16)}
    )

    # happy path: B12 present in both prior and preloaded TOA
    tr, pr = _resolve_psf_shift_reference(obs, surface, ("B12",), (16, 16), template, None)
    assert tr is not None and pr is not None

    # reference band absent from the prior -> (None, None)
    tr, pr = _resolve_psf_shift_reference(obs, surface, ("B11",), (16, 16), template, None)
    assert (tr, pr) == (None, None)

    # non-banded prior -> (None, None)
    flat_surface = replace(surface, boa=surface.boa.isel(band=0, drop=True))
    tr, pr = _resolve_psf_shift_reference(obs, flat_surface, ("B12",), (16, 16), template, None)
    assert (tr, pr) == (None, None)

    # B12 not preloaded but available via the lazy band loader
    obs_no_b12 = replace(obs, toa=obs.toa.drop_vars("B12"))
    b12 = xr.DataArray(_textured(16, 16), dims=["y", "x"])
    tr, pr = _resolve_psf_shift_reference(
        obs_no_b12, surface, ("B12",), (16, 16), template, lambda *_a, **_k: b12
    )
    assert tr is not None and pr is not None

    # loader raises -> band skipped -> (None, None)
    def _raising_loader(name, native=True):
        raise KeyError(name)

    tr, pr = _resolve_psf_shift_reference(
        obs_no_b12, surface, ("B12",), (16, 16), template, _raising_loader
    )
    assert (tr, pr) == (None, None)


def test_assemble_grids_observation_psf_disabled_is_unchanged():
    from siac.algorithms.grid.assembler import assemble_grids

    obs, atmo, surface = _msi_scene_with_b12()
    bundle = assemble_grids(
        obs,
        atmo,
        surface,
        object(),
        aerosol_resolution_m=10.0,
        solver_band_names=("B02", "B04"),
        toa_psf_config=ToaPsfConfig(enabled=False),
    )
    assert list(bundle.toa.coords["band"].values) == ["B02", "B04"]


def test_assemble_grids_psf_grid_mode():
    from siac.algorithms.grid.assembler import assemble_grids

    obs, atmo, surface = _msi_scene_with_b12()
    bundle = assemble_grids(
        obs,
        atmo,
        surface,
        object(),
        aerosol_resolution_m=10.0,
        solver_band_names=("B02", "B04"),
        toa_psf_config=ToaPsfConfig(
            enabled=True, convolve_resolution="grid", reference_bands=("B12",)
        ),
    )
    assert list(bundle.toa.coords["band"].values) == ["B02", "B04"]


def test_assemble_grids_psf_native_uses_band_loader():
    import dataclasses

    from siac.algorithms.grid.assembler import assemble_grids

    obs0, atmo, surface = _msi_scene_with_b12()
    # drop a solver band from the preloaded TOA; supply it via the lazy loader
    b04 = obs0.toa["B04"]
    new_toa = obs0.toa.drop_vars("B04")
    new_toa.attrs["_siac_toa_band_loader"] = lambda *_a, **_k: b04
    obs = dataclasses.replace(obs0, toa=new_toa)
    bundle = assemble_grids(
        obs,
        atmo,
        surface,
        object(),
        aerosol_resolution_m=10.0,
        solver_band_names=("B02", "B04"),
        toa_psf_config=ToaPsfConfig(enabled=True, reference_bands=("B12",)),
    )
    assert list(bundle.toa.coords["band"].values) == ["B02", "B04"]


def test_orchestrator_convolves_and_shifts_with_refs():
    prior = _textured(80, 90)
    dy0, dx0 = 2, -3
    toa_ref = np.roll(prior, (dy0, dx0), (0, 1))
    cube = _da(np.stack([toa_ref, np.roll(prior, (dy0, dx0), (0, 1))]), bands=["B02", "B04"])
    out, res = psf_convolve_and_align_toa(
        cube,
        _da(toa_ref),
        _da(prior),
        np.ones((80, 90), bool),
        grid_resolution_m=120.0,
        cfg=ToaPsfConfig(shift_search_radius_m=600.0),
    )
    assert res.accepted
    assert (res.dx, res.dy) == (-dx0, -dy0)
    assert out.shape == cube.shape

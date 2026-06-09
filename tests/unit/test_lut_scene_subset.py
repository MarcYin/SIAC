"""Unit tests for scene-aware spectral LUT subsetting helpers."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.lut.backend import ZarrLUTBackend
from siac.domain import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles


def _spectral_lut() -> xr.Dataset:
    shape = (3, 3, 3, 4, 4, 5, 6, 7)
    dims = ("sza", "vza", "raa", "ozone", "altitude", "aot", "tcwv", "wavelength")
    coords = {
        "sza": np.array([10.0, 20.0, 30.0], dtype=np.float32),
        "vza": np.array([0.0, 10.0, 20.0], dtype=np.float32),
        "raa": np.array([0.0, 90.0, 180.0], dtype=np.float32),
        "ozone": np.array([250.0, 280.0, 310.0, 340.0], dtype=np.float32),
        "altitude": np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32),
        "aot": np.linspace(0.01, 1.0, 5, dtype=np.float32),
        "tcwv": np.linspace(0.1, 6.0, 6, dtype=np.float32),
        "wavelength": np.linspace(400.0, 700.0, 7, dtype=np.float32),
    }
    data = np.random.default_rng(0).random(shape, dtype=np.float32)
    return xr.Dataset(
        {
            "TOA_rho1": (dims, data),
            "TOA_rho2": (dims, data + 0.01),
            "Eg_rho1": (dims, data + 0.02),
            "Eg_rho2": (dims, data + 0.03),
        },
        coords=coords,
    )


def test_subset_spectral_lut_uses_angle_means_and_range_filters():
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()

    subset = backend._subset_spectral_lut_for_scene(
        lut,
        sza=np.full((2, 2), 19.5, dtype=np.float32),
        vza=np.full((2, 2), 9.8, dtype=np.float32),
        raa=np.full((2, 2), 88.0, dtype=np.float32),
        tco3=np.array([[270.0, 300.0], [290.0, 295.0]], dtype=np.float32),
        elevation=np.array([[0.2, 1.8], [1.2, 0.7]], dtype=np.float32),
    )

    # Scene-angle means are selected up front.
    assert "sza" not in subset.dims
    assert "vza" not in subset.dims
    assert "raa" not in subset.dims
    # Ozone/elevation are constrained and then reduced.
    assert "ozone" not in subset.dims
    assert "altitude" not in subset.dims
    # Dimensions needed for retrieval remain.
    assert "aot" in subset.dims
    assert "tcwv" in subset.dims
    assert "wavelength" in subset.dims


def test_build_aot_tcwv_point_coords_clips_to_lut_axes():
    backend = ZarrLUTBackend("dummy")
    lut = xr.Dataset(
        coords={
            "aot": np.array([0.01, 0.5, 1.0], dtype=np.float32),
            "tcwv": np.array([0.1, 2.0, 6.0], dtype=np.float32),
        }
    )

    coords = backend._build_aot_tcwv_point_coords(
        lut,
        aot=np.array([-1.0, 0.2, 5.0], dtype=np.float32),
        tcwv=np.array([-3.0, 3.0, 10.0], dtype=np.float32),
    )

    assert coords["aot"].dims == ("point",)
    assert coords["tcwv"].dims == ("point",)
    np.testing.assert_allclose(coords["aot"].values, np.array([0.01, 0.2, 1.0], dtype=np.float32))
    np.testing.assert_allclose(coords["tcwv"].values, np.array([0.1, 3.0, 6.0], dtype=np.float32))


def test_build_aot_tcwv_point_coords_rejects_nonfinite_values():
    backend = ZarrLUTBackend("dummy")
    lut = xr.Dataset(
        coords={
            "aot": np.array([0.01, 0.5, 1.0], dtype=np.float32),
            "tcwv": np.array([0.1, 2.0, 6.0], dtype=np.float32),
        }
    )

    with pytest.raises(ValueError, match="aot must contain only finite values"):
        backend._build_aot_tcwv_point_coords(
            lut,
            aot=np.array([np.nan, np.inf, -np.inf], dtype=np.float32),
            tcwv=np.array([np.nan, np.inf, -np.inf], dtype=np.float32),
        )


def test_subset_wavelength_for_band_trims_lut():
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()
    band = SensorBand(
        name="B02",
        center_wavelength=490.0,
        bandwidth=65.0,
        resolution=10.0,
        band_index=1,
    )

    subset = backend._subset_wavelength_for_band(lut, band)

    assert "wavelength" in subset.dims
    assert subset.sizes["wavelength"] < lut.sizes["wavelength"]


def test_subset_spectral_lut_for_scene_handles_nonfinite_optional_scene_inputs():
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()

    subset = backend._subset_spectral_lut_for_scene(
        lut,
        sza=np.full((2, 2), 20.0, dtype=np.float32),
        vza=np.full((2, 2), 10.0, dtype=np.float32),
        raa=np.full((2, 2), 90.0, dtype=np.float32),
        tco3=np.full((2, 2), np.nan, dtype=np.float32),
        elevation=np.full((2, 2), np.nan, dtype=np.float32),
    )

    assert "sza" not in subset.dims
    assert "vza" not in subset.dims
    assert "raa" not in subset.dims
    assert "ozone" not in subset.dims
    assert "altitude" not in subset.dims
    assert "aot" in subset.dims
    assert "tcwv" in subset.dims


def test_subset_spectral_lut_for_scene_rejects_nonfinite_required_angles():
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()

    with pytest.raises(ValueError, match="sza must contain only finite values"):
        backend._subset_spectral_lut_for_scene(
            lut,
            sza=np.full((2, 2), np.nan, dtype=np.float32),
            vza=np.full((2, 2), 10.0, dtype=np.float32),
            raa=np.full((2, 2), 90.0, dtype=np.float32),
            tco3=np.full((2, 2), 300.0, dtype=np.float32),
            elevation=np.full((2, 2), 1.0, dtype=np.float32),
        )


def test_subset_wavelength_for_band_uses_tabulated_srf_support():
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut().assign_coords(
        wavelength=np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32)
    )
    band = SensorBand(
        name="B02",
        center_wavelength=460.0,
        bandwidth=80.0,
        resolution=10.0,
        band_index=1,
        rsrf_wavelengths_nm=np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
        rsrf_response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )

    subset = backend._subset_wavelength_for_band(lut, band)

    np.testing.assert_allclose(
        subset.coords["wavelength"].values,
        np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32),
    )


def test_compute_coefficients_spectral_interpolates_after_wavelength_collapse(monkeypatch):
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()
    backend._lut = lut
    backend._lut_coords = {name: coord.values for name, coord in lut.coords.items()}
    calls: list[tuple[str, ...]] = []
    original = backend._interpolate_variable_fast

    def _record_interp_dims(var: xr.DataArray, coords: dict[str, xr.DataArray], method: str):
        calls.append(var.dims)
        return original(var, coords, method)

    monkeypatch.setattr(backend, "_interpolate_variable_fast", _record_interp_dims)

    band = SensorBand(
        name="B02",
        center_wavelength=490.0,
        bandwidth=65.0,
        resolution=10.0,
        band_index=1,
    )
    aot = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32)
    tcwv = np.array([[1.0, 1.2], [1.4, 1.6]], dtype=np.float32)

    xap, xbp, xcp = backend._compute_coefficients_from_spectral_lut(
        sza=np.full_like(aot, 19.5),
        vza=np.full_like(aot, 9.8),
        raa=np.full_like(aot, 88.0),
        aot=aot,
        tcwv=tcwv,
        tco3=np.full_like(aot, 290.0),
        elevation=np.full_like(aot, 1.0),
        band=band,
    )

    assert xap.shape == aot.shape
    assert xbp.shape == aot.shape
    assert xcp.shape == aot.shape
    assert calls
    assert all("wavelength" not in dims for dims in calls)


def test_spectral_integration_weights_use_tabulated_srf_over_gaussian():
    backend = ZarrLUTBackend("dummy")
    lut = xr.Dataset(
        coords={
            "wavelength": np.array(
                [430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32
            ),
        }
    )
    band = SensorBand(
        name="B02",
        center_wavelength=460.0,
        bandwidth=200.0,
        resolution=10.0,
        band_index=1,
        rsrf_wavelengths_nm=np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
        rsrf_response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )

    weights = backend._spectral_integration_weights(band, lut)

    assert weights.dims == ("wavelength",)
    np.testing.assert_allclose(
        weights.coords["wavelength"].values,
        np.array([430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32),
    )
    assert float(weights.sel(wavelength=430.0).values) == pytest.approx(0.0)
    assert float(weights.sel(wavelength=440.0).values) == pytest.approx(0.0)
    assert float(weights.sel(wavelength=480.0).values) == pytest.approx(0.0)
    assert float(weights.sel(wavelength=490.0).values) == pytest.approx(0.0)
    assert float(weights.sel(wavelength=460.0).values) > float(weights.sel(wavelength=450.0).values)


def test_weighted_spectral_mean_zero_fills_missing_srf_support():
    data = xr.DataArray(
        np.array([10.0, 1.0, 2.0, 4.0, 2.0, 1.0, 9.0], dtype=np.float32),
        dims=("wavelength",),
        coords={
            "wavelength": np.array(
                [430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32
            ),
        },
    )
    weights = xr.DataArray(
        np.array([0.0, 0.0, 0.25, 0.5, 0.25, 0.0, 0.0], dtype=np.float32),
        dims=("wavelength",),
        coords={
            "wavelength": np.array(
                [430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0], dtype=np.float32
            ),
        },
    )

    result = ZarrLUTBackend._weighted_spectral_mean(data, weights)

    assert float(result.values) == pytest.approx(3.0)


def test_compute_coefficients_retries_transient_lut_io(monkeypatch):
    backend = ZarrLUTBackend("dummy")
    geom_arr = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=("y", "x"))
    geometry = GeometryAngles(sza=geom_arr, saa=geom_arr, vza=geom_arr, vaa=geom_arr)
    atmo = AtmosphericState(
        aot=geom_arr * 0.1,
        tcwv=geom_arr * 1.5,
        tco3=geom_arr * 0.3,
        aot_unc=geom_arr * 0.02,
        tcwv_unc=geom_arr * 0.1,
        tco3_unc=geom_arr * 0.01,
        elevation=geom_arr * 0.0,
    )
    band = SensorBand(
        name="B02",
        center_wavelength=490.0,
        bandwidth=65.0,
        resolution=10.0,
        band_index=1,
    )

    calls = {"n": 0, "reloads": 0}

    def _fake_compute(*args, **kwargs):
        del args, kwargs
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("ReferenceNotReachable: temporary disconnect")
        out = np.ones((2, 2), dtype=np.float32)
        return out, out * 2.0, out * 3.0

    monkeypatch.setattr(backend, "_supports_coefficient_lut", lambda: False)
    monkeypatch.setattr(backend, "_supports_spectral_lut", lambda: True)
    monkeypatch.setattr(backend, "_compute_coefficients_from_spectral_lut", _fake_compute)
    monkeypatch.setattr(
        backend, "_reload_lut", lambda: calls.__setitem__("reloads", calls["reloads"] + 1)
    )

    coeffs = backend.compute_coefficients(geometry, atmo, band, compute_jacobian=False)

    assert calls["n"] == 2
    assert calls["reloads"] == 1
    np.testing.assert_allclose(coeffs.xap.values, np.ones((2, 2), dtype=np.float32))


def test_preload_scene_subset_retries_transient_lut_io(monkeypatch):
    backend = ZarrLUTBackend("dummy")
    geom_arr = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=("y", "x"))
    geometry = GeometryAngles(sza=geom_arr, saa=geom_arr, vza=geom_arr, vaa=geom_arr)
    atmo = AtmosphericState(
        aot=geom_arr * 0.1,
        tcwv=geom_arr * 1.5,
        tco3=geom_arr * 0.3,
        aot_unc=geom_arr * 0.02,
        tcwv_unc=geom_arr * 0.1,
        tco3_unc=geom_arr * 0.01,
        elevation=geom_arr * 0.0,
    )
    band = SensorBand(
        name="B02",
        center_wavelength=490.0,
        bandwidth=65.0,
        resolution=10.0,
        band_index=1,
    )

    calls = {"n": 0, "reloads": 0}

    def _fake_preload(*, geometry, atmo_state, bands):
        del geometry, atmo_state, bands
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("ServerDisconnectedError: temporary disconnect")

    monkeypatch.setattr(backend, "_preload_scene_subset_once", _fake_preload)
    monkeypatch.setattr(
        backend, "_reload_lut", lambda: calls.__setitem__("reloads", calls["reloads"] + 1)
    )

    backend.preload_scene_subset(geometry, atmo, [band])

    assert calls["n"] == 2
    assert calls["reloads"] == 1


def test_scene_subset_disk_cache_reused_across_backends(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A re-run loads the materialised scene subset from disk (no re-fetch)."""
    lut = _spectral_lut()
    kw = {
        "sza": np.full((2, 2), 19.5, dtype=np.float32),
        "vza": np.full((2, 2), 9.8, dtype=np.float32),
        "raa": np.full((2, 2), 88.0, dtype=np.float32),
        "tco3": np.full((2, 2), 300.0, dtype=np.float32),
        "elevation": np.full((2, 2), 0.5, dtype=np.float32),
    }
    b1 = ZarrLUTBackend("dummy.zarr.zip", scene_cache_dir=tmp_path)
    b1._lut = lut
    key1, sub1 = b1._get_or_build_spectral_scene_subset(**kw)
    path = b1._scene_subset_cache_path(key1)
    assert path is not None and path.exists()

    # Fresh backend, same cache dir: it must load from disk, not re-subset.
    b2 = ZarrLUTBackend("dummy.zarr.zip", scene_cache_dir=tmp_path)
    b2._lut = lut
    monkeypatch.setattr(
        b2,
        "_subset_spectral_lut_for_scene",
        lambda *_a, **_k: pytest.fail("subset rebuilt despite a warm disk cache"),
    )
    key2, sub2 = b2._get_or_build_spectral_scene_subset(**kw)
    assert key2 == key1
    for var in ("TOA_rho1", "Eg_rho1"):
        np.testing.assert_allclose(sub2[var].values, sub1[var].values)


def test_scene_subset_disk_cache_disabled() -> None:
    backend = ZarrLUTBackend("dummy.zarr.zip", scene_cache_enabled=False)
    assert backend._scene_cache_dir is None
    assert backend._scene_subset_cache_path((1.0, 2.0, 3.0)) is None

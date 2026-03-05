"""Unit tests for scene-aware spectral LUT subsetting helpers."""

from __future__ import annotations

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState, GeometryAngles, SensorBand
from siac.rt.lut.backend import ZarrLUTBackend


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


def test_compute_coefficients_spectral_interpolates_after_wavelength_collapse(monkeypatch):
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()
    backend._lut = lut
    backend._lut_coords = {name: coord.values for name, coord in lut.coords.items()}
    calls: list[tuple[str, ...]] = []
    original = backend._interpolate_variable

    def _record_interp_dims(var: xr.DataArray, coords: dict[str, xr.DataArray], method: str) -> xr.DataArray:
        calls.append(var.dims)
        return original(var, coords, method)

    monkeypatch.setattr(backend, "_interpolate_variable", _record_interp_dims)

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
    monkeypatch.setattr(backend, "_reload_lut", lambda: calls.__setitem__("reloads", calls["reloads"] + 1))

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
    monkeypatch.setattr(backend, "_reload_lut", lambda: calls.__setitem__("reloads", calls["reloads"] + 1))

    backend.preload_scene_subset(geometry, atmo, [band])

    assert calls["n"] == 2
    assert calls["reloads"] == 1

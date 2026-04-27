from __future__ import annotations

from threading import Lock

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.lut.backend import ZarrLUTBackend
from siac.domain import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles


def _spectral_lut() -> xr.Dataset:
    shape = (2, 2, 2, 2, 2, 3, 3, 4)
    dims = ("sza", "vza", "raa", "ozone", "altitude", "aot", "tcwv", "wavelength")
    coords = {
        "sza": np.array([10.0, 20.0], dtype=np.float32),
        "vza": np.array([0.0, 10.0], dtype=np.float32),
        "raa": np.array([0.0, 90.0], dtype=np.float32),
        "ozone": np.array([250.0, 300.0], dtype=np.float32),
        "altitude": np.array([0.0, 1.0], dtype=np.float32),
        "aot": np.array([0.1, 0.5, 1.0], dtype=np.float32),
        "tcwv": np.array([0.5, 2.0, 4.0], dtype=np.float32),
        "wavelength": np.array([450.0, 550.0, 650.0, 750.0], dtype=np.float32),
    }
    base = np.arange(np.prod(shape), dtype=np.float32).reshape(shape) / 100.0
    return xr.Dataset(
        {
            "TOA_rho1": (dims, base),
            "TOA_rho2": (dims, base + 0.1),
            "Eg_rho1": (dims, base + 0.2),
            "Eg_rho2": (dims, base + 0.3),
            "solar_irradiance": (("wavelength",), np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)),
        },
        coords=coords,
        attrs={"rho1": 0.11, "rho2": 0.44},
    )


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    arr = xr.DataArray(np.ones(shape, dtype=np.float32), dims=("y", "x"))
    return GeometryAngles(sza=arr * 20.0, saa=arr * 10.0, vza=arr * 10.0, vaa=arr * 100.0)


def _atmo(shape: tuple[int, int]) -> AtmosphericState:
    arr = xr.DataArray(np.ones(shape, dtype=np.float32), dims=("y", "x"))
    return AtmosphericState(
        aot=arr * 0.2,
        tcwv=arr * 1.5,
        tco3=arr * 0.3,
        aot_unc=arr * 0.01,
        tcwv_unc=arr * 0.1,
        tco3_unc=arr * 0.01,
        elevation=arr * 0.0,
    )


def test_spectral_scene_subset_cache_hit_and_reset(monkeypatch: pytest.MonkeyPatch) -> None:
    backend = ZarrLUTBackend("dummy")
    lut = _spectral_lut()
    backend._lut = lut
    backend._lut_coords = {name: coord.values for name, coord in lut.coords.items()}
    backend._cache_lock = Lock()
    calls = {"count": 0}
    original = backend._subset_spectral_lut_for_scene

    def _record_subset(dataset: xr.Dataset, **kwargs) -> xr.Dataset:
        calls["count"] += 1
        return original(dataset, **kwargs)

    monkeypatch.setattr(backend, "_subset_spectral_lut_for_scene", _record_subset)

    scene_args = {
        "sza": np.full((2, 2), 20.0, dtype=np.float32),
        "vza": np.full((2, 2), 10.0, dtype=np.float32),
        "raa": np.full((2, 2), 90.0, dtype=np.float32),
        "tco3": np.full((2, 2), 0.3, dtype=np.float32),
        "elevation": np.full((2, 2), 0.0, dtype=np.float32),
    }

    key0, subset0 = backend._get_or_build_spectral_scene_subset(**scene_args)
    key1, subset1 = backend._get_or_build_spectral_scene_subset(**scene_args)

    assert key0 == key1
    assert subset0.identical(subset1)
    assert calls["count"] == 1

    backend._spectral_band_grid_cache[("stale", ("B03", 560.0, 35.0))] = (
        xr.DataArray(1.0),
        xr.DataArray(1.0),
        xr.DataArray(1.0),
        xr.DataArray(1.0),
    )
    backend._scene_subset_logged = True
    backend._get_or_build_spectral_scene_subset(
        **{**scene_args, "elevation": np.full((2, 2), 1.0, dtype=np.float32)}
    )

    assert calls["count"] == 2
    assert backend._spectral_band_grid_cache == {}
    assert backend._scene_subset_logged is False


def test_spectral_band_grid_cache_reuses_existing_grids(monkeypatch: pytest.MonkeyPatch) -> None:
    backend = ZarrLUTBackend("dummy")
    backend._cache_lock = Lock()
    lut = _spectral_lut()
    band = SensorBand("B03", 560.0, 35.0, 10.0, 0)
    calls = {"subset": 0, "weights": 0, "weighted_mean": 0}

    monkeypatch.setattr(
        backend,
        "_subset_wavelength_for_band",
        lambda dataset, _band: calls.__setitem__("subset", calls["subset"] + 1) or dataset,
    )
    monkeypatch.setattr(
        backend,
        "_spectral_integration_weights",
        lambda _band, _lut: (
            calls.__setitem__("weights", calls["weights"] + 1)
            or xr.DataArray(
                np.ones(4, dtype=np.float32),
                dims=("wavelength",),
                coords={"wavelength": lut.coords["wavelength"].values},
            )
        ),
    )

    def _record_weighted_mean(data: xr.DataArray, _weights: xr.DataArray) -> xr.DataArray:
        calls["weighted_mean"] += 1
        return data.mean("wavelength")

    monkeypatch.setattr(
        "siac.algorithms.rt.lut.backend.weighted_spectral_mean",
        _record_weighted_mean,
    )

    key = (20.0, 10.0, 90.0, 0.3, 0.3, 0.0, 0.0)
    grids0 = backend._get_or_build_spectral_band_grids(key, lut, band)
    grids1 = backend._get_or_build_spectral_band_grids(key, lut, band)

    assert grids0 is grids1
    assert calls == {"subset": 1, "weights": 1, "weighted_mean": 4}


def test_run_with_transient_lut_io_retry_failure_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    backend = ZarrLUTBackend("dummy")
    backend._SPECTRAL_IO_MAX_RETRIES = 2
    monkeypatch.setattr("siac.algorithms.rt.lut.backend.time.sleep", lambda _seconds: None)

    reloads = {"count": 0}
    monkeypatch.setattr(
        backend, "_reload_lut", lambda: reloads.__setitem__("count", reloads["count"] + 1)
    )

    def _not_transient(_exc: Exception) -> bool:
        return False

    def _transient(_exc: Exception) -> bool:
        return True

    monkeypatch.setattr(backend, "_is_transient_lut_io_error", _not_transient)
    with pytest.raises(RuntimeError, match="hard failure"):
        backend._run_with_transient_lut_io_retry(
            lambda: (_ for _ in ()).throw(RuntimeError("hard failure")), operation="read"
        )
    assert reloads["count"] == 0

    monkeypatch.setattr(backend, "_is_transient_lut_io_error", _transient)
    with pytest.raises(RuntimeError, match="still failing"):
        backend._run_with_transient_lut_io_retry(
            lambda: (_ for _ in ()).throw(RuntimeError("still failing")), operation="read"
        )
    assert reloads["count"] == 1


def test_build_point_coords_weight_helpers_and_reference_attrs() -> None:
    backend = ZarrLUTBackend("dummy")

    empty_axis_lut = xr.Dataset(
        coords={
            "aot": np.array([], dtype=np.float32),
            "tcwv": np.array([0.5, 2.0], dtype=np.float32),
        }
    )
    coords = backend._build_aot_tcwv_point_coords(
        empty_axis_lut,
        aot=np.array([0.2, 0.4], dtype=np.float32),
        tcwv=np.array([1.0, 3.0], dtype=np.float32),
    )
    np.testing.assert_allclose(coords["aot"].values, np.array([0.2, 0.4], dtype=np.float32))
    np.testing.assert_allclose(coords["tcwv"].values, np.array([1.0, 2.0], dtype=np.float32))

    band = SensorBand("B03", 560.0, 35.0, 10.0, 0)
    no_wavelength = xr.Dataset()
    with pytest.raises(ValueError, match="wavelength coordinate"):
        backend._spectral_integration_weights(band, no_wavelength)

    lut = _spectral_lut()
    weights = backend._spectral_integration_weights(band, lut)
    assert weights.dims == ("wavelength",)
    assert float(weights.sel(wavelength=550.0).values) > float(weights.sel(wavelength=450.0).values)
    assert float(weights.sel(wavelength=550.0).values) > float(weights.sel(wavelength=650.0).values)

    backend._lut = xr.Dataset()
    assert backend._spectral_reference_reflectances() == (0.15, 0.5)
    backend._lut = lut
    assert backend._spectral_reference_reflectances() == (0.11, 0.44)


def test_band_rsrf_and_weighted_mean_edge_paths() -> None:
    backend = ZarrLUTBackend("dummy")

    with pytest.raises(ValueError, match="does not define a tabulated RSRF"):
        backend._band_rsrf(SensorBand("B03", 560.0, 35.0, 10.0, 0))

    passthrough = xr.DataArray(np.array([1.0, 2.0], dtype=np.float32), dims=("point",))
    same = backend._weighted_spectral_mean(
        passthrough,
        xr.DataArray(
            np.array([1.0], dtype=np.float32), dims=("wavelength",), coords={"wavelength": [550.0]}
        ),
    )
    assert same.identical(passthrough)

    single = xr.DataArray(
        np.array([4.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": [550.0]},
    )
    weights = xr.DataArray(
        np.array([0.0], dtype=np.float32),
        dims=("wavelength",),
        coords={"wavelength": [550.0]},
    )
    result = backend._weighted_spectral_mean(single, weights)
    assert np.isfinite(float(result.values))


def test_preload_scene_subset_once_short_circuits_without_spectral_support(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _unexpected_grid_validation(
        _geometry: GeometryAngles, _atmo_state: AtmosphericState
    ) -> None:
        raise AssertionError("should not validate")

    backend = ZarrLUTBackend("dummy")
    backend._lut = _spectral_lut()
    backend._lut_coords = {name: coord.values for name, coord in backend._lut.coords.items()}
    monkeypatch.setattr(backend, "_supports_spectral_lut", lambda: False)
    monkeypatch.setattr(backend, "_require_matching_grid_shapes", _unexpected_grid_validation)

    backend._preload_scene_subset_once(
        geometry=_geometry((2, 2)),
        atmo_state=_atmo((2, 2)),
        bands=[SensorBand("B03", 560.0, 35.0, 10.0, 0)],
    )

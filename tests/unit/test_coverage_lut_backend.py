"""Coverage tests for the Zarr LUT backend."""

from __future__ import annotations

import asyncio
import json
import shutil
import struct
import sys
import zipfile
from io import BytesIO
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.lut import ZarrLUTBackend, create_lut_from_py6s
from siac.algorithms.rt.lut.http_zip_store import _HTTPRangeFileSystem, _ReadOnlyZipFileSystem
from siac.domain import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles

if TYPE_CHECKING:
    from pathlib import Path


def _geometry(shape: tuple[int, int] = (3, 3)) -> GeometryAngles:
    sza = xr.DataArray(np.full(shape, np.deg2rad(30.0), dtype=np.float32), dims=["y", "x"])
    saa = xr.DataArray(np.full(shape, np.deg2rad(150.0), dtype=np.float32), dims=["y", "x"])
    vza = xr.DataArray(np.full(shape, np.deg2rad(10.0), dtype=np.float32), dims=["y", "x"])
    vaa = xr.DataArray(np.full(shape, np.deg2rad(60.0), dtype=np.float32), dims=["y", "x"])
    return GeometryAngles(sza=sza, saa=saa, vza=vza, vaa=vaa)


def _atmo(shape: tuple[int, int] = (3, 3)) -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
    )


def _write_small_lut(path: Path, consolidated: bool = True) -> Path:
    sza = np.array([0.0, 30.0, 60.0], dtype=np.float32)
    vza = np.array([0.0, 10.0, 30.0], dtype=np.float32)
    raa = np.array([0.0, 90.0, 180.0], dtype=np.float32)
    aot = np.array([0.1, 0.5], dtype=np.float32)
    tcwv = np.array([1.0, 3.0], dtype=np.float32)
    wl = np.array([490.0, 560.0], dtype=np.float32)

    shape = (len(sza), len(vza), len(raa), len(aot), len(tcwv), len(wl))
    # Deterministic values with non-zero transmittance
    base = np.ones(shape, dtype=np.float32)
    path_ref = base * 0.02
    trans_down = base * 0.85
    trans_up = base * 0.90
    sph_alb = base * 0.03

    ds = xr.Dataset(
        {
            "path_reflectance": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], path_ref),
            "transmittance_down": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_down),
            "transmittance_up": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_up),
            "spherical_albedo": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], sph_alb),
        },
        coords={"sza": sza, "vza": vza, "raa": raa, "aot": aot, "tcwv": tcwv, "wavelength": wl},
    )
    ds.to_zarr(path, mode="w", consolidated=consolidated)
    return path


def _write_small_spectral_lut(path: Path, *, include_solar: bool = True) -> tuple[Path, dict[str, float]]:
    wavelength = np.array([480.0, 490.0, 500.0], dtype=np.float32)
    sza = np.array([30.0], dtype=np.float32)
    vza = np.array([10.0], dtype=np.float32)
    raa = np.array([90.0], dtype=np.float32)
    aot = np.array([0.2], dtype=np.float32)
    tcwv = np.array([2.0], dtype=np.float32)
    ozone = np.array([300.0], dtype=np.float32)
    altitude = np.array([0.1], dtype=np.float32)

    rho1 = 0.15
    rho2 = 0.5
    toa1 = 0.05
    toa2 = 0.29
    eg1 = 0.20
    eg2 = 0.50

    toa_shape = (
        len(wavelength),
        len(sza),
        len(vza),
        len(raa),
        len(aot),
        len(tcwv),
        len(ozone),
        len(altitude),
    )
    eg_shape = (
        len(wavelength),
        len(sza),
        len(aot),
        len(tcwv),
        len(ozone),
        len(altitude),
    )
    data_vars: dict[str, tuple[list[str], np.ndarray]] = {
        "TOA_rho1": (
            ["wavelength", "sza", "vza", "raa", "aot", "tcwv", "ozone", "altitude"],
            np.full(toa_shape, toa1, dtype=np.float32),
        ),
        "TOA_rho2": (
            ["wavelength", "sza", "vza", "raa", "aot", "tcwv", "ozone", "altitude"],
            np.full(toa_shape, toa2, dtype=np.float32),
        ),
        "Eg_rho1": (
            ["wavelength", "sza", "aot", "tcwv", "ozone", "altitude"],
            np.full(eg_shape, eg1, dtype=np.float32),
        ),
        "Eg_rho2": (
            ["wavelength", "sza", "aot", "tcwv", "ozone", "altitude"],
            np.full(eg_shape, eg2, dtype=np.float32),
        ),
    }
    if include_solar:
        data_vars["extraterrestrial_solar_irradiance"] = (
            ["wavelength"],
            np.array([1.0, 2.0, 3.0], dtype=np.float32),
        )

    ds = xr.Dataset(
        data_vars,
        coords={
            "wavelength": wavelength,
            "sza": sza,
            "vza": vza,
            "raa": raa,
            "aot": aot,
            "tcwv": tcwv,
            "ozone": ozone,
            "altitude": altitude,
        },
        attrs={"rho1": rho1, "rho2": rho2},
    )
    ds.to_zarr(path, mode="w", consolidated=True)

    denom = rho2 * eg2 - rho1 * eg1
    s_term = (eg2 - eg1) / denom
    path_ref = (toa2 * rho1 * eg1 - toa1 * rho2 * eg2) / (rho1 * eg1 - rho2 * eg2)
    t_up = (toa2 - toa1) / denom
    eg0 = eg1 * (1.0 - rho1 * s_term)
    t_total = eg0 * t_up
    expected = {
        "xap": 1.0 / t_total,
        "xbp": path_ref / t_total,
        "xcp": s_term,
        "path_ref": path_ref,
        "t_total": t_total,
    }
    return path, expected


class TestZarrLUTBackend:
    def test_invalid_interpolation_method_raises(self):
        with pytest.raises(ValueError, match="interpolation_method"):
            ZarrLUTBackend("unused.zarr", interpolation_method="cubic")

    def test_missing_lut_raises(self, tmp_path: Path):
        b = ZarrLUTBackend(tmp_path / "missing.zarr")
        with pytest.raises(FileNotFoundError):
            _ = b.lut

    def test_load_lut_falls_back_from_v3_and_consolidated(self, monkeypatch):
        import siac.algorithms.rt.lut.backend as lut_backend

        calls: list[dict[str, object]] = []
        dataset = xr.Dataset(coords={"sza": [30.0], "vza": [10.0], "aot": [0.2]})

        def _fake_open_zarr(*, store, consolidated, zarr_format=None):  # noqa: ANN001
            calls.append(
                {
                    "store": store,
                    "consolidated": consolidated,
                    "zarr_format": zarr_format,
                }
            )
            if len(calls) < 3:
                raise RuntimeError("expected fallback")
            return dataset

        monkeypatch.setattr(lut_backend, "as_local_path", lambda _path: None)
        monkeypatch.setattr(
            lut_backend,
            "build_lut_store",
            lambda _path, _options: {"zarr.json": b"{}"},
        )
        monkeypatch.setattr(lut_backend.xr, "open_zarr", _fake_open_zarr)

        backend = ZarrLUTBackend("https://example.com/lut.zarr")
        loaded = backend.lut

        assert loaded is dataset
        assert calls == [
            {
                "store": {"zarr.json": b"{}"},
                "consolidated": False,
                "zarr_format": 3,
            },
            {
                "store": {"zarr.json": b"{}"},
                "consolidated": True,
                "zarr_format": None,
            },
            {
                "store": {"zarr.json": b"{}"},
                "consolidated": False,
                "zarr_format": None,
            },
        ]
        assert np.allclose(backend._lut_coords["sza"], [30.0])

    def test_store_contains_key_handles_contains_errors(self):
        class _BrokenContains:
            def __contains__(self, item):  # noqa: ANN001, ARG002
                raise RuntimeError("broken")

        assert ZarrLUTBackend._store_contains_key(_BrokenContains(), "zarr.json") is False

    def test_build_point_coords_sanitizes_optional_axes_only(self):
        backend = ZarrLUTBackend("dummy")
        backend._lut_coords = {
            "sza": np.array([0.0, 30.0, 60.0], dtype=np.float32),
            "vza": np.array([0.0, 10.0, 20.0], dtype=np.float32),
            "raa": np.array([0.0, 90.0, 180.0], dtype=np.float32),
            "aot": np.array([0.1, 0.5], dtype=np.float32),
            "tcwv": np.array([1.0, 3.0], dtype=np.float32),
            "ozone": np.array([250.0, 300.0, 350.0], dtype=np.float32),
            "altitude": np.array([0.0, 1.0, 2.0], dtype=np.float32),
        }

        coords = backend._build_point_coords(
            sza=np.array([15.0], dtype=np.float32),
            vza=np.array([12.0], dtype=np.float32),
            raa=np.array([100.0], dtype=np.float32),
            aot=np.array([0.2], dtype=np.float32),
            tcwv=np.array([2.0], dtype=np.float32),
            tco3=np.array([np.nan], dtype=np.float32),
            elevation=np.array([np.nan], dtype=np.float32),
        )

        assert float(coords["sza"].values[0]) == pytest.approx(15.0)
        assert float(coords["vza"].values[0]) == pytest.approx(12.0)
        assert float(coords["raa"].values[0]) == pytest.approx(100.0)
        assert float(coords["aot"].values[0]) == pytest.approx(0.2)
        assert float(coords["tcwv"].values[0]) == pytest.approx(2.0)
        assert float(coords["ozone"].values[0]) == pytest.approx(300.0)
        assert float(coords["altitude"].values[0]) == pytest.approx(1.0)

    def test_build_point_coords_rejects_nonfinite_required_values(self):
        backend = ZarrLUTBackend("dummy")
        backend._lut_coords = {
            "sza": np.array([0.0, 30.0, 60.0], dtype=np.float32),
            "vza": np.array([0.0, 10.0, 20.0], dtype=np.float32),
            "raa": np.array([0.0, 90.0, 180.0], dtype=np.float32),
            "aot": np.array([0.1, 0.5], dtype=np.float32),
            "tcwv": np.array([1.0, 3.0], dtype=np.float32),
        }

        with pytest.raises(ValueError, match="sza must contain only finite values"):
            backend._build_point_coords(
                sza=np.array([np.nan], dtype=np.float32),
                vza=np.array([10.0], dtype=np.float32),
                raa=np.array([90.0], dtype=np.float32),
                aot=np.array([0.2], dtype=np.float32),
                tcwv=np.array([2.0], dtype=np.float32),
                tco3=np.array([0.3], dtype=np.float32),
                elevation=np.array([0.1], dtype=np.float32),
            )

    def test_linear_and_nearest_paths(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut.zarr")
        geom = _geometry()
        atmo = _atmo()
        band = SensorBand("B02", 490.0, 65.0, 10.0, 0)

        b_lin = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs_lin = b_lin.compute_coefficients(geom, atmo, band, compute_jacobian=False)
        assert coeffs_lin.xap.shape == (3, 3)
        assert float(coeffs_lin.xap.mean()) > 0.0
        assert b_lin.backend_name == "lut"
        assert not b_lin.supports_jacobian()
        assert b_lin.is_available_for_sensor("ANY", "ANY")
        assert b_lin.get_available_wavelengths().size == 2

        b_near = ZarrLUTBackend(lut_path, interpolation_method="nearest")
        coeffs_near = b_near.compute_coefficients(geom, atmo, band, compute_jacobian=False)
        assert coeffs_near.xbp.shape == (3, 3)
        assert np.isfinite(coeffs_near.xbp.values).all()

    def test_jacobian_numerical_path(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_jac.zarr")
        geom = _geometry((2, 2))
        atmo = _atmo((2, 2))
        band = SensorBand("B03", 560.0, 35.0, 10.0, 1)

        b = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs = b.compute_coefficients(geom, atmo, band, compute_jacobian=True)
        assert coeffs.d_xap is not None
        assert coeffs.d_xbp is not None
        assert coeffs.d_xcp is not None
        assert coeffs.d_xap.sizes["param"] == 2

    def test_compute_coefficients_rejects_mismatched_geometry_and_atmo_grids(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_bad_shape.zarr")
        geom = _geometry((2, 2))
        atmo = _atmo((3, 3))
        band = SensorBand("B02", 490.0, 65.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")

        with pytest.raises(ValueError, match="must share the same grid shape"):
            backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

    def test_load_local_zipped_zarr(self, tmp_path: Path):
        lut_dir = _write_small_lut(tmp_path / "lut_zip.zarr")
        zip_path = tmp_path / "lut_zip.zarr.zip"
        shutil.make_archive(str(zip_path)[:-4], "zip", root_dir=lut_dir)

        b = ZarrLUTBackend(zip_path, interpolation_method="nearest")
        assert "path_reflectance" in b.lut

    def test_s3_storage_options_are_normalized(self, tmp_path: Path, monkeypatch):
        import fsspec

        lut_path = _write_small_lut(tmp_path / "lut_mapper.zarr")
        captured: dict[str, object] = {}
        original_get_mapper = fsspec.get_mapper

        def _fake_get_mapper(path: str, **kwargs):
            captured["path"] = path
            captured["kwargs"] = kwargs
            return original_get_mapper(str(lut_path))

        monkeypatch.setattr(fsspec, "get_mapper", _fake_get_mapper)

        b = ZarrLUTBackend(
            "s3://example/lut.zarr",
            storage_options={
                "region": "eu-west-2",
                "endpoint_url": "https://s3.eu-west-2.amazonaws.com",
                "anon": True,
            },
        )
        _ = b.lut

        assert captured["path"] == "s3://example/lut.zarr"
        kwargs = captured["kwargs"]
        assert isinstance(kwargs, dict)
        assert kwargs.get("anon") is True
        assert kwargs.get("client_kwargs") == {
            "region_name": "eu-west-2",
            "endpoint_url": "https://s3.eu-west-2.amazonaws.com",
        }

    def test_load_non_consolidated_zarr(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_no_consolidated.zarr", consolidated=False)
        b = ZarrLUTBackend(lut_path, interpolation_method="nearest")
        assert "path_reflectance" in b.lut

    def test_coefficient_lut_with_altitude_axis_skips_legacy_elevation_correction(self, tmp_path: Path, monkeypatch):
        lut_path = tmp_path / "lut_with_altitude.zarr"
        ds = xr.Dataset(
            {
                "path_reflectance": (
                    ["sza", "vza", "raa", "aot", "tcwv", "altitude", "wavelength"],
                    np.full((1, 1, 1, 1, 1, 1, 1), 0.02, dtype=np.float32),
                ),
                "transmittance_down": (
                    ["sza", "vza", "raa", "aot", "tcwv", "altitude", "wavelength"],
                    np.full((1, 1, 1, 1, 1, 1, 1), 0.85, dtype=np.float32),
                ),
                "transmittance_up": (
                    ["sza", "vza", "raa", "aot", "tcwv", "altitude", "wavelength"],
                    np.full((1, 1, 1, 1, 1, 1, 1), 0.90, dtype=np.float32),
                ),
                "spherical_albedo": (
                    ["sza", "vza", "raa", "aot", "tcwv", "altitude", "wavelength"],
                    np.full((1, 1, 1, 1, 1, 1, 1), 0.03, dtype=np.float32),
                ),
            },
            coords={
                "sza": np.array([30.0], dtype=np.float32),
                "vza": np.array([10.0], dtype=np.float32),
                "raa": np.array([90.0], dtype=np.float32),
                "aot": np.array([0.2], dtype=np.float32),
                "tcwv": np.array([2.0], dtype=np.float32),
                "altitude": np.array([0.1], dtype=np.float32),
                "wavelength": np.array([490.0], dtype=np.float32),
            },
        )
        ds.to_zarr(lut_path, mode="w", consolidated=True)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")

        def _should_not_run(*args, **kwargs):  # noqa: ANN002, ANN003
            raise AssertionError("legacy elevation correction should be skipped")

        monkeypatch.setattr(backend, "_apply_elevation_correction", _should_not_run)

        coeffs = backend.compute_coefficients(
            _geometry((1, 1)),
            _atmo((1, 1)),
            SensorBand("B02", 490.0, 20.0, 10.0, 0),
            compute_jacobian=False,
        )

        trans_total = 0.85 * 0.90
        np.testing.assert_allclose(coeffs.xap.values, 1.0 / trans_total, rtol=1e-6)
        np.testing.assert_allclose(coeffs.xbp.values, 0.02 / trans_total, rtol=1e-6)
        np.testing.assert_allclose(coeffs.xcp.values, 0.03 / 0.90, rtol=1e-6)

    def test_spectral_lut_derives_standard_coefficients(self, tmp_path: Path):
        lut_path, expected = _write_small_spectral_lut(tmp_path / "lut_spectral.zarr")
        geom = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), np.deg2rad(30.0), dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full((2, 2), np.deg2rad(10.0), dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((2, 2), np.deg2rad(90.0), dtype=np.float32), dims=["y", "x"]),
        )
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs = backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

        np.testing.assert_allclose(coeffs.xap.values, expected["xap"], rtol=1e-5)
        np.testing.assert_allclose(coeffs.xbp.values, expected["xbp"], rtol=1e-5)
        np.testing.assert_allclose(coeffs.xcp.values, expected["xcp"], rtol=1e-5)

        toa = xr.DataArray(np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"])
        expected_boa = ((0.15 - expected["path_ref"]) / expected["t_total"])
        expected_boa = expected_boa / (1.0 + expected["xcp"] * expected_boa)
        np.testing.assert_allclose(coeffs.apply_correction(toa).values, expected_boa, rtol=1e-5)

    def test_spectral_lut_without_solar_irradiance_still_returns_finite_coefficients(self, tmp_path: Path):
        lut_path, expected = _write_small_spectral_lut(
            tmp_path / "lut_spectral_no_solar.zarr",
            include_solar=False,
        )
        geom = _geometry((1, 1))
        atmo = _atmo((1, 1))
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="nearest")
        coeffs = backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

        assert np.isfinite(coeffs.xap.values).all()
        np.testing.assert_allclose(coeffs.xap.values, expected["xap"], rtol=1e-5)

    def test_spectral_lut_compute_coefficients_allows_nonfinite_optional_scene_axes(self, tmp_path: Path):
        lut_path, expected = _write_small_spectral_lut(tmp_path / "lut_spectral_nonfinite_optional.zarr")
        geom = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), np.deg2rad(30.0), dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full((2, 2), np.deg2rad(10.0), dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((2, 2), np.deg2rad(90.0), dtype=np.float32), dims=["y", "x"]),
        )
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), np.inf, dtype=np.float32), dims=["y", "x"]),
        )
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs = backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

        assert np.isfinite(coeffs.xap.values).all()
        assert np.isfinite(coeffs.xbp.values).all()
        assert np.isfinite(coeffs.xcp.values).all()
        np.testing.assert_allclose(coeffs.xap.values, expected["xap"], rtol=1e-5)
        np.testing.assert_allclose(coeffs.xbp.values, expected["xbp"], rtol=1e-5)
        np.testing.assert_allclose(coeffs.xcp.values, expected["xcp"], rtol=1e-5)

    def test_spectral_lut_compute_coefficients_rejects_nonfinite_required_scene_inputs(self, tmp_path: Path):
        lut_path, _ = _write_small_spectral_lut(tmp_path / "lut_spectral_nonfinite_required.zarr")
        geom = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
        )
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")

        with pytest.raises(ValueError, match="sza must contain only finite values"):
            backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

    def test_spectral_scene_preload_reuses_subset_and_band_grids(self, tmp_path: Path, monkeypatch):
        lut_path, _ = _write_small_spectral_lut(tmp_path / "lut_spectral_preload.zarr")
        geom = _geometry((2, 2))
        atmo = _atmo((2, 2))
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")
        calls = {"scene_subset": 0, "band_subset": 0}

        orig_scene_subset = backend._subset_spectral_lut_for_scene
        orig_band_subset = backend._subset_wavelength_for_band

        def _count_scene_subset(*args, **kwargs):  # noqa: ANN002, ANN003
            calls["scene_subset"] += 1
            return orig_scene_subset(*args, **kwargs)

        def _count_band_subset(*args, **kwargs):  # noqa: ANN002, ANN003
            calls["band_subset"] += 1
            return orig_band_subset(*args, **kwargs)

        monkeypatch.setattr(backend, "_subset_spectral_lut_for_scene", _count_scene_subset)
        monkeypatch.setattr(backend, "_subset_wavelength_for_band", _count_band_subset)

        backend.preload_scene_subset(geom, atmo, [band])
        backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)
        atmo_2 = atmo.with_updated_aot_tcwv(
            aot=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"]),
        )
        backend.compute_coefficients(geom, atmo_2, band, compute_jacobian=False)

        assert calls["scene_subset"] == 1
        assert calls["band_subset"] == 1

    def test_spectral_scene_preload_handles_nonfinite_optional_scene_axes(self, tmp_path: Path, monkeypatch):
        lut_path, expected = _write_small_spectral_lut(tmp_path / "lut_spectral_preload_nonfinite_optional.zarr")
        geom = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), np.deg2rad(30.0), dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full((2, 2), np.deg2rad(10.0), dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((2, 2), np.deg2rad(90.0), dtype=np.float32), dims=["y", "x"]),
        )
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), -np.inf, dtype=np.float32), dims=["y", "x"]),
        )
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")
        calls = {"scene_subset": 0, "band_subset": 0}

        orig_scene_subset = backend._subset_spectral_lut_for_scene
        orig_band_subset = backend._subset_wavelength_for_band

        def _count_scene_subset(*args, **kwargs):  # noqa: ANN002, ANN003
            calls["scene_subset"] += 1
            return orig_scene_subset(*args, **kwargs)

        def _count_band_subset(*args, **kwargs):  # noqa: ANN002, ANN003
            calls["band_subset"] += 1
            return orig_band_subset(*args, **kwargs)

        monkeypatch.setattr(backend, "_subset_spectral_lut_for_scene", _count_scene_subset)
        monkeypatch.setattr(backend, "_subset_wavelength_for_band", _count_band_subset)

        backend.preload_scene_subset(geom, atmo, [band])
        coeffs = backend.compute_coefficients(geom, atmo, band, compute_jacobian=False)

        assert calls["scene_subset"] == 1
        assert calls["band_subset"] == 1
        assert np.isfinite(coeffs.xap.values).all()
        np.testing.assert_allclose(coeffs.xap.values, expected["xap"], rtol=1e-5)

    def test_spectral_scene_preload_rejects_nonfinite_required_scene_inputs(self, tmp_path: Path):
        lut_path, _ = _write_small_spectral_lut(tmp_path / "lut_spectral_preload_nonfinite_required.zarr")
        geom = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((2, 2), np.nan, dtype=np.float32), dims=["y", "x"]),
        )
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        band = SensorBand("B02", 490.0, 20.0, 10.0, 0)

        backend = ZarrLUTBackend(lut_path, interpolation_method="linear")

        with pytest.raises(ValueError, match="sza must contain only finite values"):
            backend.preload_scene_subset(geom, atmo, [band])

    def test_unsupported_lut_representation_raises(self):
        backend = ZarrLUTBackend("unused.zarr")
        backend._lut = xr.Dataset(coords={"wavelength": [490.0]})
        backend._lut_coords = {"wavelength": np.array([490.0], dtype=np.float32)}

        with pytest.raises(ValueError, match="supported RT representation"):
            backend.compute_coefficients(
                _geometry((1, 1)),
                _atmo((1, 1)),
                SensorBand("B02", 490.0, 20.0, 10.0, 0),
                compute_jacobian=False,
            )

    def test_interpolate_lut_initializes_coord_cache_on_demand(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_coord_cache.zarr")
        backend = ZarrLUTBackend(lut_path, interpolation_method="nearest")

        path_ref, trans_down, trans_up, sph_alb = backend._interpolate_lut(
            np.array([[30.0]], dtype=np.float32),
            np.array([[10.0]], dtype=np.float32),
            np.array([[90.0]], dtype=np.float32),
            np.array([[0.2]], dtype=np.float32),
            np.array([[2.0]], dtype=np.float32),
            np.array([[0.3]], dtype=np.float32),
            np.array([[0.1]], dtype=np.float32),
            490.0,
        )

        assert "wavelength" in backend._lut_coords
        assert path_ref.shape == (1, 1)
        assert trans_down.shape == (1, 1)
        assert trans_up.shape == (1, 1)
        assert sph_alb.shape == (1, 1)

    def test_interpolate_lut_requires_non_empty_wavelength_axis(self):
        backend = ZarrLUTBackend("unused.zarr")
        backend._lut = xr.Dataset(coords={"sza": [30.0]})
        backend._lut_coords = {"sza": np.array([30.0], dtype=np.float32)}

        with pytest.raises(ValueError, match="wavelength coordinate"):
            backend._interpolate_lut(
                np.array([[30.0]], dtype=np.float32),
                np.array([[10.0]], dtype=np.float32),
                np.array([[90.0]], dtype=np.float32),
                np.array([[0.2]], dtype=np.float32),
                np.array([[2.0]], dtype=np.float32),
                np.array([[0.3]], dtype=np.float32),
                np.array([[0.1]], dtype=np.float32),
                490.0,
            )

    def test_build_point_coords_handles_unscaled_ozone_and_sparse_axes(self):
        backend = ZarrLUTBackend("unused.zarr")
        backend._lut_coords = {
            "sza": np.array([0.0, 60.0], dtype=np.float32),
            "vza": np.array([], dtype=np.float32),
            "ozone": np.array([0.2, 0.4], dtype=np.float32),
        }

        coords = backend._build_point_coords(
            sza=np.array([100.0], dtype=np.float32),
            vza=np.array([15.0], dtype=np.float32),
            raa=np.array([20.0], dtype=np.float32),
            aot=np.array([0.2], dtype=np.float32),
            tcwv=np.array([2.0], dtype=np.float32),
            tco3=np.array([0.3], dtype=np.float32),
            elevation=np.array([0.1], dtype=np.float32),
        )

        assert float(coords["sza"].values[0]) == 60.0
        assert float(coords["vza"].values[0]) == 15.0
        assert float(coords["ozone"].values[0]) == pytest.approx(0.3)

    def test_interpolate_variable_returns_input_without_matching_coords(self):
        var = xr.DataArray(np.array([1.0, 2.0], dtype=np.float32), dims=["band"])
        coords = {"sza": xr.DataArray(np.array([30.0], dtype=np.float32), dims=["point"])}

        result = ZarrLUTBackend._interpolate_variable(var, coords, "linear")

        assert result is var

    def test_spectral_backend_requires_wavelength_coordinate(self):
        backend = ZarrLUTBackend("unused.zarr")
        backend._lut = xr.Dataset(
            {
                "TOA_rho1": (["point"], np.array([0.05], dtype=np.float32)),
                "TOA_rho2": (["point"], np.array([0.25], dtype=np.float32)),
                "Eg_rho1": (["point"], np.array([0.2], dtype=np.float32)),
                "Eg_rho2": (["point"], np.array([0.5], dtype=np.float32)),
            }
        )
        backend._lut_coords = {}

        with pytest.raises(ValueError, match="wavelength coordinate"):
            backend._compute_coefficients_from_spectral_lut(
                np.array([[30.0]], dtype=np.float32),
                np.array([[10.0]], dtype=np.float32),
                np.array([[90.0]], dtype=np.float32),
                np.array([[0.2]], dtype=np.float32),
                np.array([[2.0]], dtype=np.float32),
                np.array([[0.3]], dtype=np.float32),
                np.array([[0.1]], dtype=np.float32),
                SensorBand("B02", 490.0, 20.0, 10.0, 0),
            )

    def test_spectral_weight_helpers_cover_fallback_and_single_wavelength(self):
        backend = ZarrLUTBackend("unused.zarr")
        backend._lut = xr.Dataset(
            {
                "extraterrestrial_solar_irradiance": (
                    ["sample"],
                    np.array([10.0, 12.0], dtype=np.float32),
                ),
                "solar_irradiance": (
                    ["wavelength", "sample"],
                    np.array(
                        [
                            [1.0, 3.0],
                            [2.0, 4.0],
                            [3.0, 5.0],
                        ],
                        dtype=np.float32,
                    ),
                ),
            },
            coords={
                "wavelength": np.array([480.0, 490.0, 500.0], dtype=np.float32),
                "sample": np.array([0, 1], dtype=np.int32),
            },
        )
        backend._lut_coords = {
            "wavelength": np.array([480.0, 490.0, 500.0], dtype=np.float32),
        }

        weights = backend._spectral_integration_weights(
            SensorBand("B02", 490.0, 20.0, 10.0, 0)
        )
        assert weights.dims == ("wavelength",)
        assert float(weights.max()) > 0.0

        passthrough = xr.DataArray(np.array([7.0], dtype=np.float32), dims=["point"])
        assert ZarrLUTBackend._weighted_spectral_mean(passthrough, weights) is passthrough

        single = xr.DataArray(
            np.array([0.25], dtype=np.float32),
            dims=["wavelength"],
            coords={"wavelength": np.array([490.0], dtype=np.float32)},
        )
        reduced = ZarrLUTBackend._weighted_spectral_mean(
            single,
            weights.sel(wavelength=[490.0]),
        )
        assert float(reduced.values) == pytest.approx(0.25)

    def test_http_zip_store_uses_reference_mapper_for_remote_zip(self, monkeypatch):
        import siac.algorithms.rt.lut.store as lut_store

        class _FakeZipFS:
            _files = {
                "": {"children": ["lut.zarr"]},
                "lut.zarr": {"children": ["lut.zarr/.zgroup"]},
                "lut.zarr/.zgroup": {"offset": 1, "size": 2},
            }

        class _FakeMapper(dict):
            def __init__(self):
                super().__init__()
                self.fs = _FakeZipFS()
                self.root = "lut.zarr"

        sentinel_store = {"kind": "reference"}
        monkeypatch.setattr(lut_store, "build_readonly_zip_mapper", lambda _path, _opts: _FakeMapper())
        monkeypatch.setattr("fsspec.get_mapper", lambda _path, **_kwargs: sentinel_store)

        store = lut_store.build_lut_store("https://example.com/lut.zarr.zip", storage_options={})
        assert store == sentinel_store

    def test_http_zip_headers_validation(self):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zarr.zip",
                {"headers": "not-a-dict"},
            )


class _FakeResponse:
    def __init__(self, *, status_code: int, content: bytes, headers: dict[str, str]):
        self.status_code = status_code
        self.content = content
        self.headers = headers
        self.ok = 200 <= status_code < 300

    def raise_for_status(self):
        if not self.ok:
            raise RuntimeError(f"HTTP {self.status_code}")


class _FakeSession:
    def __init__(self, data: bytes, *, support_range: bool = True):
        self._data = data
        self._support_range = support_range
        self.get_calls = 0

    def head(self, url, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _FakeResponse(
            status_code=200,
            content=b"",
            headers={"Content-Length": str(len(self._data))},
        )

    def get(self, url, headers=None, timeout=None):  # noqa: ARG002
        self.get_calls += 1
        headers = headers or {}
        range_header = headers.get("Range")
        if self._support_range and range_header is not None:
            bounds = range_header.removeprefix("bytes=").split("-", 1)
            if bounds[0] == "":
                suffix = int(bounds[1])
                start = max(0, len(self._data) - suffix)
                end = len(self._data) - 1
            else:
                start = int(bounds[0])
                end = int(bounds[1]) if bounds[1] else len(self._data) - 1
            payload = self._data[start:end + 1]
            return _FakeResponse(
                status_code=206,
                content=payload,
                headers={"Content-Range": f"bytes {start}-{end}/{len(self._data)}"},
            )

        return _FakeResponse(status_code=200, content=self._data, headers={})

    def close(self):
        pass


class _BytesRangeFS:
    """In-memory async byte-range source used for ZIP filesystem tests."""

    def __init__(self, payload: bytes):
        self.payload = payload

    async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ARG002
        size = len(self.payload)
        if start is None:
            start_i = 0
        elif start < 0:
            start_i = max(0, size + start)
        else:
            start_i = start

        if end is None:
            end_i = size
        elif end < 0:
            end_i = max(0, size + end)
        else:
            end_i = end

        end_i = min(end_i, size)
        if start_i >= size or end_i <= start_i:
            return b""
        return self.payload[start_i:end_i]


def _build_eocd(
    *,
    disk: int = 0,
    cd_disk: int = 0,
    entries_disk: int = 0,
    entries_total: int = 0,
    cd_size: int = 0,
    cd_offset: int = 0,
    comment_len: int = 0,
) -> bytes:
    return struct.pack(
        "<IHHHHIIH",
        0x06054B50,
        disk,
        cd_disk,
        entries_disk,
        entries_total,
        cd_size,
        cd_offset,
        comment_len,
    )


def _build_cd_header(
    *,
    signature: int = 0x02014B50,
    flag: int = 0,
    compression: int = 0,
    comp_size: int = 1,
    uncomp_size: int = 1,
    fname_len: int = 0,
    extra_len: int = 0,
    comment_len: int = 0,
    disk_start: int = 0,
    rel_offset: int = 100,
) -> bytes:
    return struct.pack(
        "<IHHHHHHIIIHHHHHII",
        signature,
        20,
        20,
        flag,
        compression,
        0,
        0,
        0,
        comp_size,
        uncomp_size,
        fname_len,
        extra_len,
        comment_len,
        disk_start,
        0,
        0,
        rel_offset,
    )


class TestHTTPRangeHelpers:
    def test_http_range_filesystem_with_range_support(self, monkeypatch):
        import requests

        data = b"0123456789abcdef"
        fake = _FakeSession(data, support_range=True)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=10.0)
        out = asyncio.run(fs._cat_file("https://example.com/test.zip", start=2, end=7))
        assert out == b"23456"
        fs.close()
        assert fake.get_calls >= 1

    def test_http_range_filesystem_server_ignores_range(self, monkeypatch):
        import requests

        data = b"abcdefghijklmnop"
        fake = _FakeSession(data, support_range=False)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=10.0)
        out = asyncio.run(fs._cat_file("https://example.com/test.zip", start=1, end=11))
        assert out == data[1:11]
        fs.close()
        assert fake.get_calls == 1

    def test_read_only_zip_filesystem_reads_zarr_entries(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr(".zattrs", "{}")
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
            zf.writestr("arr/0", b"\x01")

        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")
        async def _exercise():
            listing = await zip_fs._ls("", detail=False)
            payload = await zip_fs._cat_file(".zgroup")
            return listing, payload

        ls_root, payload = asyncio.run(_exercise())
        assert "/.zgroup" in ls_root
        assert "/arr" in ls_root
        assert payload.startswith(b"{")

    def test_zipfs_initialize_error_branches_zip64_and_central_directory(self):
        # ZIP64 sentinel values + too-short tail for ZIP64 EOCD locator.
        tail_short = b"A" * 10 + _build_eocd(
            disk=0xFFFF,
            cd_disk=0xFFFF,
            entries_disk=0xFFFF,
            entries_total=0xFFFF,
            cd_size=0xFFFFFFFF,
            cd_offset=0xFFFFFFFF,
            comment_len=0,
        )
        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_short), "dummy.zip")
        with pytest.raises(ValueError, match="ZIP64 EOCD and locator do not fit"):
            asyncio.run(zip_fs._initialize())

        # ZIP64 markers present but ZIP64 EOCD signature missing in the searched window.
        tail_long = b"A" * 120 + _build_eocd(
            disk=0xFFFF,
            cd_disk=0xFFFF,
            entries_disk=0xFFFF,
            entries_total=0xFFFF,
            cd_size=0xFFFFFFFF,
            cd_offset=0xFFFFFFFF,
            comment_len=0,
        )
        zip_fs2 = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_long), "dummy.zip")
        with pytest.raises(ValueError, match="No ZIP64 EOCD found"):
            asyncio.run(zip_fs2._initialize())

        # Multi-disk central directory is rejected.
        tail_multidisk = b"A" * 10 + _build_eocd(disk=1, cd_disk=0, entries_disk=1, entries_total=1)
        zip_fs3 = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_multidisk), "dummy.zip")
        with pytest.raises(ValueError, match="Unsupported multi-disk"):
            asyncio.run(zip_fs3._initialize())

    def test_zipfs_initialize_error_branches_cd_entries(self):
        class _TwoPhaseFS:
            def __init__(self, tail: bytes, cd_payload: bytes):
                self.tail = tail
                self.cd_payload = cd_payload

            async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ANN001, ARG002
                if start is not None and start < 0:
                    return self.tail
                return self.cd_payload

        # cd_size > returned bytes
        tail = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=50, cd_offset=0)
        fs_short = _TwoPhaseFS(tail, b"x" * 10)
        zip_fs = _ReadOnlyZipFileSystem(fs_short, "dummy.zip")
        with pytest.raises(ValueError, match="Failed to read central directory"):
            asyncio.run(zip_fs._initialize())

        # Truncated entry (< 46 bytes)
        tail2 = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=20, cd_offset=0)
        fs_trunc_entry = _TwoPhaseFS(tail2, b"x" * 20)
        zip_fs2 = _ReadOnlyZipFileSystem(fs_trunc_entry, "dummy.zip")
        with pytest.raises(ValueError, match="Truncated central directory entry"):
            asyncio.run(zip_fs2._initialize())

        # Invalid central directory signature.
        bad_sig = _build_cd_header(signature=0x12345678)
        tail3 = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=len(bad_sig), cd_offset=0)
        fs_bad_sig = _TwoPhaseFS(tail3, bad_sig)
        zip_fs3 = _ReadOnlyZipFileSystem(fs_bad_sig, "dummy.zip")
        with pytest.raises(ValueError, match="Invalid central directory signature"):
            asyncio.run(zip_fs3._initialize())

    def test_zipfs_initialize_error_branches_filename_and_extra(self):
        class _TwoPhaseFS:
            def __init__(self, tail: bytes, cd_payload: bytes):
                self.tail = tail
                self.cd_payload = cd_payload

            async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ANN001, ARG002
                if start is not None and start < 0:
                    return self.tail
                return self.cd_payload

        def _run_with_cd(cd_payload: bytes, msg: str) -> None:
            tail = b"A" * 20 + _build_eocd(
                entries_disk=1,
                entries_total=1,
                cd_size=len(cd_payload),
                cd_offset=0,
            )
            zip_fs = _ReadOnlyZipFileSystem(_TwoPhaseFS(tail, cd_payload), "dummy.zip")
            with pytest.raises(ValueError, match=msg):
                asyncio.run(zip_fs._initialize())

        # filename bytes missing
        _run_with_cd(_build_cd_header(fname_len=5), "Truncated filename")

        # extra field body missing
        cd_extra_short = _build_cd_header(fname_len=1, extra_len=5) + b"a"
        _run_with_cd(cd_extra_short, "Truncated extra field")

        # extra field header missing (need at least 4 bytes)
        cd_extra_header = _build_cd_header(fname_len=1, extra_len=3) + b"a" + b"\x01\x02\x03"
        _run_with_cd(cd_extra_header, "Truncated extra field header")

        # extra field payload missing after valid tag/size
        cd_extra_payload = _build_cd_header(fname_len=1, extra_len=5) + b"a" + struct.pack("<HH", 1, 5) + b"x"
        _run_with_cd(cd_extra_payload, "Truncated extra field payload")

    def test_zipfs_parent_directory_construction_branch(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            # No explicit directory entries; parent dirs are inferred from file names.
            zf.writestr("nested/dir/a.bin", b"a")
            zf.writestr("nested/dir/b.bin", b"b")

        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")
        asyncio.run(zip_fs._initialize())
        assert zip_fs._files is not None
        assert "nested" in zip_fs._files
        assert "nested/dir" in zip_fs._files
        assert "nested/dir/a.bin" in zip_fs._files


class _FallbackSession:
    """Session stub that forces HTTP range fallback branches."""

    def __init__(self, payload: bytes):
        self.payload = payload
        self.calls: list[tuple[str, dict[str, str]]] = []

    def head(self, url, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _FakeResponse(status_code=403, content=b"", headers={})

    def get(self, url, headers=None, timeout=None):  # noqa: ARG002
        hdrs = dict(headers or {})
        self.calls.append((url, hdrs))
        if "Range" in hdrs:
            return _FakeResponse(status_code=403, content=b"", headers={})
        return _FakeResponse(status_code=200, content=self.payload, headers={})

    def close(self):
        pass


class TestZipStoreUtilities:
    def test_helper_path_and_slice_utilities(self):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        assert zip_store._normalize_fs_path("") == ""
        assert zip_store._normalize_fs_path("/") == ""
        assert zip_store._normalize_fs_path("a/b/") == "a/b"

        assert zip_store._slice_bounds(10, None, None) == (0, 10)
        assert zip_store._slice_bounds(10, -3, None) == (7, 10)
        assert zip_store._slice_bounds(10, 3, -2) == (3, 8)
        assert zip_store._slice_bounds(10, 20, None) == (10, 10)

    def test_http_range_filesystem_full_body_fallback(self, monkeypatch):
        import requests

        payload = b"abcdefghijklmnopqrstuvwxyz"
        fake = _FallbackSession(payload)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=5.0)
        out = asyncio.run(fs._cat_file("https://example.com/fallback.zip", start=2, end=9))
        assert out == payload[2:9]
        info = asyncio.run(fs._info("https://example.com/fallback.zip"))
        assert info["size"] == len(payload)
        listing = asyncio.run(fs._ls("https://example.com/fallback.zip", detail=False))
        assert listing == ["https://example.com/fallback.zip"]
        fs.close()

    def test_local_filesystem_wrapper(self, tmp_path: Path):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        file_path = tmp_path / "x.bin"
        file_path.write_bytes(b"0123456789")

        local_fs = zip_store._build_local_filesystem()
        info = asyncio.run(local_fs._info(str(file_path)))
        assert info["type"] == "file"
        assert info["size"] == 10

        listing = asyncio.run(local_fs._ls(str(tmp_path), detail=False))
        assert str(file_path) in listing

        out = asyncio.run(local_fs._cat_file(str(file_path), start=-4, end=None))
        assert out == b"6789"

        with pytest.raises(FileNotFoundError):
            asyncio.run(local_fs._info(str(tmp_path / "missing.bin")))

    def test_read_only_zip_filesystem_error_paths(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")

        with pytest.raises(FileNotFoundError):
            asyncio.run(zip_fs._info("missing"))
        with pytest.raises(FileNotFoundError):
            asyncio.run(zip_fs._cat_file("arr"))

    def test_s3_parse_and_import_error(self):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        assert zip_store._parse_s3_url("s3://bucket/key/path") == ("bucket", "key/path")
        with pytest.raises(ValueError):
            zip_store._parse_s3_url("http://bucket/key")
        with pytest.raises(ValueError):
            zip_store._parse_s3_url("s3://bucket-only")

        try:
            fs = zip_store._build_s3_filesystem({})
        except ImportError:
            # Expected in lightweight test envs where s3fs is not installed.
            return
        assert fs is not None

    def test_s3_builder_defaults_to_ambient_credentials(self, monkeypatch):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        captured: dict[str, object] = {}

        class _FakeS3FileSystem:
            def __init__(self, **kwargs):
                captured.update(kwargs)

        monkeypatch.setitem(sys.modules, "s3fs", SimpleNamespace(S3FileSystem=_FakeS3FileSystem))
        _ = zip_store._build_s3_filesystem({})

        assert captured.get("asynchronous") is True
        assert "anon" not in captured
        assert "key" not in captured
        assert "secret" not in captured

    def test_build_mapper_validation_and_remote_compressed_branch(self, monkeypatch):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zip",
                {"headers": "bad"},
            )
        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "/tmp/lut.zip",
                {"unexpected": 1},
            )

        monkeypatch.setattr(
            zip_store,
            "_detect_zarr_prefix",
            lambda _fs: (_ for _ in ()).throw(ValueError("not stored (uncompressed)")),
        )
        with pytest.raises(ValueError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zip",
                {"timeout": 1.0},
            )

    def test_build_mapper_http_omits_none_headers_for_fs_serialization(self, monkeypatch):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        class _FakeZipFS:
            def __init__(self, fs, path, **kwargs):  # noqa: ANN001, ARG002
                self.fs = fs
                self.path = path

        monkeypatch.setattr(zip_store, "_ReadOnlyZipFileSystem", _FakeZipFS)
        monkeypatch.setattr(zip_store, "_detect_zarr_prefix", lambda _fs: "")
        def _fake_fsmap(root, fs, check=False, create=False):  # noqa: ANN001
            _ = (check, create)
            return SimpleNamespace(root=root, fs=fs)

        monkeypatch.setattr(zip_store, "FSMap", _fake_fsmap)

        mapper = zip_store.build_readonly_zip_mapper(
            "https://example.com/lut.zip",
            {"timeout": 1.0},
        )

        assert mapper.root == ""
        assert mapper.fs.fs.to_dict(include_password=False)["timeout"] == 1.0

    def test_local_compressed_zip_fallback_detects_nested_zarr_root(self, tmp_path: Path):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        zip_path = tmp_path / "lut_nested.zarr.zip"
        with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("lut.zarr/.zgroup", '{"zarr_format":2}')
            zf.writestr("lut.zarr/.zattrs", "{}")

        mapper = zip_store.build_readonly_zip_mapper(str(zip_path), {})
        assert ".zgroup" in mapper
        assert mapper[".zgroup"].startswith(b"{")

    def test_read_only_zip_listing_detail_and_close(self):
        import siac.algorithms.rt.lut.http_zip_store as zip_store

        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
            zf.writestr("arr/0", b"\x01")

        class _ClosableBytesFS(_BytesRangeFS):
            def __init__(self, payload: bytes):
                super().__init__(payload)
                self.closed = False

            def close(self):
                self.closed = True

        fs = _ClosableBytesFS(buf.getvalue())
        zip_fs = zip_store._ReadOnlyZipFileSystem(fs, "dummy.zip")

        async def _exercise():
            root = await zip_fs._ls("", detail=True)
            one = await zip_fs._ls(".zgroup", detail=True)
            return root, one

        root_listing, file_listing = asyncio.run(_exercise())
        assert any(item["type"] == "directory" for item in root_listing)
        assert file_listing[0]["type"] == "file"
        zip_fs.close()
        assert fs.closed is True

    def test_store_helper_paths_and_local_mapper_branch(self, monkeypatch, tmp_path: Path):
        import siac.algorithms.rt.lut.store as lut_store

        local_file_uri = (tmp_path / "lut.zarr").as_uri()
        local_path = lut_store.as_local_path(local_file_uri)
        assert local_path is not None
        assert str(local_path).endswith("lut.zarr")

        opts = lut_store.normalize_storage_options(
            "s3://bucket/key",
            {
                "region": "eu-west-1",
                "endpoint_url": "https://example.invalid",
                "client_kwargs": {"region_name": "keep-me"},
            },
        )
        assert opts["client_kwargs"]["region_name"] == "keep-me"
        assert opts["client_kwargs"]["endpoint_url"] == "https://example.invalid"

        captured: dict[str, object] = {}
        monkeypatch.setattr(
            "fsspec.get_mapper",
            lambda path, **kwargs: captured.update({"path": path, "kwargs": kwargs}) or {"ok": True},
        )
        out = lut_store.build_lut_store(str(tmp_path / "x.zarr"), {"anon": True})
        assert out == {"ok": True}
        assert captured["path"] == str(tmp_path / "x.zarr")
        assert captured["kwargs"] == {"anon": True}

    def test_remote_zip_builds_and_reuses_reference_json(self, monkeypatch, tmp_path: Path):
        import siac.algorithms.rt.lut.store as lut_store

        class _FakeZipFS:
            _files = {
                "": {"children": ["lut.zarr"]},
                "lut.zarr": {"children": ["lut.zarr/.zgroup", "lut.zarr/arr"]},
                "lut.zarr/arr": {"children": ["lut.zarr/arr/0"]},
                "lut.zarr/.zgroup": {"offset": 11, "size": 17},
                "lut.zarr/arr/0": {"offset": 42, "size": 3},
            }

        class _FakeMapper(dict):
            def __init__(self):
                super().__init__()
                self.fs = _FakeZipFS()
                self.root = "lut.zarr"

        build_calls = {"n": 0}
        reference_calls: list[dict[str, object]] = []

        def _fake_build(path: str, options: dict[str, object]):
            build_calls["n"] += 1
            assert path == "https://example.com/lut.zarr.zip"
            assert options == {"timeout": 5.0}
            return _FakeMapper()

        def _fake_get_mapper(path: str, **kwargs):
            assert path == "reference://"
            reference_calls.append(kwargs)
            return {"kind": "reference", "kwargs": kwargs}

        monkeypatch.setattr(lut_store, "build_readonly_zip_mapper", _fake_build)
        monkeypatch.setattr("fsspec.get_mapper", _fake_get_mapper)

        options = {"timeout": 5.0, "reference_cache_dir": str(tmp_path)}
        out1 = lut_store.build_lut_store("https://example.com/lut.zarr.zip", options)
        out2 = lut_store.build_lut_store("https://example.com/lut.zarr.zip", options)

        assert out1["kind"] == "reference"
        assert out2["kind"] == "reference"
        assert build_calls["n"] == 1
        assert len(reference_calls) == 2
        assert reference_calls[0]["remote_protocol"] == "https"
        assert reference_calls[0]["remote_options"] == {"timeout": 5.0, "asynchronous": True}
        assert isinstance(reference_calls[0]["fo"], dict)

        reference_json = lut_store._reference_json_path(
            "https://example.com/lut.zarr.zip",
            lut_store._ReferenceOptions(
                refresh=False,
                reference_json=None,
                cache_dir=tmp_path,
            ),
        )
        assert reference_json.exists()
        payload = json.loads(reference_json.read_text(encoding="utf-8"))
        assert payload["version"] == 1
        assert payload["refs"][".zgroup"] == ["https://example.com/lut.zarr.zip", 11, 17]
        assert payload["refs"]["arr/0"] == ["https://example.com/lut.zarr.zip", 42, 3]

    def test_remote_zip_rebuilds_invalid_cached_reference_json(self, monkeypatch, tmp_path: Path):
        import siac.algorithms.rt.lut.store as lut_store

        class _FakeZipFS:
            _files = {
                "": {"children": ["lut.zarr"]},
                "lut.zarr": {"children": ["lut.zarr/.zgroup"]},
                "lut.zarr/.zgroup": {"offset": 11, "size": 17},
            }

        class _FakeMapper(dict):
            def __init__(self):
                super().__init__()
                self.fs = _FakeZipFS()
                self.root = "lut.zarr"

        build_calls = {"n": 0}
        reference_calls: list[dict[str, object]] = []

        def _fake_build(path: str, options: dict[str, object]):
            build_calls["n"] += 1
            assert path == "https://example.com/lut.zarr.zip"
            assert options == {"timeout": 5.0}
            return _FakeMapper()

        def _fake_get_mapper(path: str, **kwargs):
            assert path == "reference://"
            reference_calls.append(kwargs)
            return {"kind": "reference", "kwargs": kwargs}

        monkeypatch.setattr(lut_store, "build_readonly_zip_mapper", _fake_build)
        monkeypatch.setattr("fsspec.get_mapper", _fake_get_mapper)

        reference_json = lut_store._reference_json_path(
            "https://example.com/lut.zarr.zip",
            lut_store._ReferenceOptions(
                refresh=False,
                reference_json=None,
                cache_dir=tmp_path,
            ),
        )
        reference_json.parent.mkdir(parents=True, exist_ok=True)
        reference_json.write_text('{"version": 1, "refs": {}}', encoding="utf-8")

        out = lut_store.build_lut_store(
            "https://example.com/lut.zarr.zip",
            {"timeout": 5.0, "reference_cache_dir": str(tmp_path)},
        )

        assert out["kind"] == "reference"
        assert build_calls["n"] == 1
        assert len(reference_calls) == 1
        assert reference_json.exists()
        payload = json.loads(reference_json.read_text(encoding="utf-8"))
        assert payload["refs"][".zgroup"] == ["https://example.com/lut.zarr.zip", 11, 17]

    def test_remote_zip_reference_mapper_failure_raises(self, monkeypatch, tmp_path: Path):
        import siac.algorithms.rt.lut.store as lut_store

        class _FakeZipFS:
            _files = {
                "": {"children": ["lut.zarr"]},
                "lut.zarr": {"children": ["lut.zarr/.zgroup"]},
                "lut.zarr/.zgroup": {"offset": 4, "size": 8},
            }

        class _FakeMapper(dict):
            def __init__(self):
                super().__init__()
                self.fs = _FakeZipFS()
                self.root = "lut.zarr"

        sentinel = _FakeMapper()
        monkeypatch.setattr(lut_store, "build_readonly_zip_mapper", lambda _path, _options: sentinel)

        def _raise_reference(path: str, **kwargs):
            if path == "reference://":
                raise RuntimeError("reference backend unavailable")
            return {"unexpected": True}

        monkeypatch.setattr("fsspec.get_mapper", _raise_reference)

        with pytest.raises(RuntimeError, match="reference backend unavailable"):
            _ = lut_store.build_lut_store(
                "https://example.com/lut.zarr.zip",
                {"reference_cache_dir": str(tmp_path)},
            )


class TestCreateLUTFromPy6S:
    def test_import_error_branch(self, tmp_path: Path):
        saved = sys.modules.pop("Py6S", None)
        try:
            with pytest.raises(ImportError):
                create_lut_from_py6s(tmp_path / "x.zarr", wavelengths=[500.0])
        finally:
            if saved is not None:
                sys.modules["Py6S"] = saved

    def test_create_lut_with_fake_py6s(self, tmp_path: Path, monkeypatch):
        class _FakeGeometry:
            @staticmethod
            def User():
                return SimpleNamespace()

        class _FakeWavelength:
            def __init__(self, value):
                self.value = value

        class _FakeAtmosProfile:
            @staticmethod
            def UserWaterAndOzone(tcwv, ozone):
                return (tcwv, ozone)

        class _FakeAeroProfile:
            Continental = "continental"
            Maritime = "maritime"

        class _FakeSixS:
            def __init__(self):
                self.outputs = None
                self.geometry = None
                self.aot550 = 0.0
                self.wavelength = None

            def run(self):
                wl = float(getattr(self.wavelength, "value", 0.5))
                val = 0.01 + 0.1 * self.aot550 + 0.001 * wl
                self.outputs = SimpleNamespace(
                    atmospheric_intrinsic_reflectance=val,
                    transmittance_total_scattering=SimpleNamespace(
                        downward=0.8,
                        upward=0.85,
                    ),
                    spherical_albedo=0.03,
                )

        fake_module = SimpleNamespace(
            SixS=_FakeSixS,
            AtmosProfile=_FakeAtmosProfile,
            AeroProfile=_FakeAeroProfile,
            Geometry=_FakeGeometry,
            Wavelength=_FakeWavelength,
        )
        monkeypatch.setitem(sys.modules, "Py6S", fake_module)

        out = tmp_path / "fake_lut.zarr"
        create_lut_from_py6s(
            output_path=out,
            wavelengths=[500.0],
            sza_range=(0, 1, 1),
            vza_range=(0, 1, 1),
            raa_range=(0, 1, 1),
            aot_values=[0.1],
            tcwv_values=[1.0],
            aerosol_type="unknown",  # hits fallback branch
            ozone=0.3,
        )

        ds = xr.open_zarr(out, consolidated=True)
        assert "path_reflectance" in ds
        assert ds["path_reflectance"].shape == (1, 1, 1, 1, 1, 1)
        assert ds.attrs["aerosol_type"] == "unknown"

    def test_create_lut_defaults_maritime_and_progress(self, tmp_path: Path, monkeypatch):
        import siac.algorithms.rt.lut.create as create_mod

        class _FakeGeometry:
            @staticmethod
            def User():
                return SimpleNamespace()

        class _FakeWavelength:
            def __init__(self, value):
                self.value = value

        class _FakeAtmosProfile:
            @staticmethod
            def UserWaterAndOzone(tcwv, ozone):
                return (tcwv, ozone)

        class _FakeAeroProfile:
            Continental = "continental"
            Maritime = "maritime"

        class _FakeSixS:
            def __init__(self):
                self.outputs = None
                self.geometry = None
                self.aot550 = 0.0
                self.wavelength = None

            def run(self):
                val = 0.02 + 0.1 * self.aot550 + 0.001 * float(getattr(self.wavelength, "value", 0.5))
                self.outputs = SimpleNamespace(
                    atmospheric_intrinsic_reflectance=val,
                    transmittance_total_scattering=SimpleNamespace(downward=0.8, upward=0.85),
                    spherical_albedo=0.03,
                )

        fake_module = SimpleNamespace(
            SixS=_FakeSixS,
            AtmosProfile=_FakeAtmosProfile,
            AeroProfile=_FakeAeroProfile,
            Geometry=_FakeGeometry,
            Wavelength=_FakeWavelength,
        )
        monkeypatch.setitem(sys.modules, "Py6S", fake_module)
        monkeypatch.setattr(
            create_mod.np,
            "logspace",
            lambda *_args, **_kwargs: np.array([0.2], dtype=np.float32),
        )

        # Defaults path (aot_values/tcwv_values None) + maritime branch.
        create_lut_from_py6s(
            output_path=tmp_path / "defaults.zarr",
            wavelengths=[500.0],
            sza_range=(0, 1, 1),
            vza_range=(0, 1, 1),
            raa_range=(0, 1, 1),
            aot_values=None,
            tcwv_values=None,
            aerosol_type="maritime",
        )

        logged: list[str] = []
        monkeypatch.setattr(
            create_mod.logger,
            "info",
            lambda msg: logged.append(str(msg)),
        )
        create_lut_from_py6s(
            output_path=tmp_path / "progress.zarr",
            wavelengths=[500.0],
            sza_range=(0, 10, 1),
            vza_range=(0, 10, 1),
            raa_range=(0, 10, 1),
            aot_values=[0.1],
            tcwv_values=[1.0],
            aerosol_type="continental",
        )
        assert any("Progress: 1000/1000" in m for m in logged)

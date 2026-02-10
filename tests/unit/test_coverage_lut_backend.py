"""Coverage tests for the Zarr LUT backend."""

from __future__ import annotations

import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.core.types import AtmosphericState, GeometryAngles, SensorBand
from siac.rt.lut.zarr_lut import ZarrLUTBackend, create_lut_from_py6s


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


def _write_small_lut(path: Path) -> Path:
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
    ds.to_zarr(path, mode="w", consolidated=True)
    return path


class TestZarrLUTBackend:
    def test_missing_lut_raises(self, tmp_path: Path):
        b = ZarrLUTBackend(tmp_path / "missing.zarr")
        with pytest.raises(FileNotFoundError):
            _ = b.lut

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
        assert b_lin.supports_jacobian()
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

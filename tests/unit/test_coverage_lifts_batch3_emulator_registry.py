"""Coverage lifts for rt.emulator.two_nn branches."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

import siac.algorithms.rt.emulator.two_nn as tnn
from siac.domain import AtmosphericState, GeometryAngles, SensorBand


def _geom(shape=(2, 2)):
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


def _atmo(shape=(2, 2)):
    def da(v):
        return xr.DataArray(np.full(shape, v, dtype=np.float32), dims=["y", "x"])
    return AtmosphericState(aot=da(0.15), tcwv=da(2.0), tco3=da(0.3), aot_unc=da(0.05), tcwv_unc=da(0.3), tco3_unc=da(0.03), elevation=da(0.1))


def _weights(hidden=4):
    hidden_layers = [
        [np.ones((7, hidden), dtype=np.float32) * 0.1, np.zeros(hidden, dtype=np.float32)],
        [np.ones((hidden, hidden), dtype=np.float32) * 0.1, np.zeros(hidden, dtype=np.float32)],
    ]
    output_layers = [
        [np.ones((hidden, 1), dtype=np.float32) * 0.1, np.zeros(1, dtype=np.float32)],
        [np.ones((hidden, 1), dtype=np.float32) * 0.2, np.zeros(1, dtype=np.float32)],
        [np.ones((hidden, 1), dtype=np.float32) * 0.3, np.zeros(1, dtype=np.float32)],
    ]
    return hidden_layers, output_layers


def test_two_nn_empty_dir_and_missing_band_paths(tmp_path: Path):
    emu = tnn.TwoLayerNNEmulator(tmp_path, sensor_id="MSI", satellite_id="S2A")
    assert emu.available_bands == []
    assert emu.is_available_for_sensor("MSI", "S2A") is False

    assert emu._get_emulator_path("B02") is None
    with pytest.raises(FileNotFoundError, match="No emulator found"):
        emu._load_band_emulator("B02")


def test_two_nn_compute_coeffs_jacobian_requested_but_none(monkeypatch, tmp_path: Path):
    emu = tnn.TwoLayerNNEmulator(tmp_path, sensor_id="MSI", satellite_id="S2A")

    class _FakeBand:
        def forward(self, x, compute_jacobian=False):
            n = x.shape[0]
            out = np.column_stack([
                np.full(n, 0.9, dtype=np.float32),
                np.full(n, 0.02, dtype=np.float32),
                np.full(n, 0.1, dtype=np.float32),
            ])
            return out, None

    monkeypatch.setattr(emu, "_load_band_emulator", lambda _band_name: _FakeBand())
    coeffs = emu.compute_coefficients(_geom(), _atmo(), SensorBand("B02", 490.0, 65.0, 10.0, 0), compute_jacobian=True)
    assert coeffs.d_xap is None
    assert coeffs.xap.shape == (2, 2)


def test_two_nn_load_classmethod_default_emulator_dir(monkeypatch, tmp_path: Path):
    import siac

    fake_pkg = tmp_path / "pkg"
    fake_pkg.mkdir(parents=True, exist_ok=True)
    fake_init = fake_pkg / "__init__.py"
    fake_init.write_text("# pkg")

    captured = {}

    def _fake_init(self, emulator_dir, sensor_id, satellite_id):
        captured["emulator_dir"] = Path(emulator_dir)
        captured["sensor_id"] = sensor_id
        captured["satellite_id"] = satellite_id
        self.emulator_dir = Path(emulator_dir)
        self.sensor_id = sensor_id
        self.satellite_id = satellite_id
        self._available_bands = []

    monkeypatch.setattr(siac, "__file__", str(fake_init), raising=False)
    monkeypatch.setattr(tnn.TwoLayerNNEmulator, "__init__", _fake_init)

    obj = tnn.TwoLayerNNEmulator.load("S2A", "B02", emulator_dir=None)
    assert isinstance(obj, tnn.TwoLayerNNEmulator)
    assert captured["sensor_id"] == "MSI"
    assert captured["satellite_id"] == "S2A"
    assert captured["emulator_dir"].name == "emus"


def test_band_emulator_initializes_combined_rust_output_layer(monkeypatch):
    hidden_layers, output_layers = _weights(hidden=3)
    captured = {}

    class _FakeRustNN:
        def __init__(self, w1, b1, w2, b2, w3, b3):
            captured["w1"] = w1
            captured["b1"] = b1
            captured["w2"] = w2
            captured["b2"] = b2
            captured["w3"] = w3
            captured["b3"] = b3

        def predict(self, x, compute_jacobian):
            n = x.shape[0]
            out = np.zeros((n, 3), dtype=np.float32)
            jac = np.zeros((n, 3, x.shape[1]), dtype=np.float32) if compute_jacobian else None
            return out, jac

    monkeypatch.setattr(tnn, "_RustNN", _FakeRustNN, raising=False)

    be = tnn._BandEmulator(hidden_layers, output_layers)
    out, jac = be.forward(np.zeros((2, 7), dtype=np.float32), compute_jacobian=True)
    assert out.shape == (2, 3)
    assert jac is not None and jac.shape == (2, 3, 7)
    assert captured["w1"].shape == (7, 3)
    assert captured["w3"].shape == (3, 3)
    assert captured["b3"].shape == (3,)


def test_emulator_registry_branching(monkeypatch, tmp_path: Path):
    reg = tnn.EmulatorRegistry(tmp_path)

    # available_bands empty -> None path
    class _NoBands:
        def __init__(self, *args, **kwargs):
            self.available_bands = []

    monkeypatch.setattr(tnn, "TwoLayerNNEmulator", _NoBands)
    assert reg.get_emulator("MSI", "S2A") is None

    # loader raises -> None path
    class _Boom:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("boom")

    reg2 = tnn.EmulatorRegistry(tmp_path)
    monkeypatch.setattr(tnn, "TwoLayerNNEmulator", _Boom)
    assert reg2.get_emulator("MSI", "S2B") is None

    # list_supported_sensors branch with delegated is_sensor_supported
    reg3 = tnn.EmulatorRegistry(tmp_path)
    monkeypatch.setattr(reg3, "is_sensor_supported", lambda _sensor, sat: sat in {"S2A", "L8"})
    supported = reg3.list_supported_sensors()
    assert ("MSI", "S2A") in supported
    assert ("OLI", "L8") in supported

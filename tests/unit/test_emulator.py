"""
Unit tests for neural network emulator.
"""

import numpy as np
import pytest
import xarray as xr
from pathlib import Path

from siac.core.types import AtmosphericState, GeometryAngles, RTCoefficients, SensorBand


class TestBandEmulator:
    """Tests for _BandEmulator internal class."""

    @pytest.fixture
    def mock_weights(self):
        """Create mock neural network weights."""
        np.random.seed(42)
        hidden = 64
        input_dim = 7

        # Two hidden layers + output
        hidden_layers = [
            [np.random.randn(input_dim, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
            [np.random.randn(hidden, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
        ]

        # Three output heads (xap, xbp, xcp)
        output_layers = [
            [np.random.randn(hidden, 1).astype(np.float32) * 0.1,
             np.zeros(1, dtype=np.float32)]
            for _ in range(3)
        ]

        return hidden_layers, output_layers

    def test_forward_pass(self, mock_weights):
        """Forward pass should produce outputs."""
        from siac.rt.emulator.two_nn import _BandEmulator

        hidden_layers, output_layers = mock_weights
        emulator = _BandEmulator(hidden_layers, output_layers, use_rust=False)

        # Sample input
        x = np.random.randn(10, 7).astype(np.float32)

        outputs, jacobians = emulator.forward(x, compute_jacobian=False)

        assert outputs.shape == (10, 3)  # 3 outputs
        assert jacobians is None

    def test_forward_with_jacobian(self, mock_weights):
        """Forward pass with Jacobian computation."""
        from siac.rt.emulator.two_nn import _BandEmulator

        hidden_layers, output_layers = mock_weights
        emulator = _BandEmulator(hidden_layers, output_layers, use_rust=False)

        x = np.random.randn(5, 7).astype(np.float32)

        outputs, jacobians = emulator.forward(x, compute_jacobian=True)

        assert outputs.shape == (5, 3)
        assert jacobians is not None
        assert jacobians.shape == (5, 3, 7)

    def test_jacobian_numerical_check(self):
        """Jacobian should match numerical differentiation."""
        from siac.rt.emulator.two_nn import _BandEmulator

        # Create small deterministic weights with positive biases
        # to ensure neurons are active (avoid ReLU discontinuities)
        np.random.seed(123)
        input_dim, hidden = 7, 16

        hidden_layers = [
            [np.random.randn(input_dim, hidden).astype(np.float32) * 0.1,
             np.ones(hidden, dtype=np.float32) * 0.5],  # Positive bias
            [np.random.randn(hidden, hidden).astype(np.float32) * 0.1,
             np.ones(hidden, dtype=np.float32) * 0.5],  # Positive bias
        ]
        output_layers = [
            [np.random.randn(hidden, 1).astype(np.float32) * 0.1,
             np.zeros(1, dtype=np.float32)]
            for _ in range(3)
        ]

        emulator = _BandEmulator(hidden_layers, output_layers, use_rust=False)

        # Input that produces positive activations
        x = np.array([[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]], dtype=np.float32)

        # Analytical Jacobian
        _, jacobians = emulator.forward(x, compute_jacobian=True)

        # Numerical Jacobian
        eps = 1e-4
        numerical_jac = np.zeros((1, 3, 7), dtype=np.float32)

        for i in range(7):
            x_plus = x.copy()
            x_plus[0, i] += eps
            out_plus, _ = emulator.forward(x_plus, compute_jacobian=False)

            x_minus = x.copy()
            x_minus[0, i] -= eps
            out_minus, _ = emulator.forward(x_minus, compute_jacobian=False)

            numerical_jac[0, :, i] = (out_plus[0] - out_minus[0]) / (2 * eps)

        # Allow tolerance for float32 numerical differentiation
        # Central differences have O(eps^2) error plus O(eps^-1) roundoff
        np.testing.assert_allclose(jacobians, numerical_jac, rtol=0.01, atol=1e-4)


class TestTwoLayerNNEmulator:
    """Tests for TwoLayerNNEmulator class."""

    @pytest.fixture
    def mock_emulator_dir(self, tmp_path):
        """Create mock emulator directory with test weights."""
        np.random.seed(42)
        hidden = 64
        input_dim = 7

        hidden_layers = [
            [np.random.randn(input_dim, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
            [np.random.randn(hidden, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
        ]

        output_layers = [
            [np.random.randn(hidden, 1).astype(np.float32) * 0.1,
             np.zeros(1, dtype=np.float32)]
        ]

        # Save for multiple bands
        for band in ["B02", "B03", "B04"]:
            path = tmp_path / f"S2A_{band}_emulator.npz"
            np.savez(
                path,
                Hidden_Layers=np.array(hidden_layers, dtype=object),
                Output_Layers=np.array(output_layers, dtype=object),
            )

        return tmp_path

    def test_discover_bands(self, mock_emulator_dir):
        """Should discover available bands."""
        from siac.rt.emulator.two_nn import TwoLayerNNEmulator

        emulator = TwoLayerNNEmulator(
            emulator_dir=mock_emulator_dir,
            sensor_id="MSI",
            satellite_id="S2A",
            use_rust=False,
        )

        assert len(emulator.available_bands) >= 3
        assert "B02" in emulator.available_bands

    def test_backend_properties(self, mock_emulator_dir):
        """Backend properties should be correct."""
        from siac.rt.emulator.two_nn import TwoLayerNNEmulator

        emulator = TwoLayerNNEmulator(
            emulator_dir=mock_emulator_dir,
            sensor_id="MSI",
            satellite_id="S2A",
        )

        assert emulator.backend_name == "emulator"
        assert emulator.supports_jacobian() is True
        assert emulator.is_available_for_sensor("MSI", "S2A") is True


class TestComputeCoefficientsMulti:
    """Tests for batch coefficient computation."""

    @pytest.fixture
    def mock_emulator_dir(self, tmp_path):
        """Create mock emulator directory with test weights (3 output heads)."""
        np.random.seed(42)
        hidden = 64
        input_dim = 7

        hidden_layers = [
            [np.random.randn(input_dim, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
            [np.random.randn(hidden, hidden).astype(np.float32) * 0.1,
             np.zeros(hidden, dtype=np.float32)],
        ]

        # Three output heads (xap, xbp, xcp)
        output_layers = [
            [np.random.randn(hidden, 1).astype(np.float32) * 0.1,
             np.zeros(1, dtype=np.float32)]
            for _ in range(3)
        ]

        for band in ["B02", "B03", "B04"]:
            path = tmp_path / f"S2A_{band}_emulator.npz"
            np.savez(
                path,
                Hidden_Layers=np.array(hidden_layers, dtype=object),
                Output_Layers=np.array(output_layers, dtype=object),
            )

        return tmp_path

    def test_multi_equals_single(self, mock_emulator_dir):
        """compute_coefficients_multi([band]) should equal compute_coefficients(band)."""
        from siac.rt.emulator.two_nn import TwoLayerNNEmulator

        shape = (8, 8)
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        atmo_state = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        emulator = TwoLayerNNEmulator(
            emulator_dir=mock_emulator_dir,
            sensor_id="MSI",
            satellite_id="S2A",
            use_rust=False,
        )

        band = SensorBand("B02", 490.0, 65.0, 10.0, 0)

        # Single-band call
        single = emulator.compute_coefficients(geometry, atmo_state, band, compute_jacobian=True)

        # Multi-band call with one band
        multi = emulator.compute_coefficients_multi(geometry, atmo_state, [band], compute_jacobian=True)

        np.testing.assert_allclose(single.xap.values, multi[0].xap.values)
        np.testing.assert_allclose(single.xbp.values, multi[0].xbp.values)
        np.testing.assert_allclose(single.xcp.values, multi[0].xcp.values)

    def test_multi_returns_correct_count(self, mock_emulator_dir):
        """compute_coefficients_multi should return one result per band."""
        from siac.rt.emulator.two_nn import TwoLayerNNEmulator

        shape = (4, 4)
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        atmo_state = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        emulator = TwoLayerNNEmulator(
            emulator_dir=mock_emulator_dir,
            sensor_id="MSI",
            satellite_id="S2A",
            use_rust=False,
        )

        bands = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
            SensorBand("B04", 665.0, 30.0, 10.0, 2),
        ]

        results = emulator.compute_coefficients_multi(geometry, atmo_state, bands)
        assert len(results) == 3


class TestEmulatorRegistry:
    """Tests for EmulatorRegistry class."""

    def test_creation(self, tmp_path):
        """Registry should be creatable."""
        from siac.rt.emulator.two_nn import EmulatorRegistry

        registry = EmulatorRegistry(tmp_path)
        assert registry.emulator_dir == tmp_path

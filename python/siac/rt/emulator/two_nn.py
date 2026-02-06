"""
Two-hidden-layer neural network emulator for radiative transfer.

This module provides neural network emulators that approximate 6S radiative
transfer model outputs for fast atmospheric correction. The emulators compute
the RT coefficients (xap, xbp, xcp) used in the correction equation:

    y = xap * toa - xbp
    boa = y / (1 + xcp * y)

The emulators support analytical Jacobian computation for use in optimization.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.core.types import (
    AtmosphericState,
    GeometryAngles,
    RTCoefficients,
    SensorBand,
)

logger = logging.getLogger(__name__)

# Try to import Rust implementation
try:
    from siac._rust import TwoLayerNN as _RustNN

    _HAS_RUST = True
    logger.debug("Using Rust neural network implementation")
except ImportError:
    _HAS_RUST = False
    logger.debug("Rust NN not available, using Python implementation")


class TwoLayerNNEmulator:
    """
    Two-hidden-layer neural network emulator for RT coefficients.

    Implements the RTModelBackend protocol for pre-trained neural network
    emulators that approximate 6S radiative transfer outputs.

    The network architecture is:
        Input (7) -> Hidden1 (64, ReLU) -> Hidden2 (64, ReLU) -> Output (3)

    Input features: [cos(sza), cos(vza), cos(raa), aot, tcwv, tco3, elevation]
    Output: [xap, xbp, xcp]

    Args:
        emulator_dir: Directory containing emulator weight files
        sensor_id: Sensor identifier (e.g., "MSI", "OLI")
        satellite_id: Satellite identifier (e.g., "S2A", "L8")
        use_rust: Whether to use Rust acceleration if available
    """

    # Standard emulator input feature names
    INPUT_FEATURES = ["cos_sza", "cos_vza", "cos_raa", "aot", "tcwv", "tco3", "elevation"]

    # Standard output names
    OUTPUT_NAMES = ["xap", "xbp", "xcp"]

    def __init__(
        self,
        emulator_dir: str | Path,
        sensor_id: str,
        satellite_id: str,
        use_rust: bool = True,
    ):
        self.emulator_dir = Path(emulator_dir)
        self.sensor_id = sensor_id
        self.satellite_id = satellite_id
        self._use_rust = use_rust and _HAS_RUST

        # Loaded emulators per band
        self._band_emulators: dict[str, _BandEmulator] = {}

        # Discover available bands
        self._available_bands = self._discover_bands()

        if not self._available_bands:
            logger.warning(
                f"No emulators found for {sensor_id}/{satellite_id} "
                f"in {emulator_dir}"
            )

    def _discover_bands(self) -> list[str]:
        """Discover available band emulators in the directory."""
        available = []

        # Common patterns for emulator files
        patterns = [
            f"*{self.satellite_id}*_B*.npz",
            f"*{self.sensor_id}*_B*.npz",
            f"{self.satellite_id}_B*.npz",
        ]

        found_files = set()
        for pattern in patterns:
            found_files.update(self.emulator_dir.glob(pattern))

        # Extract band names from filenames
        for path in found_files:
            name = path.stem
            # Try to extract band name (e.g., "B02", "B1")
            for part in name.split("_"):
                if part.startswith("B") and (part[1:].isdigit() or part in ["B8A"]):
                    available.append(part)
                    break

        return sorted(set(available))

    def _get_emulator_path(self, band_name: str) -> Path | None:
        """Find emulator file for a specific band."""
        patterns = [
            f"*{self.satellite_id}*_{band_name}*.npz",
            f"*{self.satellite_id}*_{band_name.lower()}*.npz",
            f"{self.satellite_id}_{band_name}*.npz",
            f"*{self.sensor_id}*_{band_name}*.npz",
        ]

        for pattern in patterns:
            matches = list(self.emulator_dir.glob(pattern))
            if matches:
                return matches[0]

        return None

    def _load_band_emulator(self, band_name: str) -> _BandEmulator:
        """Load emulator weights for a specific band."""
        if band_name in self._band_emulators:
            return self._band_emulators[band_name]

        path = self._get_emulator_path(band_name)
        if path is None:
            raise FileNotFoundError(
                f"No emulator found for band {band_name} "
                f"(sensor: {self.sensor_id}/{self.satellite_id})"
            )

        logger.debug(f"Loading emulator from {path}")

        # Load numpy model file
        data = np.load(path, allow_pickle=True)

        hidden_layers = data["Hidden_Layers"].tolist()
        output_layers = data["Output_Layers"].tolist()

        # Normalize bias shapes
        for layer in hidden_layers:
            layer[1] = np.atleast_1d(layer[1]).ravel().astype(np.float32)

        for layer in output_layers:
            layer[1] = np.atleast_1d(layer[1]).ravel().astype(np.float32)

        emulator = _BandEmulator(
            hidden_layers=hidden_layers,
            output_layers=output_layers,
            use_rust=self._use_rust,
        )

        self._band_emulators[band_name] = emulator
        return emulator

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        """
        Compute radiative transfer coefficients for atmospheric correction.

        Args:
            geometry: Viewing geometry (sza, vza, raa in radians)
            atmo_state: Atmospheric state (aot, tcwv, tco3, elevation)
            band: Sensor band specification
            compute_jacobian: Whether to compute d(coeff)/d(aot,tcwv)

        Returns:
            RTCoefficients with xap, xbp, xcp and optional Jacobians.
        """
        # Load band emulator
        emulator = self._load_band_emulator(band.name)

        # Prepare input array
        # Input order: [cos_sza, cos_vza, cos_raa, aot, tcwv, tco3, elevation]
        cos_sza = np.cos(geometry.sza.values)
        cos_vza = np.cos(geometry.vza.values)
        cos_raa = np.cos(geometry.raa.values)

        # Get atmospheric parameters
        aot = atmo_state.aot.values
        tcwv = atmo_state.tcwv.values
        tco3 = atmo_state.tco3.values
        elevation = atmo_state.elevation.values

        # Stack inputs
        original_shape = cos_sza.shape
        n_pixels = np.prod(original_shape)

        # Flatten and stack
        inputs = np.column_stack(
            [
                cos_sza.ravel(),
                cos_vza.ravel(),
                cos_raa.ravel(),
                aot.ravel(),
                tcwv.ravel(),
                tco3.ravel(),
                elevation.ravel(),
            ]
        ).astype(np.float32)

        # Run forward pass
        outputs, jacobians = emulator.forward(inputs, compute_jacobian=compute_jacobian)

        # Parse outputs - each output is [xap, xbp, xcp]
        xap = outputs[:, 0].reshape(original_shape)
        xbp = outputs[:, 1].reshape(original_shape)
        xcp = outputs[:, 2].reshape(original_shape)

        # Create DataArrays
        template = geometry.sza
        xap_da = xr.DataArray(xap, dims=template.dims, coords=template.coords)
        xbp_da = xr.DataArray(xbp, dims=template.dims, coords=template.coords)
        xcp_da = xr.DataArray(xcp, dims=template.dims, coords=template.coords)

        # Handle Jacobians if computed
        d_xap = None
        d_xbp = None
        d_xcp = None

        if compute_jacobian and jacobians is not None:
            # Jacobians have shape (n_pixels, 3, 7) for (outputs, inputs)
            # We want derivatives w.r.t. aot (idx 3) and tcwv (idx 4)

            # Extract relevant derivatives
            d_xap_aot = jacobians[:, 0, 3].reshape(original_shape)
            d_xap_tcwv = jacobians[:, 0, 4].reshape(original_shape)
            d_xbp_aot = jacobians[:, 1, 3].reshape(original_shape)
            d_xbp_tcwv = jacobians[:, 1, 4].reshape(original_shape)
            d_xcp_aot = jacobians[:, 2, 3].reshape(original_shape)
            d_xcp_tcwv = jacobians[:, 2, 4].reshape(original_shape)

            # Create DataArrays with param dimension
            d_xap = xr.concat(
                [
                    xr.DataArray(d_xap_aot, dims=template.dims, coords=template.coords),
                    xr.DataArray(d_xap_tcwv, dims=template.dims, coords=template.coords),
                ],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])

            d_xbp = xr.concat(
                [
                    xr.DataArray(d_xbp_aot, dims=template.dims, coords=template.coords),
                    xr.DataArray(d_xbp_tcwv, dims=template.dims, coords=template.coords),
                ],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])

            d_xcp = xr.concat(
                [
                    xr.DataArray(d_xcp_aot, dims=template.dims, coords=template.coords),
                    xr.DataArray(d_xcp_tcwv, dims=template.dims, coords=template.coords),
                ],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])

        return RTCoefficients(
            xap=xap_da,
            xbp=xbp_da,
            xcp=xcp_da,
            d_xap=d_xap,
            d_xbp=d_xbp,
            d_xcp=d_xcp,
        )

    def compute_coefficients_multi(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        compute_jacobian: bool = False,
    ) -> list[RTCoefficients]:
        """
        Compute RT coefficients for multiple bands, preparing inputs once.

        Args:
            geometry: Viewing geometry
            atmo_state: Atmospheric state
            bands: List of sensor bands
            compute_jacobian: Whether to compute Jacobians

        Returns:
            List of RTCoefficients, one per band.
        """
        # Prepare input array once for all bands
        cos_sza = np.cos(geometry.sza.values)
        cos_vza = np.cos(geometry.vza.values)
        cos_raa = np.cos(geometry.raa.values)

        aot = atmo_state.aot.values
        tcwv = atmo_state.tcwv.values
        tco3 = atmo_state.tco3.values
        elevation = atmo_state.elevation.values

        original_shape = cos_sza.shape
        template = geometry.sza

        inputs = np.column_stack(
            [
                cos_sza.ravel(),
                cos_vza.ravel(),
                cos_raa.ravel(),
                aot.ravel(),
                tcwv.ravel(),
                tco3.ravel(),
                elevation.ravel(),
            ]
        ).astype(np.float32)

        results = []
        for band in bands:
            emulator = self._load_band_emulator(band.name)
            outputs, jacobians = emulator.forward(inputs, compute_jacobian=compute_jacobian)

            xap = outputs[:, 0].reshape(original_shape)
            xbp = outputs[:, 1].reshape(original_shape)
            xcp = outputs[:, 2].reshape(original_shape)

            xap_da = xr.DataArray(xap, dims=template.dims, coords=template.coords)
            xbp_da = xr.DataArray(xbp, dims=template.dims, coords=template.coords)
            xcp_da = xr.DataArray(xcp, dims=template.dims, coords=template.coords)

            d_xap = d_xbp = d_xcp = None
            if compute_jacobian and jacobians is not None:
                d_xap_aot = jacobians[:, 0, 3].reshape(original_shape)
                d_xap_tcwv = jacobians[:, 0, 4].reshape(original_shape)
                d_xbp_aot = jacobians[:, 1, 3].reshape(original_shape)
                d_xbp_tcwv = jacobians[:, 1, 4].reshape(original_shape)
                d_xcp_aot = jacobians[:, 2, 3].reshape(original_shape)
                d_xcp_tcwv = jacobians[:, 2, 4].reshape(original_shape)

                d_xap = xr.concat(
                    [
                        xr.DataArray(d_xap_aot, dims=template.dims, coords=template.coords),
                        xr.DataArray(d_xap_tcwv, dims=template.dims, coords=template.coords),
                    ],
                    dim="param",
                ).assign_coords(param=["aot", "tcwv"])

                d_xbp = xr.concat(
                    [
                        xr.DataArray(d_xbp_aot, dims=template.dims, coords=template.coords),
                        xr.DataArray(d_xbp_tcwv, dims=template.dims, coords=template.coords),
                    ],
                    dim="param",
                ).assign_coords(param=["aot", "tcwv"])

                d_xcp = xr.concat(
                    [
                        xr.DataArray(d_xcp_aot, dims=template.dims, coords=template.coords),
                        xr.DataArray(d_xcp_tcwv, dims=template.dims, coords=template.coords),
                    ],
                    dim="param",
                ).assign_coords(param=["aot", "tcwv"])

            results.append(RTCoefficients(
                xap=xap_da, xbp=xbp_da, xcp=xcp_da,
                d_xap=d_xap, d_xbp=d_xbp, d_xcp=d_xcp,
            ))

        return results

    def supports_jacobian(self) -> bool:
        """Check if this backend supports analytical Jacobian computation."""
        return True

    @property
    def backend_name(self) -> str:
        """Name of the RT backend."""
        return "emulator"

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
        """Check if this backend has models for the specified sensor."""
        return (
            sensor_id == self.sensor_id
            and satellite_id == self.satellite_id
            and len(self._available_bands) > 0
        )

    @property
    def available_bands(self) -> list[str]:
        """List of bands with available emulators."""
        return self._available_bands

    @classmethod
    def load(
        cls,
        satellite_id: str,
        band_name: str,
        emulator_dir: str | Path | None = None,
    ) -> TwoLayerNNEmulator:
        """
        Convenience method to load emulator for a specific band.

        Args:
            satellite_id: Satellite identifier (e.g., "S2A", "L8")
            band_name: Band name (e.g., "B02")
            emulator_dir: Directory containing emulators (uses default if None)

        Returns:
            Loaded emulator instance
        """
        if emulator_dir is None:
            # Default to package data directory
            import siac

            emulator_dir = Path(siac.__file__).parent / "data" / "emus"

        # Infer sensor from satellite
        sensor_id = _SATELLITE_TO_SENSOR.get(satellite_id, satellite_id)

        return cls(
            emulator_dir=emulator_dir,
            sensor_id=sensor_id,
            satellite_id=satellite_id,
        )


class _BandEmulator:
    """
    Internal class for single-band emulator.

    Handles the actual forward pass through the neural network.
    """

    def __init__(
        self,
        hidden_layers: list,
        output_layers: list,
        use_rust: bool = True,
    ):
        self.hidden_layers = hidden_layers
        self.output_layers = output_layers
        self._use_rust = use_rust and _HAS_RUST

        # Initialize Rust backend if available
        if self._use_rust:
            self._init_rust_emulator()
        else:
            self._rust_nn = None

    def _init_rust_emulator(self) -> None:
        """Initialize Rust neural network with weights."""
        # Extract weights for Rust implementation
        # Hidden layers: [[w1, b1], [w2, b2]]
        # Output layers: [[w3, b3]]

        w1 = np.asarray(self.hidden_layers[0][0], dtype=np.float32)
        b1 = np.asarray(self.hidden_layers[0][1], dtype=np.float32)
        w2 = np.asarray(self.hidden_layers[1][0], dtype=np.float32)
        b2 = np.asarray(self.hidden_layers[1][1], dtype=np.float32)

        # For multi-output, we need to handle output layers differently
        # Original model has separate output layers for each coefficient
        # For Rust, we combine them into a single output layer

        if len(self.output_layers) > 1:
            # Multiple output heads - combine into single output layer
            w3_list = [np.asarray(ol[0], dtype=np.float32) for ol in self.output_layers]
            b3_list = [np.asarray(ol[1], dtype=np.float32) for ol in self.output_layers]

            w3 = np.hstack(w3_list)
            b3 = np.concatenate(b3_list)
        else:
            w3 = np.asarray(self.output_layers[0][0], dtype=np.float32)
            b3 = np.asarray(self.output_layers[0][1], dtype=np.float32)

        try:
            self._rust_nn = _RustNN(w1, b1, w2, b2, w3, b3)
        except Exception as e:
            logger.warning(f"Failed to initialize Rust NN: {e}")
            self._rust_nn = None
            self._use_rust = False

    def forward(
        self,
        x: np.ndarray,
        compute_jacobian: bool = False,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Forward pass through the neural network.

        Args:
            x: Input array of shape (n_samples, n_inputs)
            compute_jacobian: Whether to compute Jacobian

        Returns:
            Tuple of (outputs, jacobians) where jacobians is None if not requested
        """
        if self._use_rust and self._rust_nn is not None:
            return self._forward_rust(x, compute_jacobian)
        else:
            return self._forward_python(x, compute_jacobian)

    def _forward_rust(
        self,
        x: np.ndarray,
        compute_jacobian: bool,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Forward pass using Rust implementation."""
        output, jacobian = self._rust_nn.predict(x, compute_jacobian)
        return output, jacobian

    def _forward_python(
        self,
        x: np.ndarray,
        compute_jacobian: bool,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Forward pass using Python/NumPy implementation."""
        # Hidden layer 1
        w1, b1 = self.hidden_layers[0]
        a1 = x @ w1 + b1
        h1 = np.maximum(a1, 0)  # ReLU

        # Hidden layer 2
        w2, b2 = self.hidden_layers[1]
        a2 = h1 @ w2 + b2
        h2 = np.maximum(a2, 0)  # ReLU

        # Output layer(s)
        outputs = []
        jacobians = []

        for i, (w3, b3) in enumerate(self.output_layers):
            out = h2 @ w3 + b3
            outputs.append(out)

            if compute_jacobian:
                # Backpropagation for Jacobian
                jac = self._compute_jacobian_python(x, a1, a2, w1, w2, w3, i)
                jacobians.append(jac)

        # Stack outputs
        output_array = np.column_stack(outputs)

        if compute_jacobian:
            # Stack jacobians: shape (n_samples, n_outputs, n_inputs)
            jacobian_array = np.stack(jacobians, axis=1)
            return output_array, jacobian_array
        else:
            return output_array, None

    def _compute_jacobian_python(
        self,
        x: np.ndarray,
        a1: np.ndarray,
        a2: np.ndarray,
        w1: np.ndarray,
        w2: np.ndarray,
        w3: np.ndarray,
        output_idx: int,
    ) -> np.ndarray:
        """Compute Jacobian for a single output using backpropagation."""
        # ReLU derivatives
        d_relu1 = (a1 > 0).astype(np.float32)
        d_relu2 = (a2 > 0).astype(np.float32)

        # Vectorized backpropagation (no per-sample loop)
        grad_h2 = w3.ravel()  # dout/dh2 shape (hidden2,)
        grad_a2 = grad_h2[np.newaxis, :] * d_relu2    # (n_samples, hidden2)
        grad_h1 = grad_a2 @ w2.T                       # (n_samples, hidden1)
        grad_a1 = grad_h1 * d_relu1                    # (n_samples, hidden1)
        jacobians = (grad_a1 @ w1.T).astype(np.float32)  # (n_samples, n_inputs)

        return jacobians


# Satellite to sensor mapping
_SATELLITE_TO_SENSOR = {
    "S2A": "MSI",
    "S2B": "MSI",
    "S2C": "MSI",
    "L8": "OLI",
    "L9": "OLI",
}


class EmulatorRegistry:
    """
    Registry for available emulators.

    Provides a centralized way to discover and access emulators for
    different sensors and bands.
    """

    def __init__(self, emulator_dir: str | Path):
        self.emulator_dir = Path(emulator_dir)
        self._emulators: dict[tuple[str, str], TwoLayerNNEmulator] = {}

    def get_emulator(
        self,
        sensor_id: str,
        satellite_id: str,
    ) -> TwoLayerNNEmulator | None:
        """Get emulator for a sensor/satellite combination."""
        key = (sensor_id, satellite_id)

        if key not in self._emulators:
            try:
                emulator = TwoLayerNNEmulator(
                    emulator_dir=self.emulator_dir,
                    sensor_id=sensor_id,
                    satellite_id=satellite_id,
                )
                if emulator.available_bands:
                    self._emulators[key] = emulator
                else:
                    return None
            except Exception as e:
                logger.debug(f"Failed to load emulator for {sensor_id}/{satellite_id}: {e}")
                return None

        return self._emulators.get(key)

    def is_sensor_supported(self, sensor_id: str, satellite_id: str) -> bool:
        """Check if a sensor/satellite combination is supported."""
        emulator = self.get_emulator(sensor_id, satellite_id)
        return emulator is not None and len(emulator.available_bands) > 0

    def list_supported_sensors(self) -> list[tuple[str, str]]:
        """List all supported sensor/satellite combinations."""
        supported = []

        # Scan emulator directory for common patterns
        for satellite in ["S2A", "S2B", "L8", "L9"]:
            sensor = _SATELLITE_TO_SENSOR.get(satellite, satellite)
            if self.is_sensor_supported(sensor, satellite):
                supported.append((sensor, satellite))

        return supported

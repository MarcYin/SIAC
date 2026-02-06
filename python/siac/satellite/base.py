"""
Base satellite preprocessor class.

This module provides the abstract base class for satellite-specific
preprocessors and common utilities shared across sensors.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.core.types import GeometryAngles, SensorConfig

logger = logging.getLogger(__name__)


class BaseSatellitePreprocessor(ABC):
    """
    Abstract base class for satellite preprocessors.

    Subclasses must implement the abstract methods to handle
    sensor-specific data formats and processing.
    """

    def __init__(self, config: dict[str, Any] | None = None):
        """
        Initialize preprocessor.

        Args:
            config: Optional configuration dictionary
        """
        self.config = config or {}
        self._cached_metadata: dict[str, Any] | None = None

    @property
    @abstractmethod
    def sensor_config(self) -> SensorConfig:
        """Return sensor configuration."""
        ...

    @abstractmethod
    def load_toa(self, input_path: str | Path) -> xr.Dataset:
        """
        Load top-of-atmosphere reflectance.

        Args:
            input_path: Path to satellite data

        Returns:
            Dataset with TOA reflectance for each band
        """
        ...

    @abstractmethod
    def extract_geometry(self, input_path: str | Path) -> GeometryAngles:
        """
        Extract view and sun geometry.

        Args:
            input_path: Path to satellite data

        Returns:
            GeometryAngles with sza, saa, vza, vaa in radians
        """
        ...

    @abstractmethod
    def extract_cloud_mask(self, input_path: str | Path) -> xr.DataArray:
        """
        Generate cloud mask.

        Args:
            input_path: Path to satellite data

        Returns:
            Boolean DataArray (True = cloudy/invalid)
        """
        ...

    @abstractmethod
    def get_metadata(self, input_path: str | Path) -> dict[str, Any]:
        """
        Extract observation metadata.

        Args:
            input_path: Path to satellite data

        Returns:
            Dictionary with observation metadata
        """
        ...

    def preprocess(self, input_path: str | Path) -> dict[str, Any]:
        """
        Run full preprocessing pipeline.

        Args:
            input_path: Path to satellite data

        Returns:
            Dictionary with:
                - toa: TOA reflectance Dataset
                - geometry: GeometryAngles
                - cloud_mask: Cloud mask DataArray
                - metadata: Observation metadata
        """
        input_path = Path(input_path)

        logger.info(f"Preprocessing {input_path}")

        metadata = self.get_metadata(input_path)
        logger.info(f"Observation time: {metadata.get('observation_time')}")

        toa = self.load_toa(input_path)
        logger.info(f"Loaded TOA with bands: {list(toa.data_vars)}")

        geometry = self.extract_geometry(input_path)
        logger.info(f"Extracted geometry (SZA mean: {float(geometry.sza.mean()):.2f} rad)")

        cloud_mask = self.extract_cloud_mask(input_path)
        cloud_fraction = float(cloud_mask.mean())
        logger.info(f"Cloud fraction: {cloud_fraction:.1%}")

        return {
            "toa": toa,
            "geometry": geometry,
            "cloud_mask": cloud_mask,
            "metadata": metadata,
        }


# =============================================================================
# Common Utilities
# =============================================================================


def degrees_to_radians(da: xr.DataArray) -> xr.DataArray:
    """Convert angles from degrees to radians."""
    return da * (np.pi / 180.0)


def radians_to_degrees(da: xr.DataArray) -> xr.DataArray:
    """Convert angles from radians to degrees."""
    return da * (180.0 / np.pi)


def compute_relative_azimuth(
    saa: xr.DataArray,
    vaa: xr.DataArray,
) -> xr.DataArray:
    """
    Compute relative azimuth angle.

    Args:
        saa: Solar azimuth angle
        vaa: View azimuth angle

    Returns:
        Relative azimuth angle (same units as input)
    """
    raa = vaa - saa
    # Normalize to [-pi, pi] or [-180, 180] depending on input
    return raa


def apply_scale_offset(
    data: xr.DataArray,
    scale: float,
    offset: float,
) -> xr.DataArray:
    """
    Apply scale and offset to convert DN to reflectance.

    Args:
        data: Digital number values
        scale: Scale factor
        offset: Offset value

    Returns:
        Reflectance values: reflectance = data * scale + offset
    """
    return data * scale + offset


def create_valid_mask(
    data: xr.DataArray,
    min_val: float = 0.0,
    max_val: float = 1.5,
) -> xr.DataArray:
    """
    Create mask of valid reflectance values.

    Args:
        data: Reflectance data
        min_val: Minimum valid value
        max_val: Maximum valid value

    Returns:
        Boolean mask (True = valid)
    """
    return (data >= min_val) & (data <= max_val) & np.isfinite(data)


def resample_angles_to_data(
    angles: xr.DataArray,
    target: xr.DataArray,
    method: str = "linear",
) -> xr.DataArray:
    """
    Resample angle data to match target resolution.

    Satellite angle grids are often at coarser resolution than
    the image data and need interpolation.

    Args:
        angles: Angle DataArray at coarse resolution
        target: Target DataArray defining output grid
        method: Interpolation method

    Returns:
        Resampled angles matching target grid
    """
    from siac.io.reprojection import reproject_match

    return reproject_match(angles, target, resampling=method)


# =============================================================================
# Sensor Registry
# =============================================================================


_SENSOR_REGISTRY: dict[str, type[BaseSatellitePreprocessor]] = {}


def register_preprocessor(sensor_id: str):
    """
    Decorator to register a preprocessor class.

    Args:
        sensor_id: Sensor identifier (e.g., "s2", "l8")

    Example:
        @register_preprocessor("s2")
        class Sentinel2Preprocessor(BaseSatellitePreprocessor):
            ...
    """

    def decorator(cls: type[BaseSatellitePreprocessor]):
        _SENSOR_REGISTRY[sensor_id.lower()] = cls
        return cls

    return decorator


def get_preprocessor(sensor_id: str, **kwargs) -> BaseSatellitePreprocessor:
    """
    Get preprocessor instance for a sensor.

    Args:
        sensor_id: Sensor identifier
        **kwargs: Arguments passed to preprocessor constructor

    Returns:
        Preprocessor instance

    Raises:
        KeyError: If sensor not registered
    """
    sensor_id = sensor_id.lower()
    if sensor_id not in _SENSOR_REGISTRY:
        raise KeyError(
            f"Unknown sensor: {sensor_id}. "
            f"Available: {list(_SENSOR_REGISTRY.keys())}"
        )

    return _SENSOR_REGISTRY[sensor_id](**kwargs)


def detect_sensor(input_path: str | Path) -> str:
    """
    Auto-detect sensor from input data.

    Args:
        input_path: Path to satellite data

    Returns:
        Sensor identifier

    Raises:
        ValueError: If sensor cannot be detected
    """
    input_path = Path(input_path)

    # Check for Sentinel-2 SAFE format
    if input_path.suffix == ".SAFE" or (input_path / "MTD_MSIL1C.xml").exists():
        return "s2"

    # Check for Landsat MTL
    mtl_files = list(input_path.glob("*_MTL.txt")) + list(input_path.glob("*_MTL.xml"))
    if mtl_files:
        # Check if Landsat 8/9 or older
        mtl_content = mtl_files[0].read_text()
        if "LANDSAT_8" in mtl_content or "LANDSAT_9" in mtl_content:
            return "l8"
        elif "LANDSAT_5" in mtl_content or "LANDSAT_7" in mtl_content:
            return "l5"

    # Check for AWS Sentinel-2 format
    if (input_path / "metadata.xml").exists():
        return "s2"

    raise ValueError(f"Could not detect sensor for: {input_path}")


def list_available_sensors() -> list[str]:
    """Return list of registered sensor IDs."""
    return list(_SENSOR_REGISTRY.keys())

"""
Base satellite preprocessor class.

This module provides the abstract base class for satellite-specific
preprocessors and common utilities shared across sensors.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from inspect import Parameter, signature
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    import xarray as xr

    from siac.domain import SensorConfig
    from siac.runtime import GeometryAngles, ObservationBundle

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
    def extract_cloud_mask(
        self,
        input_path: str | Path,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray:
        """
        Generate cloud mask.

        Args:
            input_path: Path to satellite data
            toa: Optional preloaded TOA dataset to avoid duplicate reads

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

        toa = self._load_toa_for_preprocess(input_path, metadata=metadata)
        logger.info(f"Loaded TOA with bands: {list(toa.data_vars)}")

        with ThreadPoolExecutor(max_workers=2) as executor:
            f_geometry = executor.submit(self.extract_geometry, input_path)
            f_cloud = executor.submit(self._extract_cloud_mask_with_toa, input_path, toa)

            geometry = f_geometry.result()
            cloud_mask = f_cloud.result()

        logger.info(f"Extracted geometry (SZA mean: {float(geometry.sza.mean()):.2f} rad)")

        cloud_fraction = float(cloud_mask.mean())
        logger.info(f"Cloud fraction: {cloud_fraction:.1%}")

        result = {
            "toa": toa,
            "geometry": geometry,
            "cloud_mask": cloud_mask,
            "metadata": metadata,
        }
        cloud_classes = getattr(self, "_last_cloud_classes", None)
        if cloud_classes is not None:
            result["cloud_classes"] = cloud_classes
        return result

    def _preprocess_toa_band_names(
        self,
        input_path: Path,
        *,
        metadata: dict[str, Any],
    ) -> tuple[str, ...] | None:
        """Return an optional TOA band subset for preprocessing."""
        del input_path, metadata
        if not isinstance(self.config, dict):
            return None
        raw = self.config.get("preload_toa_bands")
        if raw is None:
            return None
        names = (raw,) if isinstance(raw, str) else tuple(str(name) for name in raw)
        return names or None

    def _load_toa_for_preprocess(
        self,
        input_path: Path,
        *,
        metadata: dict[str, Any],
    ) -> xr.Dataset:
        """Load TOA, passing a configured band subset when supported."""
        band_names = self._preprocess_toa_band_names(input_path, metadata=metadata)
        params = signature(self.load_toa).parameters.values()
        accepts_band_names = any(
            p.name == "band_names" or p.kind == Parameter.VAR_KEYWORD for p in params
        )
        if band_names is not None and accepts_band_names:
            return self.load_toa(input_path, band_names=band_names)  # type: ignore[call-arg]
        return self.load_toa(input_path)

    def _extract_cloud_mask_with_toa(
        self,
        input_path: Path,
        toa: xr.Dataset,
    ) -> xr.DataArray:
        """Call extract_cloud_mask with TOA when the implementation supports it."""
        params = signature(self.extract_cloud_mask).parameters.values()
        accepts_toa = any(p.name == "toa" or p.kind == Parameter.VAR_KEYWORD for p in params)
        if accepts_toa:
            return self.extract_cloud_mask(input_path, toa=toa)
        return self.extract_cloud_mask(input_path)

    def to_observation_bundle(
        self,
        input_path: str | Path,
        *,
        bounds: tuple[float, float, float, float],
        crs: str = "EPSG:4326",
    ) -> ObservationBundle:
        """Run preprocessing and return an :class:`ObservationBundle`.

        This is a convenience wrapper around :meth:`preprocess` that
        directly produces the pipeline's M1 contract type.

        Args:
            input_path: Path to satellite data.
            bounds: Spatial bounds as (west, south, east, north).
            crs: Coordinate reference system string.

        Returns:
            A fully populated :class:`ObservationBundle`.
        """
        from siac.adapters.satellite.observation import observation_bundle_from_raw

        raw = self.preprocess(input_path)
        return observation_bundle_from_raw(
            raw,
            sensor_config=self.sensor_config,
            crs=crs,
            bounds=bounds,
        )


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
    from siac.geo.reprojection import reproject_match

    return reproject_match(angles, target, resampling=method)


# =============================================================================
# Sensor Registry
# =============================================================================


_SENSOR_REGISTRY: dict[str, type[BaseSatellitePreprocessor]] = {}


def register_preprocessor(
    sensor_id: str,
) -> Callable[[type[BaseSatellitePreprocessor]], type[BaseSatellitePreprocessor]]:
    """
    Decorator to register a preprocessor class.

    Args:
        sensor_id: Sensor identifier (e.g., "s2", "l8")

    Example:
        @register_preprocessor("s2")
        class Sentinel2Preprocessor(BaseSatellitePreprocessor):
            ...
    """

    def decorator(cls: type[BaseSatellitePreprocessor]) -> type[BaseSatellitePreprocessor]:
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
        raise KeyError(f"Unknown sensor: {sensor_id}. Available: {list(_SENSOR_REGISTRY.keys())}")

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

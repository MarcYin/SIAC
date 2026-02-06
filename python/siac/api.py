"""
High-level SIAC API.

This module provides the main entry point for atmospheric correction.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import xarray as xr

from siac.core.config import SIACConfig

if TYPE_CHECKING:
    from siac.core.types import AtmosphericState


@dataclass
class SIACResult:
    """
    Result of SIAC atmospheric correction.

    Attributes:
        boa: Bottom-of-atmosphere reflectance dataset
        boa_unc: BOA uncertainty dataset (if computed)
        aot: Aerosol optical thickness map
        tcwv: Total column water vapor map
        metadata: Processing metadata
    """

    boa: xr.Dataset
    boa_unc: xr.Dataset | None
    aot: xr.DataArray
    tcwv: xr.DataArray
    metadata: dict


class SIAC:
    """
    Main SIAC atmospheric correction processor.

    Example:
        >>> from siac import SIAC
        >>> from siac.core.config import SIACConfig
        >>>
        >>> config = SIACConfig.from_yaml("config.yaml")
        >>> siac = SIAC(config)
        >>> result = siac.process("/path/to/S2_SAFE/")
        >>>
        >>> # Access results
        >>> print(result.boa)
        >>> print(result.aot.mean())
    """

    def __init__(self, config: SIACConfig | None = None) -> None:
        """
        Initialize SIAC processor.

        Args:
            config: Processing configuration. If None, uses defaults.
        """
        self.config = config or SIACConfig()
        self._preprocessor = None
        self._atmo_provider = None
        self._brdf_provider = None
        self._surface_deriver = None
        self._rt_backend = None
        self._solver = None

    def process(self, input_path: str | Path) -> SIACResult:
        """
        Run atmospheric correction on satellite data.

        Args:
            input_path: Path to satellite data (S2 SAFE directory or L8 MTL file)

        Returns:
            SIACResult with corrected imagery and auxiliary products
        """
        input_path = Path(input_path)

        # TODO: Implement full processing pipeline
        # This is a placeholder that will be implemented in subsequent phases

        raise NotImplementedError(
            "Full processing pipeline not yet implemented. "
            "This is a placeholder for the SIAC v2 refactoring."
        )

    def _initialize_modules(self, input_path: Path) -> None:
        """Initialize processing modules based on configuration."""
        # TODO: Implement module initialization
        pass

    def _detect_sensor(self, input_path: Path) -> tuple[str, str]:
        """Auto-detect sensor from input data."""
        # TODO: Implement sensor detection
        pass


# Convenience functions for direct use


def process_sentinel2(
    input_path: str | Path,
    config: SIACConfig | None = None,
) -> SIACResult:
    """
    Process Sentinel-2 data with atmospheric correction.

    Args:
        input_path: Path to S2 SAFE directory
        config: Optional configuration overrides

    Returns:
        SIACResult with corrected imagery
    """
    if config is None:
        config = SIACConfig(sensor="s2")
    else:
        config = config.with_overrides(sensor="s2")

    siac = SIAC(config)
    return siac.process(input_path)


def process_landsat8(
    input_path: str | Path,
    config: SIACConfig | None = None,
) -> SIACResult:
    """
    Process Landsat 8/9 data with atmospheric correction.

    Args:
        input_path: Path to L8/L9 data directory
        config: Optional configuration overrides

    Returns:
        SIACResult with corrected imagery
    """
    if config is None:
        config = SIACConfig(sensor="l8")
    else:
        config = config.with_overrides(sensor="l8")

    siac = SIAC(config)
    return siac.process(input_path)

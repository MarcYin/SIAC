"""
CAMS atmospheric prior provider.

Fetches atmospheric parameters (AOT, TCWV, TCO3) from ECMWF CAMS
(Copernicus Atmosphere Monitoring Service) near real-time data.
"""

from __future__ import annotations
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState

logger = logging.getLogger(__name__)


class CAMSProvider:
    """
    CAMS atmospheric prior provider.

    Implements the AtmosphericPriorProvider protocol.
    """

    # CAMS variable mappings
    VARIABLE_MAP = {
        "aot": "aod550",
        "tcwv": "tcwv",
        "tco3": "gtco3",
    }

    def __init__(self, data_dir: str | Path, temporal_interp: bool = True):
        self.data_dir = Path(data_dir)
        self.temporal_interp = temporal_interp

    @property
    def source_name(self) -> str:
        return "CAMS"

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        """Retrieve atmospheric priors for a given region and time."""
        # Load CAMS data for the observation date
        cams_data = self._load_cams_data(obs_time)

        if cams_data is None:
            logger.warning("CAMS data not available, using defaults")
            return self._default_prior(bounds, crs, resolution)

        # Extract and interpolate to target grid
        aot = self._extract_variable(cams_data, "aod550", bounds, crs, resolution, obs_time)
        tcwv = self._extract_variable(cams_data, "tcwv", bounds, crs, resolution, obs_time)
        tco3 = self._extract_variable(cams_data, "gtco3", bounds, crs, resolution, obs_time)

        # Convert ozone from kg/m2 to atm-cm (DU/1000)
        tco3 = tco3 / 2.1415e-5 / 1000

        # Uncertainties (empirical estimates)
        aot_unc = np.maximum(aot * 0.5, 0.05)
        tcwv_unc = np.maximum(tcwv * 0.15, 0.2)
        tco3_unc = tco3 * 0.1

        # Elevation (placeholder - should be from DEM)
        elevation = xr.zeros_like(aot)

        return AtmosphericState(
            aot=aot, tcwv=tcwv, tco3=tco3,
            aot_unc=aot_unc, tcwv_unc=tcwv_unc, tco3_unc=tco3_unc,
            elevation=elevation
        )

    def _load_cams_data(self, obs_time: datetime) -> xr.Dataset | None:
        """Load CAMS data for the observation date."""
        date_str = obs_time.strftime("%Y%m%d")

        # Try different file patterns
        patterns = [
            f"cams_nrt_{date_str}*.nc",
            f"cams_{date_str}*.nc",
            f"*{date_str}*.nc",
        ]

        for pattern in patterns:
            files = list(self.data_dir.glob(pattern))
            if files:
                try:
                    return xr.open_mfdataset(files, combine="by_coords")
                except Exception as e:
                    logger.warning(f"Failed to load CAMS: {e}")
        return None

    def _extract_variable(
        self, data: xr.Dataset, var_name: str, bounds: tuple, crs: str,
        resolution: float, obs_time: datetime
    ) -> xr.DataArray:
        """Extract and regrid a single variable."""
        if var_name not in data:
            logger.warning(f"Variable {var_name} not in CAMS data")
            return self._create_default_array(bounds, crs, resolution, 0.15 if "aod" in var_name else 1.5)

        var = data[var_name]

        # Temporal interpolation
        if "time" in var.dims and self.temporal_interp:
            var = var.interp(time=np.datetime64(obs_time), method="linear")
        elif "time" in var.dims:
            var = var.sel(time=obs_time, method="nearest")

        # Spatial subset and regrid (simplified)
        xmin, ymin, xmax, ymax = bounds
        if "latitude" in var.dims:
            var = var.sel(latitude=slice(ymax, ymin), longitude=slice(xmin, xmax))

        return var

    def _create_default_array(self, bounds: tuple, crs: str, resolution: float, value: float) -> xr.DataArray:
        """Create default array with constant value."""
        xmin, ymin, xmax, ymax = bounds
        nx = int((xmax - xmin) / resolution)
        ny = int((ymax - ymin) / resolution)
        return xr.DataArray(np.full((ny, nx), value, dtype=np.float32), dims=["y", "x"])

    def _default_prior(self, bounds: tuple, crs: str, resolution: float) -> AtmosphericState:
        """Return default atmospheric state when CAMS unavailable."""
        aot = self._create_default_array(bounds, crs, resolution, 0.15)
        tcwv = self._create_default_array(bounds, crs, resolution, 1.5)
        tco3 = self._create_default_array(bounds, crs, resolution, 0.3)

        return AtmosphericState(
            aot=aot, tcwv=tcwv, tco3=tco3,
            aot_unc=aot * 0.5, tcwv_unc=tcwv * 0.2, tco3_unc=tco3 * 0.1,
            elevation=xr.zeros_like(aot)
        )

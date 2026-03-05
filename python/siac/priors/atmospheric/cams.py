"""
CAMS atmospheric prior provider.

Fetches atmospheric parameters (AOT, TCWV, TCO3) from ECMWF CAMS
(Copernicus Atmosphere Monitoring Service) near real-time data.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState

if TYPE_CHECKING:
    from datetime import datetime

    from siac.core.auth import CredentialManager

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

    _NETCDF_SUFFIXES = {".nc", ".nc4"}
    _TIFF_SUFFIXES = {".tif", ".tiff"}
    _REQUIRED_VARIABLES = ("aod550", "tcwv", "gtco3")
    _TIFF_VARIABLE_ALIASES = {
        "aod550": "aod550",
        "aot": "aod550",
        "tcwv": "tcwv",
        "water_vapour": "tcwv",
        "water_vapor": "tcwv",
        "gtco3": "gtco3",
        "tco3": "gtco3",
        "ozone": "gtco3",
    }

    _CDS_VARIABLES = [
        "mean_sea_level_pressure",
        "surface_pressure",
        "ammonium_aerosol_optical_depth_550nm",
        "black_carbon_aerosol_optical_depth_550nm",
        "dust_aerosol_optical_depth_550nm",
        "nitrate_aerosol_optical_depth_550nm",
        "organic_matter_aerosol_optical_depth_550nm",
        "sea_salt_aerosol_optical_depth_550nm",
        "secondary_organic_aerosol_optical_depth_550nm",
        "sulphate_aerosol_optical_depth_550nm",
        "total_aerosol_optical_depth_469nm",
        "total_aerosol_optical_depth_550nm",
        "total_aerosol_optical_depth_670nm",
        "total_aerosol_optical_depth_865nm",
        "total_aerosol_optical_depth_1240nm",
        "total_column_carbon_monoxide",
        "total_column_methane",
        "total_column_nitrogen_dioxide",
        "total_column_ozone",
        "total_column_water_vapour",
    ]
    _CDS_DATASET = "cams-global-atmospheric-composition-forecasts"
    _CDS_LEAD_TIMES = [str(hour) for hour in range(24)]

    def __init__(
        self,
        data_dir: str | Path,
        temporal_interp: bool = True,
        download_missing: bool = False,
        auth: CredentialManager | None = None,
    ):
        self.data_dir = Path(data_dir)
        self.temporal_interp = temporal_interp
        self.download_missing = download_missing
        self._auth = auth

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
        direct = self._load_from_explicit_path(self.data_dir)
        if direct is not None:
            return direct

        if not self.data_dir.exists() or not self.data_dir.is_dir():
            logger.warning(f"CAMS path does not exist or is not a directory: {self.data_dir}")
            return None

        date_str = obs_time.strftime("%Y%m%d")
        iso_date = obs_time.strftime("%Y-%m-%d")

        # Try different file patterns
        patterns = [
            f"cams_nrt_{date_str}*.nc",
            f"cams_{date_str}*.nc",
            f"*{date_str}*.nc",
            f"CAMS_{iso_date}*.nc",
            f"cams_nrt_{date_str}*.nc4",
            f"cams_{date_str}*.nc4",
            f"*{date_str}*.nc4",
            f"CAMS_{iso_date}*.nc4",
        ]

        for pattern in patterns:
            files = list(self.data_dir.glob(pattern))
            if files:
                try:
                    return xr.open_mfdataset(files, combine="by_coords")
                except Exception as e:
                    logger.warning(f"Failed to load CAMS: {e}")

        tif_dataset = self._load_cams_tif_group(date_str, iso_date)
        if tif_dataset is not None:
            return tif_dataset

        if self.download_missing:
            downloaded_file = self._download_cams_file(obs_time)
            if downloaded_file is not None:
                return self._load_from_explicit_path(downloaded_file)

        return None

    def _extract_variable(
        self, data: xr.Dataset, var_name: str, bounds: tuple, crs: str,
        resolution: float, obs_time: datetime
    ) -> xr.DataArray:
        """Extract and regrid a single variable."""
        if var_name not in data:
            logger.warning(f"Variable {var_name} not in CAMS data")
            defaults = {
                "aod550": 0.15,
                "tcwv": 1.5,
                # Native CAMS ozone unit is kg/m2; this value converts to 0.3 atm-cm.
                "gtco3": 0.0064245,
            }
            return self._create_default_array(bounds, crs, resolution, defaults.get(var_name, 0.0))

        var = data[var_name]
        if "band" in var.dims and var.sizes["band"] == 1:
            var = var.squeeze("band", drop=True)

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

    def _create_default_array(self, bounds: tuple, _crs: str, resolution: float, value: float) -> xr.DataArray:
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

    def _load_from_explicit_path(self, path: Path) -> xr.Dataset | None:
        """Load CAMS data when the configured path points directly to a file."""
        if not path.exists() or not path.is_file():
            return None

        suffix = path.suffix.lower()

        if suffix in self._NETCDF_SUFFIXES:
            try:
                return xr.open_dataset(path)
            except Exception as e:
                logger.warning(f"Failed to open CAMS NetCDF file {path}: {e}")
                return None

        if suffix in self._TIFF_SUFFIXES:
            return self._load_tif_dataset(path)

        logger.warning(f"Unsupported CAMS file extension: {path}")
        return None

    def _load_cams_tif_group(self, date_str: str, iso_date: str) -> xr.Dataset | None:
        """Load CAMS variables from GeoTIFF files in a directory."""
        patterns = [
            f"*{date_str}*.tif*",
            f"*{iso_date}*.tif*",
            "*.tif",
            "*.tiff",
        ]

        for pattern in patterns:
            files = sorted(
                p for p in self.data_dir.glob(pattern)
                if p.is_file() and p.suffix.lower() in self._TIFF_SUFFIXES
            )
            if not files:
                continue
            dataset = self._merge_tif_files(files)
            if dataset is not None:
                return dataset
        return None

    def _merge_tif_files(self, files: list[Path]) -> xr.Dataset | None:
        """Merge a set of GeoTIFF files into one CAMS dataset."""
        merged: xr.Dataset | None = None

        for path in files:
            ds = self._load_tif_dataset(path)
            if ds is None:
                continue
            merged = ds if merged is None else xr.merge([merged, ds], compat="override")

        return merged

    def _load_tif_dataset(self, path: Path) -> xr.Dataset | None:
        """Open a CAMS GeoTIFF and map bands/files to CAMS variable names."""
        try:
            da = xr.open_dataarray(path, engine="rasterio")
        except Exception as e:
            logger.warning(f"Failed to open CAMS GeoTIFF {path}: {e}")
            return None

        if "band" in da.dims and da.sizes["band"] > 1:
            return self._dataset_from_multiband_tif(da)

        variable = self._infer_variable_name(path.name)
        if variable is None:
            logger.warning(f"Could not infer CAMS variable from GeoTIFF name: {path.name}")
            return None

        if "band" in da.dims and da.sizes["band"] == 1:
            da = da.squeeze("band", drop=True)
        return xr.Dataset({variable: da})

    def _dataset_from_multiband_tif(self, da: xr.DataArray) -> xr.Dataset:
        """Map GeoTIFF bands to CAMS variables using labels then fallback ordering."""
        band_count = int(da.sizes["band"])
        band_labels = self._extract_band_labels(da)
        data_vars: dict[str, xr.DataArray] = {}

        for idx in range(1, band_count + 1):
            variable = None
            if idx - 1 < len(band_labels):
                variable = self._normalize_variable_name(band_labels[idx - 1])

            if variable is None and idx - 1 < len(self._REQUIRED_VARIABLES):
                variable = self._REQUIRED_VARIABLES[idx - 1]

            if variable is None or variable in data_vars:
                continue

            band = da.sel(band=idx, drop=True)
            data_vars[variable] = band

        return xr.Dataset(data_vars)

    def _extract_band_labels(self, da: xr.DataArray) -> list[str]:
        """Get band labels, if present, from raster metadata."""
        long_name = da.attrs.get("long_name")
        if isinstance(long_name, (list, tuple)):
            return [str(item) for item in long_name]
        if isinstance(long_name, str):
            return [long_name]
        return []

    def _infer_variable_name(self, file_name: str) -> str | None:
        """Infer CAMS variable from filename fragments."""
        return self._normalize_variable_name(file_name)

    def _normalize_variable_name(self, name: str) -> str | None:
        """Normalize a filename/label token into one of required CAMS variables."""
        lowered = name.lower()
        for key, variable in self._TIFF_VARIABLE_ALIASES.items():
            if key in lowered:
                return variable
        return None

    def _download_cams_file(self, obs_time: datetime) -> Path | None:
        """Download CAMS data using the official CDS API request pattern."""
        try:
            import cdsapi  # type: ignore[import-not-found]
        except ImportError:
            logger.warning("cdsapi is not installed; cannot auto-download CAMS data")
            return None

        output_dir = self.data_dir if self.data_dir.is_dir() else self.data_dir.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"CAMS_{obs_time:%Y-%m-%d}.nc"
        date_range = obs_time.strftime("%Y-%m-%d/%Y-%m-%d")

        request = {
            "variable": self._CDS_VARIABLES,
            "date": [date_range],
            "time": "00:00",
            "type": ["forecast"],
            "data_format": "netcdf",
            "area": [90, -180, -90, 180],
            "leadtime_hour": self._CDS_LEAD_TIMES,
        }

        client_kwargs: dict[str, str] = {}
        if self._auth is not None and self._auth.has_credentials("cds"):
            cds_cred = self._auth.get_credentials("cds")
            if cds_cred.key:
                client_kwargs["key"] = cds_cred.key

        try:
            cdsapi.Client(**client_kwargs).retrieve(self._CDS_DATASET, request).download(str(output_path))
        except Exception as e:
            logger.warning(f"Failed to download CAMS data for {obs_time:%Y-%m-%d}: {e}")
            return None
        return output_path

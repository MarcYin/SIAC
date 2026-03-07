"""
CAMS atmospheric prior provider.

Fetches atmospheric parameters (AOT, TCWV, TCO3) from ECMWF CAMS
(Copernicus Atmosphere Monitoring Service) near real-time data.
"""

from __future__ import annotations

import logging
import re
from contextlib import contextmanager
from datetime import datetime, timedelta
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import urlparse

import numpy as np
import xarray as xr

from siac.core.auth import CredentialManager
from siac.core.types import AtmosphericState

if TYPE_CHECKING:
    from collections.abc import Callable

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
    _DEFAULT_REMOTE_CACHE_DIR = Path.home() / ".cache" / "siac" / "cams"
    _JASMIN_CAMS_BASE_URL = "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"
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
    _CDSE_CAMS_FILE_PATTERN = re.compile(
        r"z_cams_c_ecmf_(?P<run>\d{14})_prod_(?P<mode>an|fc)_sfc_(?P<lead>\d{3})_(?P<var>[a-z0-9_]+)$"
    )
    _CDSE_CAMS_VARIABLE_SUFFIXES = {
        "aod550": "aod550",
        "gtco3": "gtco3",
        "tcwv": "tcwv",
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
        cache_dir: str | Path | None = None,
    ):
        self.data_dir = self._normalize_data_source(data_dir)
        self.temporal_interp = temporal_interp
        self.download_missing = download_missing
        self._auth = auth
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None

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
            return self._complete_cams_dataset(direct, obs_time)

        if self._is_remote_source(self.data_dir):
            remote = self._load_from_remote_base(str(self.data_dir), obs_time)
            if remote is not None:
                return self._complete_cams_dataset(remote, obs_time)
            return self._complete_cams_dataset(None, obs_time)

        data_dir = self._require_local_data_dir()

        if not data_dir.exists():
            if self.download_missing:
                if data_dir.suffix:
                    data_dir.parent.mkdir(parents=True, exist_ok=True)
                else:
                    data_dir.mkdir(parents=True, exist_ok=True)
            else:
                logger.warning(f"CAMS path does not exist: {data_dir}")
                return None

        if data_dir.exists() and not data_dir.is_dir():
            logger.warning(f"CAMS path is not a directory: {data_dir}")
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
            files = list(data_dir.glob(pattern))
            if files:
                try:
                    return self._complete_cams_dataset(
                        xr.open_mfdataset(files, combine="by_coords"),
                        obs_time,
                    )
                except Exception as e:
                    logger.warning(f"Failed to load CAMS: {e}")

        tif_dataset = self._load_cams_tif_group(date_str, iso_date)
        if tif_dataset is not None:
            return self._complete_cams_dataset(tif_dataset, obs_time)

        return self._complete_cams_dataset(None, obs_time)

    @classmethod
    def _missing_required_variables(cls, dataset: xr.Dataset | None) -> list[str]:
        if dataset is None:
            return list(cls._REQUIRED_VARIABLES)
        return [name for name in cls._REQUIRED_VARIABLES if name not in dataset.data_vars]

    @staticmethod
    def _reference_variable(dataset: xr.Dataset | None) -> xr.DataArray | None:
        if dataset is None:
            return None
        for name in CAMSProvider._REQUIRED_VARIABLES:
            if name in dataset.data_vars:
                return dataset[name]
        for data_var in dataset.data_vars.values():
            return data_var
        return None

    @staticmethod
    def _align_missing_variable(
        reference: xr.DataArray | None,
        variable: xr.DataArray,
    ) -> xr.DataArray:
        if reference is None:
            return variable

        aligned = variable
        for dim in ("latitude", "longitude", "time"):
            if (
                dim in aligned.dims
                and dim in reference.dims
                and dim in aligned.coords
                and dim in reference.coords
                and not aligned.coords[dim].identical(reference.coords[dim])
            ):
                try:
                    aligned = aligned.interp({dim: reference.coords[dim]}, method="linear")
                except Exception:
                    aligned = aligned.interp({dim: reference.coords[dim]}, method="nearest")
        return aligned

    @classmethod
    def _merge_missing_variables(
        cls,
        primary: xr.Dataset | None,
        fallback: xr.Dataset,
    ) -> tuple[xr.Dataset, list[str]]:
        if primary is None:
            return fallback, list(fallback.data_vars)

        added = [name for name in fallback.data_vars if name not in primary.data_vars]
        if not added:
            return primary, []
        reference = cls._reference_variable(primary)
        aligned = {
            name: cls._align_missing_variable(reference, fallback[name])
            for name in added
        }
        return xr.merge([primary, xr.Dataset(aligned)], compat="override", join="exact"), added

    def _load_cds_dataset(self, obs_time: datetime) -> xr.Dataset | None:
        if not self.download_missing:
            return None
        downloaded_file = self._download_cams_file(obs_time)
        if downloaded_file is None:
            return None
        return self._load_from_local_explicit_path(downloaded_file, source_name=downloaded_file.name)

    def _uses_jasmin_source(self) -> bool:
        if isinstance(self.data_dir, Path):
            return False
        return str(self.data_dir).rstrip("/").startswith(self._JASMIN_CAMS_BASE_URL.rstrip("/"))

    def _complete_cams_dataset(self, dataset: xr.Dataset | None, obs_time: datetime) -> xr.Dataset | None:
        missing = self._missing_required_variables(dataset)
        if not missing:
            return dataset

        candidates: list[tuple[str, Callable[[], xr.Dataset | None]]] = []
        if not self._uses_jasmin_source():
            candidates.append(
                (
                    "JASMIN CAMS mirror",
                    lambda: self._load_from_remote_base(self._JASMIN_CAMS_BASE_URL, obs_time),
                )
            )
        candidates.append(("CDS download", lambda: self._load_cds_dataset(obs_time)))

        merged = dataset
        for source_name, loader in candidates:
            if not missing:
                break
            fallback = loader()
            if fallback is None:
                continue

            fallback_vars = [name for name in missing if name in fallback.data_vars]
            if not fallback_vars:
                continue

            merged, added = self._merge_missing_variables(merged, fallback[fallback_vars])
            logger.info(
                "Supplemented CAMS variables from %s: %s",
                source_name,
                ", ".join(added),
            )
            missing = self._missing_required_variables(merged)

        return merged

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
        var = self._standardize_temporal_dims(var, obs_time)

        # Temporal interpolation
        if "time" in var.dims and self.temporal_interp:
            var = var.interp(time=np.datetime64(obs_time), method="linear")
        elif "time" in var.dims:
            var = var.sel(time=obs_time, method="nearest")

        # Spatial subset and regrid (simplified)
        xmin, ymin, xmax, ymax = bounds
        if "latitude" in var.dims:
            if crs and str(crs).upper() != "EPSG:4326":
                from siac.io.reprojection import transform_bounds

                xmin, ymin, xmax, ymax = transform_bounds(bounds, crs, "EPSG:4326")

            xmin, xmax = sorted((float(xmin), float(xmax)))
            ymin, ymax = sorted((float(ymin), float(ymax)))

            latitude_vals = var.coords["latitude"].values
            lat_descending = bool(latitude_vals[0] > latitude_vals[-1])
            lat_slice = slice(ymax, ymin) if lat_descending else slice(ymin, ymax)

            longitude_vals = var.coords["longitude"].values
            if float(np.nanmin(longitude_vals)) >= 0.0 and xmax <= 180.0:
                xmin = xmin % 360.0
                xmax = xmax % 360.0
            lon_descending = bool(longitude_vals[0] > longitude_vals[-1])
            lon_slice = slice(xmax, xmin) if lon_descending else slice(xmin, xmax)

            var = var.sel(latitude=lat_slice, longitude=lon_slice)

        return var

    def _standardize_temporal_dims(self, var: xr.DataArray, obs_time: datetime) -> xr.DataArray:
        """Normalize forecast-style CAMS time axes to a simple ``time`` axis."""
        if "forecast_reference_time" in var.dims:
            if var.sizes["forecast_reference_time"] == 1:
                var = var.squeeze("forecast_reference_time", drop=True)
            else:
                var = var.sel(forecast_reference_time=np.datetime64(obs_time), method="nearest")

        if "forecast_period" not in var.dims:
            return var

        valid_time = var.coords.get("valid_time")
        if valid_time is None:
            return var

        if "forecast_reference_time" in valid_time.dims:
            if valid_time.sizes.get("forecast_reference_time", 0) == 1:
                valid_time = valid_time.squeeze("forecast_reference_time", drop=True)
            else:
                valid_time = valid_time.sel(
                    forecast_reference_time=np.datetime64(obs_time),
                    method="nearest",
                )

        if valid_time.ndim != 1 or valid_time.dims != ("forecast_period",):
            return var

        return (
            var.assign_coords(time=("forecast_period", valid_time.values))
            .swap_dims({"forecast_period": "time"})
            .drop_vars("forecast_period", errors="ignore")
        )

    def _create_default_array(self, bounds: tuple, _crs: str, resolution: float, value: float) -> xr.DataArray:
        """Create default array with constant value."""
        xmin, ymin, xmax, ymax = bounds
        xmin, xmax = sorted((float(xmin), float(xmax)))
        ymin, ymax = sorted((float(ymin), float(ymax)))
        nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
        ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
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

    @classmethod
    def _normalize_data_source(cls, value: str | Path) -> str | Path:
        """Keep remote URLs as strings and normalize local paths."""
        if isinstance(value, Path):
            return value.expanduser()

        source = str(value).strip()
        if cls._is_remote_source(source):
            return source
        return Path(source).expanduser()

    @staticmethod
    def _is_remote_source(value: str | Path) -> bool:
        if isinstance(value, Path):
            return False
        scheme = urlparse(str(value)).scheme.lower()
        return scheme in {"http", "https", "s3"}

    @staticmethod
    def _remote_scheme(value: str | Path) -> str:
        if isinstance(value, Path):
            return ""
        return urlparse(str(value)).scheme.lower()

    @staticmethod
    def _source_suffix(value: str | Path) -> str:
        if isinstance(value, Path):
            return value.suffix.lower()
        return Path(urlparse(str(value)).path).suffix.lower()

    def _require_local_data_dir(self) -> Path:
        if isinstance(self.data_dir, Path):
            return self.data_dir
        raise TypeError("Local CAMS path requested for a remote CAMS source.")

    def _remote_cache_root(self) -> Path:
        cache_root = self.cache_dir or self._DEFAULT_REMOTE_CACHE_DIR
        cache_root.mkdir(parents=True, exist_ok=True)
        return cache_root

    def _remote_candidate_urls(self, base_url: str, obs_time: datetime) -> list[str]:
        base = base_url.rstrip("/")
        date_compact = obs_time.strftime("%Y%m%d")
        iso_date = obs_time.strftime("%Y-%m-%d")
        return [
            f"{base}/{iso_date}.nc",
            f"{base}/{iso_date}.nc4",
            f"{base}/CAMS_{iso_date}.nc",
            f"{base}/CAMS_{iso_date}.nc4",
            f"{base}/cams_{date_compact}.nc",
            f"{base}/cams_{date_compact}.nc4",
            f"{base}/cams_nrt_{date_compact}.nc",
            f"{base}/cams_nrt_{date_compact}.nc4",
        ]

    def _load_from_explicit_path(self, path: str | Path) -> xr.Dataset | None:
        """Load CAMS data when the configured source points directly to a file."""
        if self._is_remote_source(path):
            if self._source_suffix(path) not in self._NETCDF_SUFFIXES | self._TIFF_SUFFIXES:
                return None
            return self._load_from_remote_url(str(path))
        return self._load_from_local_explicit_path(Path(path))

    def _load_from_local_explicit_path(
        self,
        path: Path,
        *,
        source_name: str | None = None,
    ) -> xr.Dataset | None:
        """Load CAMS data from a local NetCDF or GeoTIFF file."""
        if not path.exists() or not path.is_file():
            return None

        suffix = self._source_suffix(source_name or path)

        if suffix in self._NETCDF_SUFFIXES:
            try:
                return xr.open_dataset(path, decode_timedelta=True)
            except Exception as e:
                origin = source_name or str(path)
                logger.warning(f"Failed to open CAMS NetCDF file {origin}: {e}")
                return None

        if suffix in self._TIFF_SUFFIXES:
            return self._load_tif_dataset(path, source_name=source_name)

        logger.warning(f"Unsupported CAMS file extension: {source_name or path}")
        return None

    def _load_from_remote_base(self, base_url: str, obs_time: datetime) -> xr.Dataset | None:
        if self._remote_scheme(base_url) == "s3":
            return self._load_from_remote_s3_base(base_url, obs_time)
        for url in self._remote_candidate_urls(base_url, obs_time):
            dataset = self._load_from_remote_url(url, missing_ok=True)
            if dataset is not None:
                return dataset
        return None

    @contextmanager
    def _remote_storage_options_context(self, url: str):
        if self._remote_scheme(url) != "s3":
            yield {}
            return

        parsed = urlparse(url)
        if (
            parsed.netloc == "eodata"
            and self._auth is not None
            and self._auth.cdse().has_credentials()
        ):
            with self._auth.cdse().temporary_s3_credentials() as creds:
                yield creds.storage_options()
            return

        if self._auth is not None and self._auth.aws().has_credentials():
            yield self._auth.aws().storage_options()
            return

        yield {}

    def _cache_remote_file(self, url: str, storage_options: dict | None = None) -> Path:
        import fsspec

        if self._remote_scheme(url) == "s3":
            parsed = urlparse(url)
            remote_path = f"{parsed.netloc}{parsed.path}"
            local_path = self._remote_cache_root() / parsed.netloc / parsed.path.lstrip("/")
            local_path.parent.mkdir(parents=True, exist_ok=True)
            if not local_path.exists():
                fs = fsspec.filesystem("s3", **(storage_options or {}))
                fs.get(remote_path, str(local_path))
            return local_path

        local_path = fsspec.open_local(
            f"simplecache::{url}",
            simplecache={"cache_storage": str(self._remote_cache_root())},
        )
        return Path(local_path)

    def _load_from_remote_url(
        self,
        url: str,
        *,
        missing_ok: bool = False,
        storage_options: dict | None = None,
    ) -> xr.Dataset | None:
        suffix = self._source_suffix(url)
        if suffix not in self._NETCDF_SUFFIXES and suffix not in self._TIFF_SUFFIXES:
            if not missing_ok:
                logger.warning(f"Unsupported CAMS remote file extension: {url}")
            return None

        if storage_options is not None:
            try:
                local_path = self._cache_remote_file(url, storage_options=storage_options)
            except FileNotFoundError:
                if not missing_ok:
                    logger.warning(f"CAMS remote file not found: {url}")
                return None
            except Exception as exc:
                log = logger.debug if missing_ok else logger.warning
                log(f"Failed to cache remote CAMS file {url}: {exc}")
                return None
        else:
            with self._remote_storage_options_context(url) as resolved_storage_options:
                try:
                    local_path = self._cache_remote_file(url, storage_options=resolved_storage_options)
                except FileNotFoundError:
                    if not missing_ok:
                        logger.warning(f"CAMS remote file not found: {url}")
                    return None
                except Exception as exc:
                    log = logger.debug if missing_ok else logger.warning
                    log(f"Failed to cache remote CAMS file {url}: {exc}")
                    return None

        return self._load_from_local_explicit_path(
            local_path,
            source_name=Path(urlparse(url).path).name,
        )

    def _load_from_remote_s3_base(self, base_url: str, obs_time: datetime) -> xr.Dataset | None:
        if not self._looks_like_cdse_cams_base(base_url):
            return None

        with self._remote_storage_options_context(base_url) as storage_options:
            selected = self._select_cdse_cams_files(base_url, obs_time, storage_options)
            if not selected:
                return None

            datasets: list[xr.Dataset] = []
            for url in selected:
                dataset = self._load_from_remote_url(
                    url,
                    missing_ok=True,
                    storage_options=storage_options,
                )
                if dataset is not None:
                    datasets.append(dataset)

            if not datasets:
                return None
            return xr.merge(datasets, compat="override")

    @staticmethod
    def _looks_like_cdse_cams_base(base_url: str) -> bool:
        parsed = urlparse(base_url)
        return parsed.scheme.lower() == "s3" and parsed.netloc == "eodata" and "/CAMS/" in parsed.path

    def _select_cdse_cams_files(
        self,
        base_url: str,
        obs_time: datetime,
        storage_options: dict,
    ) -> list[str]:
        import fsspec

        parsed = urlparse(base_url)
        day_prefix = f"{parsed.netloc}{parsed.path.rstrip('/')}/{obs_time:%Y/%m/%d}"
        fs = fsspec.filesystem("s3", **storage_options)

        try:
            entries = fs.ls(day_prefix, detail=False)
        except FileNotFoundError:
            return []

        chosen: dict[str, tuple[timedelta, bool, str]] = {}
        for entry in entries:
            name = entry.rsplit("/", 1)[-1]
            match = self._CDSE_CAMS_FILE_PATTERN.fullmatch(name)
            if match is None:
                continue

            variable = match.group("var")
            if variable not in self._CDSE_CAMS_VARIABLE_SUFFIXES.values():
                continue

            run_time = datetime.strptime(match.group("run"), "%Y%m%d%H%M%S")
            lead_hours = int(match.group("lead"))
            valid_time = run_time + timedelta(hours=lead_hours)
            delta = abs(valid_time - obs_time)
            is_forecast = match.group("mode") == "fc"
            candidate = (delta, is_forecast, entry)
            existing = chosen.get(variable)
            if existing is None or candidate < existing:
                chosen[variable] = candidate

        urls: list[str] = []
        for variable in ("aod550", "gtco3", "tcwv"):
            candidate = chosen.get(variable)
            if candidate is None:
                continue
            entry = candidate[2]
            name = entry.rsplit("/", 1)[-1]
            urls.append(f"s3://{entry}/{name}.nc")
        return urls

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
                p for p in self._require_local_data_dir().glob(pattern)
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

    def _load_tif_dataset(self, path: Path, *, source_name: str | None = None) -> xr.Dataset | None:
        """Open a CAMS GeoTIFF and map bands/files to CAMS variable names."""
        try:
            da = xr.open_dataarray(path, engine="rasterio")
        except Exception as e:
            origin = source_name or str(path)
            logger.warning(f"Failed to open CAMS GeoTIFF {origin}: {e}")
            return None

        if "band" in da.dims and da.sizes["band"] > 1:
            return self._dataset_from_multiband_tif(da)

        variable = self._infer_variable_name(source_name or path.name)
        if variable is None:
            logger.warning(f"Could not infer CAMS variable from GeoTIFF name: {source_name or path.name}")
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

        if isinstance(self.data_dir, Path) and self.data_dir.is_dir():
            output_dir = self.data_dir
        elif isinstance(self.data_dir, Path):
            output_dir = self.data_dir.parent
        else:
            output_dir = self._remote_cache_root()
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

        cds_auth = self._auth.cds() if self._auth is not None else CredentialManager().cds()
        client_kwargs = cds_auth.client_kwargs()
        has_any_credentials = cds_auth.has_any_credentials()

        if not has_any_credentials:
            logger.warning("CDS credentials are not configured; cannot auto-download CAMS data")
            return None

        try:
            cdsapi.Client(**client_kwargs).retrieve(self._CDS_DATASET, request).download(str(output_path))
        except Exception as e:
            logger.warning(f"Failed to download CAMS data for {obs_time:%Y-%m-%d}: {e}")
            return None
        return output_path

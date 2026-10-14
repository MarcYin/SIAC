"""
Landsat 8/9 OLI Collection-2 Level-1 preprocessor.

Reads an unpacked USGS Collection-2 Level-1 scene (``*_MTL.txt`` or
``*_MTL.json`` plus the per-band GeoTIFFs) into the M1 contract:

- TOA reflectance for the reflective OLI bands ``B1``-``B7``, rescaled with the
  per-band ``REFLECTANCE_MULT/ADD`` coefficients from the MTL and corrected for
  the per-pixel solar zenith angle from the Collection-2 ``SZA`` angle band
  (scene-centre ``SUN_ELEVATION`` when the angle band is absent).
- Sun and view geometry from the Collection-2 ``SZA``/``SAA``/``VZA``/``VAA``
  angle bands, block-averaged to a coarse grid (the Sentinel-2 path keeps the
  5 km metadata grid; every downstream consumer resamples to its own grid).
- Standardized cloud classes, either from the configured cloud-mask provider or,
  with ``provider = "product_qa"``, from the product's own ``QA_PIXEL`` band.
"""

from __future__ import annotations

import json
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.adapters.rsrf import load_sensor_config_with_rsrf
from siac.adapters.satellite.base import (
    BaseSatellitePreprocessor,
    degrees_to_radians,
    register_preprocessor,
)
from siac.algorithms.cloud import build_cloud_classes, classes_to_bool_mask
from siac.geo import reproject_match
from siac.runtime import GeometryAngles
from siac.storage.readers import read_raster, read_raster_window

if TYPE_CHECKING:
    from siac.domain import SensorConfig

logger = logging.getLogger(__name__)
_TOA_BAND_LOADER_ATTR = "_siac_toa_band_loader"

_SATELLITE_BY_SPACECRAFT = {"LANDSAT_8": "L8", "LANDSAT_9": "L9"}

# Collection-2 angle bands are int16 in hundredths of a degree.
_ANGLE_SCALE = 0.01
_ANGLE_SUFFIXES = {
    "sza": ("FILE_NAME_ANGLE_SOLAR_ZENITH_BAND_4", "_SZA.TIF"),
    "saa": ("FILE_NAME_ANGLE_SOLAR_AZIMUTH_BAND_4", "_SAA.TIF"),
    "vza": ("FILE_NAME_ANGLE_SENSOR_ZENITH_BAND_4", "_VZA.TIF"),
    "vaa": ("FILE_NAME_ANGLE_SENSOR_AZIMUTH_BAND_4", "_VAA.TIF"),
}
_QA_PIXEL_KEY = ("FILE_NAME_QUALITY_L1_PIXEL", "_QA_PIXEL.TIF")

# Collection-2 QA_PIXEL bit positions.
_QA_FILL = 0
_QA_DILATED_CLOUD = 1
_QA_CIRRUS = 2
_QA_CLOUD = 3
_QA_CLOUD_SHADOW = 4

_DEFAULT_ANGLE_BLOCK_PIXELS = 100  # 3 km blocks at 30 m


def parse_mtl_text(text: str) -> dict[str, Any]:
    """Parse a Landsat ``*_MTL.txt`` ODL document into nested dictionaries.

    Quoted values stay strings; unquoted values are returned as the raw token,
    leaving numeric conversion to the caller (MTL mixes dates, times, integers
    and floats under the same syntax).
    """
    root: dict[str, Any] = {}
    stack: list[dict[str, Any]] = [root]
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line == "END":
            continue
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key == "GROUP":
            group: dict[str, Any] = {}
            stack[-1][value] = group
            stack.append(group)
        elif key == "END_GROUP":
            if len(stack) > 1:
                stack.pop()
        else:
            if len(value) >= 2 and value[0] == value[-1] == '"':
                value = value[1:-1]
            stack[-1][key] = value
    return root


def flatten_mtl(tree: dict[str, Any]) -> dict[str, str]:
    """Flatten nested MTL groups to a ``{KEY: value}`` map (first occurrence wins)."""
    flat: dict[str, str] = {}

    def _walk(node: dict[str, Any]) -> None:
        for key, value in node.items():
            if isinstance(value, dict):
                _walk(value)
            elif key not in flat:
                flat[key] = str(value)

    _walk(tree)
    return flat


def qa_pixel_to_cloud_classes(qa: np.ndarray) -> np.ndarray:
    """Map Collection-2 ``QA_PIXEL`` bits to SIAC classes {0,1,2,3}.

    Fill -> 0 (missing); cloud, dilated cloud or high-confidence cirrus -> 2;
    cloud shadow -> 3; everything else (clear, water, snow) -> 1.
    """
    bits = np.asarray(qa).astype(np.uint16)

    def _bit(position: int) -> np.ndarray:
        return cast("np.ndarray", ((bits >> np.uint16(position)) & np.uint16(1)).astype(bool))

    classes: np.ndarray = np.ones(bits.shape, dtype=np.uint8)
    classes[_bit(_QA_CLOUD_SHADOW)] = 3
    cloud = _bit(_QA_CLOUD) | _bit(_QA_DILATED_CLOUD) | _bit(_QA_CIRRUS)
    classes[cloud] = 2
    classes[_bit(_QA_FILL)] = 0
    return classes


@register_preprocessor("landsat")
@register_preprocessor("l9")
@register_preprocessor("l8")
class LandsatPreprocessor(BaseSatellitePreprocessor):
    """Preprocessor for Landsat 8/9 OLI Collection-2 Level-1 scenes."""

    def __init__(self, config: dict[str, Any] | None = None):
        super().__init__(config)
        self._satellite_id: str | None = None
        self._sensor_config_cache: SensorConfig | None = None
        self._sensor_config_cache_satellite_id: str | None = None
        self._scene_dir: Path | None = None
        self._mtl: dict[str, str] | None = None
        self._resolved_input_path: Path | None = None
        self._reference_grid: xr.DataArray | None = None
        self._cos_sza: xr.DataArray | None = None
        self._last_cloud_classes: xr.DataArray | None = None
        self._subset_bounds: tuple[float, float, float, float] | None = None
        self._subset_crs: str | None = None

    # ------------------------------------------------------------------
    # Spatial subset
    # ------------------------------------------------------------------

    def set_spatial_subset(
        self,
        bounds: tuple[float, float, float, float],
        *,
        crs: str,
    ) -> None:
        """Restrict raster reads to an AOI window."""
        xmin, ymin, xmax, ymax = bounds
        self._subset_bounds = (float(xmin), float(ymin), float(xmax), float(ymax))
        self._subset_crs = str(crs)
        self._reference_grid = None
        self._cos_sza = None

    def clear_spatial_subset(self) -> None:
        """Clear any active AOI read window."""
        self._subset_bounds = None
        self._subset_crs = None
        self._reference_grid = None
        self._cos_sza = None

    # ------------------------------------------------------------------
    # Sensor configuration
    # ------------------------------------------------------------------

    @property
    def sensor_config(self) -> SensorConfig:
        """Return the RSRF-backed OLI configuration for the resolved platform."""
        satellite_id = self._satellite_id or "L8"
        if (
            self._sensor_config_cache is not None
            and self._sensor_config_cache_satellite_id == satellite_id
        ):
            return self._sensor_config_cache
        try:
            sensor_config = load_sensor_config_with_rsrf(
                "OLI",
                satellite_id,
                rsrf_root=self.config.get("rsrf_root"),
            )
        except Exception as exc:
            raise RuntimeError(
                f"Unable to load Landsat RSRF metadata for {satellite_id}. Configure "
                f"paths.rsrf_root to a catalog containing {satellite_id} OLI band responses."
            ) from exc
        self._sensor_config_cache = sensor_config
        self._sensor_config_cache_satellite_id = satellite_id
        return sensor_config

    @property
    def band_names(self) -> tuple[str, ...]:
        return tuple(band.name for band in self.sensor_config.bands)

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------

    def get_metadata(self, input_path: str | Path) -> dict[str, Any]:
        """Extract observation metadata from the scene MTL."""
        input_path = Path(input_path).expanduser()
        mtl = self._resolve_scene(input_path)

        metadata: dict[str, Any] = {
            "sensor": "OLI",
            "satellite": self._satellite_id,
            "input_path": str(input_path),
            "product_id": mtl.get("LANDSAT_PRODUCT_ID"),
            "processing_level": mtl.get("PROCESSING_LEVEL"),
            "collection_number": mtl.get("COLLECTION_NUMBER"),
            "collection_category": mtl.get("COLLECTION_CATEGORY"),
        }
        path, row = mtl.get("WRS_PATH"), mtl.get("WRS_ROW")
        if path is not None and row is not None:
            metadata["wrs_path"] = int(path)
            metadata["wrs_row"] = int(row)
            metadata["tile_id"] = f"{int(path):03d}{int(row):03d}"
        for key, name in (
            ("SUN_AZIMUTH", "sun_azimuth"),
            ("SUN_ELEVATION", "sun_elevation"),
            ("EARTH_SUN_DISTANCE", "earth_sun_distance"),
            ("CLOUD_COVER", "cloud_cover"),
        ):
            if key in mtl:
                metadata[name] = float(mtl[key])

        observation_time = self._parse_observation_time(mtl)
        if observation_time is None:
            observation_time = self._parse_observation_time_from_name(
                mtl.get("LANDSAT_PRODUCT_ID") or input_path.name
            )
        if observation_time is not None:
            metadata["observation_time"] = observation_time

        metadata["reflectance_mult"] = {}
        metadata["reflectance_add"] = {}
        for band_name in self.band_names:
            number = self._band_number(band_name)
            mult = mtl.get(f"REFLECTANCE_MULT_BAND_{number}")
            add = mtl.get(f"REFLECTANCE_ADD_BAND_{number}")
            if mult is not None and add is not None:
                metadata["reflectance_mult"][band_name] = float(mult)
                metadata["reflectance_add"][band_name] = float(add)
        return metadata

    @staticmethod
    def _parse_observation_time(mtl: dict[str, str]) -> datetime | None:
        date_text = mtl.get("DATE_ACQUIRED")
        time_text = mtl.get("SCENE_CENTER_TIME")
        if not date_text:
            return None
        if not time_text:
            return datetime.strptime(date_text, "%Y-%m-%d")
        # SCENE_CENTER_TIME carries up to 7 fractional digits; datetime takes 6.
        match = re.match(r"(\d{2}):(\d{2}):(\d{2})(?:\.(\d+))?", time_text.strip().rstrip("Z"))
        if match is None:
            logger.warning("Could not parse Landsat SCENE_CENTER_TIME %r", time_text)
            return datetime.strptime(date_text, "%Y-%m-%d")
        hour, minute, second, fraction = match.groups()
        microsecond = int((fraction or "0")[:6].ljust(6, "0"))
        day = datetime.strptime(date_text, "%Y-%m-%d")
        return day.replace(
            hour=int(hour), minute=int(minute), second=int(second), microsecond=microsecond
        )

    @staticmethod
    def _parse_observation_time_from_name(name: str) -> datetime | None:
        """Infer the acquisition date from a Collection-2 product id."""
        match = re.search(r"_\d{6}_(\d{8})_\d{8}_", name)
        if match is None:
            return None
        try:
            return datetime.strptime(match.group(1), "%Y%m%d")
        except ValueError:
            return None

    # ------------------------------------------------------------------
    # TOA reflectance
    # ------------------------------------------------------------------

    def load_toa(
        self,
        input_path: str | Path,
        *,
        band_names: tuple[str, ...] | list[str] | None = None,
    ) -> xr.Dataset:
        """Load TOA reflectance for the reflective OLI bands."""
        input_path = Path(input_path).expanduser()
        metadata = self.get_metadata(input_path)
        requested = self._resolve_requested_band_names(band_names)

        reference_name = "B4" if "B4" in self.band_names else self.band_names[0]
        reference = self._load_band(reference_name, metadata=metadata, reference_grid=None)
        self._reference_grid = reference

        loaded: dict[str, xr.DataArray] = {}
        if reference_name in requested:
            loaded[reference_name] = reference
        for band_name in requested:
            if band_name not in loaded:
                loaded[band_name] = self._load_band(
                    band_name, metadata=metadata, reference_grid=reference
                )

        ds = xr.Dataset(loaded)
        ds.attrs[_TOA_BAND_LOADER_ATTR] = self._make_band_loader(
            metadata=metadata,
            reference_grid=reference,
            loaded_band_names=tuple(loaded),
        )
        ds.attrs["loaded_toa_bands"] = tuple(loaded)
        ds.attrs["sensor"] = "OLI"
        ds.attrs["satellite"] = self._satellite_id
        obs_time = metadata.get("observation_time")
        ds.attrs["observation_time"] = obs_time.isoformat() if obs_time else ""
        return ds

    def _resolve_requested_band_names(
        self, band_names: tuple[str, ...] | list[str] | None
    ) -> tuple[str, ...]:
        available = self.band_names
        requested = tuple(str(name) for name in band_names) if band_names else available
        resolved: list[str] = []
        for name in requested:
            if name not in available:
                logger.warning("Ignoring unsupported Landsat band request %s", name)
                continue
            if self._band_file(name) is None:
                logger.warning("Band %s not found", name)
                continue
            if name not in resolved:
                resolved.append(name)
        if not resolved:
            raise RuntimeError("No requested Landsat bands were available for TOA loading.")
        return tuple(resolved)

    def _load_band(
        self,
        band_name: str,
        *,
        metadata: dict[str, Any],
        reference_grid: xr.DataArray | None,
    ) -> xr.DataArray:
        band_file = self._band_file(band_name)
        if band_file is None:
            raise KeyError(f"Band {band_name!r} is not available in the Landsat scene.")
        mult = metadata["reflectance_mult"].get(band_name)
        add = metadata["reflectance_add"].get(band_name)
        if mult is None or add is None:
            raise ValueError(
                f"Landsat MTL lacks REFLECTANCE_MULT/ADD for {band_name}; "
                "Collection-2 Level-1 metadata is required."
            )
        dn = self._read(band_file)
        if reference_grid is not None and not self._coords_match(dn, reference_grid):
            dn = reproject_match(dn, reference_grid, resampling="nearest")
        reflectance = self._calibrate(dn, mult=mult, add=add, metadata=metadata)
        reflectance.name = band_name
        return reflectance

    def _calibrate(
        self,
        dn: xr.DataArray,
        *,
        mult: float,
        add: float,
        metadata: dict[str, Any],
    ) -> xr.DataArray:
        """DN -> TOA reflectance: (M*Q + A) / cos(SZA); DN 0 is fill."""
        values = dn.astype(np.float32)
        fill = ~np.isfinite(values) | (values == 0)
        rho_prime = values * np.float32(mult) + np.float32(add)
        cos_sza = self._cos_solar_zenith(dn, metadata=metadata)
        rho = (rho_prime / cos_sza).clip(np.float32(0.0), np.float32(1.5))
        return cast("xr.DataArray", rho.where(~fill).astype(np.float32))

    def _cos_solar_zenith(
        self, template: xr.DataArray, *, metadata: dict[str, Any]
    ) -> xr.DataArray | np.float32:
        """Per-pixel cos(SZA) from the angle band, else scene-centre sun elevation."""
        cached = self._cos_sza
        if cached is not None and self._coords_match(cached, template):
            return cached
        sza_file = self._auxiliary_file(*_ANGLE_SUFFIXES["sza"])
        if sza_file is not None:
            sza = self._read(sza_file).astype(np.float32) * np.float32(_ANGLE_SCALE)
            if not self._coords_match(sza, template):
                sza = reproject_match(sza, template, resampling="bilinear")
            cos_sza = np.cos(np.deg2rad(sza)).astype(np.float32)
            cos_sza = cast("xr.DataArray", cos_sza.where(cos_sza > 0))
            self._cos_sza = cos_sza
            return cos_sza
        elevation = metadata.get("sun_elevation")
        if elevation is None:
            raise ValueError("Landsat scene has neither an SZA angle band nor SUN_ELEVATION.")
        logger.warning(
            "Landsat SZA angle band not found; using scene-centre SUN_ELEVATION for the "
            "TOA solar-zenith correction."
        )
        return np.float32(np.sin(np.deg2rad(float(elevation))).item())

    def _make_band_loader(
        self,
        *,
        metadata: dict[str, Any],
        reference_grid: xr.DataArray,
        loaded_band_names: tuple[str, ...],
    ) -> Any:
        cache: dict[tuple[str, bool], xr.DataArray] = {}
        loaded = set(loaded_band_names)

        def _load(band_name: str, *, native: bool = False) -> xr.DataArray:
            key = (band_name, native)
            if key in cache:
                return cache[key]
            if band_name in loaded and not native:
                raise KeyError(f"Band {band_name!r} is already present in the loaded TOA dataset.")
            band = self._load_band(
                band_name,
                metadata=metadata,
                reference_grid=None if native else reference_grid,
            )
            cache[key] = band
            return band

        return _load

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def extract_geometry(self, input_path: str | Path) -> GeometryAngles:
        """Block-average the Collection-2 angle bands to a coarse geometry grid."""
        input_path = Path(input_path).expanduser()
        self._resolve_scene(input_path)
        block = int(self.config.get("angle_block_pixels", _DEFAULT_ANGLE_BLOCK_PIXELS))
        if block < 1:
            raise ValueError(f"angle_block_pixels must be >= 1, got {block}")

        grids: dict[str, xr.DataArray] = {}
        for name, (mtl_key, suffix) in _ANGLE_SUFFIXES.items():
            path = self._auxiliary_file(mtl_key, suffix)
            if path is None:
                raise FileNotFoundError(
                    f"Landsat angle band {suffix} not found for {input_path}. "
                    "Collection-2 Level-1 scenes ship SZA/SAA/VZA/VAA bands."
                )
            raw = self._read(path).astype(np.float32)
            degrees = raw * np.float32(_ANGLE_SCALE)
            # Angle bands carry 0 outside the imaged footprint; the zero test is
            # safe for SZA but not for VZA/VAA, so footprint comes from QA/TOA.
            degrees = degrees.where(self._footprint(raw))
            grids[name] = self._block_reduce(degrees, block, circular=name in ("saa", "vaa"))

        return GeometryAngles(
            sza=degrees_to_radians(grids["sza"]),
            saa=degrees_to_radians(grids["saa"]),
            vza=degrees_to_radians(grids["vza"]),
            vaa=degrees_to_radians(grids["vaa"]),
        )

    def _footprint(self, template: xr.DataArray) -> xr.DataArray:
        qa_file = self._auxiliary_file(*_QA_PIXEL_KEY)
        if qa_file is not None:
            qa = self._read(qa_file)
            if not self._coords_match(qa, template):
                qa = reproject_match(qa, template, resampling="nearest")
            valid = np.isfinite(qa) & (qa_pixel_to_cloud_classes(np.nan_to_num(qa.values)) != 0)
            return xr.DataArray(np.asarray(valid), dims=template.dims, coords=template.coords)
        if self._reference_grid is not None and self._coords_match(self._reference_grid, template):
            return cast("xr.DataArray", np.isfinite(self._reference_grid))
        return cast("xr.DataArray", np.isfinite(template) & (template != 0))

    @staticmethod
    def _block_reduce(da: xr.DataArray, block: int, *, circular: bool) -> xr.DataArray:
        """NaN-aware block mean (circular for azimuths), then fill empty blocks."""
        from scipy.ndimage import distance_transform_edt

        values = np.asarray(da.values, dtype=np.float64)
        height, width = values.shape
        rows, cols = -(-height // block), -(-width // block)
        padded = np.full((rows * block, cols * block), np.nan)
        padded[:height, :width] = values
        tiles = padded.reshape(rows, block, cols, block).transpose(0, 2, 1, 3)
        tiles = tiles.reshape(rows, cols, block * block)
        finite = np.isfinite(tiles)
        count = finite.sum(axis=-1)
        if circular:
            radians = np.deg2rad(np.where(finite, tiles, 0.0))
            sine = np.where(finite, np.sin(radians), 0.0).sum(axis=-1)
            cosine = np.where(finite, np.cos(radians), 0.0).sum(axis=-1)
            reduced = np.mod(np.rad2deg(np.arctan2(sine, cosine)), 360.0)
        else:
            total = np.where(finite, tiles, 0.0).sum(axis=-1)
            reduced = np.divide(total, np.maximum(count, 1))
        reduced = np.where(count > 0, reduced, np.nan)
        if not np.isfinite(reduced).any():
            raise ValueError("Landsat angle band has no valid pixels inside the scene footprint.")
        empty = ~np.isfinite(reduced)
        if empty.any():
            _, (iy, ix) = distance_transform_edt(empty, return_indices=True)
            reduced = reduced[iy, ix]
        if circular:
            unwrapped = np.unwrap(np.unwrap(np.deg2rad(reduced), axis=1), axis=0)
            reduced = np.rad2deg(unwrapped)

        x = np.asarray(da.coords["x"].values, dtype=np.float64)
        y = np.asarray(da.coords["y"].values, dtype=np.float64)
        dx = float(x[1] - x[0]) if x.size > 1 else 30.0
        dy = float(y[1] - y[0]) if y.size > 1 else -30.0
        offset = (block - 1) / 2.0
        x_coarse = x[0] + (np.arange(cols) * block + offset) * dx
        y_coarse = y[0] + (np.arange(rows) * block + offset) * dy
        out = xr.DataArray(
            reduced.astype(np.float32),
            dims=("y", "x"),
            coords={"y": y_coarse, "x": x_coarse},
        )
        crs = da.rio.crs
        if crs is not None:
            out = out.rio.write_crs(crs)
        return out

    # ------------------------------------------------------------------
    # Cloud mask
    # ------------------------------------------------------------------

    def extract_cloud_mask(
        self,
        input_path: str | Path,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray:
        """Build standardized cloud classes and return the boolean mask."""
        input_path = Path(input_path).expanduser()
        if toa is None:
            toa = self.load_toa(input_path)
        settings = self._cloud_mask_settings(input_path)

        if settings["mode"] == "auto" and settings["provider"] == "product_qa":
            cloud_classes = self._product_qa_classes(toa)
        else:
            cloud_classes = build_cloud_classes(
                toa,
                self.sensor_config,
                mode=settings["mode"],
                provider=settings["provider"],
                class_mapping=settings["class_mapping"],
                external_mask_path=settings["external_mask_path"],
                user_callable=settings["user_callable"],
                target_resolution_m=settings["target_resolution_m"],
                resolution_policy=settings["resolution_policy"],
                allow_upsample_to_target=settings["allow_upsample_to_target"],
                unmapped_to_missing=settings["unmapped_to_missing"],
                cloud_cache_dir=settings.get("cache_dir"),
                inference_device=settings.get("inference_device"),
            )
        self._last_cloud_classes = cloud_classes
        return classes_to_bool_mask(cloud_classes)

    def _product_qa_classes(self, toa: xr.Dataset) -> xr.DataArray:
        qa_file = self._auxiliary_file(*_QA_PIXEL_KEY)
        if qa_file is None:
            raise FileNotFoundError(
                "cloud_mask.provider='product_qa' needs the Collection-2 QA_PIXEL band."
            )
        reference = next(iter(toa.data_vars.values()))
        qa = self._read(qa_file)
        if not self._coords_match(qa, reference):
            qa = reproject_match(qa, reference, resampling="nearest")
        classes = qa_pixel_to_cloud_classes(np.nan_to_num(np.asarray(qa.values), nan=1.0))
        classes[~np.isfinite(np.asarray(reference.values))] = 0
        out = xr.DataArray(
            classes, dims=reference.dims, coords=reference.coords, name="cloud_classes"
        )
        out.attrs["class_values"] = {0: "missing", 1: "clear", 2: "cloud", 3: "shadow"}
        out.attrs["source"] = "landsat_c2_qa_pixel"
        return out

    def _cloud_mask_settings(self, input_path: Path) -> dict[str, Any]:
        defaults: dict[str, Any] = {
            "mode": "auto",
            "provider": "omnicloudmask",
            "class_mapping": None,
            "external_mask_path": None,
            "user_callable": None,
            "target_resolution_m": 30.0,
            "resolution_policy": "auto",
            "allow_upsample_to_target": False,
            "unmapped_to_missing": True,
            "inference_device": "cpu",
        }
        cloud_mask_config = (
            self.config.get("cloud_mask", {}) if isinstance(self.config, dict) else {}
        )
        if isinstance(cloud_mask_config, dict):
            defaults.update(cloud_mask_config)
        for key in ("mode", "provider", "resolution_policy"):
            defaults[key] = str(getattr(defaults[key], "value", defaults[key]))
        external = defaults.get("external_mask_path")
        if isinstance(external, (str, Path)) and str(external):
            p = Path(external).expanduser()
            if not p.is_absolute():
                p = (self._scene_dir or input_path) / p
            defaults["external_mask_path"] = p
        return defaults

    # ------------------------------------------------------------------
    # Scene resolution and I/O
    # ------------------------------------------------------------------

    def _resolve_scene(self, input_path: Path) -> dict[str, str]:
        if self._resolved_input_path == input_path and self._mtl is not None:
            return self._mtl
        if input_path.is_file():
            mtl_path = input_path
            scene_dir = input_path.parent
        else:
            scene_dir = input_path
            candidates = sorted(scene_dir.glob("*_MTL.txt")) or sorted(scene_dir.glob("*_MTL.json"))
            if not candidates:
                raise FileNotFoundError(f"No Landsat *_MTL.txt or *_MTL.json under {scene_dir}")
            if len(candidates) > 1:
                raise ValueError(
                    f"Multiple Landsat MTL files under {scene_dir}: "
                    + ", ".join(path.name for path in candidates)
                )
            mtl_path = candidates[0]

        text = mtl_path.read_text(encoding="utf-8")
        tree = json.loads(text) if mtl_path.suffix.lower() == ".json" else parse_mtl_text(text)
        mtl = flatten_mtl(tree)

        spacecraft = mtl.get("SPACECRAFT_ID", "").upper()
        satellite_id = _SATELLITE_BY_SPACECRAFT.get(spacecraft)
        if satellite_id is None:
            raise ValueError(
                f"Unsupported Landsat spacecraft {spacecraft!r}; only Landsat 8/9 OLI "
                "Collection-2 Level-1 scenes are supported."
            )
        level = mtl.get("PROCESSING_LEVEL", "")
        if level and not level.startswith("L1"):
            raise ValueError(
                f"Landsat product processing level is {level!r}; SIAC expects a Level-1 "
                "(L1TP/L1GT/L1GS) scene, not surface reflectance."
            )

        self._satellite_id = satellite_id
        self._scene_dir = scene_dir
        self._mtl = mtl
        self._resolved_input_path = input_path
        self._reference_grid = None
        self._cos_sza = None
        return mtl

    def _band_file(self, band_name: str) -> Path | None:
        number = self._band_number(band_name)
        return self._auxiliary_file(f"FILE_NAME_BAND_{number}", f"_B{number}.TIF")

    def _auxiliary_file(self, mtl_key: str, suffix: str) -> Path | None:
        if self._scene_dir is None or self._mtl is None:
            raise RuntimeError("Landsat scene is not resolved; call get_metadata() first.")
        named = self._mtl.get(mtl_key)
        if named:
            candidate = self._scene_dir / named
            if candidate.exists():
                return candidate
        matches = [
            path
            for path in sorted(self._scene_dir.iterdir())
            if path.name.upper().endswith(suffix.upper())
        ]
        return matches[0] if matches else None

    @staticmethod
    def _band_number(band_name: str) -> int:
        match = re.fullmatch(r"B(\d{1,2})", band_name)
        if match is None:
            raise ValueError(f"Not a Landsat OLI band name: {band_name!r}")
        return int(match.group(1))

    def _read(self, path: Path) -> xr.DataArray:
        if self._subset_bounds is None or self._subset_crs is None:
            da = read_raster(path, masked=False)
        else:
            da = read_raster_window(
                path,
                bounds=self._subset_bounds,
                bounds_crs=self._subset_crs,
                masked=False,
            )
        if "band" in da.dims:
            da = da.squeeze("band", drop=True)
        return da

    @staticmethod
    def _coords_match(a: xr.DataArray, b: xr.DataArray) -> bool:
        if a.shape != b.shape:
            return False
        for dim in ("x", "y"):
            if dim not in a.coords or dim not in b.coords:
                return False
            if not np.allclose(a.coords[dim].values, b.coords[dim].values):
                return False
        return True

"""
Sentinel-2 MSI preprocessor.

This module handles reading and preprocessing Sentinel-2 Level-1C data
in SAFE format (ESA SciHub or AWS formats).
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
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
from siac.catalog import SENTINEL2A_CONFIG, get_sensor_config
from siac.geo import reproject_match
from siac.runtime import GeometryAngles
from siac.storage.readers import read_raster, read_raster_window

if TYPE_CHECKING:
    from siac.domain import SensorConfig

logger = logging.getLogger(__name__)
_TOA_BAND_LOADER_ATTR = "_siac_toa_band_loader"


@register_preprocessor("s2")
class Sentinel2Preprocessor(BaseSatellitePreprocessor):
    """
    Preprocessor for Sentinel-2 MSI Level-1C data.

    Supports both ESA SciHub SAFE format and AWS format.
    """

    # Band file patterns
    BAND_PATTERNS = {
        "B01": "*B01*.jp2",
        "B02": "*B02*.jp2",
        "B03": "*B03*.jp2",
        "B04": "*B04*.jp2",
        "B05": "*B05*.jp2",
        "B06": "*B06*.jp2",
        "B07": "*B07*.jp2",
        "B08": "*B08*.jp2",
        "B8A": "*B8A*.jp2",
        "B09": "*B09*.jp2",
        "B10": "*B10*.jp2",
        "B11": "*B11*.jp2",
        "B12": "*B12*.jp2",
    }
    _BAND_TOKEN_BY_NAME = {band_name: f"_{band_name}" for band_name in BAND_PATTERNS}
    _BAND_NAME_BY_OFFSET_ID = {
        0: "B01",
        1: "B02",
        2: "B03",
        3: "B04",
        4: "B05",
        5: "B06",
        6: "B07",
        7: "B08",
        8: "B8A",
        9: "B09",
        10: "B10",
        11: "B11",
        12: "B12",
    }

    def __init__(self, config: dict[str, Any] | None = None):
        super().__init__(config)
        self._satellite_id: str | None = None
        self._sensor_config_cache: SensorConfig | None = None
        self._sensor_config_cache_satellite_id: str | None = None
        self._granule_path: Path | None = None
        self._last_cloud_classes: xr.DataArray | None = None
        self._resolved_input_path: Path | None = None
        self._reference_grid: xr.DataArray | None = None
        self._subset_bounds: tuple[float, float, float, float] | None = None
        self._subset_crs: str | None = None

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

    def clear_spatial_subset(self) -> None:
        """Clear any active AOI read window."""
        self._subset_bounds = None
        self._subset_crs = None

    @property
    def sensor_config(self) -> SensorConfig:
        """Return sensor configuration based on satellite platform."""
        if (
            self._sensor_config_cache is not None
            and self._sensor_config_cache_satellite_id == self._satellite_id
        ):
            return self._sensor_config_cache
        if self._satellite_id is None:
            return SENTINEL2A_CONFIG
        try:
            sensor_config = load_sensor_config_with_rsrf(
                "MSI",
                self._satellite_id,
                rsrf_root=self.config.get("rsrf_root"),
            )
        except Exception as exc:
            logger.warning(
                "Falling back to built-in Sentinel-2 band metadata for %s because RSRF lookup failed (%s)",
                self._satellite_id,
                exc,
            )
            sensor_config = get_sensor_config("MSI", self._satellite_id)
        self._sensor_config_cache = sensor_config
        self._sensor_config_cache_satellite_id = self._satellite_id
        return sensor_config

    def load_toa(
        self,
        input_path: str | Path,
        *,
        band_names: tuple[str, ...] | list[str] | None = None,
    ) -> xr.Dataset:
        """Load TOA reflectance from Sentinel-2 SAFE directory."""
        input_path = Path(input_path).expanduser()
        self._resolve_paths(input_path)

        # Get scaling parameters
        metadata = self.get_metadata(input_path)
        quantification = self._quantification_value(metadata)
        offsets = metadata.get("radiometric_offsets", {})

        # Find and read band files
        img_data_path = self._get_img_data_path()
        band_files = self._discover_band_files(img_data_path)
        if not band_files:
            raise RuntimeError(f"No Sentinel-2 bands found under {img_data_path}")
        requested_band_names = self._resolve_requested_band_names(
            band_names=band_names,
            available_band_files=band_files,
        )

        ref_name = self._reference_band_name_for_alignment(
            available_band_files=band_files,
            requested_band_names=requested_band_names,
        )
        ref_file = band_files[ref_name]
        logger.debug(f"Reading reference band {ref_name} from {ref_file}")
        ref_da = self._load_aligned_band(
            band_name=ref_name,
            band_file=ref_file,
            reference_grid=None,
            quantification=quantification,
            offset=float(offsets.get(ref_name, 0.0)),
            subset_bounds=self._subset_bounds,
            subset_crs=self._subset_crs,
        )
        self._reference_grid = self._ensure_float32(ref_da, name=ref_name)
        aligned_vars: dict[str, xr.DataArray] = {}
        if ref_name in requested_band_names:
            aligned_vars[ref_name] = self._reference_grid

        # Align all bands to a single high-resolution reference grid before creating
        # the Dataset. Building a Dataset from mixed-resolution coordinates causes
        # xarray to align on the union grid, introducing sparse/warped coordinates.
        for band_name in requested_band_names:
            if band_name == ref_name:
                continue
            band_file = band_files[band_name]
            aligned_vars[band_name] = self._load_aligned_band(
                band_name=band_name,
                band_file=band_file,
                reference_grid=self._reference_grid,
                quantification=quantification,
                offset=float(offsets.get(band_name, 0.0)),
                subset_bounds=self._subset_bounds,
                subset_crs=self._subset_crs,
            )

        ds = xr.Dataset(aligned_vars)
        ds.attrs[_TOA_BAND_LOADER_ATTR] = self._make_band_loader(
            band_files=band_files,
            quantification=quantification,
            offsets=offsets,
            reference_grid=self._reference_grid,
            subset_bounds=self._subset_bounds,
            subset_crs=self._subset_crs,
            loaded_band_names=tuple(aligned_vars),
        )
        ds.attrs["loaded_toa_bands"] = tuple(aligned_vars)

        # Add metadata
        ds.attrs["sensor"] = "MSI"
        ds.attrs["satellite"] = self._satellite_id
        obs_time = metadata.get("observation_time")
        if isinstance(obs_time, datetime):
            ds.attrs["observation_time"] = obs_time.isoformat()
        elif isinstance(obs_time, str):
            ds.attrs["observation_time"] = obs_time
        else:
            ds.attrs["observation_time"] = ""

        return ds

    def _resolve_requested_band_names(
        self,
        *,
        band_names: tuple[str, ...] | list[str] | None,
        available_band_files: dict[str, Path],
    ) -> tuple[str, ...]:
        requested = band_names or tuple(self.BAND_PATTERNS)
        resolved: list[str] = []
        for band_name in requested:
            name = str(band_name)
            if name not in self.BAND_PATTERNS:
                logger.warning("Ignoring unsupported Sentinel-2 band request %s", name)
                continue
            if name not in available_band_files:
                logger.warning("Band %s not found", name)
                continue
            if name not in resolved:
                resolved.append(name)
        if not resolved:
            raise RuntimeError("No requested Sentinel-2 bands were available for TOA loading.")
        return tuple(resolved)

    def _reference_band_name_for_alignment(
        self,
        *,
        available_band_files: dict[str, Path],
        requested_band_names: tuple[str, ...],
    ) -> str:
        preference = {"B04": 0, "B03": 1, "B02": 2, "B08": 3}

        def _sort_key(name: str) -> tuple[float, int, str]:
            try:
                band = self.sensor_config.get_band(name)
                resolution = float(band.resolution)
            except Exception:
                resolution = float("inf")
            return (resolution, preference.get(name, 99), name)

        available_names = tuple(available_band_files)
        if available_names:
            return min(available_names, key=_sort_key)
        return requested_band_names[0]

    @classmethod
    def _discover_band_files(cls, img_data_path: Path) -> dict[str, Path]:
        discovered: dict[str, Path] = {}
        for band_file in sorted(img_data_path.glob("*.jp2")):
            name = band_file.name.upper()
            for band_name, token in cls._BAND_TOKEN_BY_NAME.items():
                if band_name in discovered:
                    continue
                if token in name:
                    discovered[band_name] = band_file
                    break
        return discovered

    @staticmethod
    def _ensure_float32(da: xr.DataArray, *, name: str | None = None) -> xr.DataArray:
        out = da if da.dtype == np.float32 else da.astype(np.float32)
        if name is not None and out.name != name:
            out = out.rename(name)
        return out

    def _calibrate_band(
        self,
        da: xr.DataArray,
        *,
        band_name: str,
        quantification: float,
        offset: float,
    ) -> xr.DataArray:
        out = da.astype(np.float32)
        if offset != 0.0:
            out = out + np.float32(offset)
        out = out / np.float32(quantification)
        out = out.clip(np.float32(0.0), np.float32(1.5))
        out.name = band_name
        return out

    def _load_aligned_band(
        self,
        *,
        band_name: str,
        band_file: Path,
        reference_grid: xr.DataArray | None,
        quantification: float,
        offset: float,
        subset_bounds: tuple[float, float, float, float] | None,
        subset_crs: str | None,
    ) -> xr.DataArray:
        logger.debug(f"Reading {band_name} from {band_file}")
        da = self._calibrate_band(
            self._read_band_with_subset(
                band_file,
                subset_bounds=subset_bounds,
                subset_crs=subset_crs,
            ),
            band_name=band_name,
            quantification=quantification,
            offset=offset,
        )

        if reference_grid is None or self._coords_match(da, reference_grid):
            return self._ensure_float32(da, name=band_name)

        try:
            aligned = reproject_match(da, reference_grid, resampling="bilinear")
        except Exception as exc:
            logger.debug(
                "Failed to reproject %s onto %s grid; falling back to interp (%s)",
                band_name,
                reference_grid.name or "reference",
                exc,
            )
            aligned = da.interp(
                y=reference_grid.coords["y"],
                x=reference_grid.coords["x"],
                method="linear",
            )

        return self._ensure_float32(aligned, name=band_name)

    def _make_band_loader(
        self,
        *,
        band_files: dict[str, Path],
        quantification: float,
        offsets: dict[str, Any],
        reference_grid: xr.DataArray,
        subset_bounds: tuple[float, float, float, float] | None,
        subset_crs: str | None,
        loaded_band_names: tuple[str, ...],
    ) -> Any:
        cache: dict[str, xr.DataArray] = {}
        loaded = set(loaded_band_names)

        def _load_band(band_name: str, *, native: bool = False) -> xr.DataArray:
            cache_key = (band_name, native)
            if cache_key in cache:
                return cache[cache_key]
            if band_name in loaded and not native:
                raise KeyError(f"Band {band_name!r} is already present in the loaded TOA dataset.")
            band_file = band_files.get(band_name)
            if band_file is None:
                raise KeyError(f"Band {band_name!r} is not available in the Sentinel-2 SAFE.")
            band = self._load_aligned_band(
                band_name=band_name,
                band_file=band_file,
                reference_grid=None if native else reference_grid,
                quantification=quantification,
                offset=float(offsets.get(band_name, 0.0)),
                subset_bounds=subset_bounds,
                subset_crs=subset_crs,
            )
            cache[cache_key] = band
            return band

        return _load_band

    def extract_geometry(self, input_path: str | Path) -> GeometryAngles:
        """Extract sun and view angles from metadata."""
        input_path = Path(input_path).expanduser()
        self._resolve_paths(input_path)

        # Parse angle grids from XML
        granule_xml = self._find_granule_xml()
        if granule_xml is None:
            raise FileNotFoundError(
                f"No granule metadata XML found under {self._require_granule_path()}"
            )
        tree = ET.parse(granule_xml)
        root = tree.getroot()

        # Extract sun angles (single grid for whole tile)
        sun_angles = self._parse_sun_angles(root)

        # Extract view angles (per-band, per-detector)
        view_angles = self._parse_view_angles(root)

        # Get reference band for grid
        ref_da = self._reference_grid
        if ref_da is None:
            img_data_path = self._get_img_data_path()
            band_files = self._discover_band_files(img_data_path)
            ref_name = next(
                (name for name in ("B04", "B03", "B02", "B08") if name in band_files),
                None,
            )
            if ref_name is None:
                raise FileNotFoundError(f"No reference band found under {img_data_path}")
            ref_da = self._read_band(band_files[ref_name])

        # Keep angle grids at their native 23×23 resolution (5 km spacing).
        # Every downstream consumer resamples to its own target grid, so
        # upsampling to 10980×10980 here would be wasted work.
        sza, saa, vza, vaa = self._georeference_angle_grids(
            [sun_angles["zenith"], sun_angles["azimuth"],
             view_angles["zenith"], view_angles["azimuth"]],
            ref_da,
        )

        # Convert to radians
        return GeometryAngles(
            sza=degrees_to_radians(sza),
            saa=degrees_to_radians(saa),
            vza=degrees_to_radians(vza),
            vaa=degrees_to_radians(vaa),
        )

    def _read_band(self, path: str | Path) -> xr.DataArray:
        return self._read_band_with_subset(
            path,
            subset_bounds=self._subset_bounds,
            subset_crs=self._subset_crs,
        )

    @staticmethod
    def _read_band_with_subset(
        path: str | Path,
        *,
        subset_bounds: tuple[float, float, float, float] | None,
        subset_crs: str | None,
    ) -> xr.DataArray:
        if subset_bounds is None or subset_crs is None:
            return read_raster(path)
        return read_raster_window(
            path,
            bounds=subset_bounds,
            bounds_crs=subset_crs,
        )

    def extract_cloud_mask(
        self,
        input_path: str | Path,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray:
        """
        Generate boolean cloud mask from standardized cloud classes.
        """
        input_path = Path(input_path).expanduser()

        if toa is None:
            toa = self.load_toa(input_path)
        settings = self._cloud_mask_settings(input_path)
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
        )
        self._last_cloud_classes = cloud_classes

        return classes_to_bool_mask(cloud_classes)

    def get_metadata(self, input_path: str | Path) -> dict[str, Any]:
        """Extract metadata from SAFE directory."""
        input_path = Path(input_path).expanduser()
        self._resolve_paths(input_path)

        metadata: dict[str, Any] = {
            "sensor": "MSI",
            "satellite": self._satellite_id,
            "input_path": str(input_path),
        }

        # Parse main metadata file
        mtd_file = self._find_product_xml(input_path)
        if mtd_file:
            tree = ET.parse(mtd_file)
            root = tree.getroot()

            # Extract key metadata
            metadata.update(self._parse_product_metadata(root))

        # Parse granule metadata
        granule_xml = self._find_granule_xml()
        if granule_xml:
            tree = ET.parse(granule_xml)
            root = tree.getroot()
            metadata.update(self._parse_granule_metadata(root))

        if metadata.get("observation_time") is None:
            inferred_time = self._parse_observation_time_from_name(input_path.name)
            if inferred_time is None and self._granule_path is not None:
                inferred_time = self._parse_observation_time_from_name(self._granule_path.name)
            if inferred_time is not None:
                metadata["observation_time"] = inferred_time

        return metadata

    # =========================================================================
    # Private Methods
    # =========================================================================

    def _resolve_paths(self, input_path: Path) -> None:
        """Resolve SAFE directory structure."""
        if self._resolved_input_path == input_path and self._granule_path is not None:
            return
        if self._resolved_input_path != input_path:
            self._granule_path = None
            self._satellite_id = None
            self._sensor_config_cache = None
            self._sensor_config_cache_satellite_id = None
            self._reference_grid = None
            self._resolved_input_path = input_path
        if not input_path.exists():
            raise FileNotFoundError(f"Sentinel-2 input path does not exist: {input_path}")

        # Find granule directory
        granule_dirs = sorted(input_path.glob("GRANULE/L1C_*"))
        if not granule_dirs:
            # Try AWS format
            if (input_path / "metadata.xml").exists():
                self._granule_path = input_path
            else:
                raise FileNotFoundError(f"No granule found in {input_path}")
        else:
            self._granule_path = granule_dirs[0]

        # Detect satellite platform.
        self._satellite_id = self._detect_satellite_id(input_path)

    def _get_img_data_path(self) -> Path:
        """Get path to IMG_DATA directory."""
        granule_path = self._require_granule_path()
        img_data = granule_path / "IMG_DATA"
        if not img_data.exists():
            # AWS format
            img_data = granule_path
        return img_data

    def _find_product_xml(self, input_path: Path) -> Path | None:
        """Find main product metadata XML."""
        candidates = [
            input_path / "MTD_MSIL1C.xml",
            input_path / "metadata.xml",
        ]
        for path in candidates:
            if path.exists():
                return path
        return None

    def _find_granule_xml(self) -> Path | None:
        """Find granule metadata XML."""
        granule_path = self._require_granule_path()
        candidates = list(granule_path.glob("MTD_TL.xml"))
        candidates += list(granule_path.glob("*MTD*.xml"))

        for path in candidates:
            if path.exists():
                return path
        return None

    def _parse_product_metadata(self, root: ET.Element) -> dict[str, Any]:
        """Parse product-level metadata."""
        metadata: dict[str, Any] = {}

        # Find namespace
        ns = self._get_namespace(root)

        # Processing baseline
        baseline_elem = self._find_descendant(root, "PROCESSING_BASELINE", ns)
        if baseline_elem is not None and baseline_elem.text is not None:
            metadata["processing_baseline"] = baseline_elem.text.strip()

        # Quantification value
        quant_elem = self._find_descendant(root, "QUANTIFICATION_VALUE", ns)
        if quant_elem is not None and quant_elem.text is not None:
            metadata["quantification_value"] = float(quant_elem.text)

        # Radiometric offsets (Processing Baseline >= 04.00)
        offsets: dict[str, float] = {}
        for offset_elem in self._findall_descendants(root, "RADIO_ADD_OFFSET", ns):
            band_id = offset_elem.get("band_id")
            if band_id and offset_elem.text is not None:
                band_name = self._band_name_for_offset_id(band_id)
                if band_name is None:
                    logger.warning("Ignoring unknown Sentinel-2 RADIO_ADD_OFFSET band_id=%r", band_id)
                    continue
                offsets[band_name] = float(offset_elem.text)
        if offsets:
            metadata["radiometric_offsets"] = offsets

        return metadata

    def _parse_granule_metadata(self, root: ET.Element) -> dict[str, Any]:
        """Parse granule-level metadata."""
        metadata: dict[str, Any] = {}
        ns = self._get_namespace(root)

        # Sensing time
        time_elem = self._find_descendant(root, "SENSING_TIME", ns)
        if time_elem is not None and time_elem.text:
            time_text = time_elem.text.strip()
            try:
                metadata["observation_time"] = datetime.strptime(time_text, "%Y-%m-%dT%H:%M:%S.%fZ")
            except ValueError:
                try:
                    metadata["observation_time"] = datetime.fromisoformat(time_text.replace("Z", "+00:00"))
                except ValueError:
                    logger.warning("Could not parse Sentinel-2 sensing time %r", time_text)

        # Tile ID
        tile_elem = self._find_descendant(root, "TILE_ID", ns)
        if tile_elem is not None and tile_elem.text is not None:
            metadata["tile_id"] = tile_elem.text.strip()

        return metadata

    @staticmethod
    def _parse_observation_time_from_name(name: str) -> datetime | None:
        """Infer sensing time from SAFE/granule names when XML metadata lacks it."""
        match = re.search(r"_(\d{8}T\d{6})_", name)
        if match is None:
            return None
        try:
            return datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")
        except ValueError:
            return None

    def _detect_satellite_id(self, input_path: Path) -> str:
        """Resolve satellite platform from SAFE and granule names with explicit fallback logging."""
        candidates = [input_path.name]
        if self._granule_path is not None:
            candidates.append(self._granule_path.name)
        for candidate in candidates:
            satellite_id = self._resolve_satellite_id(candidate)
            if satellite_id is not None:
                return satellite_id
        logger.warning(
            "Could not infer Sentinel-2 platform from %s; defaulting to S2A",
            ", ".join(repr(candidate) for candidate in candidates),
        )
        return "S2A"

    @staticmethod
    def _resolve_satellite_id(name: str) -> str | None:
        """Infer Sentinel-2 platform identifier from a SAFE or granule name."""
        for satellite_id in ("S2A", "S2B", "S2C"):
            if satellite_id in name:
                return satellite_id
        return None

    def _parse_sun_angles(self, root: ET.Element) -> dict[str, np.ndarray]:
        """Parse sun angle grids from XML."""
        ns = self._get_namespace(root)

        angles = {}

        # Find Sun_Angles_Grid
        sun_grid = root.find(f".//{ns}Sun_Angles_Grid")
        if sun_grid is None:
            # Try without namespace
            sun_grid = root.find(".//Sun_Angles_Grid")

        if sun_grid is not None:
            zenith = sun_grid.find(f".//{ns}Zenith")
            if zenith is None:
                zenith = sun_grid.find(".//Zenith")

            azimuth = sun_grid.find(f".//{ns}Azimuth")
            if azimuth is None:
                azimuth = sun_grid.find(".//Azimuth")

            if zenith is not None:
                angles["zenith"] = self._parse_angle_grid(zenith, ns)
            if azimuth is not None:
                angles["azimuth"] = self._parse_angle_grid(azimuth, ns)

        return angles

    def _parse_view_angles(self, root: ET.Element) -> dict[str, np.ndarray]:
        """Parse view angle grids from XML (mean across detectors)."""
        ns = self._get_namespace(root)

        # Get mean viewing angles
        mean_vza: float | None = None
        mean_vaa: float | None = None

        # Look for Mean_Viewing_Incidence_Angle
        for mean_elem in root.findall(f".//{ns}Mean_Viewing_Incidence_Angle"):
            zenith = mean_elem.find(f"{ns}ZENITH_ANGLE")
            azimuth = mean_elem.find(f"{ns}AZIMUTH_ANGLE")

            if zenith is not None and azimuth is not None and zenith.text is not None and azimuth.text is not None:
                if mean_vza is None:
                    mean_vza = float(zenith.text)
                    mean_vaa = float(azimuth.text)
                else:
                    assert mean_vaa is not None
                    # Average across bands
                    mean_vza = (mean_vza + float(zenith.text)) / 2
                    mean_vaa = (mean_vaa + float(azimuth.text)) / 2

        # Fallback values
        if mean_vza is None or mean_vaa is None:
            mean_vza = 5.0
            mean_vaa = 100.0

        # Create uniform grids (simplified - full implementation would parse per-detector grids)
        zenith_grid = np.full((23, 23), mean_vza)
        azimuth_grid = np.full((23, 23), mean_vaa)

        return {
            "zenith": zenith_grid,
            "azimuth": azimuth_grid,
        }

    def _parse_angle_grid(self, elem: ET.Element, ns: str) -> np.ndarray:
        """Parse angle values from grid element."""
        values_list = elem.find(f"{ns}Values_List")
        if values_list is None:
            values_list = elem.find("Values_List")

        if values_list is None:
            # Return default grid
            return np.full((23, 23), 30.0)

        rows: list[list[float]] = []
        for values in values_list.findall(f"{ns}VALUES"):
            if values is None:
                values = values_list.find("VALUES")
            if values is not None and values.text:
                row = [float(v) for v in values.text.split()]
                rows.append(row)

        if not rows:
            return np.full((23, 23), 30.0)

        return cast("np.ndarray[Any, np.dtype[np.float32]]", np.asarray(rows, dtype=np.float32))

    def _angles_to_grid(
        self,
        angles: np.ndarray,
        target: xr.DataArray,
    ) -> xr.DataArray:
        """Interpolate angle grid to image resolution."""
        # Create DataArray from angle grid
        # Angles are on 5km grid (23x23 for 10980x10980 image at 10m)
        height, width = angles.shape

        # Approximate coordinates
        bounds = target.rio.bounds()
        x = np.linspace(bounds[0], bounds[2], width)
        y = np.linspace(bounds[3], bounds[1], height)  # Note: y is inverted

        da = xr.DataArray(
            angles.astype(np.float32),
            dims=["y", "x"],
            coords={"y": y, "x": x},
        )
        da = da.rio.write_crs(target.rio.crs)

        # Resample to target grid
        return reproject_match(da, target, resampling="bilinear")

    def _angles_to_grid_batch(
        self,
        angle_grids: list[np.ndarray],
        target: xr.DataArray,
    ) -> list[xr.DataArray]:
        """Upsample multiple angle grids to image resolution using PIL resize.

        The angle grids (23×23) live on the same CRS/extent as the target — no
        actual reprojection is needed, just bilinear upsampling.  PIL's
        float32-mode resize is ~7× faster than scipy.ndimage.zoom for this
        extreme upsampling ratio.
        """
        from PIL import Image

        h_out, w_out = target.sizes["y"], target.sizes["x"]
        results: list[xr.DataArray] = []
        for angles in angle_grids:
            src = angles.astype(np.float32)
            # PIL.resize takes (width, height) — opposite of numpy (rows, cols).
            out = np.array(
                Image.fromarray(src, mode="F").resize((w_out, h_out), Image.BILINEAR),
                dtype=np.float32,
            )
            da = xr.DataArray(
                out,
                dims=["y", "x"],
                coords={"y": target.coords["y"], "x": target.coords["x"]},
                attrs=target.attrs,
            )
            if target.rio.crs is not None:
                da = da.rio.write_crs(target.rio.crs)
                da = da.rio.write_transform(target.rio.transform())
            results.append(da)
        return results

    def _georeference_angle_grids(
        self,
        angle_grids: list[np.ndarray],
        ref_da: xr.DataArray,
    ) -> list[xr.DataArray]:
        """Wrap raw angle grids (23×23) as georeferenced DataArrays.

        The grids keep their native resolution but receive y/x coordinates
        and CRS metadata derived from *ref_da* so that downstream
        resampling (``resample_field_to_template``, ``_resample_da``, etc.)
        can correctly interpolate them to any target grid.
        """
        bounds = ref_da.rio.bounds()          # (left, bottom, right, top)
        crs = ref_da.rio.crs

        results: list[xr.DataArray] = []
        for angles in angle_grids:
            src = np.asarray(angles, dtype=np.float32)
            h, w = src.shape
            # Pixel-centre coordinates spanning the tile extent.
            x = np.linspace(bounds[0], bounds[2], w, dtype=np.float64)
            y = np.linspace(bounds[3], bounds[1], h, dtype=np.float64)  # top→bottom
            da = xr.DataArray(
                src,
                dims=["y", "x"],
                coords={"y": y, "x": x},
            )
            if crs is not None:
                da = da.rio.write_crs(crs)
            results.append(da)
        return results

    def _cloud_mask_settings(self, input_path: Path) -> dict[str, Any]:
        """Resolve cloud-mask settings from preprocessor config and input path."""
        defaults: dict[str, Any] = {
            "mode": "auto",
            "provider": "omnicloudmask",
            "class_mapping": None,
            "external_mask_path": None,
            "user_callable": None,
            "target_resolution_m": 10.0,
            "resolution_policy": "auto",
            "allow_upsample_to_target": False,
            "unmapped_to_missing": True,
        }
        cloud_mask_config = self.config.get("cloud_mask", {}) if isinstance(self.config, dict) else {}
        if isinstance(cloud_mask_config, dict):
            defaults.update(cloud_mask_config)
        # Resolve external file path relative to SAFE directory when needed.
        external = defaults.get("external_mask_path")
        if isinstance(external, (str, Path)) and str(external):
            p = Path(external).expanduser()
            if not p.is_absolute():
                p = input_path / p
            defaults["external_mask_path"] = p
        return defaults

    @staticmethod
    def _quantification_value(metadata: dict[str, Any]) -> float:
        raw = metadata.get("quantification_value", 10000.0)
        quantification = float(raw)
        if not np.isfinite(quantification) or quantification <= 0.0:
            raise ValueError(
                f"Sentinel-2 quantification_value must be finite and > 0, got {raw!r}"
            )
        return quantification

    def _require_granule_path(self) -> Path:
        if self._granule_path is None:
            raise RuntimeError("Sentinel-2 granule path is not resolved; call _resolve_paths() first.")
        return self._granule_path

    @staticmethod
    def _coords_match(a: xr.DataArray, b: xr.DataArray) -> bool:
        """Return True when two arrays share the same x/y grid."""
        if "x" not in a.coords or "y" not in a.coords or "x" not in b.coords or "y" not in b.coords:
            return False
        return (
            a.sizes.get("x") == b.sizes.get("x")
            and a.sizes.get("y") == b.sizes.get("y")
            and np.array_equal(a.coords["x"].values, b.coords["x"].values)
            and np.array_equal(a.coords["y"].values, b.coords["y"].values)
        )

    def _get_namespace(self, root: ET.Element) -> str:
        """Extract XML namespace from root element."""
        match = re.match(r"\{(.+)\}", root.tag)
        if match:
            return "{" + match.group(1) + "}"
        return ""

    @staticmethod
    def _find_descendant(root: ET.Element, tag: str, ns: str) -> ET.Element | None:
        elem = root.find(f".//{ns}{tag}") if ns else None
        if elem is None:
            elem = root.find(f".//{tag}")
        return elem

    @staticmethod
    def _findall_descendants(root: ET.Element, tag: str, ns: str) -> list[ET.Element]:
        elems = root.findall(f".//{ns}{tag}") if ns else []
        if elems:
            return list(elems)
        return list(root.findall(f".//{tag}"))

    @classmethod
    def _band_name_for_offset_id(cls, band_id: str) -> str | None:
        try:
            index = int(band_id)
        except (TypeError, ValueError):
            return None
        return cls._BAND_NAME_BY_OFFSET_ID.get(index)

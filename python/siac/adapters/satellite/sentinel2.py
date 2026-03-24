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

    def __init__(self, config: dict[str, Any] | None = None):
        super().__init__(config)
        self._satellite_id: str | None = None
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
        self._subset_bounds = tuple(float(value) for value in bounds)
        self._subset_crs = str(crs)

    def clear_spatial_subset(self) -> None:
        """Clear any active AOI read window."""
        self._subset_bounds = None
        self._subset_crs = None

    @property
    def sensor_config(self) -> SensorConfig:
        """Return sensor configuration based on satellite platform."""
        if self._satellite_id is None:
            return SENTINEL2A_CONFIG
        try:
            return load_sensor_config_with_rsrf(
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
            return get_sensor_config("MSI", self._satellite_id)

    def load_toa(self, input_path: str | Path) -> xr.Dataset:
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

        ref_name = next(
            (name for name in ("B04", "B03", "B02", "B08") if name in band_files),
            next(iter(band_files)),
        )
        ref_file = band_files[ref_name]
        logger.debug(f"Reading reference band {ref_name} from {ref_file}")
        ref_da = self._calibrate_band(
            self._read_band(ref_file),
            band_name=ref_name,
            quantification=quantification,
            offset=float(offsets.get(ref_name, 0.0)),
        )
        self._reference_grid = self._ensure_float32(ref_da, name=ref_name)
        aligned_vars: dict[str, xr.DataArray] = {}
        aligned_vars[ref_name] = self._reference_grid

        # Align all bands to a single high-resolution reference grid before creating
        # the Dataset. Building a Dataset from mixed-resolution coordinates causes
        # xarray to align on the union grid, introducing sparse/warped coordinates.
        for band_name in self.BAND_PATTERNS:
            band_file = band_files.get(band_name)
            if band_file is None:
                logger.warning(f"Band {band_name} not found")
                continue
            if band_name == ref_name:
                continue

            logger.debug(f"Reading {band_name} from {band_file}")
            da = self._calibrate_band(
                self._read_band(band_file),
                band_name=band_name,
                quantification=quantification,
                offset=float(offsets.get(band_name, 0.0)),
            )

            if self._coords_match(da, self._reference_grid):
                aligned_vars[band_name] = da
                continue

            try:
                aligned = reproject_match(da, self._reference_grid, resampling="bilinear")
            except Exception as exc:
                logger.debug(
                    "Failed to reproject %s onto %s grid; falling back to interp (%s)",
                    band_name,
                    ref_name,
                    exc,
                )
                aligned = da.interp(
                    y=self._reference_grid.coords["y"],
                    x=self._reference_grid.coords["x"],
                    method="linear",
                )

            aligned.name = band_name
            aligned_vars[band_name] = self._ensure_float32(aligned, name=band_name)

        ds = xr.Dataset(aligned_vars)

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
        values = np.asarray(da.data, dtype=np.float32)
        if values is da.data and (name is None or da.name == name):
            return da
        out = da.copy(deep=False)
        out.data = values
        if name is not None:
            out.name = name
        return out

    def _calibrate_band(
        self,
        da: xr.DataArray,
        *,
        band_name: str,
        quantification: float,
        offset: float,
    ) -> xr.DataArray:
        values = np.array(da.data, dtype=np.float32, copy=True)
        if offset != 0.0:
            np.add(values, np.float32(offset), out=values)
        np.divide(values, np.float32(quantification), out=values)
        np.clip(values, np.float32(0.0), np.float32(1.5), out=values)

        out = da.copy(deep=False)
        out.data = values
        out.name = band_name
        return out

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

        # Resample angles to image grid
        sza = self._angles_to_grid(sun_angles["zenith"], ref_da)
        saa = self._angles_to_grid(sun_angles["azimuth"], ref_da)

        # Use mean view angles (average across detectors)
        vza = self._angles_to_grid(view_angles["zenith"], ref_da)
        vaa = self._angles_to_grid(view_angles["azimuth"], ref_da)

        # Convert to radians
        return GeometryAngles(
            sza=degrees_to_radians(sza),
            saa=degrees_to_radians(saa),
            vza=degrees_to_radians(vza),
            vaa=degrees_to_radians(vaa),
        )

    def _read_band(self, path: str | Path) -> xr.DataArray:
        if self._subset_bounds is None or self._subset_crs is None:
            return read_raster(path)
        return read_raster_window(
            path,
            bounds=self._subset_bounds,
            bounds_crs=self._subset_crs,
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
        try:
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
        except ValueError as exc:
            if settings["mode"] != "auto" or "Could not find any" not in str(exc):
                raise
            logger.warning(
                "Auto cloud-mask mode needs red/green/nir bands; "
                "falling back to clear mask because required bands are missing."
            )
            cloud_classes = build_cloud_classes(
                toa,
                self.sensor_config,
                mode="none",
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
        baseline_elem = root.find(f".//{ns}PROCESSING_BASELINE")
        if baseline_elem is not None and baseline_elem.text is not None:
            metadata["processing_baseline"] = baseline_elem.text.strip()

        # Quantification value
        quant_elem = root.find(f".//{ns}QUANTIFICATION_VALUE")
        if quant_elem is not None and quant_elem.text is not None:
            metadata["quantification_value"] = float(quant_elem.text)

        # Radiometric offsets (Processing Baseline >= 04.00)
        offsets: dict[str, float] = {}
        for offset_elem in root.findall(f".//{ns}RADIO_ADD_OFFSET"):
            band_id = offset_elem.get("band_id")
            if band_id and offset_elem.text is not None:
                offsets[f"B{int(band_id):02d}"] = float(offset_elem.text)
        if offsets:
            metadata["radiometric_offsets"] = offsets

        return metadata

    def _parse_granule_metadata(self, root: ET.Element) -> dict[str, Any]:
        """Parse granule-level metadata."""
        metadata: dict[str, Any] = {}
        ns = self._get_namespace(root)

        # Sensing time
        time_elem = root.find(f".//{ns}SENSING_TIME")
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
        tile_elem = root.find(f".//{ns}TILE_ID")
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

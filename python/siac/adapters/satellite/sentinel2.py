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
        """Return RSRF-backed sensor configuration based on satellite platform."""
        satellite_id = self._satellite_id or "S2A"
        if (
            self._sensor_config_cache is not None
            and self._sensor_config_cache_satellite_id == satellite_id
        ):
            return self._sensor_config_cache
        try:
            sensor_config = load_sensor_config_with_rsrf(
                "MSI",
                satellite_id,
                rsrf_root=self.config.get("rsrf_root"),
            )
        except Exception as exc:
            raise RuntimeError(
                "Unable to load Sentinel-2 RSRF metadata for "
                f"{satellite_id}. Configure paths.rsrf_root to a catalog containing "
                f"{satellite_id} MSI band responses."
            ) from exc
        self._sensor_config_cache = sensor_config
        self._sensor_config_cache_satellite_id = satellite_id
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
        cache: dict[tuple[str, bool], xr.DataArray] = {}
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
            [
                sun_angles["zenith"],
                self._unwrap_azimuth_grid_degrees(sun_angles["azimuth"]),
                view_angles["zenith"],
                self._unwrap_azimuth_grid_degrees(view_angles["azimuth"]),
            ],
            ref_da,
            metadata_root=root,
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
            cloud_cache_dir=settings.get("cache_dir"),
            inference_device=settings.get("inference_device"),
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
                    logger.warning(
                        "Ignoring unknown Sentinel-2 RADIO_ADD_OFFSET band_id=%r", band_id
                    )
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
                    metadata["observation_time"] = datetime.fromisoformat(
                        time_text.replace("Z", "+00:00")
                    )
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
        """Parse sun angle grids from XML.

        Wave 15: warn loudly when expected elements are missing. Previously
        this method silently returned an empty / partial dict; the caller
        had no signal that sun angles were missing and would fall through
        to whatever its own default was — exactly the silent-fallback
        pattern that hid the wave-14 namespace bug.
        """
        ns = self._get_namespace(root)

        angles: dict[str, np.ndarray] = {}

        sun_grid = self._find_descendant(root, "Sun_Angles_Grid", ns)
        if sun_grid is None:
            logger.warning(
                "Sentinel-2 MTD_TL.xml has no Sun_Angles_Grid element — "
                "downstream code will see an empty sun-angle dict. This is "
                "almost always a bug (every well-formed S2 SAFE has this "
                "element); check whether the MTD_TL.xml is truncated."
            )
            return angles

        zenith = self._find_descendant(sun_grid, "Zenith", ns)
        azimuth = self._find_descendant(sun_grid, "Azimuth", ns)

        if zenith is not None:
            angles["zenith"] = self._parse_angle_grid(zenith, ns)
        else:
            logger.warning(
                "Sentinel-2 Sun_Angles_Grid is missing the Zenith sub-element; "
                "downstream sun-zenith will fall back to whatever default the "
                "caller uses."
            )

        if azimuth is not None:
            angles["azimuth"] = self._parse_angle_grid(azimuth, ns)
        else:
            logger.warning(
                "Sentinel-2 Sun_Angles_Grid is missing the Azimuth sub-element; "
                "downstream sun-azimuth will fall back to whatever default the "
                "caller uses."
            )

        return angles

    def _parse_view_angles(self, root: ET.Element) -> dict[str, np.ndarray]:
        """Parse spatial view-angle grids and combine bands/detectors.

        REVIEW.md §1.1 #5: previously this used a running ``(a + b) / 2`` step
        per band, which biases toward the last-seen entry rather than producing
        the arithmetic mean. We now collect all values and average once.

        Wave 14: also uses the namespace-tolerant ``_findall_descendants`` /
        ``_find_child`` helpers because S2 MTD_TL.xml puts ``xmlns:n1`` on
        the root only, leaving every descendant in the empty namespace.
        The previous code searched with the prefix and missed every
        match, silently falling back to ``DEFAULT_S2_VZA_DEG = 5.0`` and
        ``DEFAULT_S2_VAA_DEG = 100.0`` on every real scene.
        """
        ns = self._get_namespace(root)

        zenith_values: list[float] = []
        azimuth_values: list[float] = []

        for mean_elem in self._findall_descendants(root, "Mean_Viewing_Incidence_Angle", ns):
            zenith = self._find_child(mean_elem, "ZENITH_ANGLE", ns)
            azimuth = self._find_child(mean_elem, "AZIMUTH_ANGLE", ns)
            if (
                zenith is not None
                and azimuth is not None
                and zenith.text is not None
                and azimuth.text is not None
            ):
                zenith_values.append(float(zenith.text))
                azimuth_values.append(float(azimuth.text))

        if zenith_values and azimuth_values:
            mean_vza = float(np.mean(zenith_values))
            mean_vaa = self._circular_nanmean_degrees(np.asarray(azimuth_values))
        else:
            # Fallback to defaults from ``siac.constants``; warn so the
            # operator notices the XML didn't carry the angles we needed
            # (REVIEW.md §3.3 sentinel2.py).
            from siac.constants import DEFAULT_S2_VAA_DEG, DEFAULT_S2_VZA_DEG

            logger.warning(
                "Mean_Viewing_Incidence_Angle entries not found in MTD_TL.xml; "
                "falling back to default VZA=%.1fdeg, VAA=%.1fdeg.",
                DEFAULT_S2_VZA_DEG,
                DEFAULT_S2_VAA_DEG,
            )
            mean_vza = DEFAULT_S2_VZA_DEG
            mean_vaa = DEFAULT_S2_VAA_DEG

        detector_entries = self._findall_descendants(root, "Viewing_Incidence_Angles_Grids", ns)
        by_band: dict[int, list[tuple[np.ndarray, np.ndarray]]] = {}
        for entry in detector_entries:
            band_id = entry.get("bandId")
            if band_id is None:
                continue
            try:
                zenith = self._parse_named_angle_grid(entry, "Zenith", ns)
                azimuth = self._parse_named_angle_grid(entry, "Azimuth", ns)
                by_band.setdefault(int(band_id), []).append((zenith, azimuth))
            except (TypeError, ValueError):
                logger.warning(
                    "Ignoring malformed Sentinel-2 detector angle grid bandId=%r detectorId=%r",
                    band_id,
                    entry.get("detectorId"),
                    exc_info=True,
                )

        band_zenith: list[np.ndarray] = []
        band_azimuth: list[np.ndarray] = []
        for entries in by_band.values():
            zenith_stack = np.stack([entry[0] for entry in entries])
            azimuth_stack = np.stack([entry[1] for entry in entries])
            band_zenith.append(self._nanmean(zenith_stack, axis=0))
            band_azimuth.append(self._circular_nanmean_grid_degrees(azimuth_stack, axis=0))

        if band_zenith:
            zenith_grid = self._nanmean(np.stack(band_zenith), axis=0)
            azimuth_grid = self._circular_nanmean_grid_degrees(np.stack(band_azimuth), axis=0)
            zenith_grid = np.where(np.isfinite(zenith_grid), zenith_grid, mean_vza)
            azimuth_grid = np.where(np.isfinite(azimuth_grid), azimuth_grid, mean_vaa)
        else:
            logger.warning(
                "Sentinel-2 detector viewing grids are unavailable; using tile-mean viewing "
                "geometry."
            )
            zenith_grid = np.full((23, 23), mean_vza)
            azimuth_grid = np.full((23, 23), mean_vaa)

        return {
            "zenith": np.asarray(zenith_grid, dtype=np.float32),
            "azimuth": np.asarray(azimuth_grid, dtype=np.float32),
        }

    def _parse_named_angle_grid(
        self,
        parent: ET.Element,
        name: str,
        ns: str,
    ) -> np.ndarray:
        angle = self._find_child(parent, name, ns)
        if angle is None:
            raise ValueError(f"Missing {name} angle grid")
        return self._parse_angle_grid(angle, ns)

    @staticmethod
    def _nanmean(values: np.ndarray, *, axis: int) -> np.ndarray:
        finite = np.isfinite(values)
        count = np.sum(finite, axis=axis)
        total = np.sum(np.where(finite, values, 0.0), axis=axis)
        return np.divide(
            total,
            count,
            out=np.full_like(total, np.nan, dtype=np.float64),
            where=count > 0,
        )

    @staticmethod
    def _circular_nanmean_grid_degrees(values: np.ndarray, *, axis: int) -> np.ndarray:
        radians = np.deg2rad(values)
        finite = np.isfinite(radians)
        count = np.sum(finite, axis=axis)
        sine = np.sum(np.where(finite, np.sin(radians), 0.0), axis=axis)
        cosine = np.sum(np.where(finite, np.cos(radians), 0.0), axis=axis)
        angle = np.mod(np.rad2deg(np.arctan2(sine, cosine)), 360.0)
        return np.where(count > 0, angle, np.nan)

    @classmethod
    def _circular_nanmean_degrees(cls, values: np.ndarray) -> float:
        return float(cls._circular_nanmean_grid_degrees(values, axis=0))

    @staticmethod
    def _unwrap_azimuth_grid_degrees(values: np.ndarray) -> np.ndarray:
        radians = np.deg2rad(np.asarray(values, dtype=np.float64))
        unwrapped = np.unwrap(np.unwrap(radians, axis=1), axis=0)
        return np.asarray(np.rad2deg(unwrapped), dtype=np.float32)

    def _parse_angle_grid(self, elem: ET.Element, ns: str) -> np.ndarray:
        """Parse angle values from grid element.

        Wave 14: uses the namespace-tolerant ``_find_child`` /
        ``_findall_children`` helpers. The previous code searched
        ``f"{ns}Values_List"`` first then fell back to ``"Values_List"``,
        but the ``VALUES`` rows inside used the prefix-only form with no
        fallback — so when the elements were unprefixed (the normal
        S2 case), ``rows`` came back empty and the helper returned a
        magic 30.0° uniform grid.
        """
        from siac.constants import DEFAULT_S2_ANGLE_GRID_DEG

        values_list = self._find_child(elem, "Values_List", ns)

        if values_list is None:
            logger.warning(
                "Sentinel-2 angle grid Values_List missing; falling back to uniform %.1fdeg grid.",
                DEFAULT_S2_ANGLE_GRID_DEG,
            )
            return np.full((23, 23), DEFAULT_S2_ANGLE_GRID_DEG)

        rows: list[list[float]] = []
        for values in self._findall_children(values_list, "VALUES", ns):
            if values.text:
                row = [float(v) for v in values.text.split()]
                rows.append(row)

        if not rows:
            logger.warning(
                "Sentinel-2 angle grid Values_List empty; falling back to uniform %.1fdeg grid.",
                DEFAULT_S2_ANGLE_GRID_DEG,
            )
            return np.full((23, 23), DEFAULT_S2_ANGLE_GRID_DEG)

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
            src: np.ndarray = np.asarray(angles, dtype=np.float32)
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
        *,
        metadata_root: ET.Element | None = None,
    ) -> list[xr.DataArray]:
        """Wrap raw angle grids (23×23) as georeferenced DataArrays.

        The grids keep their native resolution but receive y/x coordinates
        and CRS metadata derived from *ref_da* so that downstream
        resampling (``resample_field_to_template``, ``_resample_da``, etc.)
        can correctly interpolate them to any target grid.
        """
        bounds = ref_da.rio.bounds()  # (left, bottom, right, top)
        crs = ref_da.rio.crs
        origin_and_step = self._angle_grid_origin_and_step(metadata_root)
        if origin_and_step is not None:
            ulx, uly, col_step, row_step, metadata_crs = origin_and_step
            crs = metadata_crs or crs

        results: list[xr.DataArray] = []
        for angles in angle_grids:
            src = np.asarray(angles, dtype=np.float32)
            h, w = src.shape
            if origin_and_step is None:
                x = np.linspace(bounds[0], bounds[2], w, dtype=np.float64)
                y = np.linspace(bounds[3], bounds[1], h, dtype=np.float64)
            else:
                x = ulx + np.arange(w, dtype=np.float64) * col_step
                y = uly - np.arange(h, dtype=np.float64) * row_step
            da = xr.DataArray(
                src,
                dims=["y", "x"],
                coords={"y": y, "x": x},
            )
            if crs is not None:
                da = da.rio.write_crs(crs)
            results.append(da)
        return results

    def _angle_grid_origin_and_step(
        self,
        root: ET.Element | None,
    ) -> tuple[float, float, float, float, str | None] | None:
        if root is None:
            return None
        ns = self._get_namespace(root)
        geopositions = self._findall_descendants(root, "Geoposition", ns)
        geoposition = next(
            (item for item in geopositions if item.get("resolution") == "10"),
            geopositions[0] if geopositions else None,
        )
        sun_grid = self._find_descendant(root, "Sun_Angles_Grid", ns)
        if geoposition is None or sun_grid is None:
            return None
        zenith = self._find_child(sun_grid, "Zenith", ns)
        if zenith is None:
            return None
        try:
            ulx = float(self._required_child_text(geoposition, "ULX", ns))
            uly = float(self._required_child_text(geoposition, "ULY", ns))
            col_step = float(self._required_child_text(zenith, "COL_STEP", ns))
            row_step = float(self._required_child_text(zenith, "ROW_STEP", ns))
        except (TypeError, ValueError):
            logger.warning(
                "Could not parse Sentinel-2 angle-grid georeferencing; deriving it from the "
                "loaded raster extent.",
                exc_info=True,
            )
            return None
        crs_element = self._find_descendant(root, "HORIZONTAL_CS_CODE", ns)
        metadata_crs = None if crs_element is None else crs_element.text
        return ulx, uly, col_step, row_step, metadata_crs

    def _required_child_text(self, parent: ET.Element, name: str, ns: str) -> str:
        child = self._find_child(parent, name, ns)
        if child is None or child.text is None:
            raise ValueError(f"Missing {name} below {parent.tag}")
        return child.text

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
            # Wave 19a: default to CPU for bit-deterministic cross-process
            # OmniCloudMask inference. Override via the config field
            # ``algorithms.cloud_mask.inference_device`` (e.g. ``"auto"``,
            # ``"mps"``, ``"cuda"``) when determinism isn't required.
            "inference_device": "cpu",
        }
        cloud_mask_config = (
            self.config.get("cloud_mask", {}) if isinstance(self.config, dict) else {}
        )
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
        # Wave 15: previously this silently substituted the standard
        # ``10000.0`` whenever the metadata dict lacked ``quantification_value``.
        # That value IS the conventional S2 default, but a missing value in
        # metadata indicates a real parse failure in
        # ``_parse_product_metadata`` (the XML's QUANTIFICATION_VALUE element).
        # Warn so the operator notices, but keep the fallback so production
        # runs against truncated metadata don't hard-fail.
        if "quantification_value" not in metadata:
            logger.warning(
                "Sentinel-2 product metadata is missing QUANTIFICATION_VALUE; "
                "falling back to the standard 10000.0. The MTD_MSIL1C.xml is "
                "almost certainly truncated or non-standard — verify the SAFE "
                "integrity."
            )
        raw = metadata.get("quantification_value", 10000.0)
        quantification = float(raw)
        if not np.isfinite(quantification) or quantification <= 0.0:
            raise ValueError(f"Sentinel-2 quantification_value must be finite and > 0, got {raw!r}")
        return quantification

    def _require_granule_path(self) -> Path:
        if self._granule_path is None:
            raise RuntimeError(
                "Sentinel-2 granule path is not resolved; call _resolve_paths() first."
            )
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

    @staticmethod
    def _find_child(parent: ET.Element, tag: str, ns: str) -> ET.Element | None:
        """Find a direct child by tag, trying namespaced first then unprefixed.

        Critical for S2 MTD_TL.xml: the root carries ``xmlns:n1=...`` but
        the child elements (``Mean_Viewing_Incidence_Angle``, ``ZENITH_ANGLE``,
        ``Values_List``, etc.) are unprefixed and live in *no* namespace.
        A naïve ``parent.find(f"{ns}TAG")`` returns ``None`` for those,
        which was silently making ``_parse_view_angles`` /
        ``_parse_angle_grid`` fall back to the magic
        ``DEFAULT_S2_*`` constants on every real scene.
        """
        elem = parent.find(f"{ns}{tag}") if ns else None
        if elem is None:
            elem = parent.find(tag)
        return elem

    @staticmethod
    def _findall_children(parent: ET.Element, tag: str, ns: str) -> list[ET.Element]:
        """findall direct children by tag, trying namespaced first then unprefixed.

        See :meth:`_find_child` for the rationale.
        """
        elems = parent.findall(f"{ns}{tag}") if ns else []
        if elems:
            return list(elems)
        return list(parent.findall(tag))

    @classmethod
    def _band_name_for_offset_id(cls, band_id: str) -> str | None:
        try:
            index = int(band_id)
        except (TypeError, ValueError):
            return None
        return cls._BAND_NAME_BY_OFFSET_ID.get(index)

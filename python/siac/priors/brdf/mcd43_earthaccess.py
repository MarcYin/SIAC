"""Earthaccess-backed BRDF providers for MCD43, VNP43, and MCD19."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from siac.core.types import BRDFKernelWeights, SensorBand
from siac.io.earthaccess_catalog import EarthAccessCatalog
from siac.io.earthaccess_source import EarthAccessSource
from siac.priors.earthdata_common import (
    ProductBandDefinition,
    apply_scale_and_mask,
    granule_intersects_bounds,
    make_native_grid_dataarray,
    parse_granule_date,
    parse_tile_indices,
    read_hdf4_dataset,
    read_hdf5_dataset,
    reproject_native_to_target,
)

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import datetime

logger = logging.getLogger(__name__)


class _EarthAccessBRDFProvider:
    """Shared logic for Earthaccess BRDF product wrappers."""

    product_key: str = ""
    _source_name: str = ""
    _product_bands: tuple[ProductBandDefinition, ...] = ()
    _legacy_band_map: dict[int, str] = {}

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        probe_earthdata: bool = True,
        max_granules: int = 8,
    ) -> None:
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
        self.source = source or EarthAccessSource(provider=provider)
        self.catalog = catalog or EarthAccessCatalog(source=self.source)
        self.short_name = short_name
        self.provider = provider
        self.probe_earthdata = probe_earthdata
        self.max_granules = max(1, int(max_granules))

    @property
    def source_name(self) -> str:
        return self._source_name

    @property
    def _bands_by_label(self) -> dict[str, ProductBandDefinition]:
        return {band.label: band for band in self._product_bands}

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int] | Sequence[SensorBand],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights:
        """Return BRDF kernel weights on the requested grid."""
        if target_resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {target_resolution}")
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        requested = self._resolve_requested_bands(list(bands))
        if self.probe_earthdata:
            paths = self._download_granules(bounds, crs, obs_time, temporal_window)
            if paths:
                try:
                    return self._load_from_granules(
                        paths,
                        requested=requested,
                        bounds=bounds,
                        crs=crs,
                        target_resolution=target_resolution,
                    )
                except Exception as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s BRDF granule parsing failed; using defaults (%s)",
                        self._source_name,
                        exc,
                    )

        return self._default_weights(
            bounds,
            target_resolution,
            [coord for coord, _band in requested],
        )

    def _resolve_requested_bands(
        self,
        bands: list[int] | list[SensorBand],
    ) -> list[tuple[int | str, ProductBandDefinition]]:
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        product_bands = self._product_bands
        if isinstance(bands[0], SensorBand):
            resolved: list[tuple[int | str, ProductBandDefinition]] = []
            for sensor_band in bands:
                match = min(
                    product_bands,
                    key=lambda candidate: abs(candidate.wavelength_nm - sensor_band.center_wavelength),
                )
                resolved.append((sensor_band.name, match))
            return resolved

        resolved = []
        for band_id in bands:
            label = self._legacy_band_map.get(int(band_id))
            if label is None:
                raise KeyError(f"Band {band_id!r} is not available in {self._source_name}")
            resolved.append((int(band_id), self._bands_by_label[label]))
        return resolved

    def _download_granules(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        temporal_window: int,
    ) -> list[Path]:
        short_name = self.short_name or self.catalog.resolve_short_name(self.product_key)
        temporal = EarthAccessSource.temporal_window(obs_time, temporal_window)
        granules = self.source.search_granules(
            short_name=short_name,
            bounds=bounds,
            crs=crs,
            temporal=temporal,
            provider=self.provider,
            count=self.max_granules,
        )
        if not granules:
            logger.warning(
                "%s granule probe returned no results for AOI/time window",
                self._source_name,
            )
            return []

        dest = self.cache_dir or Path.home() / ".cache" / "siac" / "earthdata" / short_name
        downloaded = self.source.download_granules(granules, dest)
        return self._select_candidate_paths(downloaded, obs_time, bounds, crs)

    @staticmethod
    def _select_candidate_paths(
        paths: list[Path],
        obs_time: datetime,
        bounds: tuple[float, float, float, float],
        crs: str,
    ) -> list[Path]:
        if not paths:
            return []

        selected: list[tuple[tuple[int, int], float, str, Path]] = []
        for path in paths:
            try:
                delta = abs((parse_granule_date(path) - obs_time).total_seconds())
                intersects = granule_intersects_bounds(path, bounds=bounds, crs=crs)
            except Exception:
                return paths

            if intersects:
                try:
                    tile = parse_tile_indices(path)
                except Exception:
                    tile = (999, 999)
                selected.append((tile, delta, Path(path).name, path))

        if not selected:
            return []
        return [item[3] for item in sorted(selected, key=lambda value: (value[0], value[1], value[2]))]

    def _load_from_granules(
        self,
        paths: list[Path],
        *,
        requested: list[tuple[int | str, ProductBandDefinition]],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
    ) -> BRDFKernelWeights:
        f0_list: list[xr.DataArray] = []
        f1_list: list[xr.DataArray] = []
        f2_list: list[xr.DataArray] = []
        f0_unc_list: list[xr.DataArray] = []
        f1_unc_list: list[xr.DataArray] = []
        f2_unc_list: list[xr.DataArray] = []

        for band_coord, product_band in requested:
            param_tiles: list[tuple[xr.DataArray, xr.DataArray, xr.DataArray]] = []
            qa_tiles: list[xr.DataArray] = []
            for path in paths:
                params, qa = self._load_native_band_stack(path, product_band)
                qa_tiles.append(qa)
                param_tiles.append(params)

            f0 = self._merge_reprojected_tiles(
                [params[0] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            f1 = self._merge_reprojected_tiles(
                [params[1] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            f2 = self._merge_reprojected_tiles(
                [params[2] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            qa = self._merge_reprojected_tiles(
                qa_tiles,
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
            )
            f0 = f0.fillna(0.20)
            f1 = f1.fillna(0.05)
            f2 = f2.fillna(0.02)
            unc = self._qa_to_uncertainty(qa).fillna(0.08)

            f0_list.append(f0.expand_dims(band=[band_coord]))
            f1_list.append(f1.expand_dims(band=[band_coord]))
            f2_list.append(f2.expand_dims(band=[band_coord]))
            f0_unc_list.append(unc.expand_dims(band=[band_coord]))
            f1_unc_list.append((unc * 1.1).expand_dims(band=[band_coord]))
            f2_unc_list.append((unc * 1.1).expand_dims(band=[band_coord]))

        return BRDFKernelWeights(
            f0=xr.concat(f0_list, dim="band").transpose("band", "y", "x"),
            f1=xr.concat(f1_list, dim="band").transpose("band", "y", "x"),
            f2=xr.concat(f2_list, dim="band").transpose("band", "y", "x"),
            f0_unc=xr.concat(f0_unc_list, dim="band").transpose("band", "y", "x"),
            f1_unc=xr.concat(f1_unc_list, dim="band").transpose("band", "y", "x"),
            f2_unc=xr.concat(f2_unc_list, dim="band").transpose("band", "y", "x"),
        )

    @staticmethod
    def _merge_reprojected_tiles(
        arrays: list[xr.DataArray],
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        resampling: Resampling,
        nodata: float | None,
    ) -> xr.DataArray:
        if not arrays:
            raise ValueError("Expected at least one array to merge")

        reprojected = [
            reproject_native_to_target(
                arr,
                target_bounds=bounds,
                target_crs=crs,
                target_resolution=target_resolution,
                resampling=resampling,
                nodata=nodata,
            )
            for arr in arrays
        ]

        merged = reprojected[0]
        for arr in reprojected[1:]:
            merged = merged.combine_first(arr)
        return merged.astype(np.float32)

    @staticmethod
    def _qa_to_uncertainty(qa: xr.DataArray) -> xr.DataArray:
        qa_values = qa.values.astype(np.float32)
        unc = np.full(qa_values.shape, 0.08, dtype=np.float32)
        unc = np.where(qa_values == 0, 0.03, unc)
        unc = np.where(qa_values == 1, 0.05, unc)
        unc = np.where(np.isfinite(qa_values), unc, np.nan)
        return xr.DataArray(unc, dims=qa.dims, coords=qa.coords)

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        raise NotImplementedError

    @staticmethod
    def _grid(bounds: tuple[float, float, float, float], resolution: float) -> tuple[np.ndarray, np.ndarray]:
        xmin, ymin, xmax, ymax = bounds
        if resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {resolution}")

        nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
        ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
        x = xmin + (np.arange(nx, dtype=np.float32) + 0.5) * resolution
        y = ymax - (np.arange(ny, dtype=np.float32) + 0.5) * resolution
        return y, x

    def _constant_band_array(
        self,
        bands: list[int | str],
        bounds: tuple[float, float, float, float],
        resolution: float,
        value: float,
    ) -> xr.DataArray:
        y, x = self._grid(bounds, resolution)
        arr = np.full((len(bands), y.size, x.size), value, dtype=np.float32)
        return xr.DataArray(
            arr,
            dims=["band", "y", "x"],
            coords={"band": np.array(bands, dtype=object), "y": y, "x": x},
        )

    def _default_weights(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        bands: list[int | str],
    ) -> BRDFKernelWeights:
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        f0 = self._constant_band_array(bands, bounds, resolution, 0.20)
        f1 = self._constant_band_array(bands, bounds, resolution, 0.05)
        f2 = self._constant_band_array(bands, bounds, resolution, 0.02)

        return BRDFKernelWeights(
            f0=f0,
            f1=f1,
            f2=f2,
            f0_unc=xr.full_like(f0, 0.03),
            f1_unc=xr.full_like(f1, 0.02),
            f2_unc=xr.full_like(f2, 0.02),
        )


class _StackParameterProvider(_EarthAccessBRDFProvider):
    """BRDF provider where one dataset stores all three kernel parameters."""

    _read_dataset = staticmethod(read_hdf4_dataset)

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        params_raw, params_attrs = self._read_dataset(path, product_band.parameter_dataset)
        qa_raw, qa_attrs = self._read_dataset(path, product_band.qa_dataset or "")

        params = apply_scale_and_mask(params_raw, params_attrs)
        if params.ndim != 3 or params.shape[-1] != 3:
            raise ValueError(
                f"Expected a 3-parameter BRDF stack for {product_band.parameter_dataset}, "
                f"got shape={params.shape}"
            )

        qa = apply_scale_and_mask(qa_raw, qa_attrs)
        param_da = make_native_grid_dataarray(
            np.moveaxis(params, -1, 0),
            granule_path=path,
            dims=("parameter", "y", "x"),
            coords={"parameter": ["f0", "f1", "f2"]},
        )
        qa_da = make_native_grid_dataarray(qa, granule_path=path)
        return (
            param_da.sel(parameter="f0", drop=True),
            param_da.sel(parameter="f1", drop=True),
            param_da.sel(parameter="f2", drop=True),
        ), qa_da


class MCD43EarthAccessProvider(_StackParameterProvider):
    """MODIS MCD43 BRDF provider using real MCD43A1 kernel parameters."""

    product_key = "mcd43_brdf"
    _source_name = "MCD43"
    _product_bands = (
        ProductBandDefinition("Band1", 645.0, "BRDF_Albedo_Parameters_Band1", "BRDF_Albedo_Band_Mandatory_Quality_Band1"),
        ProductBandDefinition("Band2", 858.5, "BRDF_Albedo_Parameters_Band2", "BRDF_Albedo_Band_Mandatory_Quality_Band2"),
        ProductBandDefinition("Band3", 469.0, "BRDF_Albedo_Parameters_Band3", "BRDF_Albedo_Band_Mandatory_Quality_Band3"),
        ProductBandDefinition("Band4", 555.0, "BRDF_Albedo_Parameters_Band4", "BRDF_Albedo_Band_Mandatory_Quality_Band4"),
        ProductBandDefinition("Band5", 1240.0, "BRDF_Albedo_Parameters_Band5", "BRDF_Albedo_Band_Mandatory_Quality_Band5"),
        ProductBandDefinition("Band6", 1640.0, "BRDF_Albedo_Parameters_Band6", "BRDF_Albedo_Band_Mandatory_Quality_Band6"),
        ProductBandDefinition("Band7", 2130.0, "BRDF_Albedo_Parameters_Band7", "BRDF_Albedo_Band_Mandatory_Quality_Band7"),
    )
    _legacy_band_map = {index + 1: band.label for index, band in enumerate(_product_bands)}


class VNP43EarthAccessProvider(_StackParameterProvider):
    """VIIRS VNP43 BRDF provider using real VNP43MA1 kernel parameters."""

    product_key = "vnp43_brdf"
    _source_name = "VNP43"
    _read_dataset = staticmethod(read_hdf5_dataset)
    _product_bands = (
        ProductBandDefinition(
            "M1",
            412.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M1",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M1",
        ),
        ProductBandDefinition(
            "M2",
            445.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M2",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M2",
        ),
        ProductBandDefinition(
            "M3",
            488.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M3",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M3",
        ),
        ProductBandDefinition(
            "M4",
            555.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M4",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M4",
        ),
        ProductBandDefinition(
            "M5",
            672.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M5",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M5",
        ),
        ProductBandDefinition(
            "M7",
            865.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M7",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M7",
        ),
        ProductBandDefinition(
            "M8",
            1240.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M8",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M8",
        ),
        ProductBandDefinition(
            "M10",
            1610.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M10",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M10",
        ),
        ProductBandDefinition(
            "M11",
            2250.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M11",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M11",
        ),
    )
    _legacy_band_map = {
        1: "M5",
        2: "M7",
        3: "M3",
        4: "M4",
        5: "M8",
        6: "M10",
        7: "M11",
    }


class MCD19EarthAccessProvider(_EarthAccessBRDFProvider):
    """MODIS MAIAC BRDF provider using real MCD19A3 kernel parameters."""

    product_key = "mcd19_brdf"
    _source_name = "MCD19"
    _product_bands = (
        ProductBandDefinition("Band1", 645.0, "Kiso_Band1", "Status_QA"),
        ProductBandDefinition("Band2", 858.5, "Kiso_Band2", "Status_QA"),
        ProductBandDefinition("Band3", 469.0, "Kiso_Band3", "Status_QA"),
        ProductBandDefinition("Band4", 555.0, "Kiso_Band4", "Status_QA"),
        ProductBandDefinition("Band5", 1240.0, "Kiso_Band5", "Status_QA"),
        ProductBandDefinition("Band6", 1640.0, "Kiso_Band6", "Status_QA"),
        ProductBandDefinition("Band7", 2130.0, "Kiso_Band7", "Status_QA"),
        ProductBandDefinition("Band8", 412.0, "Kiso_Band8", "Status_QA"),
    )
    _legacy_band_map = {index + 1: band.label for index, band in enumerate(_product_bands)}

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        label = product_band.label.replace("Band", "")
        f0_raw, f0_attrs = read_hdf4_dataset(path, f"Kiso_Band{label}")
        f1_raw, f1_attrs = read_hdf4_dataset(path, f"Kvol_Band{label}")
        f2_raw, f2_attrs = read_hdf4_dataset(path, f"Kgeo_Band{label}")
        qa_raw, qa_attrs = read_hdf4_dataset(path, "Status_QA")

        f0 = make_native_grid_dataarray(apply_scale_and_mask(f0_raw, f0_attrs), granule_path=path)
        f1 = make_native_grid_dataarray(apply_scale_and_mask(f1_raw, f1_attrs), granule_path=path)
        f2 = make_native_grid_dataarray(apply_scale_and_mask(f2_raw, f2_attrs), granule_path=path)
        qa = make_native_grid_dataarray(apply_scale_and_mask(qa_raw, qa_attrs), granule_path=path)
        return (f0, f1, f2), qa

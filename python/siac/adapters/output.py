"""Configured output writing for SIAC correction products."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.storage.product_writers import write_dataset, write_rgb_quicklook
from siac.storage.raster_writers import write_cog, write_netcdf, write_raster, write_zarr

if TYPE_CHECKING:
    from siac.config.schema import OutputDefaultsConfig
    from siac.runtime import CorrectionResult


def _uint16_encode(data: xr.DataArray, *, scale: float, nodata: float) -> xr.DataArray:
    scaled = data * scale
    filled = scaled.where(np.isfinite(scaled), other=nodata)
    clipped = filled.clip(min=0, max=np.iinfo(np.uint16).max)
    return clipped.round().astype(np.uint16)


def _cast_dataset(
    dataset: xr.Dataset,
    *,
    dtype: str,
    scale: float,
    nodata: float,
) -> xr.Dataset:
    if dtype == "uint16":
        return xr.Dataset(
            {name: _uint16_encode(dataset[name], scale=scale, nodata=nodata) for name in dataset.data_vars},
            coords=dataset.coords,
            attrs=dataset.attrs,
        )
    return dataset.astype(dtype)


@dataclass(frozen=True)
class ConfiguredOutputWriter:
    """Write correction products using the resolved output configuration."""

    defaults: OutputDefaultsConfig

    def write(
        self,
        result: CorrectionResult,
        output_dir: str | Path,
    ) -> dict[str, Path]:
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)

        if self.defaults.format in {"geotiff", "cog"}:
            return self._write_raster_products(result, destination)
        if self.defaults.format == "netcdf":
            return self._write_netcdf_products(result, destination)
        if self.defaults.format == "zarr":
            return self._write_zarr_products(result, destination)
        raise ValueError(f"Unsupported output format: {self.defaults.format!r}")

    def _write_raster_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = {}
        raster_dir = output_dir / "boa"
        prepared_boa = _cast_dataset(
            result.boa,
            dtype=self.defaults.boa_dtype,
            scale=self.defaults.boa_scale,
            nodata=self.defaults.boa_nodata,
        )
        nodata = int(self.defaults.boa_nodata) if self.defaults.boa_dtype == "uint16" else None
        boa_paths = write_dataset(
            prepared_boa,
            raster_dir,
            compression=self.defaults.compression,
            dtype=self.defaults.boa_dtype,
            as_cog=self.defaults.format == "cog",
            nodata=nodata,
        )
        artifacts.update({f"boa.{name}": path for name, path in boa_paths.items()})

        if self.defaults.include_uncertainty and result.boa_unc is not None:
            unc_dir = output_dir / "boa_unc"
            prepared_unc = _cast_dataset(
                result.boa_unc,
                dtype=self.defaults.boa_dtype,
                scale=self.defaults.boa_scale,
                nodata=self.defaults.boa_nodata,
            )
            unc_paths = write_dataset(
                prepared_unc,
                unc_dir,
                compression=self.defaults.compression,
                dtype=self.defaults.boa_dtype,
                as_cog=self.defaults.format == "cog",
                nodata=nodata,
            )
            artifacts.update({f"boa_unc.{name}": path for name, path in unc_paths.items()})

        if self.defaults.include_auxiliary:
            write_fn = write_cog if self.defaults.format == "cog" else write_raster
            aux_dir = output_dir / "auxiliary"
            aux_dir.mkdir(parents=True, exist_ok=True)
            artifacts["auxiliary.aot"] = write_fn(result.aot.astype(np.float32), aux_dir / "aot.tif")
            artifacts["auxiliary.tcwv"] = write_fn(result.tcwv.astype(np.float32), aux_dir / "tcwv.tif")
            cloud_mask = result.cloud_mask.astype(np.uint8)
            artifacts["auxiliary.cloud_mask"] = write_fn(
                cloud_mask,
                aux_dir / "cloud_mask.tif",
                compression="lzw",
                dtype="uint8",
                nodata=255,
            )

        quicklook = self._write_rgb_if_available(result, output_dir)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        return artifacts

    def _write_netcdf_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = {}
        prepared_boa = _cast_dataset(
            result.boa,
            dtype="float64" if self.defaults.boa_dtype == "float64" else "float32",
            scale=self.defaults.boa_scale,
            nodata=self.defaults.boa_nodata,
        )
        artifacts["boa"] = write_netcdf(prepared_boa, output_dir / "boa.nc")
        if self.defaults.include_uncertainty and result.boa_unc is not None:
            artifacts["boa_unc"] = write_netcdf(result.boa_unc.astype(np.float32), output_dir / "boa_unc.nc")
        if self.defaults.include_auxiliary:
            aux_ds = xr.Dataset(
                {
                    "aot": result.aot.astype(np.float32),
                    "tcwv": result.tcwv.astype(np.float32),
                    "cloud_mask": result.cloud_mask.astype(np.uint8),
                }
            )
            artifacts["auxiliary"] = write_netcdf(aux_ds, output_dir / "auxiliary.nc")
        quicklook = self._write_rgb_if_available(result, output_dir)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        return artifacts

    def _write_zarr_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = {}
        prepared_boa = _cast_dataset(
            result.boa,
            dtype="float64" if self.defaults.boa_dtype == "float64" else "float32",
            scale=self.defaults.boa_scale,
            nodata=self.defaults.boa_nodata,
        )
        artifacts["boa"] = write_zarr(prepared_boa, output_dir / "boa.zarr")
        if self.defaults.include_uncertainty and result.boa_unc is not None:
            artifacts["boa_unc"] = write_zarr(result.boa_unc.astype(np.float32), output_dir / "boa_unc.zarr")
        if self.defaults.include_auxiliary:
            aux_ds = xr.Dataset(
                {
                    "aot": result.aot.astype(np.float32),
                    "tcwv": result.tcwv.astype(np.float32),
                    "cloud_mask": result.cloud_mask.astype(np.uint8),
                }
            )
            artifacts["auxiliary"] = write_zarr(aux_ds, output_dir / "auxiliary.zarr")
        quicklook = self._write_rgb_if_available(result, output_dir)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        return artifacts

    def _write_rgb_if_available(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> Path | None:
        if not self.defaults.include_rgb:
            return None
        if not {"B04", "B03", "B02"} <= set(result.boa.data_vars):
            return None
        return cast("Path", write_rgb_quicklook(result.boa, output_dir / "quicklook.tif"))


__all__ = ["ConfiguredOutputWriter"]

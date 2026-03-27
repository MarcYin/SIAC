"""Configured output writing for SIAC correction products."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr
from scipy.ndimage import maximum_filter, zoom

from siac.runtime.models import copy_spatial_metadata_like
from siac.storage.product_writers import write_aot_scatter_plot, write_dataset, write_rgb_quicklook
from siac.storage.raster_writers import write_cog, write_raster, write_zarr
from siac.storage.writers import write_netcdf

if TYPE_CHECKING:
    from siac.config.schema import OutputDefaultsConfig
    from siac.runtime import CorrectionResult

logger = logging.getLogger(__name__)


def _shares_grid(field: xr.DataArray, template: xr.DataArray) -> bool:
    if tuple(field.shape) != tuple(template.shape):
        return False
    if tuple(field.dims) != tuple(template.dims):
        return False
    for dim in template.dims:
        if dim not in field.coords or dim not in template.coords:
            return False
        if not np.array_equal(np.asarray(field.coords[dim].values), np.asarray(template.coords[dim].values), equal_nan=True):
            return False
    return True


def _axis_resolution(values: np.ndarray) -> float | None:
    if values.size <= 1:
        return None
    diffs = np.abs(np.diff(np.asarray(values, dtype=np.float64)))
    finite = diffs[np.isfinite(diffs) & (diffs > 0.0)]
    if finite.size == 0:
        return None
    return float(np.median(finite))


def _cloud_mask_on_template_grid(mask: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    if _shares_grid(mask, template):
        return mask.astype(np.uint8)

    if (
        len(mask.dims) == 2
        and mask.dims == template.dims
        and all(dim in mask.coords for dim in mask.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        factor_y = 1
        factor_x = 1
        src_y_res = _axis_resolution(np.asarray(mask.coords[mask.dims[0]].values))
        src_x_res = _axis_resolution(np.asarray(mask.coords[mask.dims[1]].values))
        dst_y_res = _axis_resolution(np.asarray(template.coords[template.dims[0]].values))
        dst_x_res = _axis_resolution(np.asarray(template.coords[template.dims[1]].values))
        if src_y_res is not None and dst_y_res is not None and dst_y_res > src_y_res:
            factor_y = max(1, int(np.ceil(dst_y_res / src_y_res)))
        if src_x_res is not None and dst_x_res is not None and dst_x_res > src_x_res:
            factor_x = max(1, int(np.ceil(dst_x_res / src_x_res)))

        source = np.asarray(mask.values, dtype=np.float32)
        if factor_y > 1 or factor_x > 1:
            source = maximum_filter(source, size=(factor_y, factor_x), mode="nearest")
        remapped = xr.DataArray(
            source,
            dims=mask.dims,
            coords={dim: mask.coords[dim] for dim in mask.dims},
            attrs=mask.attrs,
        )
        try:
            interpolated = remapped.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="nearest",
            )
            out = np.asarray(interpolated.values, dtype=np.float32)
            out = np.where(np.isfinite(out), out, 1.0)
        except Exception:
            out = np.asarray(mask.values, dtype=np.float32)
    else:
        out = np.asarray(mask.values, dtype=np.float32)

    h_out = int(template.sizes[template.dims[0]])
    w_out = int(template.sizes[template.dims[1]])
    if out.shape != (h_out, w_out):
        factor_y = max(1, out.shape[0] // h_out)
        factor_x = max(1, out.shape[1] // w_out)
        dilated = maximum_filter(out, size=(factor_y, factor_x), mode="nearest")
        out = zoom(dilated, (h_out / dilated.shape[0], w_out / dilated.shape[1]), order=0)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded = np.ones((h_out, w_out), dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    resampled = xr.DataArray(
        (out > 0.5).astype(np.uint8),
        dims=template.dims,
        coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
        attrs=mask.attrs,
    )
    return copy_spatial_metadata_like(resampled, template)


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


def _cast_monthly_composite_reflectance(
    dataset: xr.Dataset,
    *,
    dtype: str,
    scale: float,
    nodata: float,
) -> xr.Dataset:
    return _cast_dataset(dataset, dtype=dtype, scale=scale, nodata=nodata)


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
            as_cog = self.defaults.format == "cog"
            return self._write_raster_products(result, destination, as_cog=as_cog)
        if self.defaults.format == "netcdf":
            return self._write_netcdf_products(result, destination)
        if self.defaults.format == "zarr":
            return self._write_zarr_products(result, destination)
        raise ValueError(f"Unsupported output format: {self.defaults.format!r}")

    def _write_raster_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
        *,
        as_cog: bool,
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
            as_cog=as_cog,
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
                as_cog=as_cog,
                nodata=nodata,
            )
            artifacts.update({f"boa_unc.{name}": path for name, path in unc_paths.items()})

        if result.surface_prior is not None:
            prior_dir = output_dir / "surface_prior"
            prepared_prior = _cast_dataset(
                result.surface_prior,
                dtype=self.defaults.boa_dtype,
                scale=self.defaults.boa_scale,
                nodata=self.defaults.boa_nodata,
            )
            prior_paths = write_dataset(
                prepared_prior,
                prior_dir,
                compression=self.defaults.compression,
                dtype=self.defaults.boa_dtype,
                as_cog=as_cog,
                nodata=nodata,
            )
            artifacts.update({f"surface_prior.{name}": path for name, path in prior_paths.items()})

        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            prior_unc_dir = output_dir / "surface_prior_unc"
            prepared_prior_unc = _cast_dataset(
                result.surface_prior_unc,
                dtype=self.defaults.boa_dtype,
                scale=self.defaults.boa_scale,
                nodata=self.defaults.boa_nodata,
            )
            prior_unc_paths = write_dataset(
                prepared_prior_unc,
                prior_unc_dir,
                compression=self.defaults.compression,
                dtype=self.defaults.boa_dtype,
                as_cog=as_cog,
                nodata=nodata,
            )
            artifacts.update({f"surface_prior_unc.{name}": path for name, path in prior_unc_paths.items()})

        if result.monthly_composites:
            write_fn = write_cog if as_cog else write_raster
            monthly_root = output_dir / "monthly_composites"
            monthly_root.mkdir(parents=True, exist_ok=True)
            for label, composite in sorted(result.monthly_composites.items()):
                composite_dir = monthly_root / label
                prepared_reflectance = _cast_monthly_composite_reflectance(
                    composite.reflectance,
                    dtype=self.defaults.boa_dtype,
                    scale=self.defaults.boa_scale,
                    nodata=self.defaults.boa_nodata,
                )
                composite_paths = write_dataset(
                    prepared_reflectance,
                    composite_dir,
                    compression=self.defaults.compression,
                    dtype=self.defaults.boa_dtype,
                    as_cog=as_cog,
                    nodata=nodata,
                )
                artifacts.update({f"monthly_composites.{label}.{name}": path for name, path in composite_paths.items()})
                artifacts[f"monthly_composites.{label}.quality"] = write_fn(
                    composite.quality.astype(np.float32),
                    composite_dir / "quality.tif",
                )
                artifacts[f"monthly_composites.{label}.sample_index"] = write_fn(
                    composite.sample_index.astype(np.int16),
                    composite_dir / "sample_index.tif",
                    compression="lzw",
                    dtype="int16",
                    nodata=-1,
                )

        if self.defaults.include_auxiliary:
            write_fn = write_cog if as_cog else write_raster
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
        artifacts.update(self._write_scatter_diagnostics_if_available(result, output_dir))
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
        if result.surface_prior is not None:
            artifacts["surface_prior"] = write_netcdf(result.surface_prior.astype(np.float32), output_dir / "surface_prior.nc")
        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            artifacts["surface_prior_unc"] = write_netcdf(
                result.surface_prior_unc.astype(np.float32),
                output_dir / "surface_prior_unc.nc",
            )
        if result.monthly_composites:
            monthly_root = output_dir / "monthly_composites"
            monthly_root.mkdir(parents=True, exist_ok=True)
            for label, composite in sorted(result.monthly_composites.items()):
                composite_ds = composite.reflectance.copy()
                composite_ds["quality"] = composite.quality.astype(np.float32)
                composite_ds["sample_index"] = composite.sample_index.astype(np.int16)
                artifacts[f"monthly_composites.{label}"] = write_netcdf(
                    composite_ds,
                    monthly_root / f"{label}.nc",
                )
        if self.defaults.include_auxiliary:
            aux_template = result.aot
            aux_ds = xr.Dataset(
                {
                    "aot": aux_template.astype(np.float32),
                    "tcwv": result.tcwv.astype(np.float32),
                    "cloud_mask": _cloud_mask_on_template_grid(result.cloud_mask, aux_template),
                }
            )
            artifacts["auxiliary"] = write_netcdf(aux_ds, output_dir / "auxiliary.nc")
        quicklook = self._write_rgb_if_available(result, output_dir)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        artifacts.update(self._write_scatter_diagnostics_if_available(result, output_dir))
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
        if result.surface_prior is not None:
            artifacts["surface_prior"] = write_zarr(result.surface_prior.astype(np.float32), output_dir / "surface_prior.zarr")
        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            artifacts["surface_prior_unc"] = write_zarr(
                result.surface_prior_unc.astype(np.float32),
                output_dir / "surface_prior_unc.zarr",
            )
        if result.monthly_composites:
            monthly_root = output_dir / "monthly_composites"
            monthly_root.mkdir(parents=True, exist_ok=True)
            for label, composite in sorted(result.monthly_composites.items()):
                composite_ds = composite.reflectance.copy()
                composite_ds["quality"] = composite.quality.astype(np.float32)
                composite_ds["sample_index"] = composite.sample_index.astype(np.int16)
                artifacts[f"monthly_composites.{label}"] = write_zarr(
                    composite_ds,
                    monthly_root / f"{label}.zarr",
                )
        if self.defaults.include_auxiliary:
            aux_template = result.aot
            aux_ds = xr.Dataset(
                {
                    "aot": aux_template.astype(np.float32),
                    "tcwv": result.tcwv.astype(np.float32),
                    "cloud_mask": _cloud_mask_on_template_grid(result.cloud_mask, aux_template),
                }
            )
            artifacts["auxiliary"] = write_zarr(aux_ds, output_dir / "auxiliary.zarr")
        quicklook = self._write_rgb_if_available(result, output_dir)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        artifacts.update(self._write_scatter_diagnostics_if_available(result, output_dir))
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

    def _write_scatter_diagnostics_if_available(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> dict[str, Path]:
        if not self.defaults.include_auxiliary:
            return {}
        plots = tuple(getattr(result.diagnostics, "aot_scatter_plots", ()))
        if not plots:
            return {}

        diagnostics_dir = output_dir / "diagnostics"
        artifacts: dict[str, Path] = {}
        for plot in plots:
            path = diagnostics_dir / f"aot_scatter_{plot.band_name}.png"
            artifacts[f"diagnostics.scatter.{plot.band_name}"] = cast(
                "Path",
                write_aot_scatter_plot(plot, path),
            )
        return artifacts


__all__ = ["ConfiguredOutputWriter"]

"""Configured output writing for SIAC correction products.

Produces satellite-standard L2A output with naming convention::

    {satellite}_L2A_{YYYYMMDDTHHMMSS}_{product}[_{band}].ext

and a STAC Item JSON with eo:bands, view angles, and projection metadata.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.geo.resample import resample_mask_to_template
from siac.storage.product_writers import (
    write_aot_scatter_plot,
    write_cloud_mask_preview,
    write_dataset,
    write_false_colour_preview,
    write_field_preview,
    write_rgb_quicklook,
)
from siac.storage.raster_writers import write_cog, write_raster, write_zarr
from siac.storage.stac import build_stac_item_from_result
from siac.storage.writers import ensure_writable_directory, write_netcdf

if TYPE_CHECKING:
    from pathlib import Path

    from siac.config.system import OutputDefaultsConfig
    from siac.runtime import CorrectionResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _auxiliary_dataset(result: CorrectionResult) -> xr.Dataset:
    aux_template = result.aot
    resampled_mask = resample_mask_to_template(result.cloud_mask, aux_template)
    aux_vars: dict[str, xr.DataArray] = {
        "aot": aux_template.astype(np.float32),
        "tcwv": result.tcwv.astype(np.float32),
        "cloud_mask": resampled_mask.astype(np.uint8),
    }
    if result.solver_qa is not None:
        for name, field in sorted(result.solver_qa.data_vars.items()):
            if np.issubdtype(field.dtype, np.floating):
                aux_vars[name] = field.astype(np.float32)
            else:
                aux_vars[name] = resample_mask_to_template(field.astype(bool), aux_template).astype(
                    np.uint8
                )
    return xr.Dataset(aux_vars)


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
    return xr.Dataset(
        {
            name: _cast_dataarray(dataset[name], dtype=dtype, scale=scale, nodata=nodata)
            for name in dataset.data_vars
        },
        coords=dataset.coords,
        attrs=dataset.attrs,
    )


def _cast_dataarray(
    data: xr.DataArray,
    *,
    dtype: str,
    scale: float,
    nodata: float,
) -> xr.DataArray:
    if dtype == "uint16":
        return _uint16_encode(data, scale=scale, nodata=nodata)
    return data.astype(dtype)


# ---------------------------------------------------------------------------
# Scene prefix derivation
# ---------------------------------------------------------------------------


def _derive_scene_prefix(metadata: dict[str, Any]) -> str:
    """Build ``S2A_L2A_20240115T103045`` style prefix from metadata."""
    satellite = metadata.get("satellite", "")
    observation_time = metadata.get("observation_time")

    # Satellite token (e.g. "S2A", "S2B", "L8")
    sat_token = str(satellite).upper() if satellite else "SAT"

    # Timestamp token
    if isinstance(observation_time, datetime):
        ts_token = observation_time.strftime("%Y%m%dT%H%M%S")
    else:
        ts_token = "00000000T000000"

    return "_".join([sat_token, "L2A", ts_token])


# ---------------------------------------------------------------------------
# ConfiguredOutputWriter
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConfiguredOutputWriter:
    """Write correction products with satellite-standard naming and STAC metadata.

    Output layout (raster/cog format)::

        output_dir/
        ├── S2A_L2A_20240115T103045_BOA_B02.tif
        ├── S2A_L2A_20240115T103045_BOA_B03.tif
        ├── ...
        ├── S2A_L2A_20240115T103045_BOA_UNC_B02.tif   (if include_uncertainty)
        ├── S2A_L2A_20240115T103045_SURF_B02.tif      (surface prior)
        ├── S2A_L2A_20240115T103045_SURF_UNC_B02.tif
        ├── S2A_L2A_20240115T103045_AOT.tif
        ├── S2A_L2A_20240115T103045_TCWV.tif
        ├── S2A_L2A_20240115T103045_CLOUD.tif
        ├── S2A_L2A_20240115T103045_QA_{flag}.tif
        ├── S2A_L2A_20240115T103045_RGB.tif
        ├── S2A_L2A_20240115T103045.json               (STAC Item)
        └── preview/
            ├── false_colour.png        (NIR-R-G composite)
            ├── aot.png                 (colour-mapped AOT field)
            ├── tcwv.png                (colour-mapped TCWV field)
            ├── cloud_mask.png          (RGB + cloud mask overlay)
            └── aot_scatter_B02.png     (solver scatter plots)
    """

    defaults: OutputDefaultsConfig

    def write(
        self,
        result: CorrectionResult,
        output_dir: str | Path,
    ) -> dict[str, Path]:
        destination = ensure_writable_directory(output_dir)

        prefix = _derive_scene_prefix(result.metadata)

        if self.defaults.format in {"geotiff", "cog"}:
            as_cog = self.defaults.format == "cog"
            artifacts = self._write_raster_products(result, destination, prefix, as_cog=as_cog)
        elif self.defaults.format == "netcdf":
            artifacts = self._write_netcdf_products(result, destination, prefix)
        elif self.defaults.format == "zarr":
            artifacts = self._write_zarr_products(result, destination, prefix)
        else:
            raise ValueError(f"Unsupported output format: {self.defaults.format!r}")

        return self._write_stac_item(result, destination, prefix, artifacts)

    def open_correction_boa_stream(
        self,
        output_dir: str | Path,
        *,
        metadata: dict[str, Any],
    ) -> _RasterCorrectionBoaStream | None:
        """Open an inline BOA writer for the correction stage when supported."""
        if self.defaults.format not in {"geotiff", "cog"}:
            return None
        if self.defaults.boa_dtype not in {"float32", "float64"}:
            return None

        destination = ensure_writable_directory(output_dir)
        return _RasterCorrectionBoaStream(
            writer=self,
            output_dir=destination,
            prefix=_derive_scene_prefix(metadata),
            as_cog=self.defaults.format == "cog",
        )

    def _write_stac_item(
        self,
        result: CorrectionResult,
        output_dir: Path,
        prefix: str,
        artifacts: dict[str, Path],
    ) -> dict[str, Path]:
        stac_path = output_dir / f"{prefix}.json"
        stac_item = build_stac_item_from_result(result, output_dir=output_dir, artifacts=artifacts)
        stac_path.write_text(json.dumps(stac_item, indent=2), encoding="utf-8")
        artifacts["stac_item"] = stac_path
        return artifacts

    # -----------------------------------------------------------------
    # Raster (GeoTIFF / COG)
    # -----------------------------------------------------------------

    def _write_raster_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
        prefix: str,
        *,
        as_cog: bool,
        prewritten_artifacts: dict[str, Path] | None = None,
        skip_boa: bool = False,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = dict(prewritten_artifacts or {})
        write_fn = write_cog if as_cog else write_raster
        nodata = int(self.defaults.boa_nodata) if self.defaults.boa_dtype == "uint16" else None

        def _write_bands(
            dataset: xr.Dataset,
            product_tag: str,
            artifact_prefix: str,
        ) -> None:
            for band_name in dataset.data_vars:
                prepared = _cast_dataarray(
                    dataset[band_name],
                    dtype=self.defaults.boa_dtype,
                    scale=self.defaults.boa_scale,
                    nodata=self.defaults.boa_nodata,
                )
                path = output_dir / f"{prefix}_{product_tag}_{band_name}.tif"
                write_fn(prepared, path, compression=self.defaults.compression, nodata=nodata)
                artifacts[f"{artifact_prefix}.{band_name}"] = path

        if not skip_boa and result.boa.data_vars:
            _write_bands(result.boa, "BOA", "boa")
        if self.defaults.include_uncertainty and result.boa_unc is not None:
            _write_bands(result.boa_unc, "BOA_UNC", "boa_unc")
        if result.surface_prior is not None:
            _write_bands(result.surface_prior, "SURF", "surface_prior")
        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            _write_bands(result.surface_prior_unc, "SURF_UNC", "surface_prior_unc")

        # --- Monthly composites ---
        if result.monthly_composites:
            for label, composite in sorted(result.monthly_composites.items()):
                prepared_reflectance = _cast_dataset(
                    composite.reflectance,
                    dtype=self.defaults.boa_dtype,
                    scale=self.defaults.boa_scale,
                    nodata=self.defaults.boa_nodata,
                )
                composite_dir = output_dir / "monthly_composites" / label
                composite_paths = write_dataset(
                    prepared_reflectance,
                    composite_dir,
                    compression=self.defaults.compression,
                    dtype=self.defaults.boa_dtype,
                    as_cog=as_cog,
                    nodata=nodata,
                )
                artifacts.update(
                    {
                        f"monthly_composites.{label}.{name}": path
                        for name, path in composite_paths.items()
                    }
                )
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

        # --- Auxiliary (AOT, TCWV, cloud mask, QA) ---
        if self.defaults.include_auxiliary:
            aux_ds = _auxiliary_dataset(result)
            artifacts["auxiliary.aot"] = write_fn(
                aux_ds["aot"],
                output_dir / f"{prefix}_AOT.tif",
            )
            artifacts["auxiliary.tcwv"] = write_fn(
                aux_ds["tcwv"],
                output_dir / f"{prefix}_TCWV.tif",
            )
            artifacts["auxiliary.cloud_mask"] = write_fn(
                aux_ds["cloud_mask"],
                output_dir / f"{prefix}_CLOUD.tif",
                compression="lzw",
                dtype="uint8",
                nodata=255,
            )
            for name, field in aux_ds.data_vars.items():
                if name in {"aot", "tcwv", "cloud_mask"}:
                    continue
                is_float_qa = np.issubdtype(field.dtype, np.floating)
                artifacts[f"auxiliary.qa.{name}"] = write_fn(
                    field,
                    output_dir / f"{prefix}_QA_{name}.tif",
                    compression="lzw",
                    **(
                        {
                            "dtype": "float32",
                            "nodata": float("nan"),
                        }
                        if is_float_qa
                        else {
                            "dtype": "uint8",
                            "nodata": 255,
                        }
                    ),
                )

        # --- RGB quicklook ---
        quicklook = self._write_rgb_if_available(result, output_dir, prefix)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook

        # --- Preview PNGs ---
        artifacts.update(self._write_previews(result, output_dir))

        return artifacts

    def _finish_correction_boa_stream(
        self,
        result: CorrectionResult,
        stream: _RasterCorrectionBoaStream,
    ) -> dict[str, Path]:
        artifacts = self._write_raster_products(
            result,
            stream.output_dir,
            stream.prefix,
            as_cog=stream.as_cog,
            prewritten_artifacts=stream.artifacts,
            skip_boa=True,
        )
        return self._write_stac_item(result, stream.output_dir, stream.prefix, artifacts)

    # -----------------------------------------------------------------
    # NetCDF
    # -----------------------------------------------------------------

    def _write_netcdf_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
        prefix: str,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = {}
        if result.boa.data_vars:
            prepared_boa = _cast_dataset(
                result.boa,
                dtype="float64" if self.defaults.boa_dtype == "float64" else "float32",
                scale=self.defaults.boa_scale,
                nodata=self.defaults.boa_nodata,
            )
            artifacts["boa"] = write_netcdf(prepared_boa, output_dir / f"{prefix}_BOA.nc")
        if self.defaults.include_uncertainty and result.boa_unc is not None:
            artifacts["boa_unc"] = write_netcdf(
                result.boa_unc.astype(np.float32), output_dir / f"{prefix}_BOA_UNC.nc"
            )
        if result.surface_prior is not None:
            artifacts["surface_prior"] = write_netcdf(
                result.surface_prior.astype(np.float32), output_dir / f"{prefix}_SURF.nc"
            )
        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            artifacts["surface_prior_unc"] = write_netcdf(
                result.surface_prior_unc.astype(np.float32),
                output_dir / f"{prefix}_SURF_UNC.nc",
            )
        if result.monthly_composites:
            monthly_root = ensure_writable_directory(output_dir / "monthly_composites")
            for label, composite in sorted(result.monthly_composites.items()):
                composite_ds = composite.reflectance.copy()
                composite_ds["quality"] = composite.quality.astype(np.float32)
                composite_ds["sample_index"] = composite.sample_index.astype(np.int16)
                artifacts[f"monthly_composites.{label}"] = write_netcdf(
                    composite_ds,
                    monthly_root / f"{label}.nc",
                )
        if self.defaults.include_auxiliary:
            aux_ds = _auxiliary_dataset(result)
            artifacts["auxiliary"] = write_netcdf(aux_ds, output_dir / f"{prefix}_AUX.nc")
        quicklook = self._write_rgb_if_available(result, output_dir, prefix)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        artifacts.update(self._write_previews(result, output_dir))
        return artifacts

    # -----------------------------------------------------------------
    # Zarr
    # -----------------------------------------------------------------

    def _write_zarr_products(
        self,
        result: CorrectionResult,
        output_dir: Path,
        prefix: str,
    ) -> dict[str, Path]:
        artifacts: dict[str, Path] = {}
        if result.boa.data_vars:
            prepared_boa = _cast_dataset(
                result.boa,
                dtype="float64" if self.defaults.boa_dtype == "float64" else "float32",
                scale=self.defaults.boa_scale,
                nodata=self.defaults.boa_nodata,
            )
            artifacts["boa"] = write_zarr(prepared_boa, output_dir / f"{prefix}_BOA.zarr")
        if self.defaults.include_uncertainty and result.boa_unc is not None:
            artifacts["boa_unc"] = write_zarr(
                result.boa_unc.astype(np.float32), output_dir / f"{prefix}_BOA_UNC.zarr"
            )
        if result.surface_prior is not None:
            artifacts["surface_prior"] = write_zarr(
                result.surface_prior.astype(np.float32), output_dir / f"{prefix}_SURF.zarr"
            )
        if self.defaults.include_uncertainty and result.surface_prior_unc is not None:
            artifacts["surface_prior_unc"] = write_zarr(
                result.surface_prior_unc.astype(np.float32),
                output_dir / f"{prefix}_SURF_UNC.zarr",
            )
        if result.monthly_composites:
            monthly_root = ensure_writable_directory(output_dir / "monthly_composites")
            for label, composite in sorted(result.monthly_composites.items()):
                composite_ds = composite.reflectance.copy()
                composite_ds["quality"] = composite.quality.astype(np.float32)
                composite_ds["sample_index"] = composite.sample_index.astype(np.int16)
                artifacts[f"monthly_composites.{label}"] = write_zarr(
                    composite_ds,
                    monthly_root / f"{label}.zarr",
                )
        if self.defaults.include_auxiliary:
            aux_ds = _auxiliary_dataset(result)
            artifacts["auxiliary"] = write_zarr(aux_ds, output_dir / f"{prefix}_AUX.zarr")
        quicklook = self._write_rgb_if_available(result, output_dir, prefix)
        if quicklook is not None:
            artifacts["quicklook.rgb"] = quicklook
        artifacts.update(self._write_previews(result, output_dir))
        return artifacts

    # -----------------------------------------------------------------
    # Shared helpers
    # -----------------------------------------------------------------

    def _write_rgb_if_available(
        self,
        result: CorrectionResult,
        output_dir: Path,
        prefix: str,
    ) -> Path | None:
        if not self.defaults.include_rgb:
            return None
        if not {"B04", "B03", "B02"} <= set(result.boa.data_vars):
            return None
        return cast("Path", write_rgb_quicklook(result.boa, output_dir / f"{prefix}_RGB.tif"))

    def _write_previews(
        self,
        result: CorrectionResult,
        output_dir: Path,
    ) -> dict[str, Path]:
        """Write all preview PNGs: false colour, AOT/TCWV maps, cloud overlay, scatter."""
        if not self.defaults.include_rgb:
            return {}

        preview_dir = output_dir / "preview"
        artifacts: dict[str, Path] = {}

        # False-colour composite (NIR-Red-Green)
        try:
            fc_path = write_false_colour_preview(
                result.boa,
                preview_dir / "false_colour.png",
            )
            if fc_path is not None:
                artifacts["preview.false_colour"] = fc_path
        except Exception:
            logger.debug("Skipped false-colour preview.", exc_info=True)

        # AOT colour map
        try:
            aot_path = write_field_preview(
                result.aot,
                preview_dir / "aot.png",
                vmin=0.0,
                vmax=1.0,
                palette="magma",
                title="AOT 550 nm",
                cloud_mask=result.cloud_mask,
            )
            artifacts["preview.aot"] = aot_path
        except Exception:
            logger.debug("Skipped AOT preview.", exc_info=True)

        # TCWV colour map
        try:
            tcwv_path = write_field_preview(
                result.tcwv,
                preview_dir / "tcwv.png",
                vmin=0.0,
                palette="viridis",
                title="TCWV (cm)",
                unit="cm",
                cloud_mask=result.cloud_mask,
            )
            artifacts["preview.tcwv"] = tcwv_path
        except Exception:
            logger.debug("Skipped TCWV preview.", exc_info=True)

        # Cloud mask overlay on true-colour
        try:
            cloud_path = write_cloud_mask_preview(
                result.boa,
                result.cloud_mask,
                preview_dir / "cloud_mask.png",
            )
            if cloud_path is not None:
                artifacts["preview.cloud_mask"] = cloud_path
        except Exception:
            logger.debug("Skipped cloud mask preview.", exc_info=True)

        # Scatter diagnostics
        plots = tuple(getattr(result.diagnostics, "aot_scatter_plots", ()))
        for plot in plots:
            try:
                path = preview_dir / f"aot_scatter_{plot.band_name}.png"
                artifacts[f"preview.scatter.{plot.band_name}"] = cast(
                    "Path",
                    write_aot_scatter_plot(plot, path),
                )
            except Exception:
                logger.debug("Skipped scatter plot for %s.", plot.band_name, exc_info=True)

        return artifacts


@dataclass
class _RasterCorrectionBoaStream:
    """Write BOA COG/GeoTIFF bands during correction and finish the product later."""

    writer: ConfiguredOutputWriter
    output_dir: Path
    prefix: str
    as_cog: bool
    artifacts: dict[str, Path] = field(default_factory=dict)

    @property
    def has_written(self) -> bool:
        return bool(self.artifacts)

    def write_boa_band(self, band_name: str, data: xr.DataArray) -> xr.DataArray:
        write_fn = write_cog if self.as_cog else write_raster
        nodata = (
            int(self.writer.defaults.boa_nodata)
            if self.writer.defaults.boa_dtype == "uint16"
            else None
        )
        prepared = _cast_dataarray(
            data,
            dtype=self.writer.defaults.boa_dtype,
            scale=self.writer.defaults.boa_scale,
            nodata=self.writer.defaults.boa_nodata,
        )
        path = self.output_dir / f"{self.prefix}_BOA_{band_name}.tif"
        write_fn(prepared, path, compression=self.writer.defaults.compression, nodata=nodata)
        self.artifacts[f"boa.{band_name}"] = path

        if not self.writer.defaults.reopen_streamed_boa:
            prepared.name = band_name
            prepared.attrs.update(data.attrs)
            return prepared

        from siac.storage.readers import read_raster

        reopened = read_raster(path, masked=True)
        if self.writer.defaults.boa_dtype in {"float32", "float64"}:
            reopened = reopened.astype(self.writer.defaults.boa_dtype)
        reopened.name = band_name
        reopened.attrs.update(data.attrs)
        return reopened

    def finish(self, result: CorrectionResult) -> dict[str, Path]:
        return self.writer._finish_correction_boa_stream(result, self)

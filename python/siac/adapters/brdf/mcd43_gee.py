"""BRDF provider sourcing MCD43A1 kernel weights from Google Earth Engine.

This provider downloads MODIS ``MCD43A1`` BRDF/Albedo model parameters from the
``MODIS/061/MCD43A1`` Earth Engine collection via the ``edown`` tool, then
composes the same :class:`BRDFKernelWeights` payload the Earthaccess provider
produces. Only the *source* of the kernels changes; the downstream surface-prior
composition (kernel_model / whittaker / monthly_database) is untouched, so it
stays valid for cross-approach comparison.

``edown`` writes one AOI-mosaicked GeoTIFF per date with bands named by the
Earth Engine band id (``BRDF_Albedo_Parameters_Band{N}_{iso,vol,geo}`` and
``BRDF_Albedo_Band_Mandatory_Quality_Band{N}``), int16 scaled by 0.001, nodata
32767, in the MODIS sinusoidal CRS. Because edown already mosaics the AOI, this
provider skips the multi-tile sinusoidal merge the HDF path needs and reprojects
each date directly onto the requested target grid.
"""

from __future__ import annotations

import logging
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
from shlex import quote
from shutil import which
from typing import TYPE_CHECKING, Any

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr

from siac.adapters.brdf._mcd43_qa import qa_to_uncertainty
from siac.adapters.brdf._product_specs import MCD43_PRODUCT_BANDS
from siac.adapters.earthdata_common import ProductBandDefinition, build_target_template
from siac.domain.sensors import SensorBand
from siac.geo.reprojection import reproject_match
from siac.runtime import BRDFKernelWeights

if TYPE_CHECKING:
    from collections.abc import Sequence

logger = logging.getLogger(__name__)

#: Default path to the ``edown`` executable. edown ships in a separate
#: environment from SIAC (it pulls in earthengine-api), so the provider shells
#: out rather than importing it. Override via ``providers.brdf.data_path`` ->
#: ``edown_executable`` in config or the constructor.
DEFAULT_EDOWN_EXECUTABLE = "/home/users/marcyin/.pixi/envs/base/bin/edown"
REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_EDOWN_RUNTIME = REPO_ROOT / "tools" / "edown_runtime"

GEE_COLLECTION_ID = "MODIS/061/MCD43A1"
#: MCD43A1 stores BRDF parameters as int16 scaled by 1000.
_PARAMETER_SCALE = 0.001
_PARAMETER_FILL = 32767
_KERNEL_SUFFIXES = ("iso", "vol", "geo")


def _parameter_band_id(label: str, kernel: str) -> str:
    return f"BRDF_Albedo_Parameters_{label}_{kernel}"


def _quality_band_id(label: str) -> str:
    return f"BRDF_Albedo_Band_Mandatory_Quality_{label}"


class MCD43GEEProvider:
    """MCD43A1 BRDF provider backed by Earth Engine downloads (via edown)."""

    _product_bands: tuple[ProductBandDefinition, ...] = MCD43_PRODUCT_BANDS
    _source_name = "MCD43 (GEE)"
    # Must match the Earthaccess MCD43 provider so the spectral-library RSRF
    # lookup resolves the source sensor identically.
    _rsrf_sensor_unit_id = "terra_modis"
    _rsrf_representation_variant = "band_average"

    def __init__(
        self,
        *,
        cache_dir: str | Path,
        edown_executable: str | Path = DEFAULT_EDOWN_EXECUTABLE,
        download_timeout_s: float = 1800.0,
    ):
        self.cache_dir = Path(cache_dir).expanduser()
        self.edown_executable = str(edown_executable)
        self.download_timeout_s = download_timeout_s

    @property
    def edown_executable(self) -> str:
        return self._edown_executable

    @edown_executable.setter
    def edown_executable(self, value: str) -> None:
        self._edown_executable = str(value)
        self._edown_executables = self._resolve_edown_executables(self._edown_executable)

    @property
    def source_name(self) -> str:
        return self._source_name

    @property
    def source_bands(self) -> tuple[SensorBand, ...]:
        return tuple(
            SensorBand(
                name=band.label,
                center_wavelength=band.wavelength_nm,
                bandwidth=band.bandwidth_nm,
                resolution=500.0,
                band_index=index,
                rsrf_sensor_unit_id=self._rsrf_sensor_unit_id,
                rsrf_representation_variant=self._rsrf_representation_variant,
                rsrf_band_id=band.rsrf_band_id or band.label,
            )
            for index, band in enumerate(self._product_bands)
        )

    # ── public BRDF contract ──────────────────────────────────────────

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[SensorBand],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights:
        """Single best-quality BRDF composite over the temporal window."""
        requested = self._resolve_requested_bands(list(bands))
        template = build_target_template(bounds, crs, target_resolution)
        per_date = self._download_and_read(
            bounds, crs, obs_time, temporal_window, requested, template
        )
        if not per_date:
            return self._default_weights(requested, template)
        composite = self._composite_nearest_valid(per_date, obs_time)
        return self._weights_from_band_layers(composite, requested, template)

    def get_temporal_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[SensorBand],
        temporal_window: int = 16,
        sample_dates: Sequence[Any] | None = None,
    ) -> BRDFKernelWeights:
        """Temporal BRDF stack with dims ``(time, band, y, x)``."""
        requested = self._resolve_requested_bands(list(bands))
        template = build_target_template(bounds, crs, target_resolution)
        per_date = self._download_and_read(
            bounds, crs, obs_time, temporal_window, requested, template
        )
        time_axis = (
            np.asarray([np.datetime64(d, "D") for d in sample_dates])
            if sample_dates is not None
            else self._time_axis(obs_time, temporal_window)
        )
        return self._temporal_weights_from_layers(per_date, requested, template, time_axis)

    def get_temporal_brdf_parameters_batch(
        self,
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_times: Sequence[datetime],
        target_resolution: float,
        bands: Sequence[SensorBand],
        temporal_windows: Sequence[int],
        sample_date_sets: Sequence[Sequence[Any] | None],
    ) -> list[BRDFKernelWeights]:
        """Multiple temporal BRDF stacks (one per request).

        The Earthaccess provider batches granule downloads here to dedupe; the
        GEE provider's per-window edown downloads are content-addressed and
        cached, so a straightforward per-request loop reuses overlapping data
        without a bespoke batch path.
        """
        if target_resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {target_resolution}")
        if not bands:
            raise ValueError("bands must be a non-empty sequence")
        if not (len(obs_times) == len(temporal_windows) == len(sample_date_sets)):
            raise ValueError(
                "obs_times, temporal_windows, and sample_date_sets must have the same length"
            )
        return [
            self.get_temporal_brdf_parameters(
                bounds,
                crs,
                obs_time,
                target_resolution,
                bands,
                temporal_window=temporal_window,
                sample_dates=sample_dates,
            )
            for obs_time, temporal_window, sample_dates in zip(
                obs_times, temporal_windows, sample_date_sets, strict=True
            )
        ]

    # ── band resolution ───────────────────────────────────────────────

    def _resolve_requested_bands(
        self, bands: Sequence[SensorBand]
    ) -> list[tuple[str, ProductBandDefinition]]:
        if not bands:
            raise ValueError("bands must be a non-empty sequence")
        resolved: list[tuple[str, ProductBandDefinition]] = []
        for band in bands:
            if not isinstance(band, SensorBand):
                raise TypeError(f"Bands must be SensorBand instances, got {type(band).__name__}.")
            match = min(
                self._product_bands,
                key=lambda candidate: abs(candidate.wavelength_nm - band.center_wavelength),
            )
            resolved.append((band.name, match))
        return resolved

    # ── edown download + read ─────────────────────────────────────────

    def _download_and_read(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        temporal_window: int,
        requested: Sequence[tuple[str, ProductBandDefinition]],
        template: xr.DataArray,
    ) -> dict[np.datetime64, dict[str, xr.DataArray]]:
        """Download the window and return ``{date: {band_label: layers}}``.

        ``layers`` is a DataArray with a ``parameter`` dim (f0/f1/f2) plus a
        ``unc`` companion stored under ``.attrs['unc']``, both on the target
        grid.
        """
        from pyproj import Transformer

        labels = [definition.label for _coord, definition in requested]
        band_ids: list[str] = []
        for label in labels:
            band_ids.extend(_parameter_band_id(label, k) for k in _KERNEL_SUFFIXES)
            band_ids.append(_quality_band_id(label))

        # edown expects a WGS84 bbox; transform the (possibly UTM) target bounds.
        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        xs, ys = transformer.transform(
            [bounds[0], bounds[2], bounds[0], bounds[2]],
            [bounds[1], bounds[1], bounds[3], bounds[3]],
        )
        wgs84_bbox = (min(xs), min(ys), max(xs), max(ys))

        start = (obs_time - timedelta(days=temporal_window)).date()
        end = (obs_time + timedelta(days=temporal_window)).date()
        out_root = self._download_window(wgs84_bbox, start, end, band_ids)
        if out_root is None:
            return {}

        per_date: dict[np.datetime64, dict[str, xr.DataArray]] = {}
        tiffs = sorted(out_root.glob("images/*/*.tif"))
        for tiff in tiffs:
            date = self._date_from_tiff_name(tiff.name)
            if date is None:
                continue
            try:
                layers = self._read_tiff_to_band_layers(tiff, requested, template)
            except Exception as exc:  # noqa: BLE001 - skip a corrupt date, keep the rest
                logger.warning("Skipping unreadable BRDF GeoTIFF %s: %s", tiff.name, exc)
                continue
            if layers:
                per_date[date] = layers
        logger.info(
            "%s: read %d/%d MCD43A1 dates over window %s..%s",
            self._source_name,
            len(per_date),
            len(tiffs),
            start,
            end,
        )
        return per_date

    def _download_window(
        self,
        wgs84_bbox: tuple[float, float, float, float],
        start: Any,
        end: Any,
        band_ids: Sequence[str],
    ) -> Path | None:
        """Run edown for the window into a content-addressed cache dir."""
        key = "_".join(
            [f"{c:.4f}".replace("-", "m") for c in wgs84_bbox]
            + [str(start), str(end), str(len(band_ids))]
        )
        out_root = self.cache_dir / key
        # edown resumes within output-root; treat a populated images/ dir as done.
        if any(out_root.glob("images/*/*.tif")):
            return out_root
        out_root.mkdir(parents=True, exist_ok=True)
        command_core = [
            "download",
            "--collection-id",
            GEE_COLLECTION_ID,
            "--start-date",
            str(start),
            "--end-date",
            str(end),
            "--bbox",
            *[str(c) for c in wgs84_bbox],
            "--band",
            ",".join(band_ids),
            "--output-root",
            str(out_root),
        ]
        for candidate in self._edown_executables:
            try:
                result = subprocess.run(
                    [candidate, *command_core],
                    capture_output=True,
                    text=True,
                    timeout=self.download_timeout_s,
                    check=False,
                )
            except OSError as exc:
                logger.error("edown invocation failed (%s): %s", candidate, exc)
                if any(out_root.glob("images/*/*.tif")):
                    return out_root
                if candidate != self._edown_executables[-1]:
                    logger.warning("Attempting fallback edown executable for %s", GEE_COLLECTION_ID)
                    continue
                return None
            except subprocess.TimeoutExpired as exc:
                logger.error("edown invocation timed out (%s): %s", candidate, exc)
                return None
            if result.returncode == 0 or any(out_root.glob("images/*/*.tif")):
                return out_root
            stderr = result.stderr or result.stdout or ""
            if (
                ("No module named 'edown'" in stderr)
                and candidate != self._edown_executables[-1]
            ):
                logger.warning(
                    "edown command missing module in %s; trying fallback executable: %s",
                    quote(candidate),
                    self._edown_executables[self._edown_executables.index(candidate) + 1],
                )
                continue
            logger.error(
                "edown download failed (rc=%s) for %s: %s",
                result.returncode,
                GEE_COLLECTION_ID,
                stderr.strip()[-500:],
            )
            return None
        return out_root

    @staticmethod
    def _resolve_edown_executables(edown_executable: str | Path) -> tuple[str, ...]:
        """Resolve candidate edown executable paths, with fallbacks.

        The configured path can be either the executable path, a legacy runtime
        directory, or an explicit PATH lookup string. We try the requested path
        first, then fallback to known-good repository/runtime locations.
        """
        raw = Path(edown_executable).expanduser()
        candidates: list[Path] = []

        def _add_candidate(path: Path | None) -> None:
            if path is None:
                return
            path = Path(path).expanduser()
            if not path.exists() and (found := which(str(path))) is not None:
                path = Path(found)
            if path.exists() and path.is_dir():
                path = path / "bin" / "edown"
            path = path.expanduser()
            if path.exists() and str(path) not in seen:
                seen.add(str(path))
                candidates.append(path)

        seen: set[str] = set()
        _add_candidate(raw)
        _add_candidate(LOCAL_EDOWN_RUNTIME)
        _add_candidate(LOCAL_EDOWN_RUNTIME / "bin" / "edown")
        _add_candidate(Path(DEFAULT_EDOWN_EXECUTABLE))

        if not candidates:
            # Keep the explicit path as a final attempt for transparent error handling.
            candidates.append(raw)
        return tuple(str(path) for path in candidates)

    def _read_tiff_to_band_layers(
        self,
        tiff: Path,
        requested: Sequence[tuple[str, ProductBandDefinition]],
        template: xr.DataArray,
    ) -> dict[str, xr.DataArray]:
        """Read one date's GeoTIFF into per-band f0/f1/f2 + unc on the target grid."""
        raster = rioxarray.open_rasterio(tiff, masked=False)
        if not isinstance(raster, xr.DataArray):
            raise TypeError(f"Expected a single multiband raster from {tiff}, got {type(raster)}")
        descriptions = list(raster.attrs.get("long_name", ())) or [
            raster.attrs.get(f"band_{i + 1}_name", "") for i in range(raster.sizes["band"])
        ]
        by_name = {name: idx for idx, name in enumerate(descriptions) if name}

        layers: dict[str, xr.DataArray] = {}
        for band_coord, definition in requested:
            label = definition.label
            kernel_arrays = []
            for kernel in _KERNEL_SUFFIXES:
                band_id = _parameter_band_id(label, kernel)
                if band_id not in by_name:
                    break
                raw = raster.isel(band=by_name[band_id], drop=True)
                # Mask the int16 fill, scale to reflectance, and set NaN nodata so
                # AOI areas outside the MODIS tile reproject to NaN, not 32767.
                scaled = (raw.where(raw != _PARAMETER_FILL) * _PARAMETER_SCALE).rio.write_nodata(
                    np.nan
                )
                kernel_arrays.append(
                    reproject_match(scaled, template, resampling="bilinear", nodata=np.nan)
                )
            if len(kernel_arrays) != len(_KERNEL_SUFFIXES):
                continue
            params = xr.concat(kernel_arrays, dim="parameter")
            params = params.assign_coords(parameter=["f0", "f1", "f2"])

            qa_id = _quality_band_id(label)
            if qa_id in by_name:
                qa_raw = raster.isel(band=by_name[qa_id], drop=True)
                qa_masked = qa_raw.where(qa_raw <= 4).rio.write_nodata(np.nan)
                qa = reproject_match(qa_masked, template, resampling="nearest", nodata=np.nan)
                # Scale the QA uncertainty by the band's nadir reflectance (f0)
                # so dark targets get a tight prior the aerosol solver can use,
                # instead of a flat ~0.015 absolute (~50% relative over dark land).
                unc = qa_to_uncertainty(qa, reflectance=params.sel(parameter="f0", drop=True))
            else:
                unc = xr.full_like(template, np.nan)
            params.attrs["unc"] = unc
            layers[band_coord] = params
        return layers

    # ── compositing / assembly ────────────────────────────────────────

    def _composite_nearest_valid(
        self,
        per_date: dict[np.datetime64, dict[str, xr.DataArray]],
        obs_time: datetime,
    ) -> dict[str, xr.DataArray]:
        """Per pixel, pick the temporally nearest valid observation per band."""
        obs = np.datetime64(obs_time.date(), "D")
        ordered = sorted(
            per_date, key=lambda d: abs((d - obs).astype("timedelta64[D]").astype(int))
        )
        band_coords = {coord for date in per_date for coord in per_date[date]}
        composite: dict[str, xr.DataArray] = {}
        for coord in band_coords:
            chosen_params: xr.DataArray | None = None
            chosen_unc: xr.DataArray | None = None
            for date in ordered:
                layer = per_date[date].get(coord)
                if layer is None:
                    continue
                if chosen_params is None:
                    chosen_params = layer.copy()
                    chosen_unc = layer.attrs["unc"].copy()
                    continue
                assert chosen_unc is not None  # set together with chosen_params
                # ``fill`` is (y, x): pixels still missing in the running
                # composite, to be backfilled from this (further) date.
                fill = chosen_params.isel(parameter=0, drop=True).isnull()
                chosen_params = chosen_params.where(~fill, layer)
                chosen_unc = chosen_unc.where(~fill, layer.attrs["unc"])
            if chosen_params is not None:
                chosen_params.attrs["unc"] = chosen_unc
                composite[coord] = chosen_params
        return composite

    def _weights_from_band_layers(
        self,
        composite: dict[str, xr.DataArray],
        requested: Sequence[tuple[str, ProductBandDefinition]],
        template: xr.DataArray,
    ) -> BRDFKernelWeights:
        band_coords = [coord for coord, _definition in requested]
        params_stack = []
        unc_stack = []
        for coord in band_coords:
            layer = composite.get(coord)
            if layer is None:
                nan = xr.full_like(template, np.nan)
                params_stack.append(
                    xr.concat([nan, nan, nan], dim="parameter").assign_coords(
                        parameter=["f0", "f1", "f2"]
                    )
                )
                unc_stack.append(nan)
            else:
                params_stack.append(layer)
                unc_stack.append(layer.attrs["unc"])
        params = xr.concat(params_stack, dim="band").assign_coords(band=band_coords)
        unc = xr.concat(unc_stack, dim="band").assign_coords(band=band_coords)
        params = params.transpose("parameter", "band", "y", "x")
        unc = unc.transpose("band", "y", "x")
        scaled_unc = unc * np.float32(1.1)
        return BRDFKernelWeights(
            f0=params.sel(parameter="f0", drop=True),
            f1=params.sel(parameter="f1", drop=True),
            f2=params.sel(parameter="f2", drop=True),
            f0_unc=unc,
            f1_unc=scaled_unc,
            f2_unc=scaled_unc,
            reflectance_unc=unc,
        )

    def _temporal_weights_from_layers(
        self,
        per_date: dict[np.datetime64, dict[str, xr.DataArray]],
        requested: Sequence[tuple[str, ProductBandDefinition]],
        template: xr.DataArray,
        time_axis: np.ndarray,
    ) -> BRDFKernelWeights:
        band_coords = [coord for coord, _definition in requested]
        nan2d = xr.full_like(template, np.nan)
        nan_params = xr.concat([nan2d, nan2d, nan2d], dim="parameter").assign_coords(
            parameter=["f0", "f1", "f2"]
        )
        time_slabs_params = []
        time_slabs_unc = []
        for day in time_axis:
            day_d = np.datetime64(day, "D")
            date_layers = per_date.get(day_d, {})
            params_stack = []
            unc_stack = []
            for coord in band_coords:
                layer = date_layers.get(coord)
                params_stack.append(layer if layer is not None else nan_params)
                unc_stack.append(layer.attrs["unc"] if layer is not None else nan2d)
            time_slabs_params.append(
                xr.concat(params_stack, dim="band").assign_coords(band=band_coords)
            )
            time_slabs_unc.append(xr.concat(unc_stack, dim="band").assign_coords(band=band_coords))
        params = xr.concat(time_slabs_params, dim="time").assign_coords(time=time_axis)
        unc = xr.concat(time_slabs_unc, dim="time").assign_coords(time=time_axis)
        params = params.transpose("time", "parameter", "band", "y", "x")
        unc = unc.transpose("time", "band", "y", "x")
        scaled_unc = unc * np.float32(1.1)
        return BRDFKernelWeights(
            f0=params.sel(parameter="f0", drop=True),
            f1=params.sel(parameter="f1", drop=True),
            f2=params.sel(parameter="f2", drop=True),
            f0_unc=unc,
            f1_unc=scaled_unc,
            f2_unc=scaled_unc,
            reflectance_unc=unc,
        )

    def _default_weights(
        self,
        requested: Sequence[tuple[str, ProductBandDefinition]],
        template: xr.DataArray,
    ) -> BRDFKernelWeights:
        logger.warning("%s: no MCD43A1 data downloaded; returning NaN weights", self._source_name)
        return self._weights_from_band_layers({}, requested, template)

    # ── helpers ───────────────────────────────────────────────────────

    @staticmethod
    def _time_axis(obs_time: datetime, temporal_window: int) -> np.ndarray:
        center = np.datetime64(obs_time.date(), "D")
        offsets = np.arange(-temporal_window, temporal_window + 1)
        return center + offsets.astype("timedelta64[D]")

    @staticmethod
    def _date_from_tiff_name(name: str) -> np.datetime64 | None:
        # MODIS_061_MCD43A1_2024_02_16.tif -> 2024-02-16
        stem = name.rsplit(".", 1)[0]
        parts = stem.split("_")
        for i in range(len(parts) - 2):
            chunk = parts[i : i + 3]
            if len(chunk[0]) == 4 and all(p.isdigit() for p in chunk):
                try:
                    return np.datetime64(f"{chunk[0]}-{chunk[1]}-{chunk[2]}", "D")
                except ValueError:
                    continue
        return None

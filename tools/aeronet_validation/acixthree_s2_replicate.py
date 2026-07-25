#!/usr/bin/env python
"""Replicate the acixThree S2 AOD path on a SIAC/AERONET matchup.

This stages a Sentinel-2 L1C matchup into the NetCDF contract expected by
``/home/users/marcyin/acixThree`` and then calls the acixThree aerosol-species
AOD solver from ``compute_gas_transmittance.py``. Optionally it runs the
``hyper_siac.py`` final atmospheric-correction function on the solved AOD.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable

import numpy as np
import rioxarray  # noqa: F401 - registers the xarray .rio accessor
import xarray as xr

from tools.aeronet_validation.common import (
    DEFAULT_DATA_ROOT,
    ExperimentPaths,
    require_slurm_execution,
)
from tools.aeronet_validation.run_retrieval import RunSettings, build_config_payload


HARNESS = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
ACIX = Path("/home/users/marcyin/acixThree")
DEFAULT_OUT = HARNESS / "acixthree_s2"

S2_TO_ACIX = {
    "B01": "B1",
    "B02": "B2",
    "B03": "B3",
    "B04": "B4",
    "B05": "B5",
    "B06": "B6",
    "B07": "B7",
    "B08": "B8",
    "B8A": "B8A",
    "B09": "B9",
    "B10": "B10",
    "B11": "B11",
    "B12": "B12",
}
ACIX_S2_BANDS = ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"]
ACIX_AUX_BANDS = ["vza", "vaa", "DEM", "TCWV", "cloud_shadow"]
HYPER_CW_NM = np.array(
    [442.7, 492.7, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4],
    dtype=np.float32,
)
HYPER_FWHM_NM = np.array(
    [21.0, 66.0, 36.0, 31.0, 15.0, 15.0, 20.0, 106.0, 21.0, 20.0, 31.0, 91.0, 175.0],
    dtype=np.float32,
)
SPECIES_COMP_COL = {"B1": 0, "B2": 1, "B3": 2, "B4": 3, "B8A": 4, "B11": 5, "B12": 6}
ACIX_AOD_LUT = np.concatenate(
    [
        np.arange(0.01, 0.2, 0.01),
        np.arange(0.2, 0.5, 0.025),
        np.arange(0.5, 1.5, 0.05),
        np.arange(1.5, 2.5, 0.1),
    ]
).astype(np.float32)


@contextmanager
def pushd(path: Path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def _rows() -> dict[str, dict[str, str]]:
    return {r["matchup_id"]: r for r in csv.DictReader(open(HARNESS / "matchups" / "matchups.csv"))}


def _finite_mean(arr: xr.DataArray, default: float) -> float:
    vals = np.asarray(arr.values, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if vals.size else float(default)


def _write_crs(da: xr.DataArray, crs: object) -> xr.DataArray:
    out = da.rio.write_crs(crs)
    out = out.rio.set_spatial_dims("x", "y")
    return out.rio.write_coordinate_system()


def _coarsen_to_resolution(da: xr.DataArray, target_m: float) -> xr.DataArray:
    x = np.asarray(da.coords["x"].values, dtype=float)
    if x.size < 2:
        return da
    native = abs(float(np.nanmedian(np.diff(x))))
    factor = max(1, int(round(target_m / native)))
    if factor <= 1:
        return da
    out = da.coarsen(x=factor, y=factor, boundary="trim").mean()
    return _write_crs(out, da.rio.crs)


def _reproject_like(da: xr.DataArray, template: xr.DataArray, *, resampling: str = "bilinear") -> xr.DataArray:
    if da.rio.crs is None and template.rio.crs is not None:
        da = da.rio.write_crs(template.rio.crs)
    try:
        return da.rio.reproject_match(template, resampling=resampling)
    except Exception:
        return da.interp(x=template.x, y=template.y, method="linear")


def _s2_config():
    paths = ExperimentPaths(root=DEFAULT_DATA_ROOT)
    settings = RunSettings(0.05, 1500.0, 120.0, "gcs", 4, False, False)
    payload = build_config_payload(paths, "kernel_model", settings)
    payload.setdefault("providers", {}).setdefault("s2", {})["processing_level"] = "L1C"
    return paths, payload


def _resolve_s2(matchup: dict[str, str], half_width_deg: float):
    from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor
    from siac.api import resolve_s2_input
    from siac.config import SIACConfig

    _paths, payload = _s2_config()
    cfg = SIACConfig.model_validate(payload)
    pre = Sentinel2Preprocessor(config={})
    lon = float(matchup["longitude"])
    lat = float(matchup["latitude"])
    pre.set_spatial_subset((lon - half_width_deg, lat - half_width_deg, lon + half_width_deg, lat + half_width_deg), crs="EPSG:4326")
    safe = resolve_s2_input(matchup["product_id"], cfg)
    return cfg, pre, safe


def _atmo_for_template(cfg, template: xr.DataArray, when) -> object:
    from siac.app._assembly_providers import resolve_atmo_provider

    bounds = tuple(float(v) for v in template.rio.bounds())
    provider = resolve_atmo_provider(cfg)
    return provider(bounds, str(template.rio.crs), when, 120.0)


def _maybe_maiac_prior(mid: str) -> float | None:
    path = HARNESS / "maiac_qa" / f"{mid}.json"
    if not path.exists():
        return None
    try:
        rec = json.loads(path.read_text())
        val = rec.get("aot")
        return float(val) if val is not None and np.isfinite(float(val)) else None
    except Exception:
        return None


def _cams_prior(mid: str) -> float | None:
    path = HARNESS / "prior_quality" / f"{mid}.json"
    if not path.exists():
        return None
    try:
        rec = json.loads(path.read_text())
        val = rec.get("cams_aot")
        return float(val) if val is not None and np.isfinite(float(val)) else None
    except Exception:
        return None


def _build_mask_dataset(s2a: xr.DataArray, hyper_band: xr.DataArray) -> xr.Dataset:
    def no_band_coord(da: xr.DataArray) -> xr.DataArray:
        return da.reset_coords("band", drop=True) if "band" in da.coords and "band" not in da.dims else da

    valid_s2 = np.isfinite(s2a.sel(band=["B1", "B2", "B3", "B4", "B8"])).all(dim="band")
    cloud_shadow = no_band_coord(s2a.sel(band="cloud_shadow") != 0)
    cirrus = no_band_coord(s2a.sel(band="B10") >= 0.015)
    ndwi = (s2a.sel(band="B3") - s2a.sel(band="B8")) / (s2a.sel(band="B3") + s2a.sel(band="B8"))
    water = no_band_coord(ndwi >= 0.1)
    aod_mask = (~valid_s2) | cloud_shadow | cirrus | water
    aod_mask = no_band_coord(aod_mask)
    outside = no_band_coord(~np.isfinite(s2a.sel(band=ACIX_S2_BANDS)).any(dim="band"))
    aod_mask = aod_mask | outside
    large = xr.concat(
        [xr.zeros_like(aod_mask, dtype=bool), xr.zeros_like(aod_mask, dtype=bool), xr.zeros_like(aod_mask, dtype=bool)],
        dim=xr.DataArray(["large_invalid_data_mask", "large_water_mask", "large_cloud_mask"], dims="mask", name="mask"),
    )
    tcwv_mask = no_band_coord(water | cirrus | no_band_coord(s2a.sel(band="B8A") <= 0.05) | outside)
    toa_ref_mask = ~(np.isfinite(hyper_band) & (hyper_band > 0))
    ds = xr.Dataset(
        {
            "toa_ref_mask": toa_ref_mask,
            "cirrus_mask": cirrus,
            "aod_mask": aod_mask,
            "cloud_shadow_mask": cloud_shadow,
            "outside_obs_extend_mask": outside,
            "large_aera_invalid_data_mask": large,
            "tcwv_mask": tcwv_mask,
        }
    )
    ds = ds.rio.write_crs(s2a.rio.crs).rio.set_spatial_dims("x", "y").rio.write_coordinate_system()
    return ds


def _set_geotransform(ds: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    x = np.asarray(ds.coords["x"].values, dtype=float)
    y = np.asarray(ds.coords["y"].values, dtype=float)
    if x.size < 2 or y.size < 2:
        return ds
    xres = float(np.nanmedian(np.diff(x)))
    yres = float(np.nanmedian(np.diff(y)))
    gt = (x[0] - xres / 2.0, xres, 0.0, y[0] - yres / 2.0, 0.0, yres)
    if "spatial_ref" in ds.coords:
        ds.coords["spatial_ref"].attrs["GeoTransform"] = " ".join(f"{v:.12g}" for v in gt)
    return ds


def stage_s2_bundle(
    mid: str,
    out_dir: Path,
    *,
    half_width_deg: float,
    input_resolution_m: float,
    aod_prior: str,
    overwrite: bool,
) -> Path:
    s2_toa_file = out_dir / "s2a_toa.nc"
    if s2_toa_file.exists() and not overwrite:
        return s2_toa_file
    hyper_file = Path(str(out_dir) + ".nc")
    mask_file = out_dir / "invalid_data_masks.nc"
    if overwrite:
        for path in (s2_toa_file, hyper_file, mask_file):
            if path.exists():
                path.unlink()

    matchup = _rows()[mid]
    cfg, pre, safe = _resolve_s2(matchup, half_width_deg)
    out_dir.mkdir(parents=True, exist_ok=True)

    toa = pre.load_toa(safe, band_names=tuple(S2_TO_ACIX))
    geometry = pre.extract_geometry(safe)
    metadata = pre.get_metadata(safe)
    template = _coarsen_to_resolution(toa["B02"], input_resolution_m)
    crs = template.rio.crs

    arrays: list[xr.DataArray] = []
    for siac_name, acix_name in S2_TO_ACIX.items():
        band = _coarsen_to_resolution(toa[siac_name], input_resolution_m)
        if not np.array_equal(band.x, template.x) or not np.array_equal(band.y, template.y):
            band = _reproject_like(band, template)
        arrays.append(band.rename(acix_name))

    def geom_deg(field: xr.DataArray) -> xr.DataArray:
        da = _reproject_like(field, template)
        vals = np.degrees(np.asarray(da.values, dtype=np.float32))
        vals[~np.isfinite(vals)] = np.nanmean(vals)
        return xr.DataArray(vals, dims=template.dims, coords=template.coords)

    atmo = _atmo_for_template(cfg, template, metadata["observation_time"])
    tcwv = _reproject_like(atmo.tcwv, template)
    dem_m = _reproject_like(atmo.elevation, template) * 1000.0
    tcwv_vals = np.asarray(tcwv.values, dtype=np.float32)
    dem_vals = np.asarray(dem_m.values, dtype=np.float32)
    tcwv_vals[~np.isfinite(tcwv_vals)] = np.nanmean(tcwv_vals)
    dem_vals[~np.isfinite(dem_vals)] = np.nanmean(dem_vals)

    aux = {
        "vza": geom_deg(geometry.vza),
        "vaa": geom_deg(geometry.vaa),
        "DEM": xr.DataArray(dem_vals, dims=template.dims, coords=template.coords),
        "TCWV": xr.DataArray(tcwv_vals, dims=template.dims, coords=template.coords),
        "cloud_shadow": xr.zeros_like(template, dtype=np.float32),
    }
    arrays.extend(da.rename(name) for name, da in aux.items())
    s2a = xr.concat(arrays, dim=xr.DataArray(ACIX_S2_BANDS + ACIX_AUX_BANDS, dims="band", name="band"))
    s2a = _write_crs(s2a.astype(np.float32), crs)

    prior_val = None
    if aod_prior == "maiac":
        prior_val = _maybe_maiac_prior(mid)
    if prior_val is None:
        prior_val = _cams_prior(mid)
    if prior_val is None:
        prior_val = _finite_mean(atmo.aot, 0.2)
    sza = geom_deg(geometry.sza)
    saa = geom_deg(geometry.saa)
    s2a.attrs.update(
        {
            "sza": float(np.nanmean(sza.values)),
            "saa": float(np.nanmean(saa.values)),
            "aod": float(prior_val),
            "ozone": float(_finite_mean(atmo.tco3, 0.30) * 1000.0),
            "sensing_time": metadata["observation_time"].strftime("%Y-%m-%d %H:%M:%S"),
            "matchup_id": mid,
            "product_id": matchup["product_id"],
        }
    )
    s2a = _set_geotransform(s2a)
    s2a.to_netcdf(s2_toa_file)

    hyper = s2a.sel(band=ACIX_S2_BANDS).assign_coords(band=HYPER_CW_NM).rename("reflectance")
    hyper.attrs.update({"cw": HYPER_CW_NM, "fwhm": HYPER_FWHM_NM, "sensing_time": s2a.attrs["sensing_time"]})
    hds = xr.Dataset({"reflectance": hyper})
    hds = hds.rio.write_crs(crs).rio.set_spatial_dims("x", "y").rio.write_coordinate_system()
    hds = _set_geotransform(hds)
    hds.to_netcdf(hyper_file)

    masks = _build_mask_dataset(s2a, hyper)
    masks = _set_geotransform(masks)
    mask_encoding = {}
    for name in masks.data_vars:
        masks[name].attrs = {k: v for k, v in masks[name].attrs.items() if k != "_FillValue"}
        if masks[name].dtype == bool:
            masks[name] = masks[name].astype("uint8")
            mask_encoding[name] = {"dtype": "uint8", "_FillValue": np.uint8(255)}
    masks.to_netcdf(mask_file, encoding=mask_encoding)
    return s2_toa_file


def _npz_band_da(npz, band: str) -> tuple[xr.DataArray, xr.DataArray]:
    comp = np.asarray(npz["comp"], dtype=np.float32)
    col = SPECIES_COMP_COL[band]
    stack = comp[:, col]
    med = np.nanmedian(stack, axis=0)
    rmse = np.sqrt(np.nanmean((stack - med) ** 2, axis=0))
    epsg = int(npz["epsg"])
    a, b, xoff, d, e, yoff = [float(v) for v in npz["transform"]]
    h, w = med.shape
    x = xoff + (np.arange(w) + 0.5) * a
    y = yoff + (np.arange(h) + 0.5) * e
    med_da = xr.DataArray(med, dims=("y", "x"), coords={"y": y, "x": x}, name=band).rio.write_crs(f"EPSG:{epsg}")
    rmse_da = xr.DataArray(rmse, dims=("y", "x"), coords={"y": y, "x": x}, name=band).rio.write_crs(f"EPSG:{epsg}")
    return med_da, rmse_da


def _read_acix_s2a_toa(s2_toa_file: Path) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    s2_toa = xr.open_dataarray(s2_toa_file, decode_coords="all")
    masks = xr.open_dataset(s2_toa_file.with_name("invalid_data_masks.nc"), decode_coords="all")
    aod_mask = masks.aod_mask.astype(bool)
    s2_toa = s2_toa.where(~aod_mask)
    s2_bands = [
        "B1",
        "B2",
        "B3",
        "B4",
        "B5",
        "B6",
        "B7",
        "B8",
        "B8A",
        "B9",
        "B10",
        "B11",
        "B12",
        "vza",
        "vaa",
        "DEM",
        "TCWV",
    ]
    out = s2_toa.sel(band=s2_bands).coarsen(x=4, y=4, boundary="pad").mean(keep_attrs=True)
    out.attrs.update(s2_toa.attrs)
    out = _write_crs(out, s2_toa.rio.crs)
    return (
        out,
        masks.outside_obs_extend_mask.astype(bool),
        masks.large_aera_invalid_data_mask.astype(bool),
    )


def _torch_nanmin(tensor, dim=None, keepdim=False):
    import torch

    max_value = torch.finfo(tensor.dtype).max
    return tensor.nan_to_num(max_value).min(dim=dim, keepdim=keepdim)


def _get_med_diffs(diffs, sum_size: int = 10):
    import torch.nn as nn

    x_shape, y_shape = diffs.shape[1:]
    half = sum_size // 2
    pad = nn.ConstantPad2d((half, sum_size + sum_size, half, sum_size + sum_size), np.nan)
    padded = pad(diffs)
    unfolded = padded.unfold(1, sum_size, 1).unfold(2, sum_size, 1)
    med = unfolded.reshape(unfolded.shape[:3] + (-1,)).nanmedian(dim=-1)[0]
    return med[:, :x_shape, :y_shape]


def _acix_aerosol_retrieval(
    s2a_toa: xr.DataArray,
    band_names: list[str],
    rsrs: list[object],
    prior: xr.DataArray,
    rmse: xr.DataArray,
    unique_aerosol_inds: np.ndarray,
) -> tuple[xr.DataArray, xr.DataArray]:
    import pandas as pd
    import torch
    import aerosol_species_rt.atmo_trans_utils as atmo_trans_utils
    from aerosol_species_rt.atmo_solver_utils import (
        compute_aerosol_coefs,
        prepare_emu_inputs,
        prepare_wavelengths,
    )

    orig_unique_rows = atmo_trans_utils.unique_rows

    def unique_rows_flat(ar):
        rows, inverse = orig_unique_rows(ar)
        return rows, np.asarray(inverse).ravel()

    atmo_trans_utils.unique_rows = unique_rows_flat
    band_index = {name: band_names.index(name) for name in band_names}
    emu_inps, tg_bands, atmo_solver_aux = prepare_emu_inputs(s2a_toa, band_names, rsrs)
    wavelengths_config = prepare_wavelengths(atmo_solver_aux)
    device = atmo_solver_aux.scales.device

    def to_torch_da(data: xr.DataArray):
        return torch.tensor(data.values.astype(np.float32)).to(device)

    with torch.no_grad():
        tgp1 = to_torch_da(tg_bands.sel(tg="tgp1"))
        tgp2 = to_torch_da(tg_bands.sel(tg="tgp2"))
        tgtot = to_torch_da(tg_bands.sel(tg="tgtot"))
        toa_ref = to_torch_da(s2a_toa.sel(band=band_names))
        boa_ref = to_torch_da(prior.sel(band=band_names))
        boa_rmse = to_torch_da(rmse.sel(band=band_names))
        aod_lut_torch = torch.tensor(ACIX_AOD_LUT).to(device)
        mcd19_aod = np.atleast_1d(float(s2a_toa.attrs["aod"]))
        mcd19_aod[mcd19_aod <= 0.15] = mcd19_aod[mcd19_aod <= 0.15] - 0.033
        mcd19_aod = float(np.maximum(mcd19_aod, 0.0)[0])
        aod_unc = float(np.max([mcd19_aod * 0.5, 0.01]))
        costs = []
        for mix in unique_aerosol_inds:
            mix = int(mix)
            emu_arr_model = [
                atmo_solver_aux.t_down_up_arrModel[mix],
                atmo_solver_aux.romix_arrModel[mix],
                atmo_solver_aux.roray_arrModel[mix],
                atmo_solver_aux.sast_arrModel[mix],
            ]
            emu_scales = atmo_solver_aux.scales[mix]
            diffs = []
            for aod in aod_lut_torch:
                emu_inps[:, :, 4] = atmo_solver_aux.transform_dict["AOT"](aod)
                outs = compute_aerosol_coefs(
                    emu_inps,
                    emu_arr_model,
                    emu_scales,
                    wavelengths_config,
                    return_aerosol_coefs=False,
                )
                ratm2 = (outs[:, 1] - outs[:, 2]) * tgp2 + outs[:, 2] * tgp1
                lut_boa = (toa_ref - ratm2) / tgtot / outs[:, 0]
                lut_boa = lut_boa / (1 + lut_boa * outs[:, 3])
                diff = (
                    (lut_boa[band_index["B4"]] - boa_ref[band_index["B4"]]) ** 2
                    / boa_rmse[band_index["B4"]] ** 2
                    + (lut_boa[band_index["B2"]] - boa_ref[band_index["B2"]]) ** 2
                    / boa_rmse[band_index["B2"]] ** 2
                )
                diff += (mcd19_aod - aod) ** 2 / aod_unc**2
                diffs.append(diff)
            diffs = _get_med_diffs(torch.stack(diffs), sum_size=10)
            diff_min, diff_min_inds = _torch_nanmin(diffs, dim=0)
            invalid = ~torch.isfinite(diffs[0])
            diff_min[invalid] = float("nan")
            aod_map = aod_lut_torch[diff_min_inds]
            aod_map[invalid] = float("nan")
            costs.append(torch.stack([aod_map, diff_min]))

    aod_costs = torch.stack(costs).cpu().numpy()
    aod_costs = xr.DataArray(
        aod_costs,
        coords={
            "aerosol_species": unique_aerosol_inds.astype(int),
            "band": ["aod", "cost"],
            "y": s2a_toa.y,
            "x": s2a_toa.x,
        },
        dims=["aerosol_species", "band", "y", "x"],
    )
    total_cost = aod_costs.sel(band="cost").sum(dim=["x", "y"])
    best_order = np.argsort(total_cost.values)
    best = aod_costs.sel(aerosol_species=int(unique_aerosol_inds[best_order[0]]))
    aerosol_model = pd.read_csv(ACIX / "data" / "aerosol_type_lut.csv")
    best_models = aerosol_model.iloc[unique_aerosol_inds[best_order[:5]].astype(int)]
    best.rio.write_crs(s2a_toa.rio.crs, inplace=True)
    aod_costs.rio.write_crs(s2a_toa.rio.crs, inplace=True)
    best.attrs["best_aerosol_models"] = json.dumps(np.asarray(best_models.values).tolist())
    return best, aod_costs


def _gap_fill_aod(aod_da: xr.DataArray) -> xr.DataArray:
    from smoothn import smoothn

    aod_da = aod_da.copy()
    mean = float(aod_da.mean(skipna=True).values)
    for coarsen_level in [64, 32, 16, 8, 4, 2, 1]:
        coarse = aod_da.coarsen(x=coarsen_level, y=coarsen_level, boundary="pad").mean()
        invalid = ~np.isfinite(coarse.values)
        if coarsen_level == 64:
            coarse.values[invalid] = mean
        else:
            background = background.rio.reproject_match(coarse, resampling=1, nodata=np.nan)
            coarse.values[invalid] = background.values[invalid]
        background = coarse.copy()
    filled = background
    filled.values[:] = smoothn(filled.values, s=1)
    return filled


def _apply_acix_masks(
    aod_da: xr.DataArray,
    large_area_invalid_data_mask: xr.DataArray,
    outside_obs_extend_mask: xr.DataArray,
) -> xr.DataArray:
    aod_da = aod_da.rio.reproject_match(outside_obs_extend_mask, resampling=1, nodata=np.nan)
    large = large_area_invalid_data_mask.any(dim="mask").values
    vals = np.asarray(aod_da.values)
    mean_aod = float(np.nanmean(vals[~large])) if np.isfinite(vals[~large]).any() else float(np.nanmean(vals))
    vals[~np.isfinite(vals)] = mean_aod
    aod_da.values[:] = vals
    aod_da = aod_da.where(~large, mean_aod)
    return aod_da.where(~outside_obs_extend_mask, np.nan)


def write_species_prior(mid: str, s2_toa_file: Path, surface_dir: Path, *, overwrite: bool) -> tuple[Path, Path]:
    prior_file = s2_toa_file.with_name("s2a_boa_prediction_species.nc")
    rmse_file = s2_toa_file.with_name("s2a_boa_prediction_rmse_species.nc")
    if prior_file.exists() and rmse_file.exists() and not overwrite:
        return prior_file, rmse_file
    if overwrite:
        for path in (prior_file, rmse_file):
            if path.exists():
                path.unlink()
    species_file = surface_dir / f"{mid}.npz"
    if not species_file.exists():
        raise FileNotFoundError(f"Missing L1C species surface prior: {species_file}")

    s2a_120, _outside, _large = _read_acix_s2a_toa(s2_toa_file)
    template = s2a_120.sel(band="B1")
    bands = ["B1", "B2", "B3", "B4"]
    priors = []
    rmses = []
    with np.load(species_file) as z:
        for band in bands:
            med_da, rmse_da = _npz_band_da(z, band)
            med = _reproject_like(med_da, template)
            rmse = _reproject_like(rmse_da, template)
            floor = 0.006
            priors.append(med.rename(band))
            rmses.append(xr.where(np.isfinite(rmse), np.maximum(rmse, floor), floor).rename(band))
    prior = xr.concat(priors, dim=xr.DataArray(bands, dims="band", name="band")).astype(np.float32)
    rmse = xr.concat(rmses, dim=xr.DataArray(bands, dims="band", name="band")).astype(np.float32)
    prior = _write_crs(prior, template.rio.crs)
    rmse = _write_crs(rmse, template.rio.crs)
    prior.to_netcdf(prior_file)
    rmse.to_netcdf(rmse_file)
    return prior_file, rmse_file


def _site_mean(da: xr.DataArray, lon: float, lat: float, radius_m: float) -> float:
    from pyproj import Transformer

    tx = Transformer.from_crs("EPSG:4326", da.rio.crs, always_xy=True)
    x0, y0 = tx.transform(lon, lat)
    x = np.asarray(da.x.values)
    y = np.asarray(da.y.values)
    mask = (np.abs(x[None, :] - x0) <= radius_m) & (np.abs(y[:, None] - y0) <= radius_m)
    vals = np.asarray(da.values)[mask]
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if vals.size else float("nan")


def run_acix_aod(mid: str, s2_toa_file: Path, *, overwrite: bool) -> Path:
    aod_file = s2_toa_file.with_name("aod_map_areosol_species.nc")
    if aod_file.exists() and not overwrite:
        return aod_file
    if overwrite:
        for path in (
            aod_file,
            s2_toa_file.with_name("aod_cost_best_model.nc"),
            s2_toa_file.with_name("aod_costs_all.nc"),
        ):
            if path.exists():
                path.unlink()
    sys.path.insert(0, str(ACIX))
    with pushd(ACIX):
        from aerosol_species_rt.atmo_solver_utils import read_aerosol_climatology
        from aerosol_species_rt.atmo_trans_utils import create_s2_rsrs

        s2a_toa, outside_obs_extend_mask, large_aera_invalid_data_mask = _read_acix_s2a_toa(s2_toa_file)
        prior = xr.open_dataarray(s2_toa_file.with_name("s2a_boa_prediction_species.nc"), decode_coords="all")
        rmse = xr.open_dataarray(s2_toa_file.with_name("s2a_boa_prediction_rmse_species.nc"), decode_coords="all") + 1e-6
        rmse = rmse.rio.reproject_match(s2a_toa, resampling=5, nodata=np.nan)
        aerosol_inds, _types, _mix = read_aerosol_climatology(s2a_toa)
        unique_aerosol_inds = np.unique(aerosol_inds)
        band_names, rsrs = create_s2_rsrs()
        aerosol_bands = ["B1", "B2", "B3", "B4"]
        inds = [band_names.index(b) for b in aerosol_bands]
        best, costs = _acix_aerosol_retrieval(
            s2a_toa,
            aerosol_bands,
            [rsrs[i] for i in inds],
            prior,
            rmse,
            unique_aerosol_inds,
        )
        aod_da = best.sel(band="aod")
        aod_da = _gap_fill_aod(aod_da)
        aod_da = _apply_acix_masks(aod_da, large_aera_invalid_data_mask, outside_obs_extend_mask)
        species = int(np.asarray(best.coords["aerosol_species"].values).item())
        aod_da = aod_da.assign_coords(aerosol_species=species)
        aod_da.attrs["best_aerosol_models"] = best.attrs.get("best_aerosol_models", "[]")
        aod_da.to_netcdf(aod_file)
        best.to_netcdf(s2_toa_file.with_name("aod_cost_best_model.nc"))
        costs.to_netcdf(s2_toa_file.with_name("aod_costs_all.nc"))
    return aod_file


def run_hyper(s2_toa_file: Path, backend: str) -> str:
    if backend == "none":
        return "skipped"
    sys.path.insert(0, str(ACIX))
    with pushd(ACIX):
        import hyper_siac
        import preprocessing.ac_utils as ac_utils

        orig_open_dataset = hyper_siac.xr.open_dataset

        def open_dataset_cast_masks(*args, **kwargs):
            ds = orig_open_dataset(*args, **kwargs)
            try:
                path = Path(str(args[0]))
            except Exception:
                return ds
            if path.name == "invalid_data_masks.nc":
                for name in (
                    "toa_ref_mask",
                    "cirrus_mask",
                    "aod_mask",
                    "cloud_shadow_mask",
                    "outside_obs_extend_mask",
                    "large_aera_invalid_data_mask",
                    "tcwv_mask",
                ):
                    if name in ds:
                        ds[name] = ds[name].fillna(0).astype(bool)
            return ds

        hyper_siac.xr.open_dataset = open_dataset_cast_masks
        hyper_siac.cal_rsr_esun = lambda cw, sigma: ac_utils.cal_rsr_esun(
            cw, sigma, solar_file=str(ACIX / "data" / "6S_LUTs.npz")
        )
        if backend == "6s":
            hyper_siac.ac_cor_6S(str(s2_toa_file))
            return "ok:6s"
        if backend == "lut":
            lut = Path(os.environ.get("ACIXTHREE_LRT_LUT", "/work/scratch-pw3/marc/libradtran_continental_average_lut_1nm.zarr"))
            if not lut.exists():
                raise FileNotFoundError(
                    f"Missing hyper_siac libRadtran LUT {lut}; set ACIXTHREE_LRT_LUT or restore the Zarr."
                )
            hyper_siac.ac_using_libradtran_lut(str(s2_toa_file))
            return "ok:lut"
    raise ValueError(f"Unknown hyper backend: {backend}")


def write_summary(mid: str, out_dir: Path, aod_file: Path, hyper_status: str, radius_m: float) -> Path:
    row = _rows()[mid]
    aod = xr.open_dataarray(aod_file, decode_coords="all")
    site_aod = _site_mean(aod, float(row["longitude"]), float(row["latitude"]), radius_m)
    truth = float(row["aeronet_aod550_mean"])
    ee = 0.05 + 0.15 * truth
    rec = {
        "matchup_id": mid,
        "site": row["site"],
        "truth_aod550": truth,
        "acixthree_s2_aod550_window_mean": site_aod,
        "error": site_aod - truth if np.isfinite(site_aod) else None,
        "expected_error": ee,
        "within_expected_error": bool(np.isfinite(site_aod) and abs(site_aod - truth) <= ee),
        "aod_file": str(aod_file),
        "hyper_status": hyper_status,
    }
    path = out_dir / "summary.json"
    path.write_text(json.dumps(rec, indent=2, default=float))
    print(json.dumps(rec, indent=2, default=float), flush=True)
    return path


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("matchup_id")
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUT)
    p.add_argument("--surface-dir", type=Path, default=HARNESS / "l1c_seasonal_species_lut")
    p.add_argument("--half-width-deg", type=float, default=0.05)
    p.add_argument("--input-resolution-m", type=float, default=30.0)
    p.add_argument("--aod-prior", choices=["maiac", "cams"], default="maiac")
    p.add_argument("--hyper-backend", choices=["none", "6s", "lut"], default="6s")
    p.add_argument("--extract-radius-m", type=float, default=1500.0)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--allow-login",
        action="store_true",
        help="Run locally on an interactive node (not recommended).",
    )
    return p.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    require_slurm_execution(
        "acixThree S2 replicate",
        allow_login=args.allow_login,
        suggestion=(
            "Submit through an SBATCH script such as "
            "`campaign250_acixthree_sd60_submit.sbatch` from the "
            "aeronet_validation folder."
        ),
    )
    out_dir = args.output_root / args.matchup_id
    s2_toa = stage_s2_bundle(
        args.matchup_id,
        out_dir,
        half_width_deg=args.half_width_deg,
        input_resolution_m=args.input_resolution_m,
        aod_prior=args.aod_prior,
        overwrite=args.overwrite,
    )
    write_species_prior(args.matchup_id, s2_toa, args.surface_dir, overwrite=args.overwrite)
    aod = run_acix_aod(args.matchup_id, s2_toa, overwrite=args.overwrite)
    hyper_status = run_hyper(s2_toa, args.hyper_backend)
    write_summary(args.matchup_id, out_dir, aod, hyper_status, args.extract_radius_m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

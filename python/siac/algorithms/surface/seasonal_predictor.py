"""Seasonal scene-conditioned surface-prior prediction.

This module ports the compact ExtraTree prior used by the AERONET seasonal
research harness into a reusable production helper.  It replaces selected
visible BOA prior bands with predictions from the scene's aerosol-robust NIR/SWIR
anchor plus a local climatological surface fingerprint.  The scene anchor bands
are corrected TOA->BOA through the configured RT backend — the same RT model the
solver uses — so the prediction and the AOD solve share one radiative-transfer
space (validated 2026-07-02: the real-RT anchor makes the AERONET-learned
prior debias unnecessary).
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from siac.runtime import ObservationBundle, SurfacePrior


ANCHOR_BANDS: tuple[str, ...] = ("B8A", "B11", "B12")
ANCHOR_COLUMNS: tuple[int, ...] = (4, 5, 6)
VISIBLE_COLUMNS: tuple[int, ...] = (0, 1, 2, 3)
DEFAULT_BAND_COLUMNS: dict[str, int] = {
    "B01": 0,
    "B02": 1,
    "B03": 2,
    "B04": 3,
    "B8A": 4,
    "B11": 5,
    "B12": 6,
}
_MAD_TO_STD = 1.4826
logger = logging.getLogger(__name__)


def _as_template_grid(da: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    if "band" in da.dims:
        da = da.isel(band=0, drop=True)
    return cast("xr.DataArray", da.rio.reproject_match(template))


def _mean_angle_degrees(da: xr.DataArray) -> float:
    values = np.asarray(da.values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return float(np.degrees(np.mean(finite))) if finite.size else 0.0


def _finite_mean(da: xr.DataArray, default: float) -> float:
    values = np.asarray(da.values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return float(np.mean(finite)) if finite.size else float(default)


def _full_template(template: xr.DataArray, value: float) -> xr.DataArray:
    return cast("xr.DataArray", xr.full_like(template, float(value)).astype(np.float32))


def _correct_anchor_reflectance(
    observation: ObservationBundle,
    *,
    atmo_prior: Any,
    rt_model: Any,
    template: xr.DataArray,
    anchor_grids: Mapping[str, xr.DataArray],
    valid: np.ndarray,
    anchor_aot: float,
) -> np.ndarray:
    """Correct scene anchor bands TOA->BOA with the configured RT backend.

    The anchor must be corrected in the same RT space the solver uses; there is
    deliberately no analytic fallback.
    """
    from siac.runtime import AtmosphericState, GeometryAngles

    geometry = observation.geometry
    sza_deg = _mean_angle_degrees(geometry.sza)
    vza_deg = _mean_angle_degrees(geometry.vza)
    raa_deg = _mean_angle_degrees(geometry.raa)
    atmo = AtmosphericState(
        aot=_full_template(template, float(anchor_aot)),
        tcwv=_full_template(template, _finite_mean(atmo_prior.tcwv, 2.0)),
        tco3=_full_template(template, _finite_mean(atmo_prior.tco3, 0.3)),
        aot_unc=_full_template(template, 0.1),
        tcwv_unc=_full_template(template, _finite_mean(atmo_prior.tcwv_unc, 0.5)),
        tco3_unc=_full_template(template, _finite_mean(atmo_prior.tco3_unc, 0.05)),
        elevation=_full_template(template, _finite_mean(atmo_prior.elevation, 0.0)),
    )
    geom = GeometryAngles(
        sza=_full_template(template, np.radians(sza_deg)),
        saa=_full_template(template, 0.0),
        vza=_full_template(template, np.radians(vza_deg)),
        vaa=_full_template(template, np.radians(raa_deg)),
    )
    corrected = []
    for band_name in ANCHOR_BANDS:
        band = observation.sensor_config.get_band(band_name)
        coeffs = rt_model.compute_coefficients(geom, atmo, band)
        corrected.append(coeffs.apply_correction(anchor_grids[band_name]).values.ravel()[valid])
    return np.column_stack(corrected)


def _robust_clip_composites(comp: np.ndarray, clip: float) -> np.ndarray:
    if clip <= 0.0 or comp.shape[0] < 3:
        return comp
    with np.errstate(invalid="ignore"):
        median = np.nanmedian(comp, axis=0)
        mad = np.nanmedian(np.abs(comp - median[np.newaxis]), axis=0) * _MAD_TO_STD
        keep = np.abs(comp - median[np.newaxis]) <= (
            float(clip) * np.maximum(mad, 1.0e-4)[np.newaxis]
        )
    return np.where(keep, comp, np.nan)


def _composite_reference_on_scene(
    comp: np.ndarray,
    *,
    band_column: int,
    epsg: int,
    x: np.ndarray,
    y: np.ndarray,
    template: xr.DataArray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    band_stack = comp[:, int(band_column)]
    with np.errstate(invalid="ignore"):
        median = np.nanmedian(band_stack, axis=0)
        rmse = np.sqrt(np.nanmean((band_stack - median[np.newaxis]) ** 2, axis=0))
    median_da = xr.DataArray(median, dims=("y", "x"), coords={"y": y, "x": x}).rio.write_crs(
        f"EPSG:{int(epsg)}"
    )
    rmse_da = xr.DataArray(rmse, dims=("y", "x"), coords={"y": y, "x": x}).rio.write_crs(
        f"EPSG:{int(epsg)}"
    )
    return (
        _as_template_grid(median_da, template).values.ravel()[valid],
        _as_template_grid(rmse_da, template).values.ravel()[valid],
    )


def seasonal_extra_tree_prior(
    prior: SurfacePrior,
    observation: ObservationBundle,
    *,
    seasonal_composites: np.ndarray,
    epsg: int,
    transform: Sequence[float],
    anchor_aot: float,
    atmo_prior: Any,
    rt_model: Any,
    target_band_columns: Mapping[str, int] | None = None,
    debias: Mapping[str, tuple[float, float]] | None = None,
    debias_scale: float = 1.0,
    uncertainty_floor: float = 0.006,
    b01_uncertainty_floor: float | None = None,
    min_samples_leaf: int = 5,
    random_state: int = 0,
    robust_clip: float = 0.0,
    composite_blend_weight: float = 0.0,
) -> SurfacePrior:
    """Replace visible prior bands with a seasonal ExtraTree ensemble prediction.

    Parameters mirror the validated research harness:

    - correct the scene ``B8A/B11/B12`` anchor TOA->BOA at ``anchor_aot``
      through ``rt_model`` (the same RT backend the solver uses);
    - train one ``ExtraTreeRegressor`` per seasonal realization;
    - features are corrected ``B8A/B11/B12`` plus the per-pixel mean visible
      climatology ``B01/B02/B03/B04``;
    - prediction is median across realizations;
    - uncertainty is MAD*1.4826 across realization predictions, floored.
    """
    from sklearn.tree import ExtraTreeRegressor

    if rt_model is None or atmo_prior is None:
        raise ValueError(
            "seasonal_extra_tree_prior requires the configured RT backend and an "
            "atmospheric prior to correct the scene anchor bands."
        )

    boa = prior.boa
    if "band" not in boa.dims:
        return prior
    present = {str(name) for name in boa["band"].values}
    targets = dict(target_band_columns or {"B02": 1, "B04": 3})
    targets = {name: col for name, col in targets.items() if name in present}
    if not targets:
        return prior
    if any(name not in observation.toa.data_vars for name in ANCHOR_BANDS):
        return prior

    comp = np.asarray(seasonal_composites, dtype=np.float64)
    if comp.ndim != 4 or comp.shape[1] < 7:
        raise ValueError(
            f"seasonal_composites must have shape (n_realizations, 7, y, x); got {comp.shape}"
        )
    if comp.shape[0] == 0:
        return prior
    comp = _robust_clip_composites(comp, float(robust_clip))

    template = boa.isel(band=0, drop=True)
    anchor_grids = {
        name: _as_template_grid(observation.toa[name], template).astype(np.float64)
        for name in ANCHOR_BANDS
    }
    anchor = np.stack(
        [anchor_grids[name].values.astype(np.float64) for name in ANCHOR_BANDS],
        axis=-1,
    )
    flat_anchor = anchor.reshape(-1, len(ANCHOR_BANDS))
    valid = np.all(np.isfinite(flat_anchor) & (flat_anchor > 0.0) & (flat_anchor < 1.2), axis=1)
    if not np.any(valid):
        return prior

    corrected_anchor = _correct_anchor_reflectance(
        observation,
        atmo_prior=atmo_prior,
        rt_model=rt_model,
        template=template,
        anchor_grids=anchor_grids,
        valid=valid,
        anchor_aot=float(anchor_aot),
    )

    n_real, _, height, width = comp.shape
    tr = [float(value) for value in transform]
    x = tr[2] + (np.arange(width) + 0.5) * tr[0]
    y = tr[5] + (np.arange(height) + 0.5) * tr[4]
    mean_visible = np.nanmean(comp[:, list(VISIBLE_COLUMNS)], axis=0)
    localizer_comp = mean_visible.reshape(len(VISIBLE_COLUMNS), -1).T
    localizer_scene_parts = []
    for index in range(len(VISIBLE_COLUMNS)):
        da = xr.DataArray(
            mean_visible[index],
            dims=("y", "x"),
            coords={"y": y, "x": x},
        ).rio.write_crs(f"EPSG:{int(epsg)}")
        localizer_scene_parts.append(_as_template_grid(da, template).values.ravel()[valid])
    localizer_scene = np.column_stack(localizer_scene_parts)
    scene_features = np.column_stack([corrected_anchor, localizer_scene])

    target_names = list(targets)
    target_cols = [int(targets[name]) for name in target_names]
    predictions: dict[str, list[np.ndarray]] = {name: [] for name in target_names}
    used_realizations = 0
    for index in range(n_real):
        composite = comp[index].reshape(comp.shape[1], -1).T
        good = np.all(np.isfinite(composite), axis=1)
        if int(np.count_nonzero(good)) < 200:
            continue
        train_x = np.column_stack([composite[good][:, list(ANCHOR_COLUMNS)], localizer_comp[good]])
        train_y = composite[good][:, target_cols]
        tree = ExtraTreeRegressor(
            random_state=int(random_state),
            min_samples_leaf=int(min_samples_leaf),
        ).fit(train_x, train_y)
        pred = tree.predict(scene_features)
        if pred.ndim == 1:
            pred = pred[:, np.newaxis]
        for out_index, band_name in enumerate(target_names):
            predictions[band_name].append(np.asarray(pred[:, out_index], dtype=np.float64))
        used_realizations += 1
    if used_realizations == 0:
        return prior

    boa_new = boa.copy()
    unc_new = prior.boa_unc.copy()
    debias = debias or {}
    b01_floor = float(
        b01_uncertainty_floor if b01_uncertainty_floor is not None else uncertainty_floor
    )
    blend_weight = min(max(float(composite_blend_weight), 0.0), 1.0)
    composite_reference = (
        {
            band_name: _composite_reference_on_scene(
                comp,
                band_column=int(targets[band_name]),
                epsg=int(epsg),
                x=x,
                y=y,
                template=template,
                valid=valid,
            )
            for band_name in target_names
        }
        if blend_weight > 0.0
        else {}
    )

    for band_name in target_names:
        stack = np.stack(predictions[band_name], axis=0)
        ensemble = np.median(stack, axis=0)
        spread = np.median(np.abs(stack - ensemble[np.newaxis]), axis=0) * _MAD_TO_STD
        intercept, slope = debias.get(band_name, (0.0, 0.0))
        ensemble = ensemble + float(debias_scale) * (
            float(intercept) + float(slope) * float(anchor_aot)
        )
        floor = b01_floor if band_name == "B01" else float(uncertainty_floor)
        uncertainty = np.maximum(spread, floor)
        if blend_weight > 0.0:
            reference, reference_unc = composite_reference[band_name]
            finite_reference = np.isfinite(reference) & np.isfinite(reference_unc)
            ensemble = ensemble.copy()
            uncertainty = uncertainty.copy()
            ensemble[finite_reference] = (1.0 - blend_weight) * ensemble[
                finite_reference
            ] + blend_weight * reference[finite_reference]
            uncertainty[finite_reference] = np.maximum(
                np.sqrt(
                    ((1.0 - blend_weight) * uncertainty[finite_reference]) ** 2
                    + (blend_weight * reference_unc[finite_reference]) ** 2
                ),
                floor,
            )
        predictions[band_name] = [ensemble, uncertainty]

    for band_name in target_names:
        ensemble, uncertainty = predictions[band_name]
        image = np.full(flat_anchor.shape[0], np.nan, dtype=np.float64)
        sigma = np.full(flat_anchor.shape[0], np.nan, dtype=np.float64)
        image[valid] = ensemble
        sigma[valid] = uncertainty
        image_2d = image.reshape(anchor.shape[:2])
        sigma_2d = sigma.reshape(anchor.shape[:2])
        ok = np.isfinite(image_2d) & (image_2d > 0.001) & (image_2d < 0.6)
        current = boa.sel(band=band_name).values.copy()
        current_unc = unc_new.sel(band=band_name).values.copy()
        current[ok] = image_2d[ok]
        current_unc[ok] = sigma_2d[ok]
        boa_new.loc[{"band": band_name}] = xr.DataArray(
            current,
            dims=template.dims,
            coords=template.coords,
        )
        unc_new.loc[{"band": band_name}] = xr.DataArray(
            current_unc,
            dims=template.dims,
            coords=template.coords,
        )

    return replace(prior, boa=boa_new, boa_unc=unc_new)

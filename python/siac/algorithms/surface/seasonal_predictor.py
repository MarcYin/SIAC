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


def _full_template(template: xr.DataArray, value: float) -> xr.DataArray:
    return cast("xr.DataArray", xr.full_like(template, float(value)).astype(np.float32))


def _same_template_grid(source: xr.DataArray, template: xr.DataArray) -> bool:
    """Whether two spatial fields already share exactly the same grid."""
    if source.dims != template.dims or source.shape != template.shape:
        return False
    return all(
        name in source.coords
        and name in template.coords
        and np.array_equal(source.coords[name].values, template.coords[name].values)
        for name in template.dims
    )


def _field_on_template(
    source: xr.DataArray | None,
    template: xr.DataArray,
    *,
    fallback: float,
    label: str,
) -> xr.DataArray:
    """Align one physical field to the anchor grid without losing its structure.

    Atmospheric providers and Sentinel-2 angle grids can have different native
    resolutions.  A physical anchor correction must preserve those fields up to
    the predictor grid instead of first collapsing them to a scene statistic.
    A scalar fallback is retained only for malformed or wholly missing inputs.
    """
    fallback_value = float(fallback)
    if source is None:
        return _full_template(template, fallback_value)
    if "band" in source.dims:
        source = source.isel(band=0, drop=True)
    source = source.astype(np.float32)
    # Some same-scene solver fields are constructed from x/y coordinates but
    # lose their rioxarray CRS before an anchor-iteration pass.  Their grid is
    # still in the target scene's coordinate system, so retain that spatial
    # structure instead of silently reducing the field to a scalar fallback.
    if source.rio.crs is None and template.rio.crs is not None:
        source = source.rio.write_crs(template.rio.crs)
    values = np.asarray(source.values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size:
        fallback_value = float(np.median(finite))
    try:
        aligned = (
            source
            if _same_template_grid(source, template)
            else _as_template_grid(source, template).astype(np.float32)
        )
        aligned_values = np.asarray(aligned.values, dtype=np.float32)
    except Exception:  # noqa: BLE001 - preserve a usable anchor for sparse auxiliaries
        logger.warning(
            "Anchor correction could not align %s; using its scalar median fallback.",
            label,
            exc_info=True,
        )
        return _full_template(template, fallback_value)
    if not np.isfinite(aligned_values).any():
        logger.warning(
            "Anchor correction has no finite %s support on the predictor grid; using fallback.",
            label,
        )
        return _full_template(template, fallback_value)
    aligned_values = np.where(np.isfinite(aligned_values), aligned_values, fallback_value)
    return xr.DataArray(
        aligned_values.astype(np.float32),
        dims=template.dims,
        coords=template.coords,
    )


def _circular_angle_on_template(
    source: xr.DataArray,
    template: xr.DataArray,
    *,
    label: str,
) -> xr.DataArray:
    """Resample an azimuth field through sine/cosine to avoid 0/2pi wrapping."""
    sin_angle = _field_on_template(
        xr.apply_ufunc(np.sin, source), template, fallback=0.0, label=f"sin({label})"
    )
    cos_angle = _field_on_template(
        xr.apply_ufunc(np.cos, source), template, fallback=1.0, label=f"cos({label})"
    )
    values = np.mod(np.arctan2(sin_angle.values, cos_angle.values), 2.0 * np.pi)
    return xr.DataArray(values.astype(np.float32), dims=template.dims, coords=template.coords)


def _scene_mean_angle(field: xr.DataArray, *, circular: bool) -> xr.DataArray:
    """Broadcast one image-level angle while retaining the target grid."""
    values = np.asarray(field.values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        mean = 0.0
    elif circular:
        mean = float(np.mod(np.arctan2(np.mean(np.sin(finite)), np.mean(np.cos(finite))), 2.0 * np.pi))
    else:
        mean = float(np.mean(finite))
    return xr.full_like(field, np.float32(mean))


def _correct_anchor_reflectance(
    observation: ObservationBundle,
    *,
    atmo_prior: Any,
    rt_model: Any,
    template: xr.DataArray,
    anchor_grids: Mapping[str, xr.DataArray],
    valid: np.ndarray,
    anchor_aot: float,
    anchor_aot_field: xr.DataArray | None = None,
    scene_mean_geometry: bool = False,
) -> np.ndarray:
    """Correct scene anchor bands TOA->BOA with the configured RT backend.

    The anchor must be corrected in the same RT space the solver uses; there is
    deliberately no analytic fallback.
    """
    from siac.runtime import AtmosphericState, GeometryAngles

    geometry = observation.geometry
    atmo = AtmosphericState(
        aot=_field_on_template(
            anchor_aot_field, template, fallback=float(anchor_aot), label="anchor AOD"
        ),
        tcwv=_field_on_template(atmo_prior.tcwv, template, fallback=2.0, label="TCWV"),
        tco3=_field_on_template(atmo_prior.tco3, template, fallback=0.3, label="TCO3"),
        aot_unc=_field_on_template(
            atmo_prior.aot_unc, template, fallback=0.1, label="AOD uncertainty"
        ),
        tcwv_unc=_field_on_template(
            atmo_prior.tcwv_unc, template, fallback=0.5, label="TCWV uncertainty"
        ),
        tco3_unc=_field_on_template(
            atmo_prior.tco3_unc, template, fallback=0.05, label="TCO3 uncertainty"
        ),
        elevation=_field_on_template(
            atmo_prior.elevation, template, fallback=0.0, label="elevation"
        ),
    )
    if hasattr(geometry, "saa") and hasattr(geometry, "vaa"):
        saa = _circular_angle_on_template(geometry.saa, template, label="solar azimuth")
        vaa = _circular_angle_on_template(geometry.vaa, template, label="view azimuth")
    else:
        # ObservationBundle always carries absolute azimuths.  Retain this
        # compatibility branch for small injected test/legacy geometries that
        # only expose their relative azimuth.
        saa = _full_template(template, 0.0)
        vaa = _field_on_template(geometry.raa, template, fallback=0.0, label="relative azimuth")
    sza = _field_on_template(geometry.sza, template, fallback=0.0, label="solar zenith")
    vza = _field_on_template(geometry.vza, template, fallback=0.0, label="view zenith")
    if scene_mean_geometry:
        sza = _scene_mean_angle(sza, circular=False)
        saa = _scene_mean_angle(saa, circular=True)
        vza = _scene_mean_angle(vza, circular=False)
        vaa = _scene_mean_angle(vaa, circular=True)
    geom = GeometryAngles(sza=sza, saa=saa, vza=vza, vaa=vaa)
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


def _anchor_match_weights(
    historical_anchor: np.ndarray,
    target_anchor: np.ndarray,
    *,
    scale: float,
) -> np.ndarray:
    """Return stable per-pixel realization weights from NIR/SWIR distance."""
    if scale <= 0.0:
        raise ValueError("anchor_match_scale must be positive")
    if historical_anchor.ndim != 3 or target_anchor.ndim != 2:
        raise ValueError(
            "anchor arrays must have shapes (realization, pixel, band) and (pixel, band)"
        )
    if historical_anchor.shape[1:] != target_anchor.shape:
        raise ValueError("historical and target anchor shapes are incompatible")

    with np.errstate(invalid="ignore", over="ignore"):
        distance2 = np.mean(
            np.square((historical_anchor - target_anchor[np.newaxis]) / float(scale)),
            axis=2,
        )
        finite = np.isfinite(distance2)
        minimum = np.min(np.where(finite, distance2, np.inf), axis=0)
        weights = np.exp(-0.5 * (distance2 - minimum[np.newaxis]))
    return np.where(finite, weights, 0.0)


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted median over realizations for every flattened scene pixel."""
    if values.ndim != 2 or values.shape != weights.shape:
        raise ValueError("values and weights must have matching (realization, pixel) shapes")
    usable = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    effective_weights = np.where(usable, weights, 0.0)
    missing = np.sum(effective_weights, axis=0) <= 0.0
    effective_weights[:, missing] = np.where(np.isfinite(values[:, missing]), 1.0, 0.0)
    order = np.argsort(np.where(np.isfinite(values), values, np.inf), axis=0)
    ordered_values = np.take_along_axis(values, order, axis=0)
    ordered_weights = np.take_along_axis(effective_weights, order, axis=0)
    cumulative = np.cumsum(ordered_weights, axis=0)
    threshold = 0.5 * np.sum(ordered_weights, axis=0)
    index = np.argmax(cumulative >= threshold[np.newaxis], axis=0)
    result = np.take_along_axis(ordered_values, index[np.newaxis], axis=0)[0]
    return np.where(threshold > 0.0, result, np.nan)


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
    anchor_aot_field: xr.DataArray | None = None,
    target_band_columns: Mapping[str, int] | None = None,
    debias: Mapping[str, tuple[float, float]] | None = None,
    debias_scale: float = 1.0,
    uncertainty_floor: float = 0.006,
    b01_uncertainty_floor: float | None = None,
    min_samples_leaf: int = 5,
    random_state: int = 0,
    predictor_model: str = "extra_tree",
    robust_clip: float = 0.0,
    composite_blend_weight: float = 0.0,
    attach_tau_predictor: bool = False,
    ensemble_aggregation: str = "median",
    anchor_match_scale: float = 0.05,
    scene_mean_geometry: bool = False,
) -> SurfacePrior:
    """Replace visible prior bands with a seasonal ExtraTree ensemble prediction.

    Parameters mirror the validated research harness:

    - correct the scene ``B8A/B11/B12`` anchor TOA->BOA through ``rt_model``
      using the MAIAC/anchor AOD field plus aligned WVP, ozone, elevation and
      Sentinel-2 geometry (the same RT backend the solver uses). When
      ``scene_mean_geometry`` is enabled, one image-level geometry is used
      while atmospheric and elevation fields remain spatial;
    - train one ``ExtraTreeRegressor`` per seasonal realization;
    - features are corrected ``B8A/B11/B12`` plus the per-pixel mean visible
      climatology ``B01/B02/B03/B04``;
    - prediction is a median across realizations, optionally weighted by each
      realization's local NIR/SWIR similarity to the target scene;
    - uncertainty is MAD*1.4826 across realization predictions, floored.
    """
    from sklearn.tree import ExtraTreeRegressor

    if predictor_model == "extra_trees_20":
        from sklearn.ensemble import ExtraTreesRegressor
    elif predictor_model != "extra_tree":
        raise ValueError(f"Unknown predictor_model {predictor_model!r}")

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
    aggregation = str(ensemble_aggregation).strip().lower()
    if aggregation not in {"median", "anchor_weighted"}:
        raise ValueError(f"unsupported seasonal predictor aggregation {ensemble_aggregation!r}")
    if float(anchor_match_scale) <= 0.0:
        raise ValueError("anchor_match_scale must be positive")

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
        anchor_aot_field=anchor_aot_field,
        scene_mean_geometry=bool(scene_mean_geometry),
    )
    finite_corrected = np.all(np.isfinite(corrected_anchor), axis=1)
    if not np.any(finite_corrected):
        return prior
    if not np.all(finite_corrected):
        valid_indices = np.flatnonzero(valid)
        valid[valid_indices[~finite_corrected]] = False
        corrected_anchor = corrected_anchor[finite_corrected]

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
    fitted_trees: list[Any] = []
    realization_anchors: list[np.ndarray] = []
    used_realizations = 0
    for index in range(n_real):
        composite = comp[index].reshape(comp.shape[1], -1).T
        good = np.all(np.isfinite(composite), axis=1)
        if int(np.count_nonzero(good)) < 200:
            continue
        train_x = np.column_stack([composite[good][:, list(ANCHOR_COLUMNS)], localizer_comp[good]])
        train_y = composite[good][:, target_cols]
        if predictor_model == "extra_trees_20":
            # A 20-tree ensemble per realization: individually smoother
            # predictions shrink the across-realization spread that becomes
            # the prior sigma, sharpening the solver's cost curve.
            tree = ExtraTreesRegressor(
                n_estimators=20,
                min_samples_leaf=int(min_samples_leaf),
                random_state=int(random_state),
                n_jobs=1,
            ).fit(train_x, train_y)
        else:
            tree = ExtraTreeRegressor(
                random_state=int(random_state),
                min_samples_leaf=int(min_samples_leaf),
            ).fit(train_x, train_y)
        fitted_trees.append(tree)
        pred = tree.predict(scene_features)
        if pred.ndim == 1:
            pred = pred[:, np.newaxis]
        for out_index, band_name in enumerate(target_names):
            predictions[band_name].append(np.asarray(pred[:, out_index], dtype=np.float64))
        if aggregation == "anchor_weighted":
            anchor_parts = []
            for column in ANCHOR_COLUMNS:
                da = xr.DataArray(
                    comp[index, int(column)],
                    dims=("y", "x"),
                    coords={"y": y, "x": x},
                ).rio.write_crs(f"EPSG:{int(epsg)}")
                anchor_parts.append(_as_template_grid(da, template).values.ravel()[valid])
            realization_anchors.append(np.column_stack(anchor_parts))
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

    anchor_weights = None
    if aggregation == "anchor_weighted":
        anchor_weights = _anchor_match_weights(
            np.stack(realization_anchors, axis=0),
            corrected_anchor,
            scale=float(anchor_match_scale),
        )

    for band_name in target_names:
        stack = np.stack(predictions[band_name], axis=0)
        if anchor_weights is None:
            ensemble = np.median(stack, axis=0)
            spread = np.median(np.abs(stack - ensemble[np.newaxis]), axis=0) * _MAD_TO_STD
        else:
            ensemble = _weighted_median(stack, anchor_weights)
            spread = (
                _weighted_median(np.abs(stack - ensemble[np.newaxis]), anchor_weights) * _MAD_TO_STD
            )
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

    tau_payload = None
    if attach_tau_predictor:
        # Full localizer planes on the prior template (banded, so M4 can
        # resample with the standard machinery), plus the fitted trees. M5 uses
        # these to re-predict the visible prior at EACH candidate AOD.
        localizer_planes = []
        for index in range(len(VISIBLE_COLUMNS)):
            da = xr.DataArray(
                mean_visible[index],
                dims=("y", "x"),
                coords={"y": y, "x": x},
            ).rio.write_crs(f"EPSG:{int(epsg)}")
            localizer_planes.append(_as_template_grid(da, template))
        localizer_da = xr.concat(localizer_planes, dim="band").assign_coords(
            band=[f"loc{i}" for i in range(len(VISIBLE_COLUMNS))]
        )
        tau_payload = {
            "trees": fitted_trees,
            "localizer": localizer_da,
            "anchor_bands": tuple(ANCHOR_BANDS),
            "target_bands": tuple(target_names),
            "debias": {name: tuple(debias[name]) for name in debias},
            "debias_scale": float(debias_scale),
            "anchor_aot": float(anchor_aot),
        }
    return replace(prior, boa=boa_new, boa_unc=unc_new, tau_predictor=tau_payload)

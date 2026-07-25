"""Reimplementation of the acixThree surface-prior approach for the lowcloud152 tests.

Faithful to the FINAL (active) code paths of ``~/acixThree`` — the original code
is not modified or imported. Sources: ``aux_data/get_candidate_images.py``
(GEE best-pixel monthly mosaics) and ``get_surfacereflectance_model.py``
(kNN predictor + self-prediction RMSE), config ``config/config.yaml``
(QA_BAND=cs_cdf, CLEAR_THRESHOLD=0.60, k_neighbours=10).

Replicated design:
- S2_SR_HARMONIZED (Sen2Cor L2A BOA, as-is — no re-AC) linked with Cloud
  Score+ ``cs_cdf``; date windows = +/-45 days around the 15th of the obs
  month, years 2018-2023.
- Daily clean mask = cs_cdf >= 0.60 (500 m focal-min erosion) AND
  SCL not in {1, 3, 8, 10}; candidate days = clean coverage >= month mean.
- Per-day weight = sigmoid(same-day MCD19 AOD, threshold 60th pct, raw int
  units) + sigmoid(clean coverage, center 0.5); per-pixel mosaic weight =
  (cs_cdf - 0.6)/0.4 + day weight + AOI coverage ratio;
  ``qualityMosaic('weight')`` per calendar month of candidate days.
- Model: kNN (k=10, L2) over all pixel x month mosaic samples; features =
  aerosol-transparent bands + per-pixel climatology (temporal mean of the
  visible targets); prediction = mean of the 10 neighbours' absolute
  reflectances; sigma = per-pixel per-band temporal RMSE of mosaic
  self-prediction. The solve-time cost in acixThree used ABSOLUTE reflectance
  (ratios only build the prior + sigma), replicated accordingly.

Declared adaptations (necessary for the current tests, not behaviour changes):
- Scene-day features use B8A instead of B6/B8 because the Phase-D calibration
  dumps only carry B8A/B11/B12 TOA; the library uses the same bands so the
  feature space is train/inference-consistent.
- faiss GPU -> sklearn NearestNeighbors (pure infrastructure substitution).
- Scene-day NIR/SWIR BOA "AC'd at the prior AOD" is reconstructed from the
  calibration dump's TOA + RT coefficient curves at the stored anchor AOD
  (same role as acixThree's NN AC at prior AOD).
- Daily clean coverage and MCD19 AOD are reduced server-side to per-day
  scalars (coverage fraction / best-QA median) instead of downloading daily
  mask rasters; the weights consume identical quantities.
- Visible targets include B01 (coastal) alongside B02/B03/B04.

Stages: ``mosaics`` (GEE fetch, cached npz), ``predict`` (kNN + truth-referenced
evaluation vs the current et20 prior), ``report`` (aggregate CSV/figures).

L1C variant (``mosaics_l1c``/``predict_l1c``) — the FAITHFUL library space:
acixThree mosaics L1C TOA (S2_HARMONIZED) and applies its OWN atmospheric
correction at the per-pixel winning day's MCD19 AOD (Sen2Cor L2A is only used
for the SCL clean mask). Replicated here by carrying per-image provenance
bands (day MCD19 AOD, day CAMS ozone, scene mean angles, linked L2A WVP)
through ``qualityMosaic`` and correcting each realization offline with the
production ``ZarrLUTBackend`` — the same RT space as the solve and the truth
reference, so this variant has no Sen2Cor-space handicap. Adaptation: the
continental-average LUT stands in for acixThree's species-aware emulators
(their extra axis), and Sentinel-2A RSRFs are used for all mosaic pixels.
"""

from __future__ import annotations

import argparse
import calendar
import datetime as dt
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_DUMPS = ROOT / "calib_dumps_crossrt_dt_20260719"
DEFAULT_MOSAICS = ROOT / "acix3_style_mosaics_20260719"
DEFAULT_OUTPUT = ROOT / "analysis/acix3_style_prior_20260719"
DEFAULT_RESULTS = (
    ROOT
    / "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_dt_20260718"
)

QA_BAND = "cs_cdf"
CLEAR_THRESHOLD = 0.60
K_NEIGHBOURS = 10
LIBRARY_YEARS = range(2018, 2024)
WINDOW_HALF_DAYS = 45
TARGET_BANDS = ("B1", "B2", "B3", "B4")  # GEE names; == B01, B02, B03, B04
FEATURE_BANDS = ("B8A", "B11", "B12")
ALL_BANDS = TARGET_BANDS + FEATURE_BANDS
DUMP_BAND_BY_GEE = {"B1": "B01", "B2": "B02", "B3": "B03", "B4": "B04"}


def _ee():
    import ee

    key = os.environ["GEE_SERVICE_ACCOUNT_KEY"]
    account = os.environ["GEE_SERVICE_ACCOUNT"]
    credentials = ee.ServiceAccountCredentials(account, key)
    ee.Initialize(
        credentials, opt_url="https://earthengine-highvolume.googleapis.com"
    )
    return ee


def _window_filters(ee, obs_date: dt.date):
    filters = []
    for year in LIBRARY_YEARS:
        middle = dt.datetime(year, obs_date.month, 15)
        start = middle - dt.timedelta(days=WINDOW_HALF_DAYS)
        end = middle + dt.timedelta(days=WINDOW_HALF_DAYS)
        filters.append(ee.Filter.date(start.isoformat(), end.isoformat()))
    return ee.Filter.Or(*filters)


def _clean_mask(ee, image):
    clean = image.select(QA_BAND).gte(CLEAR_THRESHOLD)
    clean = clean.focal_min(500, "circle", "meters")
    scl = image.select("SCL")
    scl_mask = scl.neq(1).And(scl.neq(3)).And(scl.neq(8)).And(scl.neq(10))
    return clean.And(scl_mask)


def _day_stats(ee, geometry, day: str):
    """Per-day clean-coverage fraction and best-QA MCD19 AOD median (raw ints)."""
    date = ee.Date(day)
    cs_plus = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    s2_sur = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geometry)
        .filterDate(date, date.advance(1, "day"))
        .linkCollection(cs_plus, [QA_BAND])
    )
    # Days present in L1C but absent from L2A (no SCL) yield an empty
    # collection; they get coverage 0 and can never become candidates —
    # matching the original code's dependence on the L2A SCL mask.
    clean = ee.Image(
        ee.Algorithms.If(
            s2_sur.size().gt(0),
            s2_sur.map(lambda img: _clean_mask(ee, img)).median().unmask(0).rename("m"),
            ee.Image.constant(0).clip(geometry).rename("m"),
        )
    )
    coverage = clean.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=geometry, scale=120, maxPixels=1e9
    ).get("m")

    def best_aot(image):
        qa = image.select("AOD_QA")
        qa_mask = (1 << 8) | (1 << 9) | (1 << 10) | (1 << 11)
        aod_qa = qa.bitwiseAnd(qa_mask).rightShift(8)
        return image.updateMask(aod_qa.neq(5))

    noon = date.advance(12, "hour")
    mcd19 = (
        ee.ImageCollection("MODIS/061/MCD19A2_GRANULES")
        .filterBounds(geometry)
        .filter(ee.Filter.date(noon.advance(-16, "hour"), noon.advance(16, "hour")))
        .map(best_aot)
        .select("Optical_Depth_055")
        .median()
    )
    aod = mcd19.reduceRegion(
        reducer=ee.Reducer.median(), geometry=geometry, scale=1000, maxPixels=1e9
    ).get("Optical_Depth_055", None)
    packed = ee.Dictionary({"coverage": coverage, "aod": aod})
    values = packed.getInfo()
    coverage_value = float(values.get("coverage") or 0.0)
    aod_value = values.get("aod")
    return coverage_value, (float(aod_value) if aod_value is not None else np.nan)


def _aod_weight(daily_aod: np.ndarray, threshold_percentile: float = 60.0) -> np.ndarray:
    aod_min = np.nanmin(daily_aod)
    aod_max = np.nanpercentile(daily_aod, threshold_percentile)
    if aod_max == aod_min:
        weight = np.isfinite(daily_aod).astype("float32")
    else:
        weight = 1 - (
            1
            / (
                1
                + np.exp(
                    -0.2 * (np.minimum(daily_aod, aod_max) - (aod_max - aod_min) / 2)
                )
            )
        )
    weight[~np.isfinite(daily_aod)] = np.nan
    return weight


def _coverage_weight(coverage: np.ndarray, minimum: float = 0.5) -> np.ndarray:
    return 1 / (1 + np.exp(-20 * (coverage - minimum)))


def _grid_request(dump: dict[str, np.ndarray]) -> tuple[dict, tuple[int, int]]:
    height, width = (int(v) for v in dump["template_shape"])
    a, b, c, d, e, f = (float(v) for v in dump["template_transform"])
    crs = str(dump["template_crs"])
    grid = {
        "dimensions": {"width": width, "height": height},
        "affineTransform": {
            "scaleX": a,
            "shearX": b,
            "translateX": c,
            "shearY": d,
            "scaleY": e,
            "translateY": f,
        },
        "crsCode": crs,
    }
    return grid, (height, width)


def build_mosaics(matchup_id: str, dumps_dir: Path, output_dir: Path) -> Path:
    import ee as ee_module  # noqa: F401  (import check before auth)

    ee = _ee()
    with np.load(dumps_dir / f"{matchup_id}.npz") as handle:
        dump = {k: handle[k] for k in ("template_shape", "template_transform", "template_crs")}
    grid, (height, width) = _grid_request(dump)

    day_token = matchup_id.split("_")[-1][:8]
    obs_date = dt.date(int(day_token[:4]), int(day_token[4:6]), int(day_token[6:8]))

    transform = grid["affineTransform"]
    xs = transform["translateX"] + transform["scaleX"] * np.array([0, width])
    ys = transform["translateY"] + transform["scaleY"] * np.array([0, height])
    rect = ee.Geometry.Rectangle(
        [float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())],
        proj=grid["crsCode"],
        evenOdd=False,
    )
    geometry = rect.transform("EPSG:4326", 1)

    cs_plus = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geometry)
        .filter(_window_filters(ee, obs_date))
        .linkCollection(cs_plus, [QA_BAND])
    )
    times = s2.aggregate_array("system:time_start").getInfo()
    days = sorted({dt.datetime.utcfromtimestamp(t / 1000).strftime("%Y-%m-%d") for t in times})
    if not days:
        raise SystemExit(f"no S2_SR days for {matchup_id}")

    with ThreadPoolExecutor(max_workers=8) as pool:
        stats = list(pool.map(lambda d: _day_stats(ee, geometry, d), days))
    coverage = np.array([s[0] for s in stats])
    daily_aod = np.array([s[1] for s in stats])

    # Candidate days: clean coverage at or above the month's mean (and > 0).
    months = np.array([d[:7] for d in days])
    candidate = np.zeros(len(days), dtype=bool)
    for month in np.unique(months):
        month_mask = months == month
        month_mean = np.nanmean(coverage[month_mask])
        candidate[month_mask] = (coverage[month_mask] >= month_mean) & (
            coverage[month_mask] > 0
        )
    day_weight = _aod_weight(daily_aod) + _coverage_weight(coverage)
    day_weight[~np.isfinite(day_weight)] = 0.0

    mosaics, mosaic_months = [], []
    for month in np.unique(months[candidate]):
        month_days = [d for d, c, m in zip(days, candidate, months) if c and m == month]
        weights = {
            d: float(w)
            for d, w, c, m in zip(days, day_weight, candidate, months)
            if c and m == month
        }
        weight_dict = ee.Dictionary(weights)
        date_filter = ee.Filter.Or(
            *[ee.Filter.date(d, ee.Date(d).advance(1, "day")) for d in month_days]
        )
        month_collection = s2.filter(date_filter)

        def add_weight(image):
            key = ee.Date(image.get("system:time_start")).format("YYYY-MM-dd")
            day_w = ee.Number(weight_dict.get(key, 0))
            weight = (
                image.select(QA_BAND).subtract(CLEAR_THRESHOLD).divide(1 - CLEAR_THRESHOLD)
            ).add(ee.Image.constant(day_w).toFloat())
            intersection = image.geometry().intersection(geometry, maxError=1)
            ratio = intersection.area(maxError=1).divide(geometry.area(maxError=1))
            weight = weight.add(ee.Image.constant(ratio).toFloat()).rename("weight")
            masked = image.updateMask(_clean_mask(ee, image))
            return masked.addBands(weight)

        mosaic = (
            month_collection.map(add_weight)
            .qualityMosaic("weight")
            .select(list(ALL_BANDS))
            .toFloat()
        )
        pixels = ee.data.computePixels(
            {
                "expression": mosaic,
                "fileFormat": "NUMPY_NDARRAY",
                "grid": grid,
            }
        )
        stack = np.stack([pixels[band].astype(np.float32) for band in ALL_BANDS])
        # computePixels renders masked pixels as 0 or +/-inf depending on band.
        stack[~np.isfinite(stack) | (stack == 0)] = np.nan
        stack /= 10000.0
        mosaics.append(stack)
        mosaic_months.append(month)

    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{matchup_id}.npz"
    np.savez_compressed(
        out,
        mosaics=np.stack(mosaics),
        months=np.array(mosaic_months),
        band_names=np.array(ALL_BANDS),
        days=np.array(days),
        coverage=coverage,
        daily_aod=daily_aod,
        candidate=candidate,
        day_weight=day_weight,
    )
    print(f"{matchup_id}: {len(mosaic_months)} monthly mosaics from {int(candidate.sum())} candidate days")
    return out


def _boa_at_aod(dump: dict[str, np.ndarray], band: str, aod: float) -> np.ndarray:
    grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    xap = float(np.interp(aod, grid, np.asarray(dump[f"xap_{band}"], dtype=np.float64)))
    xbp = float(np.interp(aod, grid, np.asarray(dump[f"xbp_{band}"], dtype=np.float64)))
    xcp = float(np.interp(aod, grid, np.asarray(dump[f"xcp_{band}"], dtype=np.float64)))
    toa = np.asarray(dump[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    return y / (1.0 + xcp * y)


def predict_case(
    matchup_id: str,
    mosaics_dir: Path,
    dumps_dir: Path,
    results_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    with np.load(mosaics_dir / f"{matchup_id}.npz") as handle:
        mosaics = handle["mosaics"].astype(np.float64)  # (T, B, H, W)
        band_names = [str(b) for b in handle["band_names"]]
    with np.load(dumps_dir / f"{matchup_id}.npz") as handle:
        dump = {k: handle[k] for k in handle.files}
    result = json.loads((results_dir / f"{matchup_id}.json").read_text())
    return _knn_evaluate(
        matchup_id, mosaics, band_names, dump, result, output_dir, tag=""
    )


def _knn_evaluate(
    matchup_id: str,
    mosaics: np.ndarray,
    band_names: list[str],
    dump: dict[str, np.ndarray],
    result: dict[str, object],
    output_dir: Path,
    *,
    tag: str,
) -> dict[str, object]:
    from sklearn.neighbors import NearestNeighbors

    truth_aod = float(result["truth"])
    anchor_aod = float(dump["anchor_aot"])

    # Guard against +/-inf sentinels in cached mosaics (nanmean propagates inf).
    mosaics = np.where(np.isfinite(mosaics), mosaics, np.nan)
    n_time, n_band, height, width = mosaics.shape
    target_idx = [band_names.index(b) for b in TARGET_BANDS]
    feature_idx = [band_names.index(b) for b in FEATURE_BANDS]

    # Climatology: per-pixel temporal mean of the visible targets.
    climatology = np.nanmean(mosaics[:, target_idx], axis=0)  # (4, H, W)

    # Library samples: features (aerosol-transparent bands + climatology).
    features = np.concatenate(
        [
            mosaics[:, feature_idx].transpose(0, 2, 3, 1).reshape(-1, len(feature_idx)),
            np.broadcast_to(
                climatology.transpose(1, 2, 0), (n_time, height, width, len(target_idx))
            ).reshape(-1, len(target_idx)),
        ],
        axis=1,
    )
    targets = mosaics.transpose(0, 2, 3, 1).reshape(-1, n_band)
    valid = np.isfinite(features).all(axis=1) & np.isfinite(targets).all(axis=1)
    library_x, library_y = features[valid], targets[valid]
    if library_x.shape[0] < 100:
        raise SystemExit(f"library too small for {matchup_id}: {library_x.shape[0]}")

    index = NearestNeighbors(n_neighbors=K_NEIGHBOURS).fit(library_x)

    # Self-prediction over the mosaic stack -> per-pixel per-band temporal RMSE.
    _, neighbour_idx = index.kneighbors(library_x)
    self_pred = library_y[neighbour_idx].mean(axis=1)
    residual_sq = np.full((features.shape[0], n_band), np.nan)
    residual_sq[valid] = (self_pred - library_y) ** 2
    rmse = np.sqrt(
        np.nanmean(residual_sq.reshape(n_time, height, width, n_band), axis=0)
    ).transpose(2, 0, 1)  # (B, H, W)

    # Scene day: anchors AC'd at the prior AOD (mirrors acixThree AC-at-prior).
    scene_features = np.stack(
        [_boa_at_aod(dump, band, anchor_aod) for band in FEATURE_BANDS]
        + [climatology[i] for i in range(len(target_idx))]
    ).transpose(1, 2, 0)
    scene_valid = np.isfinite(scene_features).all(axis=2)
    prediction = np.full((n_band, height, width), np.nan)
    if scene_valid.any():
        _, scene_idx = index.kneighbors(scene_features[scene_valid])
        prediction[:, scene_valid] = library_y[scene_idx].mean(axis=1).T

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / f"{matchup_id}_prior{tag}.npz",
        prediction=prediction.astype(np.float32),
        rmse=rmse.astype(np.float32),
        band_names=np.array(band_names),
        anchor_aod=anchor_aod,
    )

    # Truth-referenced comparison vs the current et20 prior in the same dump.
    row: dict[str, object] = {
        "matchup_id": matchup_id,
        "site": matchup_id.split("__")[0],
        "truth_aod": round(truth_aod, 4),
        "anchor_aod": round(anchor_aod, 4),
        "library_samples": int(library_x.shape[0]),
        "n_mosaics": n_time,
    }
    for gee_band in TARGET_BANDS:
        dump_band = DUMP_BAND_BY_GEE[gee_band]
        band_i = band_names.index(gee_band)
        acix_pred = prediction[band_i]
        truth_key = f"toa_{dump_band}"
        if truth_key not in dump:
            continue
        target = _boa_at_aod(dump, dump_band, truth_aod)
        current_key = f"pred_{dump_band}"
        current_pred = dump[current_key] if current_key in dump else None
        mask = np.isfinite(acix_pred) & np.isfinite(target) & (target > -0.05) & (target < 1.0)
        if current_pred is not None:
            mask &= np.isfinite(current_pred)
        if mask.sum() < 50:
            continue
        diff_acix = acix_pred[mask] - target[mask]
        row[f"{dump_band}_acix_bias"] = round(float(np.median(diff_acix)), 5)
        row[f"{dump_band}_acix_mae"] = round(float(np.mean(np.abs(diff_acix))), 5)
        row[f"{dump_band}_acix_rmse"] = round(float(np.sqrt(np.mean(diff_acix**2))), 5)
        row[f"{dump_band}_acix_sigma_median"] = round(float(np.nanmedian(rmse[band_i])), 5)
        if current_pred is not None:
            diff_cur = np.asarray(current_pred, dtype=np.float64)[mask] - target[mask]
            row[f"{dump_band}_current_bias"] = round(float(np.median(diff_cur)), 5)
            row[f"{dump_band}_current_mae"] = round(float(np.mean(np.abs(diff_cur))), 5)
            row[f"{dump_band}_current_rmse"] = round(float(np.sqrt(np.mean(diff_cur**2))), 5)
    suffix = f"_metrics{tag}.json" if tag else "_metrics.json"
    (output_dir / f"{matchup_id}{suffix}").write_text(json.dumps(row, indent=2))
    print(json.dumps(row))
    return row


KG_M2_PER_DU = 2.1415e-5
CAMS_DIR = Path("/gws/ssde/j25b/nceo_ard/public/cams")


def _cams_tco3_atm_cm(day: str, lat: float, lon: float) -> float:
    """Per-day CAMS total-ozone column (atm-cm) at the site's overpass hour."""
    import xarray as xr

    path = CAMS_DIR / f"{day}.nc"
    if not path.exists():
        return 0.30
    with xr.open_dataset(path) as data:
        field = data[["gtco3"]]
        if "forecast_reference_time" in field.dims:
            field = field.squeeze("forecast_reference_time", drop=True)
        hour = int(round(10.5 - lon / 15.0)) % 24
        if "forecast_period" in field.dims:
            selector: dict = {"forecast_period": float(hour)}
        else:
            selector = {"time": np.datetime64(f"{day}T{hour:02d}:00")}
        value = float(
            field.sel(latitude=lat, longitude=lon, **selector, method="nearest")["gtco3"].item()
        )
    tco3 = value / KG_M2_PER_DU / 1000.0
    return tco3 if 0.15 <= tco3 <= 0.55 else 0.30


L1C_PROVENANCE = ("WVP", "AOD", "TO3", "sza", "saa", "vza", "vaa")


def build_mosaics_l1c(
    matchup_id: str, dumps_dir: Path, results_dir: Path, output_dir: Path
) -> Path:
    """Faithful acixThree library: L1C TOA qualityMosaic with provenance bands."""
    ee = _ee()
    with np.load(dumps_dir / f"{matchup_id}.npz") as handle:
        dump = {k: handle[k] for k in ("template_shape", "template_transform", "template_crs")}
    grid, (height, width) = _grid_request(dump)
    result = json.loads((results_dir / f"{matchup_id}.json").read_text())
    lat, lon = float(result["lat"]), float(result["lon"])

    day_token = matchup_id.split("_")[-1][:8]
    obs_date = dt.date(int(day_token[:4]), int(day_token[4:6]), int(day_token[6:8]))

    transform = grid["affineTransform"]
    xs = transform["translateX"] + transform["scaleX"] * np.array([0, width])
    ys = transform["translateY"] + transform["scaleY"] * np.array([0, height])
    rect = ee.Geometry.Rectangle(
        [float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())],
        proj=grid["crsCode"],
        evenOdd=False,
    )
    geometry = rect.transform("EPSG:4326", 1)

    cs_plus = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    l2a = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
        .filterBounds(geometry)
        .filter(_window_filters(ee, obs_date))
        .linkCollection(cs_plus, [QA_BAND])
        .linkCollection(l2a, ["SCL", "WVP"])
    )

    # Some processing baselines store linked SCL/WVP with different pixel
    # types; a mixed collection fails computePixels' homogeneity check.
    # Casting is lossless (SCL is a class code, WVP already scaled); images
    # without an L2A match keep their band-less state.
    def _harmonize_linked_types(image):
        image = ee.Image(image)
        names = image.bandNames()
        image = ee.Image(
            ee.Algorithms.If(
                names.contains("SCL"),
                image.addBands(image.select("SCL").toUint8(), None, True),
                image,
            )
        )
        return ee.Image(
            ee.Algorithms.If(
                names.contains("WVP"),
                image.addBands(image.select("WVP").toFloat(), None, True),
                image,
            )
        )

    s2 = s2.map(_harmonize_linked_types)
    times = s2.aggregate_array("system:time_start").getInfo()
    days = sorted({dt.datetime.utcfromtimestamp(t / 1000).strftime("%Y-%m-%d") for t in times})
    if not days:
        raise SystemExit(f"no L1C days for {matchup_id}")

    with ThreadPoolExecutor(max_workers=8) as pool:
        stats = list(pool.map(lambda d: _day_stats(ee, geometry, d), days))
    coverage = np.array([s[0] for s in stats])
    daily_aod = np.array([s[1] for s in stats])
    daily_tco3 = np.array([_cams_tco3_atm_cm(d, lat, lon) for d in days])

    months = np.array([d[:7] for d in days])
    candidate = np.zeros(len(days), dtype=bool)
    for month in np.unique(months):
        month_mask = months == month
        month_mean = np.nanmean(coverage[month_mask])
        candidate[month_mask] = (coverage[month_mask] >= month_mean) & (
            coverage[month_mask] > 0
        )
    # The AC needs a finite per-day AOD; a candidate day without MCD19 keeps a
    # weight of zero and the regional low-AOD fallback 0.1 for correction.
    day_weight = _aod_weight(daily_aod) + _coverage_weight(coverage)
    day_weight[~np.isfinite(day_weight)] = 0.0
    aod_for_ac = np.where(np.isfinite(daily_aod), daily_aod / 1000.0, 0.1)

    fetch_bands = list(ALL_BANDS) + list(L1C_PROVENANCE)
    mosaics, mosaic_months = [], []
    for month in np.unique(months[candidate]):
        month_days = [d for d, c, m in zip(days, candidate, months) if c and m == month]
        info = {
            d: (float(w), float(a), float(o))
            for d, w, a, o, c, m in zip(
                days, day_weight, aod_for_ac, daily_tco3, candidate, months
            )
            if c and m == month
        }
        weight_dict = ee.Dictionary({d: v[0] for d, v in info.items()})
        aod_dict = ee.Dictionary({d: v[1] for d, v in info.items()})
        tco3_dict = ee.Dictionary({d: v[2] for d, v in info.items()})
        date_filter = ee.Filter.Or(
            *[ee.Filter.date(d, ee.Date(d).advance(1, "day")) for d in month_days]
        )

        def add_weight(image):
            key = ee.Date(image.get("system:time_start")).format("YYYY-MM-dd")
            day_w = ee.Number(weight_dict.get(key, 0))
            weight = (
                image.select(QA_BAND).subtract(CLEAR_THRESHOLD).divide(1 - CLEAR_THRESHOLD)
            ).add(ee.Image.constant(day_w).toFloat())
            intersection = image.geometry().intersection(geometry, maxError=1)
            ratio = intersection.area(maxError=1).divide(geometry.area(maxError=1))
            weight = weight.add(ee.Image.constant(ratio).toFloat()).rename("weight")
            provenance = ee.Image.cat(
                ee.Image.constant(ee.Number(aod_dict.get(key, 0.1))).toFloat().rename("AOD"),
                ee.Image.constant(ee.Number(tco3_dict.get(key, 0.3))).toFloat().rename("TO3"),
                ee.Image.constant(
                    ee.Number(image.get("MEAN_SOLAR_ZENITH_ANGLE"))
                ).toFloat().rename("sza"),
                ee.Image.constant(
                    ee.Number(image.get("MEAN_SOLAR_AZIMUTH_ANGLE"))
                ).toFloat().rename("saa"),
                ee.Image.constant(
                    ee.Number(image.get("MEAN_INCIDENCE_ZENITH_ANGLE_B2"))
                ).toFloat().rename("vza"),
                ee.Image.constant(
                    ee.Number(image.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B2"))
                ).toFloat().rename("vaa"),
            )
            masked = image.updateMask(_clean_mask(ee, image))
            return masked.addBands(provenance).addBands(weight)

        mosaic = (
            s2.filter(date_filter)
            .map(add_weight)
            .qualityMosaic("weight")
            .select(fetch_bands)
            .toFloat()
        )
        pixels = ee.data.computePixels(
            {"expression": mosaic, "fileFormat": "NUMPY_NDARRAY", "grid": grid}
        )
        stack = np.stack([pixels[band].astype(np.float32) for band in fetch_bands])
        # computePixels renders masked pixels as 0 or +/-inf depending on band.
        stack[~np.isfinite(stack)] = np.nan
        reflectance = stack[: len(ALL_BANDS)]
        reflectance[reflectance == 0] = np.nan
        reflectance /= 10000.0
        stack[len(ALL_BANDS)] /= 1000.0  # WVP scale factor -> cm
        mosaics.append(stack)
        mosaic_months.append(month)

    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{matchup_id}_l1c.npz"
    np.savez_compressed(
        out,
        mosaics=np.stack(mosaics),
        months=np.array(mosaic_months),
        band_names=np.array(fetch_bands),
        days=np.array(days),
        daily_aod=daily_aod,
        candidate=candidate,
    )
    print(
        f"{matchup_id}: {len(mosaic_months)} L1C monthly mosaics from "
        f"{int(candidate.sum())} candidate days"
    )
    return out


def _lut_correct_realization(
    rt_backend,
    sensor,
    stack: np.ndarray,
    band_names: list[str],
    elevation_km: float,
) -> np.ndarray:
    """Own-AC of one TOA mosaic realization at per-pixel provenance state."""
    import xarray as xr

    from siac.runtime import AtmosphericState, GeometryAngles

    index = {name: i for i, name in enumerate(band_names)}

    def field(name: str, fallback: float) -> xr.DataArray:
        values = np.asarray(stack[index[name]], dtype=np.float32)
        finite = np.isfinite(values) & (values != 0)
        fill = float(np.nanmedian(values[finite])) if finite.any() else fallback
        return xr.DataArray(
            np.where(finite, values, fill).astype(np.float32), dims=("y", "x")
        )

    template = field("AOD", 0.1)
    geometry = GeometryAngles.from_degrees(
        field("sza", 30.0), field("saa", 0.0), field("vza", 5.0), field("vaa", 0.0)
    )
    atmo = AtmosphericState(
        aot=field("AOD", 0.1),
        tcwv=field("WVP", 2.0),
        tco3=field("TO3", 0.3),
        aot_unc=xr.full_like(template, 0.1),
        tcwv_unc=xr.full_like(template, 0.5),
        tco3_unc=xr.full_like(template, 0.05),
        elevation=xr.full_like(template, np.float32(elevation_km)),
    )
    gee_to_s2 = {"B1": "B01", "B2": "B02", "B3": "B03", "B4": "B04"}
    corrected = np.full((len(ALL_BANDS), *stack.shape[1:]), np.nan, dtype=np.float32)
    sensor_bands = [
        sensor.get_band(gee_to_s2.get(name, name)) for name in ALL_BANDS
    ]
    rt_backend.preload_scene_subset(geometry, atmo, sensor_bands)
    for i, band in enumerate(sensor_bands):
        toa = xr.DataArray(
            np.asarray(stack[index[ALL_BANDS[i]]], dtype=np.float32), dims=("y", "x")
        )
        corrected[i] = np.asarray(
            rt_backend.compute_coefficients(geometry, atmo, band)
            .apply_correction(toa)
            .values,
            dtype=np.float32,
        )
    return corrected


def predict_l1c(
    matchup_id: str,
    mosaics_dir: Path,
    dumps_dir: Path,
    results_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    from siac.adapters.rsrf import load_sensor_config_with_rsrf
    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    with np.load(mosaics_dir / f"{matchup_id}_l1c.npz") as handle:
        stacks = handle["mosaics"].astype(np.float32)
        band_names = [str(b) for b in handle["band_names"]]
    with np.load(dumps_dir / f"{matchup_id}.npz") as handle:
        dump = {k: handle[k] for k in handle.files}
    result = json.loads((results_dir / f"{matchup_id}.json").read_text())

    pairs = ROOT / "analysis/l2a_l1c_physical_pairs_v3_deterrain_lowcloud152_20260718"
    elevation_km = 0.0
    pair_npz = pairs / f"{matchup_id}.npz"
    if pair_npz.exists():
        with np.load(pair_npz) as handle:
            if "terrain_elevation_km" in handle.files:
                values = handle["terrain_elevation_km"]
                if np.isfinite(values).any():
                    elevation_km = float(np.nanmedian(values))

    rt_backend = ZarrLUTBackend(
        Path(os.environ["PHASE_D_LUT_PATH"])
        if os.environ.get("PHASE_D_LUT_PATH")
        else ROOT.parent / "libradtran_continental_average_lut_1nm.zarr.zip"
    )
    sensor = load_sensor_config_with_rsrf("MSI", "S2A")
    corrected = np.stack(
        [
            _lut_correct_realization(rt_backend, sensor, stack, band_names, elevation_km)
            for stack in stacks
        ]
    )
    return _knn_evaluate(
        matchup_id,
        corrected.astype(np.float64),
        list(ALL_BANDS),
        dump,
        result,
        output_dir,
        tag="_l1c",
    )


def report(output_dir: Path, variant: str = "") -> None:
    pattern = f"*_metrics{variant}.json" if variant else "*_metrics.json"
    rows = []
    for path in sorted(output_dir.glob(pattern)):
        rows.append(json.loads(path.read_text()))
    if not rows:
        raise SystemExit("no metrics to aggregate")
    import csv

    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with (output_dir / f"acix3{variant}_vs_current.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        f"# acixThree-style prior{variant or ''} vs current et20 prior (truth-referenced)",
        "",
    ]
    lines.append("| band | acix3 med bias | acix3 med MAE | current med bias | current med MAE | acix3 med sigma |")
    lines.append("|---|---|---|---|---|---|")
    for band in ("B01", "B02", "B03", "B04"):
        acix_bias = [r[f"{band}_acix_bias"] for r in rows if f"{band}_acix_bias" in r]
        if not acix_bias:
            continue
        med = lambda vals: float(np.median(vals))  # noqa: E731
        acix_mae = [r[f"{band}_acix_mae"] for r in rows if f"{band}_acix_mae" in r]
        cur_bias = [r[f"{band}_current_bias"] for r in rows if f"{band}_current_bias" in r]
        cur_mae = [r[f"{band}_current_mae"] for r in rows if f"{band}_current_mae" in r]
        sigma = [r[f"{band}_acix_sigma_median"] for r in rows if f"{band}_acix_sigma_median" in r]
        current_cells = (
            f"{med(cur_bias):+.4f} | {med(cur_mae):.4f}" if cur_bias else "n/a | n/a"
        )
        lines.append(
            f"| {band} | {med(acix_bias):+.4f} | {med(acix_mae):.4f} "
            f"| {current_cells} | {med(sigma):.4f} |"
        )
    (output_dir / f"summary{variant}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=("mosaics", "predict", "mosaics_l1c", "predict_l1c", "report")
    )
    parser.add_argument("matchup_id", nargs="?")
    parser.add_argument("--dumps", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--mosaics-dir", type=Path, default=DEFAULT_MOSAICS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--variant", default="")
    args = parser.parse_args()
    if args.stage == "mosaics":
        build_mosaics(args.matchup_id, args.dumps, args.mosaics_dir)
    elif args.stage == "predict":
        predict_case(args.matchup_id, args.mosaics_dir, args.dumps, args.results, args.output)
    elif args.stage == "mosaics_l1c":
        build_mosaics_l1c(args.matchup_id, args.dumps, args.results, args.mosaics_dir)
    elif args.stage == "predict_l1c":
        predict_l1c(args.matchup_id, args.mosaics_dir, args.dumps, args.results, args.output)
    else:
        report(args.output, args.variant)


if __name__ == "__main__":
    main()

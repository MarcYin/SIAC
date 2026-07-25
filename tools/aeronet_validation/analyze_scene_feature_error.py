"""Locate the scene-day transfer loss: feature-oracle ladder for the kNN prior.

Hypothesis under test (user): library self-prediction is excellent, so if the
scene-day NIR/SWIR features were correct, the kNN should match it — the
transfer loss must therefore live in the INPUT features. Three feature
variants per site, identical kNN library and identical gated evaluation
pixels (dilated NDWI water and worst-5% library-RMSE pixels rejected,
acixThree-style):

  F0  anchors BOA at the ANCHOR AOD  — what inference actually uses today
  F1  anchors BOA at the TRUTH AOD   — oracle atmospheric correction
  F2  library temporal-median anchors — noiseless stable-surface oracle

Transfer RMSE ladder F0 -> F1 -> F2 vs the library self-prediction sigma
attributes the loss: (F0-F1) = anchor-AOD AC error, (F1-F2) = residual
scene-day anchor noise (reference RT noise + genuine change), (F2 vs self) =
model transfer floor.  Also reports the raw feature deltas per band.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DUMPS = ROOT / "calib_dumps_crossrt_dt_20260719"
PRIORS = ROOT / "analysis/acix3_style_prior_20260719"
MOSAICS = ROOT / "acix3_style_mosaics_20260719"
RESULTS = (
    ROOT
    / "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_dt_20260718"
)
OUT = PRIORS / "feature_error"
MIDS = ROOT / "calib_dump_examples_mids.txt"
ANCHORS = ("B8A", "B11", "B12")
K_NEIGHBOURS = 10


def boa_at_aod(dump, band, aod):
    grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    xap = float(np.interp(aod, grid, np.asarray(dump[f"xap_{band}"], dtype=np.float64)))
    xbp = float(np.interp(aod, grid, np.asarray(dump[f"xbp_{band}"], dtype=np.float64)))
    xcp = float(np.interp(aod, grid, np.asarray(dump[f"xcp_{band}"], dtype=np.float64)))
    toa = np.asarray(dump[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    return y / (1.0 + xcp * y)


def analyze(mid: str):
    from sklearn.neighbors import NearestNeighbors

    mosaic_path = MOSAICS / f"{mid}_l1c.npz"
    if not mosaic_path.exists():
        return None
    with np.load(DUMPS / f"{mid}.npz") as handle:
        dump = {k: handle[k] for k in handle.files}
    # The L1C npz stores RAW TOA + provenance planes; the library must be the
    # own-AC'd surface (predict_l1c corrects in memory and does not persist).
    from siac.adapters.rsrf import load_sensor_config_with_rsrf
    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    from tools.aeronet_validation.acix3_surface_prior import _lut_correct_realization

    with np.load(mosaic_path) as handle:
        stacks = handle["mosaics"].astype(np.float32)
        band_names = [str(b) for b in handle["band_names"]]
    pairs = ROOT / "analysis/l2a_l1c_physical_pairs_v3_deterrain_lowcloud152_20260718"
    elevation_km = 0.0
    pair_npz = pairs / f"{mid}.npz"
    if pair_npz.exists():
        with np.load(pair_npz) as handle:
            if "terrain_elevation_km" in handle.files:
                values = handle["terrain_elevation_km"]
                if np.isfinite(values).any():
                    elevation_km = float(np.nanmedian(values))
    rt_backend = ZarrLUTBackend(
        ROOT.parent / "libradtran_continental_average_lut_1nm.zarr.zip"
    )
    sensor = load_sensor_config_with_rsrf("MSI", "S2A")
    mosaics = np.stack(
        [
            _lut_correct_realization(rt_backend, sensor, stack, band_names, elevation_km)
            for stack in stacks
        ]
    ).astype(np.float64)
    mosaics = np.where(np.isfinite(mosaics), mosaics, np.nan)
    result = json.loads((RESULTS / f"{mid}.json").read_text())
    truth_aod = float(result["truth"])
    anchor_aod = float(dump["anchor_aot"])

    n_time, n_band, height, width = mosaics.shape
    clim = np.nanmean(mosaics[:, :4], axis=0)
    lib_features = np.concatenate(
        [
            mosaics[:, 4:7].transpose(0, 2, 3, 1).reshape(-1, 3),
            np.broadcast_to(
                clim.transpose(1, 2, 0), (n_time, height, width, 4)
            ).reshape(-1, 4),
        ],
        axis=1,
    )
    targets = mosaics.transpose(0, 2, 3, 1).reshape(-1, n_band)
    valid_lib = np.isfinite(lib_features).all(axis=1) & np.isfinite(targets).all(axis=1)
    x, y = lib_features[valid_lib], targets[valid_lib]
    if x.shape[0] < 100:
        return None
    index = NearestNeighbors(n_neighbors=K_NEIGHBOURS).fit(x)

    _, self_idx = index.kneighbors(x)
    self_pred = y[self_idx].mean(axis=1)
    self_sigma = float(np.sqrt(np.mean((self_pred[:, 1] - y[:, 1]) ** 2)))
    rmse_flat = np.full((lib_features.shape[0],), np.nan)
    rmse_flat[valid_lib] = (self_pred[:, 1] - y[:, 1]) ** 2
    rmse_b2 = np.sqrt(
        np.nanmean(rmse_flat.reshape(n_time, height, width), axis=0)
    )

    truth_b02 = boa_at_aod(dump, "B02", truth_aod)
    toa_g = np.asarray(dump["toa_B03"], dtype=np.float64)
    toa_nir = np.asarray(dump["toa_B8A"], dtype=np.float64)
    water = binary_dilation((toa_g - toa_nir) / (toa_g + toa_nir + 1e-9) >= 0.1, iterations=10)
    finite_rmse = rmse_b2[np.isfinite(rmse_b2)]
    unreliable = rmse_b2 > np.percentile(finite_rmse, 95)

    anchors_anchor = np.stack([boa_at_aod(dump, b, anchor_aod) for b in ANCHORS])
    anchors_truth = np.stack([boa_at_aod(dump, b, truth_aod) for b in ANCHORS])
    lib_median_anchor = np.nanmedian(mosaics[:, 4:7], axis=0)

    def predict(anchor_planes):
        features = np.concatenate(
            [anchor_planes.transpose(1, 2, 0), clim.transpose(1, 2, 0)], axis=2
        )
        mask = np.isfinite(features).all(axis=2)
        pred = np.full((height, width), np.nan)
        if mask.any():
            _, idx = index.kneighbors(features[mask])
            pred[mask] = y[idx].mean(axis=1)[:, 1]
        return pred

    preds = {
        "F0_anchor_aod": predict(anchors_anchor),
        "F1_truth_aod": predict(anchors_truth),
        "F2_library_median": predict(lib_median_anchor),
    }
    gate = (
        np.isfinite(truth_b02)
        & (truth_b02 > -0.05)
        & (truth_b02 < 1.0)
        & ~water
        & ~unreliable
    )
    row = dict(
        site=mid.split("__")[0],
        truth_aod=round(truth_aod, 3),
        anchor_aod=round(anchor_aod, 3),
        library_self_rmse=round(self_sigma, 5),
        gated_pixels=int(gate.sum()),
    )
    for name, plane in preds.items():
        mask = gate & np.isfinite(plane)
        if mask.sum() < 100:
            continue
        diff = plane[mask] - truth_b02[mask]
        row[f"{name}_bias"] = round(float(np.median(diff)), 5)
        row[f"{name}_rmse"] = round(float(np.sqrt(np.mean(diff**2))), 5)
    for i, band in enumerate(ANCHORS):
        mask = gate & np.isfinite(anchors_anchor[i]) & np.isfinite(anchors_truth[i])
        lib_mask = mask & np.isfinite(lib_median_anchor[i])
        row[f"{band}_ac_delta"] = round(
            float(np.median(np.abs((anchors_anchor[i] - anchors_truth[i])[mask]))), 5
        )
        row[f"{band}_scene_vs_lib"] = round(
            float(np.median((anchors_truth[i] - lib_median_anchor[i])[lib_mask])), 5
        )
    return row


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for mid in [m for m in MIDS.read_text().split() if m]:
        row = analyze(mid)
        if row is None:
            print(f"skip {mid}")
            continue
        rows.append(row)
        print(json.dumps(row))
    import csv

    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with (OUT / "feature_error_ladder.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    med = lambda key: float(np.median([r[key] for r in rows if key in r]))  # noqa: E731
    print("\n=== medians over sites ===")
    print(f"library self RMSE          {med('library_self_rmse'):.4f}")
    for name in ("F0_anchor_aod", "F1_truth_aod", "F2_library_median"):
        print(
            f"{name:26s} bias {med(name + '_bias'):+.4f}  RMSE {med(name + '_rmse'):.4f}"
        )
    for band in ANCHORS:
        print(
            f"{band}: |AC delta anchor-vs-truth| {med(band + '_ac_delta'):.4f}   "
            f"scene@truth - libmedian {med(band + '_scene_vs_lib'):+.4f}"
        )


if __name__ == "__main__":
    main()

"""Outlier forensics for the three-way surface-prior comparison.

For each example case: identify scatter outliers (robust MAD rule on the
B02 residual against the truth-AOD surface), attribute each to a cause using
acixThree's own diagnostics — NDWI water (their apply_water_mask rule,
dilated), scene-day bright/dark contamination (scene surface vs library
median), unreliable-library pixels (their prior_surface_reflectance_mask
95th-percentile RMSE gate) — and measure the accuracy ladder:

  library self-prediction (in-domain model skill, acixThree's remembered
  performance plot) -> scene-day transfer, all pixels -> scene-day transfer
  after the acixThree-style gates.

Outputs per-site demonstration figures, a summary table, and a self-contained
HTML report.
"""

from __future__ import annotations

import base64
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
OUT = PRIORS / "outliers"
MIDS = ROOT / "calib_dump_examples_mids.txt"
DEEP_DIVE = (
    "AAQ12_Kx_Xiaoliuqiu",
    "Gozo",
    "Karachi",
    "Rhyl_MO",
    "CUIABA-MIRANDA",
    "NEON_NOGP",
    "Bangkok",
    "MetObs_Lindenberg",
)
GEE_BANDS = ("B1", "B2", "B3", "B4", "B8A", "B11", "B12")
CLASS_COLORS = {
    "inlier": "#d9d9d9",
    "water": "#2a78d6",
    "scene_bright": "#e8a33d",
    "scene_dark": "#7d54c9",
    "unreliable": "#898781",
    "other": "#e34948",
}
K_NEIGHBOURS = 10


def boa_at_aod(dump, band, aod):
    grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    xap = float(np.interp(aod, grid, np.asarray(dump[f"xap_{band}"], dtype=np.float64)))
    xbp = float(np.interp(aod, grid, np.asarray(dump[f"xbp_{band}"], dtype=np.float64)))
    xcp = float(np.interp(aod, grid, np.asarray(dump[f"xcp_{band}"], dtype=np.float64)))
    toa = np.asarray(dump[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    return y / (1.0 + xcp * y)


def library_self_prediction(mosaics: np.ndarray) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    """acixThree train_model diagnostic: self-kNN prediction over the library."""
    from sklearn.neighbors import NearestNeighbors

    n_time, n_band, height, width = mosaics.shape
    ti = list(range(4))
    fi = [4, 5, 6]
    clim = np.nanmean(mosaics[:, ti], axis=0)
    features = np.concatenate(
        [
            mosaics[:, fi].transpose(0, 2, 3, 1).reshape(-1, 3),
            np.broadcast_to(
                clim.transpose(1, 2, 0), (n_time, height, width, 4)
            ).reshape(-1, 4),
        ],
        axis=1,
    )
    targets = mosaics.transpose(0, 2, 3, 1).reshape(-1, n_band)
    valid = np.isfinite(features).all(axis=1) & np.isfinite(targets).all(axis=1)
    x, y = features[valid], targets[valid]
    index = NearestNeighbors(n_neighbors=K_NEIGHBOURS).fit(x)
    _, idx = index.kneighbors(x)
    pred = y[idx].mean(axis=1)
    metrics = {}
    for i, band in enumerate(GEE_BANDS[:4]):
        diff = pred[:, i] - y[:, i]
        metrics[band] = {
            "rmse": float(np.sqrt(np.mean(diff**2))),
            "mad_sigma": float(1.4826 * np.median(np.abs(diff - np.median(diff)))),
        }
    return metrics, y[:, 1], pred[:, 1]  # B2 observed / predicted for the hist


def analyze_site(mid: str):
    site = mid.split("__")[0]
    dump_path = DUMPS / f"{mid}.npz"
    prior_path = PRIORS / f"{mid}_prior_l1c.npz"
    mosaic_path = MOSAICS / f"{mid}_l1c.npz"
    variant = "L1C own-AC"
    if not prior_path.exists():
        prior_path = PRIORS / f"{mid}_prior.npz"
        mosaic_path = MOSAICS / f"{mid}.npz"
        variant = "Sen2Cor"
    if not (dump_path.exists() and prior_path.exists() and mosaic_path.exists()):
        return None
    with np.load(dump_path) as handle:
        dump = {k: handle[k] for k in handle.files}
    with np.load(prior_path) as handle:
        names = [str(b) for b in handle["band_names"]]
        pred = np.where(np.isfinite(handle["prediction"]), handle["prediction"], np.nan)
        rmse_planes = np.where(np.isfinite(handle["rmse"]), handle["rmse"], np.nan)
    with np.load(mosaic_path) as handle:
        mosaics = handle["mosaics"].astype(np.float64)[:, : len(GEE_BANDS)]
    mosaics = np.where(np.isfinite(mosaics), mosaics, np.nan)
    result = json.loads((RESULTS / f"{mid}.json").read_text())
    truth_aod = float(result["truth"])

    b2 = names.index("B2")
    truth = boa_at_aod(dump, "B02", truth_aod)
    prior = pred[b2]
    current = (
        np.asarray(dump["pred_B02"], dtype=np.float64) if "pred_B02" in dump else None
    )
    valid = np.isfinite(prior) & np.isfinite(truth) & (truth > -0.05) & (truth < 1.0)
    residual = np.where(valid, prior - truth, np.nan)

    # --- classification, using acixThree's own rules -----------------------
    toa_g = np.asarray(dump["toa_B03"], dtype=np.float64)
    toa_nir = np.asarray(dump["toa_B8A"], dtype=np.float64)
    ndwi = (toa_g - toa_nir) / (toa_g + toa_nir + 1e-9)
    water = binary_dilation(ndwi >= 0.1, iterations=10)

    library_median = np.nanmedian(mosaics[:, 1], axis=0)  # B2
    scene_delta = truth - library_median
    scene_bright = scene_delta > 0.04
    scene_dark = scene_delta < -0.04

    rmse_b2 = rmse_planes[b2]
    finite_rmse = rmse_b2[np.isfinite(rmse_b2)]
    unreliable = (
        rmse_b2 > np.percentile(finite_rmse, 95) if finite_rmse.size else np.zeros_like(water)
    )

    dev = residual - np.nanmedian(residual)
    mad = 1.4826 * np.nanmedian(np.abs(dev))
    is_outlier = valid & (np.abs(dev) > max(0.015, 4 * mad))

    classes = np.full(truth.shape, "", dtype=object)
    classes[valid] = "inlier"
    for name, mask in (
        ("other", is_outlier),
        ("unreliable", is_outlier & unreliable),
        ("scene_dark", is_outlier & scene_dark),
        ("scene_bright", is_outlier & scene_bright),
        ("water", is_outlier & water),
    ):
        classes[valid & mask] = name

    # acixThree-style gate: water out, worst-5% library RMSE out, scene-day
    # bright/dark contamination out.
    gated = valid & ~water & ~unreliable & ~scene_bright & ~scene_dark

    def stats(mask):
        if mask.sum() < 30:
            return np.nan, np.nan
        d = (prior - truth)[mask]
        return float(np.median(d)), float(np.sqrt(np.mean(d**2)))

    bias_all, rmse_all = stats(valid)
    bias_gated, rmse_gated = stats(gated)
    cur_rmse = np.nan
    if current is not None:
        mask = valid & np.isfinite(current)
        if mask.sum() >= 30:
            cur_rmse = float(np.sqrt(np.mean(((current - truth)[mask]) ** 2)))

    self_metrics, self_obs, self_pred = library_self_prediction(mosaics)
    counts = {
        name: int((classes == name).sum())
        for name in ("water", "scene_bright", "scene_dark", "unreliable", "other")
    }
    row = dict(
        site=site,
        variant=variant,
        truth_aod=round(truth_aod, 3),
        n_valid=int(valid.sum()),
        outlier_fraction=round(float(is_outlier.sum() / max(valid.sum(), 1)), 4),
        **{f"outlier_{k}": v for k, v in counts.items()},
        library_self_rmse_B2=round(self_metrics["B2"]["rmse"], 5),
        library_self_mad_sigma_B2=round(self_metrics["B2"]["mad_sigma"], 5),
        transfer_bias_all=round(bias_all, 5),
        transfer_rmse_all=round(rmse_all, 5),
        transfer_bias_gated=round(bias_gated, 5),
        transfer_rmse_gated=round(rmse_gated, 5),
        current_rmse_all=round(cur_rmse, 5) if np.isfinite(cur_rmse) else None,
        gated_fraction=round(float(gated.sum() / max(valid.sum(), 1)), 4),
    )
    arrays = dict(
        truth=truth, prior=prior, residual=residual, classes=classes, valid=valid,
        gated=gated, self_obs=self_obs, self_pred=self_pred, variant=variant,
        truth_aod=truth_aod, site=site,
    )
    return row, arrays


def render_site(arrays, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcdefaults()
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    site = arrays["site"]
    truth, prior = arrays["truth"], arrays["prior"]
    residual, classes = arrays["residual"], arrays["classes"]
    valid, gated = arrays["valid"], arrays["gated"]

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), constrained_layout=True)
    vmin = np.nanpercentile(truth[valid], 1)
    vmax = np.nanpercentile(truth[valid], 99)
    im0 = axes[0, 0].imshow(np.where(valid, truth, np.nan), vmin=vmin, vmax=vmax, cmap="viridis")
    axes[0, 0].set_title("surface @ truth AOD (B02)", fontsize=9)
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.8)
    im1 = axes[0, 1].imshow(np.where(valid, prior, np.nan), vmin=vmin, vmax=vmax, cmap="viridis")
    axes[0, 1].set_title("acix3 predicted prior (B02)", fontsize=9)
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.8)
    im2 = axes[0, 2].imshow(residual, vmin=-0.03, vmax=0.03, cmap="RdBu_r")
    axes[0, 2].set_title("residual (prior − truth)", fontsize=9)
    fig.colorbar(im2, ax=axes[0, 2], shrink=0.8)

    order = ["inlier", "water", "scene_bright", "scene_dark", "unreliable", "other"]
    class_index = np.full(truth.shape, np.nan)
    for i, name in enumerate(order):
        class_index[classes == name] = i
    cmap = ListedColormap([CLASS_COLORS[name] for name in order])
    axes[0, 3].imshow(class_index, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5,
                      interpolation="nearest")
    axes[0, 3].set_title("outlier classes", fontsize=9)
    handles = [
        plt.Line2D([], [], marker="s", linestyle="", color=CLASS_COLORS[name], label=name)
        for name in order
    ]
    axes[0, 3].legend(handles=handles, fontsize=6, loc="lower right", framealpha=0.9)

    lim = (vmin - 0.01, vmax + 0.01)
    ax = axes[1, 0]
    ax.plot(lim, lim, color="#898781", lw=0.8)
    for name in order:
        mask = (classes == name) & valid
        if not mask.any():
            continue
        step = max(1, int(mask.sum()) // 3000)
        ax.scatter(truth[mask][::step], prior[mask][::step], s=3,
                   color=CLASS_COLORS[name], alpha=0.4 if name == "inlier" else 0.8,
                   linewidths=0, label=name)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title("scatter by class", fontsize=9)
    ax.set_xlabel("surface @ truth AOD", fontsize=8)
    ax.set_ylabel("predicted prior", fontsize=8)

    ax = axes[1, 1]
    ax.plot(lim, lim, color="#898781", lw=0.8)
    step = max(1, int(gated.sum()) // 4000)
    ax.scatter(truth[gated][::step], prior[gated][::step], s=3, color="#2a78d6",
               alpha=0.3, linewidths=0)
    d = (prior - truth)[gated]
    ax.text(0.03, 0.97, f"bias {np.median(d):+.4f}\nRMSE {np.sqrt(np.mean(d**2)):.4f}\n"
            f"kept {100 * gated.sum() / max(valid.sum(), 1):.0f}%",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title("after acixThree-style gates", fontsize=9)
    ax.set_xlabel("surface @ truth AOD", fontsize=8)

    ax = axes[1, 2]
    obs, predicted = arrays["self_obs"], arrays["self_pred"]
    hi = np.nanpercentile(obs, 99.5)
    ax.hist2d(obs, predicted, bins=150, range=[[0, hi], [0, hi]], cmin=3, cmap="Blues")
    ax.plot([0, hi], [0, hi], color="#e34948", lw=0.8, ls="--")
    d = predicted - obs
    ax.text(0.03, 0.97, f"library self-prediction\nRMSE {np.sqrt(np.mean(d**2)):.4f}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_title("library self-prediction (B2)", fontsize=9)
    ax.set_xlabel("library observed", fontsize=8)
    ax.set_ylabel("kNN predicted", fontsize=8)

    ax = axes[1, 3]
    bins = np.linspace(-0.06, 0.06, 80)
    ax.hist(residual[valid], bins=bins, color="#9ec5f4", label="all pixels")
    ax.hist((prior - truth)[gated], bins=bins, color="#2a78d6", label="gated",
            histtype="step", lw=1.5)
    ax.legend(fontsize=7)
    ax.set_title("residual histogram", fontsize=9)
    ax.set_xlabel("prior − truth", fontsize=8)

    for row_axes in axes[:1]:
        for ax in row_axes:
            ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(
        f"{site} — {arrays['variant']} library — truth AOD {arrays['truth_aod']:.2f}",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    figures = []
    for mid in [m for m in MIDS.read_text().split() if m]:
        out = analyze_site(mid)
        if out is None:
            print(f"skip {mid}")
            continue
        row, arrays = out
        rows.append(row)
        if row["site"] in DEEP_DIVE:
            fig_path = OUT / f"{row['site']}_outliers.png"
            render_site(arrays, fig_path)
            figures.append((row["site"], fig_path.name))
        print(json.dumps(row))

    import csv

    with (OUT / "outlier_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    def img(name):
        data = base64.b64encode((OUT / name).read_bytes()).decode()
        return f'<img src="data:image/png;base64,{data}" alt="{name}">'

    header_cells = (
        "site,variant,truth AOD,outliers %,water,bright,dark,unreliable,other,"
        "lib self-RMSE,transfer RMSE,gated RMSE,kept %"
    ).split(",")
    table_rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{v}</td>"
            for v in (
                r["site"], r["variant"], r["truth_aod"],
                f"{100 * r['outlier_fraction']:.1f}", r["outlier_water"],
                r["outlier_scene_bright"], r["outlier_scene_dark"],
                r["outlier_unreliable"], r["outlier_other"],
                f"{r['library_self_rmse_B2']:.4f}", f"{r['transfer_rmse_all']:.4f}",
                f"{r['transfer_rmse_gated']:.4f}", f"{100 * r['gated_fraction']:.0f}",
            )
        )
        + "</tr>"
        for r in sorted(rows, key=lambda r: -r["outlier_fraction"])
    )
    sections = "".join(
        f"<h2>{site}</h2>{img(name)}" for site, name in figures
    )
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Surface-prior outlier forensics</title>
<style>body {{ font-family: system-ui, sans-serif; margin: 1.5rem auto; max-width: 1200px; }}
img {{ max-width: 100%; border: 1px solid #e5e5e5; }}
table {{ border-collapse: collapse; }} td, th {{ border: 1px solid #ccc;
padding: 0.3rem 0.55rem; font-size: 0.82rem; }}</style></head><body>
<h1>Surface-prior outlier forensics — acix3 library vs truth-AOD surface (B02)</h1>
<p>Outlier = residual deviating &gt; max(0.015, 4·MAD) from the site median.
Classes use acixThree's own rules: NDWI ≥ 0.1 water (dilated 10), scene-day
bright/dark = scene surface vs library median beyond ±0.04, unreliable =
worst-5% library prediction RMSE (their prior_surface_reflectance_mask gate).
"Library self-prediction" is acixThree's train_model diagnostic: kNN
self-prediction over the mosaic stack — the in-domain model skill.</p>
<table><tr>{"".join(f"<th>{c}</th>" for c in header_cells)}</tr>{table_rows}</table>
{sections}
</body></html>"""
    (OUT / "outlier_report.html").write_text(html)
    print(f"report: {OUT}/outlier_report.html")


if __name__ == "__main__":
    main()

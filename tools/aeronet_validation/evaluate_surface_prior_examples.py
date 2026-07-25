"""Truth-referenced surface-prior accuracy on representative example cases.

Reads the calibration dumps produced by ``cross_rt_calib_dump.sbatch`` (per-pixel
predicted visible prior, scene TOA, and scene-mean RT coefficient curves over an
80-node AOD grid), reconstructs the scene-day surface reflectance at the AERONET
truth AOD, and scores the prior prediction against it per band and per case.

The truth surface uses the standard coefficient inversion
``y = xap * TOA - xbp; BOA = y / (1 + xcp * y)`` with the coefficient curves
interpolated at the truth AOD. RT coefficients are scene-mean (median WVP/O3,
scene-mean geometry), so per-pixel WVP variation is not represented — adequate
for scene-level prior assessment.

Outputs (under ``analysis/surface_prior_examples_<date>/``):
- ``surface_prior_examples.csv`` — one row per case x band with prior median,
  truth-surface median, bias, MAE, RMSE, correlation, and retrieval context.
- ``scatter_<band>.png`` — 20-panel per-pixel scatter grids, prior vs truth
  surface.
- ``summary.md`` — the example list with the headline numbers.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_DUMPS = ROOT / "calib_dumps_crossrt_dt_20260719"
DEFAULT_RESULTS = (
    ROOT / "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_dt_20260718"
)
DEFAULT_OUTPUT = ROOT / "analysis/surface_prior_examples_20260719"
BANDS = ("B02", "B03", "B04")


def boa_at_aod(dump: dict[str, np.ndarray], band: str, aod: float) -> np.ndarray:
    """Per-pixel surface reflectance from stored TOA and coefficient curves."""
    grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    xap = float(np.interp(aod, grid, np.asarray(dump[f"xap_{band}"], dtype=np.float64)))
    xbp = float(np.interp(aod, grid, np.asarray(dump[f"xbp_{band}"], dtype=np.float64)))
    xcp = float(np.interp(aod, grid, np.asarray(dump[f"xcp_{band}"], dtype=np.float64)))
    toa = np.asarray(dump[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    return y / (1.0 + xcp * y)


def case_metrics(dump: dict[str, np.ndarray], band: str, truth_aod: float) -> dict[str, float]:
    pred = np.asarray(dump[f"pred_{band}"], dtype=np.float64)
    target = boa_at_aod(dump, band, truth_aod)
    valid = (
        np.isfinite(pred)
        & np.isfinite(target)
        & (pred > 0.0)
        & (target > -0.05)
        & (target < 1.0)
    )
    if valid.sum() < 50:
        return {}
    p, t = pred[valid], target[valid]
    diff = p - t
    corr = float(np.corrcoef(p, t)[0, 1]) if p.size > 2 else float("nan")
    return {
        "n_pixels": int(valid.sum()),
        "prior_median": float(np.median(p)),
        "truth_surface_median": float(np.median(t)),
        "bias_median": float(np.median(diff)),
        "bias_mean": float(np.mean(diff)),
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "corr": corr,
        "pred": p,
        "target": t,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dumps", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcdefaults()
    import matplotlib.pyplot as plt

    args.output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    scatter_data: dict[str, list[tuple[str, np.ndarray, np.ndarray, float]]] = {
        band: [] for band in BANDS
    }
    for dump_path in sorted(args.dumps.glob("*.npz")):
        mid = dump_path.stem
        result_path = args.results / f"{mid}.json"
        if not result_path.exists():
            continue
        result = json.loads(result_path.read_text())
        truth = float(result["truth"])
        with np.load(dump_path) as handle:
            dump = {key: handle[key] for key in handle.files}
        for band in BANDS:
            metrics = case_metrics(dump, band, truth)
            if not metrics:
                continue
            pred = metrics.pop("pred")
            target = metrics.pop("target")
            if band == "B02":
                scatter_data[band].append((mid, pred, target, truth))
            else:
                scatter_data[band].append((mid, pred, target, truth))
            rows.append(
                {
                    "matchup_id": mid,
                    "site": mid.split("__")[0],
                    "band": band,
                    "truth_aod": round(truth, 4),
                    "retrieved_aod": round(float(result["retrieved"]), 4),
                    "maiac_aod": round(
                        float((result.get("atmo_prior") or {}).get("aot_median") or np.nan), 4
                    ),
                    "within_ee": bool(result["within_ee"]),
                    "prior_unc": round(float((result.get("prior_unc") or {}).get(band, np.nan)), 5),
                    **{k: round(v, 5) if isinstance(v, float) else v for k, v in metrics.items()},
                }
            )

    if not rows:
        raise SystemExit(f"no evaluable dumps under {args.dumps}")

    fieldnames = list(rows[0].keys())
    with (args.output / "surface_prior_examples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    for band in BANDS:
        cases = scatter_data[band]
        if not cases:
            continue
        ncols = 5
        nrows = int(np.ceil(len(cases) / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), constrained_layout=True
        )
        for ax in np.ravel(axes):
            ax.set_visible(False)
        for ax, (mid, pred, target, truth) in zip(np.ravel(axes), cases):
            ax.set_visible(True)
            lo = float(min(np.percentile(target, 1), np.percentile(pred, 1)))
            hi = float(max(np.percentile(target, 99), np.percentile(pred, 99)))
            pad = 0.05 * (hi - lo + 1e-6)
            lim = (lo - pad, hi + pad)
            ax.plot(lim, lim, color="#898781", lw=0.8, zorder=1)
            step = max(1, pred.size // 4000)
            ax.scatter(
                target[::step], pred[::step], s=2, alpha=0.25, color="#2a78d6", linewidths=0
            )
            bias = float(np.median(pred - target))
            rmse = float(np.sqrt(np.mean((pred - target) ** 2)))
            site = mid.split("__")[0][:18]
            ax.set_title(f"{site}\nAOD {truth:.2f}  bias {bias:+.4f}  RMSE {rmse:.4f}", fontsize=7)
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.tick_params(labelsize=6)
        fig.suptitle(
            f"{band}: predicted prior (y) vs surface at AERONET-truth AOD (x)", fontsize=11
        )
        fig.savefig(args.output / f"scatter_{band}.png", dpi=150)
        plt.close(fig)

    by_band: dict[str, list[dict[str, object]]] = {band: [] for band in BANDS}
    for row in rows:
        by_band[str(row["band"])].append(row)
    lines = [
        "# Surface-prior accuracy on representative cases",
        "",
        "Prior = et20 seasonal predictor over the de-terrained CAMS-O3 harmonized",
        "histories (the 71.1% baseline configuration). Reference = per-pixel surface",
        "reflectance from the scene TOA inverted at the AERONET-truth AOD.",
        "",
        "| band | median bias | median MAE | median RMSE | median corr | cases |",
        "|---|---|---|---|---|---|",
    ]
    for band in BANDS:
        band_rows = by_band[band]
        if not band_rows:
            continue
        med = lambda key: float(np.median([float(r[key]) for r in band_rows]))  # noqa: E731
        lines.append(
            f"| {band} | {med('bias_median'):+.4f} | {med('mae'):.4f} "
            f"| {med('rmse'):.4f} | {med('corr'):.3f} | {len(band_rows)} |"
        )
    lines += ["", "Per-case detail: `surface_prior_examples.csv`; scatters: `scatter_<band>.png`."]
    (args.output / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(rows)} rows for {len({r['matchup_id'] for r in rows})} cases")
    for line in lines[6:10]:
        print(line)


if __name__ == "__main__":
    main()

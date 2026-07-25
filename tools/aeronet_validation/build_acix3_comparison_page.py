"""Three-way surface-prior scatter comparison webpage.

For every example case, renders per-pixel scatters of predicted prior (y)
against the surface at AERONET-truth AOD (x) for the three prior variants —
current et20 (harmonized de-terrained L2A histories), acixThree-style on the
Sen2Cor library, and the faithful acixThree L1C + own-AC library — one figure
per site and band, assembled into ``report/index.html`` (relative image paths,
served from the GWS public path) and ``report/self_contained.html`` (base64
images, portable single file).
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DUMPS = ROOT / "calib_dumps_crossrt_dt_20260719"
PRIORS = ROOT / "analysis/acix3_style_prior_20260719"
RESULTS = (
    ROOT
    / "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_dt_20260718"
)
import os

REPORT = PRIORS / ("report_gated" if os.environ.get("ACIX3_PAGE_GATED") == "1" else "report")
GATED = os.environ.get("ACIX3_PAGE_GATED") == "1"
MIDS = ROOT / "calib_dump_examples_mids.txt"
BANDS = ("B02", "B03", "B04")
GEE_BY_DUMP = {"B02": "B2", "B03": "B3", "B04": "B4"}
BLUE = "#2a78d6"
MUTED = "#898781"
VARIANTS = (
    ("current", "current et20 (harmonized L2A)"),
    ("s2c", "acix3-style, Sen2Cor library"),
    ("l1c", "acix3 faithful, L1C + own-AC"),
)


def boa_at_aod(dump: dict[str, np.ndarray], band: str, aod: float) -> np.ndarray:
    grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    xap = float(np.interp(aod, grid, np.asarray(dump[f"xap_{band}"], dtype=np.float64)))
    xbp = float(np.interp(aod, grid, np.asarray(dump[f"xbp_{band}"], dtype=np.float64)))
    xcp = float(np.interp(aod, grid, np.asarray(dump[f"xcp_{band}"], dtype=np.float64)))
    toa = np.asarray(dump[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    return y / (1.0 + xcp * y)


def acix_plane(path: Path, band: str) -> np.ndarray | None:
    if not path.exists():
        return None
    with np.load(path) as handle:
        names = [str(b) for b in handle["band_names"]]
        gee = GEE_BY_DUMP[band]
        if gee not in names:
            return None
        plane = handle["prediction"][names.index(gee)].astype(np.float64)
    return np.where(np.isfinite(plane), plane, np.nan)


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcdefaults()
    import matplotlib.pyplot as plt

    REPORT.mkdir(parents=True, exist_ok=True)
    mids = [m for m in MIDS.read_text().split() if m]
    cases = []
    for mid in mids:
        dump_path = DUMPS / f"{mid}.npz"
        result_path = RESULTS / f"{mid}.json"
        if not dump_path.exists() or not result_path.exists():
            continue
        result = json.loads(result_path.read_text())
        with np.load(dump_path) as handle:
            dump = {k: handle[k] for k in handle.files}
        truth_aod = float(result["truth"])
        site = mid.split("__")[0]
        figures: dict[str, str] = {}
        current_bias = np.nan
        for band in BANDS:
            target = boa_at_aod(dump, band, truth_aod)
            planes = {
                "current": (
                    np.asarray(dump[f"pred_{band}"], dtype=np.float64)
                    if f"pred_{band}" in dump
                    else None
                ),
                "s2c": acix_plane(PRIORS / f"{mid}_prior.npz", band),
                "l1c": acix_plane(PRIORS / f"{mid}_prior_l1c.npz", band),
            }
            base_mask = np.isfinite(target) & (target > -0.05) & (target < 1.0)
            if GATED:
                from scipy.ndimage import binary_dilation

                toa_g = np.asarray(dump["toa_B03"], dtype=np.float64)
                toa_n = np.asarray(dump["toa_B8A"], dtype=np.float64)
                water = binary_dilation(
                    (toa_g - toa_n) / (toa_g + toa_n + 1e-9) >= 0.1, iterations=10
                )
                base_mask &= ~water
                for suffix in ("_prior_l1c.npz", "_prior.npz"):
                    rmse_path = PRIORS / f"{mid}{suffix}"
                    if rmse_path.exists():
                        with np.load(rmse_path) as handle:
                            names = [str(b) for b in handle["band_names"]]
                            plane = handle["rmse"][names.index(GEE_BY_DUMP[band])]
                        plane = np.where(np.isfinite(plane), plane, np.nan)
                        finite = plane[np.isfinite(plane)]
                        if finite.size:
                            base_mask &= ~(plane > np.percentile(finite, 95))
                        break
            fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5), constrained_layout=True)
            finite_all = [
                plane[base_mask & np.isfinite(plane)]
                for plane in planes.values()
                if plane is not None
            ]
            t_all = target[base_mask]
            lo = min(
                [np.percentile(t_all, 1)] + [np.percentile(v, 1) for v in finite_all if v.size]
            )
            hi = max(
                [np.percentile(t_all, 99)] + [np.percentile(v, 99) for v in finite_all if v.size]
            )
            pad = 0.05 * (hi - lo + 1e-6)
            lim = (lo - pad, hi + pad)
            for ax, (key, label) in zip(axes, VARIANTS):
                plane = planes[key]
                ax.plot(lim, lim, color=MUTED, lw=0.8, zorder=1)
                if plane is None:
                    ax.text(
                        0.5, 0.5, "not available", transform=ax.transAxes,
                        ha="center", va="center", color=MUTED, fontsize=9,
                    )
                else:
                    mask = base_mask & np.isfinite(plane)
                    p, t = plane[mask], target[mask]
                    step = max(1, p.size // 4000)
                    ax.scatter(
                        t[::step], p[::step], s=2.5, alpha=0.25, color=BLUE, linewidths=0
                    )
                    bias = float(np.median(p - t))
                    rmse = float(np.sqrt(np.mean((p - t) ** 2)))
                    if key == "current" and band == "B02":
                        current_bias = bias
                    ax.text(
                        0.03, 0.97, f"bias {bias:+.4f}\nRMSE {rmse:.4f}",
                        transform=ax.transAxes, ha="left", va="top", fontsize=8,
                    )
                ax.set_title(label, fontsize=9)
                ax.set_xlim(lim)
                ax.set_ylim(lim)
                ax.tick_params(labelsize=7)
                ax.set_xlabel("surface @ truth AOD", fontsize=8)
            axes[0].set_ylabel("predicted prior", fontsize=8)
            fig.suptitle(f"{site} — {band} — truth AOD {truth_aod:.2f}", fontsize=11)
            name = f"{site}_{band}.png".replace("/", "_")
            fig.savefig(REPORT / name, dpi=110)
            plt.close(fig)
            figures[band] = name
        cases.append(
            dict(site=site, truth_aod=truth_aod, figures=figures, current_bias=current_bias)
        )
        print(f"rendered {site}")

    cases.sort(key=lambda c: c["current_bias"] if np.isfinite(c["current_bias"]) else 0.0)
    summary_rows = [
        ("current et20 (harmonized L2A)", "-0.0132", "0.0165", "-0.0192", "-0.0171"),
        ("acix3-style, Sen2Cor library", "-0.0167", "0.0236", "-0.0198", "-0.0187"),
        ("acix3 faithful, L1C + own-AC", "-0.0097", "0.0163", "-0.0128", "-0.0119"),
    ]

    def page(embed: bool) -> str:
        def img(name: str) -> str:
            if not embed:
                return f'<img src="{name}" loading="lazy" alt="{name}">'
            data = base64.b64encode((REPORT / name).read_bytes()).decode()
            return f'<img src="data:image/png;base64,{data}" alt="{name}">'

        sections = []
        for band in BANDS:
            blocks = "\n".join(
                f'<div class="case"><h3>{c["site"]} <span class="aod">truth AOD '
                f'{c["truth_aod"]:.2f}</span></h3>{img(c["figures"][band])}</div>'
                for c in cases
                if band in c["figures"]
            )
            sections.append(
                f'<section class="band" id="{band}">'
                f"<h2>{band}: predicted prior vs surface at truth AOD</h2>{blocks}</section>"
            )
        table = "".join(
            f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td></tr>"
            for r in summary_rows
        )
        return f"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Surface prior three-way comparison — lowcloud152 examples</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 1.5rem auto; max-width: 1100px;
       color: #1a1a1a; }}
h1 {{ font-size: 1.4rem; }} h2 {{ font-size: 1.15rem; margin-top: 2rem; }}
h3 {{ font-size: 1rem; margin: 1.2rem 0 0.3rem; }} .aod {{ color: {MUTED}; font-weight: normal; }}
img {{ max-width: 100%; height: auto; border: 1px solid #e5e5e5; }}
table {{ border-collapse: collapse; margin: 0.8rem 0; }}
td, th {{ border: 1px solid #ccc; padding: 0.35rem 0.7rem; font-size: 0.9rem; }}
nav a {{ margin-right: 1rem; }}
.note {{ color: {MUTED}; font-size: 0.9rem; }}
</style></head><body>
<h1>Surface prior three-way comparison — truth-referenced, lowcloud152 examples{" (gated)" if GATED else ""}</h1>
{"<p class='note'>Problematic pixels rejected acixThree-style before plotting: dilated NDWI water mask + worst-5% library prediction-RMSE pixels.</p>" if GATED else ""}
<p class="note">Per pixel, x = surface reflectance from the scene TOA inverted at the
AERONET-truth AOD (continental libRadtran LUT, scene-mean RT); y = predicted prior.
Sites ordered by current-prior B02 bias (worst first). Generated 2026-07-19.</p>
<table><tr><th>prior</th><th>B02 med bias</th><th>B02 med MAE</th><th>B03 med bias</th>
<th>B04 med bias</th></tr>{table}</table>
<nav>{"".join(f'<a href="#{b}">{b}</a>' for b in BANDS)}</nav>
{"".join(sections)}
</body></html>"""

    (REPORT / "index.html").write_text(page(embed=False))
    (REPORT / "self_contained.html").write_text(page(embed=True))
    print(f"report: {REPORT}/index.html ({len(cases)} sites)")


if __name__ == "__main__":
    main()

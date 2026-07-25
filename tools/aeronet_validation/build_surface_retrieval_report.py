#!/usr/bin/env python3
"""Build the campaign-250 surface-retrieval diagnosis and experiment report.

The report is intentionally generated from saved experiment artifacts. It does
not query AERONET, Earth Engine, or any aerosol product at build time, and it
does not alter retrieval outputs. Re-run it as campaign backfills complete to
refresh the coverage-sensitive acixThree/sdspec figures.
"""

from __future__ import annotations

import html
import json
import math
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable

REPO = Path(__file__).resolve().parents[2]
ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
OUT = REPO / "notes" / "reports" / "surface-driven-c250-diagnosis-20260709"
SHELL = Path(
    "/home/users/marcyin/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.6-d37358633e00/assets/html-report-shell.html"
)
EMBED = Path(
    "/home/users/marcyin/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.6-d37358633e00/skills/build-report/scripts/"
    "embed_html_report_runtime.py"
)

R2_DIR = ROOT / "phaseD_results_campaign250_R2_full"
R2_LOCAL_DIR = ROOT / "phaseD_results_campaign250_R2_full_localdiag_20260705"
R2_COST_DIR = ROOT / "phaseD_cost_cubes_campaign250_R2_full_localdiag_20260705"
PQ_DIR = ROOT / "prior_quality"
SDSPEC_DIR = ROOT / "seasonal_val"
SDSPEC_PATTERN = "*__seas_maiac_sdspec_c250_20260708.json"
FINDINGS = REPO / "notes" / "l2a-surface-driven-campaign-findings.md"
EXPERIMENT_SUMMARY = ROOT / "proposed_surface_experiments_final_20260710.json"
PRIOR_QUALITY_FINAL = ROOT / "l1c_prior_quality_final_20260710.json"
BAND_CONSENSUS_FINAL = ROOT / "band_consensus_knn_final_20260710.json"
AGGREGATION_CONSENSUS = ROOT / "aggregation_consensus_r2_full_20260710.json"

SURFACE_DIRS = (
    "phaseD_results_campaign250_K2",
    "phaseD_results_campaign250_N1_f006",
    "phaseD_results_campaign250_N2_f006_loose",
    "phaseD_results_campaign250_N3_f006_iter",
    "phaseD_results_campaign250_O3_et20_iter",
    "phaseD_results_campaign250_O4_debias_iter",
    "phaseD_results_campaign250_O5_shape",
    "phaseD_results_campaign250_O6_b11b12",
    "phaseD_results_campaign250_P1_tau",
    "phaseD_results_campaign250_Q1_static_ext",
    "phaseD_results_campaign250_Q2_tau_ext",
    "phaseD_results_campaign250_R1_gated",
    "phaseD_results_campaign250_R2_full",
    "phaseD_results_campaign250_l2a_monthly_median3_scene_mean",
    "phaseD_results_campaign250_l2a_pc_production_mean3_scene_mean",
)

AOD_BINS = (
    (-math.inf, 0.1, "<0.1"),
    (0.1, 0.2, "0.1-0.2"),
    (0.2, 0.4, "0.2-0.4"),
    (0.4, 0.6, "0.4-0.6"),
    (0.6, 1.0, "0.6-1.0"),
    (1.0, 1.5, "1.0-1.5"),
    (1.5, math.inf, ">=1.5"),
)


def _load_jsons(directory: Path, pattern: str = "*.json") -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob(pattern)):
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
            mid = str(row.get("matchup_id") or "")
            if mid:
                rows[mid] = row
        except (OSError, ValueError, TypeError):
            continue
    return rows


def _load_json_file(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}


def _read_ids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def _is_scored(row: dict[str, Any]) -> bool:
    if "within_ee" in row:
        return row.get("within_ee") is not None
    return all(row.get(key) is not None for key in ("err", "ee"))


def _hit(row: dict[str, Any]) -> bool:
    if "within_ee" in row:
        return bool(row["within_ee"])
    return abs(float(row["err"])) <= float(row["ee"])


def _truth(row: dict[str, Any]) -> float:
    return float(row["truth"])


def _retrieved(row: dict[str, Any]) -> float:
    return float(row.get("retrieved", row.get("aod")))


def _score(rows: Iterable[dict[str, Any]]) -> tuple[int, int, float]:
    scored = [row for row in rows if _is_scored(row)]
    hits = sum(_hit(row) for row in scored)
    return hits, len(scored), 100.0 * hits / len(scored) if scored else math.nan


def _cohort_score(
    rows: dict[str, dict[str, Any]], ids: set[str]
) -> tuple[int, int, float, float, float]:
    cohort = [rows[mid] for mid in ids if mid in rows and _is_scored(rows[mid])]
    hits, count, pct = _score(cohort)
    truths = sorted(_truth(row) for row in cohort)
    median = _median(truths)
    p90 = truths[int(0.9 * (len(truths) - 1))] if truths else math.nan
    return hits, count, pct, median, p90


def _median(values: Iterable[float]) -> float:
    vals = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not vals:
        return math.nan
    middle = len(vals) // 2
    return vals[middle] if len(vals) % 2 else 0.5 * (vals[middle - 1] + vals[middle])


def _pearson(xs: Iterable[float], ys: Iterable[float]) -> float:
    pairs = [
        (float(x), float(y))
        for x, y in zip(xs, ys)
        if math.isfinite(float(x)) and math.isfinite(float(y))
    ]
    if len(pairs) < 2:
        return math.nan
    xbar = sum(x for x, _ in pairs) / len(pairs)
    ybar = sum(y for _, y in pairs) / len(pairs)
    top = sum((x - xbar) * (y - ybar) for x, y in pairs)
    sx = math.sqrt(sum((x - xbar) ** 2 for x, _ in pairs))
    sy = math.sqrt(sum((y - ybar) ** 2 for _, y in pairs))
    return top / (sx * sy) if sx > 0.0 and sy > 0.0 else math.nan


def _aod_bin_rows(rows: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for low, high, label in AOD_BINS:
        selected = [row for row in rows.values() if _is_scored(row) and low <= _truth(row) < high]
        hits, count, pct = _score(selected)
        out.append({"aod_bin": label, "hits": hits, "count": count, "pct": pct})
    return out


def _prior_calibration(r2: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    import numpy as np
    from pyproj import Transformer

    values: dict[str, list[tuple[float, float, bool]]] = {"B02": [], "B04": []}
    for mid, row in r2.items():
        if row.get("status") != "OK":
            continue
        cube_path = R2_COST_DIR / f"{mid}.npz"
        if not cube_path.exists():
            continue
        try:
            with np.load(cube_path, allow_pickle=False) as cube:
                band_names = [str(value) for value in cube["band_names"]]
                aot_axis = np.asarray(cube["aot_axis"], dtype=float)
                aot_index = int(np.argmin(np.abs(aot_axis - float(row["truth"]))))
                x = np.asarray(cube["x"], dtype=float)
                y = np.asarray(cube["y"], dtype=float)
                site_x, site_y = Transformer.from_crs(
                    "EPSG:4326", str(row["scene_crs"]), always_xy=True
                ).transform(float(row["lon"]), float(row["lat"]))
                x_index = int(np.argmin(np.abs(x - site_x)))
                y_index = int(np.argmin(np.abs(y - site_y)))
                residual = cube["band_residual_cube"]
                uncertainty = cube["boa_unc"]
                for band in values:
                    band_index = band_names.index(band)
                    bias = float(residual[band_index, aot_index, y_index, x_index])
                    sigma = max(float(uncertainty[band_index, y_index, x_index]), 1.0e-6)
                    if math.isfinite(bias) and math.isfinite(sigma):
                        values[band].append((bias, sigma, bool(row["within_ee"])))
        except (KeyError, TypeError, ValueError, OSError, IndexError):
            continue
    out = []
    for band, band_values in values.items():
        abs_bias = [item[0] for item in band_values]
        sigma = [item[1] for item in band_values]
        z = [bias / unc for bias, unc, _ in band_values]
        out.append(
            {
                "band": band,
                "count": len(band_values),
                "median_abs_bias": _median(abs_bias),
                "median_sigma": _median(sigma),
                "median_z": _median(z),
                "one_sigma_coverage_pct": 100.0 * sum(value <= 1.0 for value in z) / len(z),
                "two_sigma_coverage_pct": 100.0 * sum(value <= 2.0 for value in z) / len(z),
            }
        )
    return out


def _surface_oracle(campaign_ids: set[str]) -> tuple[int, int, float]:
    arms = [_load_jsons(ROOT / directory) for directory in SURFACE_DIRS]
    hits = 0
    for mid in campaign_ids:
        if any(mid in arm and _is_scored(arm[mid]) and _hit(arm[mid]) for arm in arms):
            hits += 1
    return hits, len(campaign_ids), 100.0 * hits / len(campaign_ids)


def _l1c_sidecar_profile(campaign_ids: set[str]) -> dict[str, Any]:
    work = ROOT / "l1c_test" / "seasonal"
    q60_values: list[float] = []
    implied_selected_aod: list[float] = []
    wvp_values: list[float] = []
    clean_count: dict[str, int] = dict.fromkeys(campaign_ids, 0)
    total_count: dict[str, int] = dict.fromkeys(campaign_ids, 0)
    for mid in campaign_ids:
        for path in (work / mid).glob("*_meta.json"):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError, TypeError):
                continue
            aod_values: list[float] = []
            for scene in payload.get("selected", []):
                try:
                    aod = float(scene["maiac"])
                    wvp = float(scene["wvp"])
                except (KeyError, TypeError, ValueError):
                    continue
                if math.isfinite(aod):
                    aod_values.append(aod)
                    total_count[mid] += 1
                    clean_count[mid] += int(aod <= 0.15)
                if math.isfinite(wvp):
                    wvp_values.append(wvp)
            if not aod_values:
                continue
            ranked = sorted(aod_values)
            keep = max(1, min(len(ranked), round(0.6 * len(ranked))))
            implied_selected_aod.extend(ranked[:keep])
            q_index = min(len(ranked) - 1, max(0, round(0.6 * (len(ranked) - 1))))
            q60_values.append(ranked[q_index])
    with_meta = [mid for mid in campaign_ids if total_count[mid] > 0]
    outside_wvp = [value for value in wvp_values if value < 0.5 or value > 4.5]
    return {
        "sites_with_metadata": len(with_meta),
        "windows": len(q60_values),
        "scenes": sum(total_count.values()),
        "median_window_q60_aod": _median(q60_values),
        "windows_q60_above_0p15_pct": 100.0
        * sum(value > 0.15 for value in q60_values)
        / len(q60_values),
        "implied_selected_above_0p15_pct": 100.0
        * sum(value > 0.15 for value in implied_selected_aod)
        / len(implied_selected_aod),
        "sites_with_metadata_zero_aod_le_0p15": sum(clean_count[mid] == 0 for mid in with_meta),
        "sites_with_metadata_fewer_than_5_aod_le_0p15": sum(
            clean_count[mid] < 5 for mid in with_meta
        ),
        "wvp_count": len(wvp_values),
        "wvp_outside_0p5_4p5_count": len(outside_wvp),
        "wvp_outside_0p5_4p5_pct": 100.0 * len(outside_wvp) / len(wvp_values),
    }


def _esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


class Tooltips:
    def __init__(self) -> None:
        self.index = 0

    def wrap(self, value: str, dataset: str, *, source: str = "local experiment artifacts") -> str:
        self.index += 1
        tip_id = f"source-value-{self.index}"
        return (
            f'<span class="source-tooltip" tabindex="0" aria-describedby="{tip_id}">'
            f"{_esc(value)}"
            f'<span class="source-tooltip-content" id="{tip_id}" role="tooltip">'
            f"Source: {_esc(source)}<br>Dataset: {_esc(dataset)}</span></span>"
        )


def _grouped_bar_svg(
    rows: list[dict[str, Any]],
    *,
    category: str,
    series: tuple[tuple[str, str, str], ...],
    aria_label: str,
) -> str:
    width, height = 960, 420
    left, right, top, bottom = 68, 24, 28, 72
    plot_w, plot_h = width - left - right, height - top - bottom
    group_w = plot_w / max(len(rows), 1)
    bar_gap = 5.0
    bar_w = min(34.0, (group_w - 18.0) / max(len(series), 1) - bar_gap)
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(aria_label)}">',
        "<style>.axis{fill:var(--secondary);font:12px var(--ds-font)}.value{fill:var(--text);font:11px var(--ds-font)}.grid{stroke:var(--grid);stroke-width:1}</style>",
    ]
    for tick in range(0, 101, 20):
        y = top + plot_h * (1.0 - tick / 100.0)
        parts.append(
            f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}"/>'
        )
        parts.append(
            f'<text class="axis" x="{left - 10}" y="{y + 4:.1f}" text-anchor="end">{tick}%</text>'
        )
    for index, row in enumerate(rows):
        center = left + group_w * (index + 0.5)
        total_w = len(series) * bar_w + (len(series) - 1) * bar_gap
        start = center - total_w / 2.0
        for si, (field, _label, color) in enumerate(series):
            value = max(0.0, min(100.0, float(row[field])))
            x = start + si * (bar_w + bar_gap)
            y = top + plot_h * (1.0 - value / 100.0)
            h = top + plot_h - y
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" '
                f'fill="{color}" rx="2"/>'
            )
            parts.append(
                f'<text class="value" x="{x + bar_w / 2:.1f}" y="{max(y - 5, 12):.1f}" '
                f'text-anchor="middle">{value:.1f}</text>'
            )
        parts.append(
            f'<text class="axis" x="{center:.1f}" y="{height - bottom + 24}" text-anchor="middle">'
            f"{_esc(row[category])}</text>"
        )
    legend_x = left
    for _field, label, color in series:
        parts.append(
            f'<rect x="{legend_x}" y="{height - 24}" width="12" height="12" fill="{color}" rx="2"/>'
        )
        parts.append(
            f'<text class="axis" x="{legend_x + 18}" y="{height - 14}">{_esc(label)}</text>'
        )
        legend_x += 180
    parts.append("</svg>")
    return "".join(parts)


def _chart_figure(
    *,
    chart_id: str,
    title: str,
    subtitle: str,
    fallback_svg: str,
    note: str,
    source_id: str,
    datasets: str,
) -> str:
    return f"""
      <div class="wide">
        <figure class="card source-figure">
          <div class="card-head"><h3>{_esc(title)}</h3><p>{_esc(subtitle)}</p></div>
          <div class="chart-wrap">
            <div data-recharts-chart="{_esc(chart_id)}">
              <div class="chart-fallback" data-recharts-fallback>{fallback_svg}</div>
              <div data-recharts-live aria-hidden="true"></div>
            </div>
          </div>
          <figcaption class="chart-note">{_esc(note)}</figcaption>
          <button type="button" class="source-tooltip" aria-describedby="{_esc(source_id)}">Source
            <span class="source-tooltip-content" id="{_esc(source_id)}" role="tooltip">Source: local experiment artifacts<br>Dataset: {_esc(datasets)}</span>
          </button>
        </figure>
      </div>
    """


def _payload_chart(
    *, chart_id: str, title: str, data: list[dict[str, Any]], category: str
) -> dict[str, Any]:
    long_rows = []
    for row in data:
        for field, label in (("r2_pct", "R2"), ("sdspec_pct", "acixThree sdspec")):
            long_rows.append(
                {category: row[category], "series": label, "within_ee": float(row[field]) / 100.0}
            )
    return {
        "id": chart_id,
        "height": 360,
        "type": "bar",
        "dataset": {
            "id": chart_id,
            "title": title,
            "data": long_rows,
            "chart_spec": {
                "id": chart_id,
                "dataset": chart_id,
                "title": title,
                "type": "bar",
                "encodings": {
                    "x": {"field": category, "type": "nominal"},
                    "y": {"field": "within_ee", "label": "Within EE", "type": "quantitative"},
                    "color": {"field": "series", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "Within expected error",
                "valueFormat": "percent",
                "settings": {"groupMode": "grouped"},
            },
        },
    }


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all62 = _read_ids(ROOT / "all62_mids.txt")
    campaign = _read_ids(ROOT / "campaign250_mids.txt")
    r2 = _load_jsons(R2_DIR)
    sdspec = _load_jsons(SDSPEC_DIR, SDSPEC_PATTERN)
    experiment_payload = _load_json_file(EXPERIMENT_SUMMARY)
    experiment_rows = {str(row.get("tag")): row for row in experiment_payload.get("summaries", [])}
    paired_rows = {str(row.get("candidate")): row for row in experiment_payload.get("paired", [])}
    band_payload = _load_json_file(BAND_CONSENSUS_FINAL)
    band_rows = {str(row.get("strategy")): row for row in band_payload.get("summaries", [])}
    aggregation_payload = _load_json_file(AGGREGATION_CONSENSUS)
    prior_quality_payload = _load_json_file(PRIOR_QUALITY_FINAL)
    acix9 = _load_jsons(
        SDSPEC_DIR,
        "*__seas_maiac_e5_acix9_fixed_knn_l2awvp_noback_smoke.json",
    )
    acix9_gate1 = _load_jsons(
        SDSPEC_DIR,
        "*__seas_maiac_e5_acix9_knn_l2awvp_noback_gate1_smoke.json",
    )
    no_back_records = _load_jsons(
        SDSPEC_DIR,
        "*__seas_maiac_e3_knn_l2awvp_nobackstop_axis4.json",
    )

    r2_hits, r2_valid, r2_valid_pct = _score(r2.values())
    sd_hits, sd_valid, sd_valid_pct = _score(sdspec.values())
    acix9_hits, acix9_valid, acix9_pct = _score(acix9.values())
    acix9_common = sorted(set(acix9) & set(no_back_records))
    acix9_gains = sum(_hit(acix9[mid]) and not _hit(no_back_records[mid]) for mid in acix9_common)
    acix9_losses = sum(_hit(no_back_records[mid]) and not _hit(acix9[mid]) for mid in acix9_common)
    acix9_gate1_hits, acix9_gate1_valid, acix9_gate1_pct = _score(acix9_gate1.values())
    acix9_gate1_common = sorted(set(acix9_gate1) & set(no_back_records))
    acix9_gate1_gains = sum(
        _hit(acix9_gate1[mid]) and not _hit(no_back_records[mid]) for mid in acix9_gate1_common
    )
    acix9_gate1_losses = sum(
        _hit(no_back_records[mid]) and not _hit(acix9_gate1[mid]) for mid in acix9_gate1_common
    )
    surface_oracle_hits, _, surface_oracle_pct = _surface_oracle(campaign)
    target_hits = 201

    r2_cohorts = {
        "all62": _cohort_score(r2, all62),
        "campaign_rest": _cohort_score(r2, campaign - all62),
    }
    sd_cohorts = {
        "all62": _cohort_score(sdspec, all62),
        "campaign_rest": _cohort_score(sdspec, campaign - all62),
    }
    cohort_chart = [
        {
            "cohort": "Selected all-62",
            "r2_pct": r2_cohorts["all62"][2],
            "sdspec_pct": sd_cohorts["all62"][2],
            "r2_n": r2_cohorts["all62"][1],
            "sdspec_n": sd_cohorts["all62"][1],
        },
        {
            "cohort": "Campaign remainder",
            "r2_pct": r2_cohorts["campaign_rest"][2],
            "sdspec_pct": sd_cohorts["campaign_rest"][2],
            "r2_n": r2_cohorts["campaign_rest"][1],
            "sdspec_n": sd_cohorts["campaign_rest"][1],
        },
    ]

    r2_bins = {row["aod_bin"]: row for row in _aod_bin_rows(r2)}
    sd_bins = {row["aod_bin"]: row for row in _aod_bin_rows(sdspec)}
    bin_chart = [
        {
            "aod_bin": label,
            "r2_pct": r2_bins[label]["pct"],
            "sdspec_pct": sd_bins[label]["pct"],
            "r2_n": r2_bins[label]["count"],
            "sdspec_n": sd_bins[label]["count"],
        }
        for _, _, label in AOD_BINS
    ]

    prior_calibration = _prior_calibration(r2)
    sidecar_profile = _l1c_sidecar_profile(campaign)
    r2_ok = [row for row in r2.values() if row.get("status") == "OK" and _is_scored(row)]
    cloud_corr = _pearson(
        (float(row.get("cloud_frac", math.nan)) for row in r2_ok),
        (abs(float(row["err"])) for row in r2_ok),
    )
    cost_corr = _pearson(
        (float(row.get("solver", {}).get("cost_final_per_band", math.nan)) for row in r2_ok),
        (abs(float(row["err"])) for row in r2_ok),
    )

    sd_maiac_rows = [
        row for row in sdspec.values() if _is_scored(row) and row.get("maiac_aot") is not None
    ]
    sd_near_maiac_005 = sum(
        abs(_retrieved(row) - float(row["maiac_aot"])) <= 0.05 for row in sd_maiac_rows
    )
    sd_near_maiac_pct = (
        100.0 * sd_near_maiac_005 / len(sd_maiac_rows) if sd_maiac_rows else math.nan
    )

    direction = Counter(
        "under" if float(row["err"]) < 0 else "over" for row in r2_ok if not _hit(row)
    )
    generated = datetime.now().astimezone().strftime("%d %B %Y %H:%M %Z")

    analysis = {
        "generated": generated,
        "metric": "abs(retrieved AOD550 - AERONET AOD550) <= 0.05 + 0.15*AERONET",
        "campaign_denominator": 250,
        "strict_above_80_target_hits": target_hits,
        "r2": {
            "hits": r2_hits,
            "valid": r2_valid,
            "raw_pct": 100.0 * r2_hits / 250.0,
            "valid_pct": r2_valid_pct,
            "under_failures": direction["under"],
            "over_failures": direction["over"],
        },
        "sdspec": {
            "hits": sd_hits,
            "available": sd_valid,
            "available_pct": sd_valid_pct,
            "raw_lower_bound_pct": 100.0 * sd_hits / 250.0,
            "raw_upper_bound_pct": 100.0 * (sd_hits + 250 - sd_valid) / 250.0,
            "within_0p05_of_maiac_pct": sd_near_maiac_pct,
            "maiac_comparison_n": len(sd_maiac_rows),
        },
        "surface_oracle": {
            "hits": surface_oracle_hits,
            "raw_pct": surface_oracle_pct,
        },
        "all62_vs_remainder": {"r2": r2_cohorts, "sdspec": sd_cohorts},
        "aod_bins": {"r2": list(r2_bins.values()), "sdspec": list(sd_bins.values())},
        "prior_calibration": prior_calibration,
        "l1c_sidecar_profile": sidecar_profile,
        "correlations": {
            "cost_vs_absolute_error_pearson": cost_corr,
            "cloud_fraction_vs_absolute_error_pearson": cloud_corr,
        },
        "coverage": {
            "l1c_priors": sum(
                (ROOT / "l1c_seasonal_species_lut" / f"{mid}.npz").exists() for mid in campaign
            ),
            "prior_quality": sum((PQ_DIR / f"{mid}.json").exists() for mid in campaign),
            "sdspec_results": sd_valid,
        },
        "staged_experiments": {
            "summary": experiment_payload,
            "band_consensus": band_payload,
            "aggregation_consensus": aggregation_payload.get("aggregation_consensus"),
            "prior_quality": prior_quality_payload,
            "acix9_fixed": {
                "hits": acix9_hits,
                "available": acix9_valid,
                "available_pct": acix9_pct,
                "paired_common": len(acix9_common),
                "paired_gains": acix9_gains,
                "paired_losses": acix9_losses,
            },
            "acix9_gate1": {
                "hits": acix9_gate1_hits,
                "available": acix9_gate1_valid,
                "available_pct": acix9_gate1_pct,
                "paired_common": len(acix9_gate1_common),
                "paired_gains": acix9_gate1_gains,
                "paired_losses": acix9_gate1_losses,
            },
        },
    }
    (OUT / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")

    tooltip = Tooltips()
    r2_dataset = "phaseD_results_campaign250_R2_full/*.json"
    sd_dataset = "seasonal_val/*__seas_maiac_sdspec_c250_20260708.json"
    prior_dataset = (
        "phaseD_results_campaign250_R2_full/*.json + "
        "phaseD_cost_cubes_campaign250_R2_full_localdiag_20260705/*.npz"
    )
    local_dataset = (
        "phaseD_results_campaign250_R2_full_localdiag_20260705 + "
        "phaseD_cost_cubes_campaign250_R2_full_localdiag_20260705"
    )

    all62_median = r2_cohorts["all62"][3]
    rest_median = r2_cohorts["campaign_rest"][3]
    all62_median_label = f"{all62_median:.3f}"
    rest_median_label = f"{rest_median:.3f}"
    r2_low_pct = float(r2_bins["<0.1"]["pct"])
    r2_mid_pct = float(r2_bins["0.4-0.6"]["pct"])
    r2_high_pct = float(r2_bins["1.0-1.5"]["pct"])
    r2_low_label = f"{r2_low_pct:.1f}%"
    r2_mid_label = f"{r2_mid_pct:.1f}%"
    r2_high_label = f"{r2_high_pct:.1f}%"
    sidecar_windows = int(sidecar_profile["windows"])
    sidecar_q60 = float(sidecar_profile["median_window_q60_aod"])
    selected_above_pct = float(sidecar_profile["implied_selected_above_0p15_pct"])
    zero_clean_sites = int(sidecar_profile["sites_with_metadata_zero_aod_le_0p15"])
    wvp_outside_pct = float(sidecar_profile["wvp_outside_0p5_4p5_pct"])
    sidecar_q60_label = f"{sidecar_q60:.3f}"
    selected_above_label = f"{selected_above_pct:.1f}%"
    wvp_outside_label = f"{wvp_outside_pct:.1f}%"
    l1c_prior_count = int(analysis["coverage"]["l1c_priors"])
    l1c_coverage_label = f"{l1c_prior_count}/250"

    def _experiment(tag: str) -> dict[str, Any]:
        return experiment_rows.get(tag, {})

    def _paired(tag: str) -> dict[str, Any]:
        return paired_rows.get(tag, {})

    no_back = _experiment("e3_knn_l2awvp_nobackstop_axis4")
    clean_knn = _experiment("e3_clean015min3_knn_model013_nobackstop")
    band_consensus = band_rows.get("spread030_median", {})
    aggregation_consensus = aggregation_payload.get("aggregation_consensus", {})
    experiment_dataset = EXPERIMENT_SUMMARY.name
    acix_dataset = "seasonal_val/*__seas_maiac_e5_acix9*_smoke.json"

    no_back_label = (
        f"{int(no_back.get('hits', 0))}/{int(no_back.get('available', 0))} = "
        f"{float(no_back.get('available_pct', math.nan)):.1f}%"
    )
    clean_knn_label = (
        f"{int(clean_knn.get('hits', 0))}/{int(clean_knn.get('available', 0))} = "
        f"{float(clean_knn.get('available_pct', math.nan)):.1f}%"
    )
    acix9_label = f"{acix9_hits}/{acix9_valid} = {acix9_pct:.1f}%"
    acix9_gate1_label = f"{acix9_gate1_hits}/{acix9_gate1_valid} = {acix9_gate1_pct:.1f}%"
    band_consensus_label = (
        f"{int(band_consensus.get('hits', 0))}/{int(band_consensus.get('count', 0))}"
    )
    aggregation_label = (
        f"{len(aggregation_consensus.get('gains', []))} gains / "
        f"{len(aggregation_consensus.get('losses', []))} losses"
    )

    title = "Surface-driven AOD: staged tests and remaining failure"
    summary_html = f"""
      <p><strong>No tested surface-only configuration qualifies for a campaign-250 run.</strong>
      The no-backstop 6S kNN retrieval scores {tooltip.wrap(no_back_label, experiment_dataset)}
      on its matched development outputs. Its best conservative band-disagreement correction reaches only
      {tooltip.wrap(band_consensus_label, BAND_CONSENSUS_FINAL.name)}, far below the promotion budget.</p>
      <p><strong>Clean-day filtering and exact WVP do not improve paired retrieval.</strong>
      Clean-day plus kNN reaches {tooltip.wrap(clean_knn_label, experiment_dataset)},
      but has {tooltip.wrap("0 gains / 1 loss", experiment_dataset)} against the original kNN result on those same sites.
      Exact historical WVP, absolute cost, shape cost, and the ACIX reliability mask are neutral or worse.</p>
      <p><strong>The corrected four-anchor ACIX replication is also below the gate.</strong>
      After fixing omitted B06 loading and the seven-band reshape, it scores
      {tooltip.wrap(acix9_label, acix_dataset)} without the reliability mask and
      {tooltip.wrap(acix9_gate1_label, acix_dataset)} with it. The independent full-campaign aggregation check gives
      {tooltip.wrap(aggregation_label, AGGREGATION_CONSENSUS.name)}. E5 was therefore not run: a full 250 sweep would only measure a configuration already falsified by its matched controls.</p>
    """

    cohort_svg = _grouped_bar_svg(
        cohort_chart,
        category="cohort",
        series=(("r2_pct", "R2", "#0169cc"), ("sdspec_pct", "acixThree sdspec", "#e25507")),
        aria_label="Within expected error for the selected all-62 cohort and the campaign remainder",
    )
    bins_svg = _grouped_bar_svg(
        bin_chart,
        category="aod_bin",
        series=(("r2_pct", "R2", "#0169cc"), ("sdspec_pct", "acixThree sdspec", "#e25507")),
        aria_label="Within expected error by AERONET AOD bin",
    )

    prior_row_parts = []
    for row in prior_calibration:
        band = str(row["band"])
        median_abs_bias = float(row["median_abs_bias"])
        median_sigma = float(row["median_sigma"])
        median_z = float(row["median_z"])
        one_sigma_pct = float(row["one_sigma_coverage_pct"])
        two_sigma_pct = float(row["two_sigma_coverage_pct"])
        bias_label = f"{median_abs_bias:.4f}"
        sigma_label = f"{median_sigma:.4f}"
        z_label = f"{median_z:.2f}x"
        one_sigma_label = f"{one_sigma_pct:.1f}%"
        two_sigma_label = f"{two_sigma_pct:.1f}%"
        prior_row_parts.append(
            f"""
            <tr><td>{_esc(band)}</td>
              <td>{tooltip.wrap(bias_label, prior_dataset)}</td>
              <td>{tooltip.wrap(sigma_label, prior_dataset)}</td>
              <td>{tooltip.wrap(z_label, prior_dataset)}</td>
              <td>{tooltip.wrap(one_sigma_label, prior_dataset)}</td>
              <td>{tooltip.wrap(two_sigma_label, prior_dataset)}</td></tr>
            """
        )
    prior_rows = "".join(prior_row_parts)

    prior_paired = {
        str(row.get("candidate")): row for row in prior_quality_payload.get("paired", [])
    }
    prior_quality_specs = (
        ("Exact historical L2A WVP", "exact"),
        ("Absolute clean-day with fallback", "clean015"),
        ("Whole-window clean-day", "global3"),
        ("ACIX nine-band composite", "acix9"),
    )

    def _prior_band_cell(bands: dict[str, Any], band: str) -> tuple[str, str]:
        values = bands.get(band, {})
        errors = (
            f"{float(values.get('baseline_median_abs_error', math.nan)):.4f} -> "
            f"{float(values.get('candidate_median_abs_error', math.nan)):.4f}"
        )
        counts = f"{int(values.get('improved', 0))} improved / {int(values.get('worse', 0))} worse"
        return errors, counts

    prior_quality_row_parts = []
    for label, key in prior_quality_specs:
        row = prior_paired.get(key, {})
        bands = row.get("bands", {})
        b02_errors, b02_counts = _prior_band_cell(bands, "B02")
        b04_errors, b04_counts = _prior_band_cell(bands, "B04")
        prior_quality_row_parts.append(
            "<tr>"
            f"<td>{_esc(label)}</td>"
            f"<td>{tooltip.wrap(str(int(row.get('common', 0))), PRIOR_QUALITY_FINAL.name)}</td>"
            f"<td>{tooltip.wrap(b02_errors, PRIOR_QUALITY_FINAL.name)}</td>"
            f"<td>{tooltip.wrap(b02_counts, PRIOR_QUALITY_FINAL.name)}</td>"
            f"<td>{tooltip.wrap(b04_errors, PRIOR_QUALITY_FINAL.name)}</td>"
            f"<td>{tooltip.wrap(b04_counts, PRIOR_QUALITY_FINAL.name)}</td>"
            "</tr>"
        )
    prior_quality_rows = "".join(prior_quality_row_parts)

    failure_groups = (
        (
            "Moderate-AOD under-read",
            14,
            "Surface prior and uncertainty are wrong or overconfident.",
        ),
        ("Clean/moderate under-read", 12, "Bright prior or backstop snaps the solution low."),
        (
            "High-cost under-read",
            10,
            "Aerosol optical model is incompatible, especially thick smoke.",
        ),
        ("Over-read", 9, "Dark/contaminated prior, adjacency, or residual cloud."),
        (
            "Low-cost high-AOD under-read",
            8,
            "Wrong aerosol optics still produces an apparently good fit.",
        ),
        ("Other hard failures", 8, "Mixed edge, missing-support, or extreme-AOD behavior."),
    )
    failure_rows = "".join(
        f"<tr><td>{_esc(name)}</td><td>{tooltip.wrap(str(count), 'campaign250_surface_unsolved61_mids.txt + ' + FINDINGS.name)}</td><td>{_esc(meaning)}</td></tr>"
        for name, count, meaning in failure_groups
    )

    experiment_specs = (
        (
            "Pure surface, auto2 cost",
            "e3_knn_l2awvp_nobackstop_axis4",
            "Reference surface-only arm; below the stage gate.",
        ),
        (
            "ACIX reliability mask",
            "e3_knn_l2awvp_nobackstop_gate1_axis4",
            "No paired AOD changes; masking the outer uncertainty tails is insufficient.",
        ),
        (
            "Absolute cost only",
            "e3_knn_l2awvp_nobackstop_abs_axis4",
            "Paired result is negative.",
        ),
        (
            "Shape cost only",
            "e3_knn_l2awvp_nobackstop_shape_axis4",
            "Paired result is negative.",
        ),
        (
            "Exact historical WVP + kNN",
            "e3_exacthistwvp_rel60_knn_l2awvp_nobackstop",
            "Removes clamping but does not improve retrieval.",
        ),
        (
            "Clean-day prior + kNN",
            "e3_clean015min3_knn_model013_nobackstop",
            "Small-cohort percentage is not a gain: paired result is net negative.",
        ),
        (
            "Global clean-window prior + kNN",
            "e3_clean015global3_knn_model013_nobackstop_gate1",
            "Filtering whole windows removes support and remains net negative.",
        ),
    )
    experiment_row_parts = []
    for name, tag, decision in experiment_specs:
        row = _experiment(tag)
        paired = _paired(tag)
        score = (
            f"{int(row.get('hits', 0))}/{int(row.get('available', 0))} = "
            f"{float(row.get('available_pct', math.nan)):.1f}%"
        )
        paired_label = (
            "reference"
            if tag == "e3_knn_l2awvp_nobackstop_axis4"
            else (
                f"+{int(paired.get('gains', 0))} / -{int(paired.get('losses', 0))} "
                f"(n={int(paired.get('common', 0))})"
            )
        )
        experiment_row_parts.append(
            "<tr>"
            f"<td>{_esc(name)}</td>"
            f"<td>{tooltip.wrap(score, experiment_dataset)}</td>"
            f"<td>{tooltip.wrap(paired_label, experiment_dataset)}</td>"
            f"<td>{_esc(decision)}</td>"
            "</tr>"
        )
    experiment_row_parts.append(
        "<tr>"
        "<td>Corrected ACIX four-anchor kNN, no mask</td>"
        f"<td>{tooltip.wrap(acix9_label, acix_dataset)}</td>"
        f"<td>{tooltip.wrap(f'+{acix9_gains} / -{acix9_losses} (n={len(acix9_common)})', acix_dataset)}</td>"
        "<td>B06/9-band defects fixed; the valid four-anchor run remains below the gate.</td>"
        "</tr>"
    )
    experiment_row_parts.append(
        "<tr>"
        "<td>Corrected ACIX four-anchor kNN, reliability mask</td>"
        f"<td>{tooltip.wrap(acix9_gate1_label, acix_dataset)}</td>"
        f"<td>{tooltip.wrap(f'+{acix9_gate1_gains} / -{acix9_gate1_losses} (n={len(acix9_gate1_common)})', acix_dataset)}</td>"
        "<td>Tests whether ACIX uncertainty-tail masking changes the corrected four-anchor result.</td>"
        "</tr>"
    )
    experiment_table_rows = "".join(experiment_row_parts)

    main = f"""
    <main data-report-audience="technical">
      <article class="reading">
        <div class="kicker">Technical diagnosis and experiment design</div>
        <header data-contract-section="title"><h1>{_esc(title)}</h1></header>
        <p class="deck">Campaign-250 AOD550, evaluated against AERONET expected error. Missing or failed cases count against the raw 250 denominator.</p>
        <section class="summary" data-contract-section="technical-summary">
          <div class="summary-label">Technical Summary</div>
          <div class="summary-body">{summary_html}</div>
        </section>
        <section class="metrics">
          <div class="metric"><div class="metric-label">Best single surface solve</div><div class="metric-value">{tooltip.wrap(f"{r2_hits}/250", r2_dataset)}</div><div class="metric-note">{tooltip.wrap(f"{100 * r2_hits / 250:.1f}%", r2_dataset)} raw within EE</div></div>
          <div class="metric"><div class="metric-label">Exact sdspec available score</div><div class="metric-value">{tooltip.wrap(f"{sd_hits}/{sd_valid}", sd_dataset)}</div><div class="metric-note">{tooltip.wrap(f"{sd_valid_pct:.1f}%", sd_dataset)} on completed outputs; partial coverage</div></div>
          <div class="metric"><div class="metric-label">Surface-only saved-arm oracle</div><div class="metric-value">{tooltip.wrap(f"{surface_oracle_hits}/250", "phaseD_results_campaign250_K2 through R2 plus fixed L2A ensembles")}</div><div class="metric-note">Still below a usable single generic method</div></div>
          <div class="metric"><div class="metric-label">Strictly above 80% target</div><div class="metric-value">{_esc(str(target_hits))}/250</div><div class="metric-note">Requires +{target_hits - r2_hits} hits over R2</div></div>
        </section>
        <section class="narrative" data-contract-section="key-findings">
          <h2>The 87% result was measured on an easier, low-AOD cohort</h2>
          <p>The selected all-62 cohort has median truth AOD {tooltip.wrap(all62_median_label, r2_dataset)}, versus
          {tooltip.wrap(rest_median_label, r2_dataset)} for the campaign remainder. The same implementation loses about thirty percentage points outside that cohort, so neither the all-62 R2 nor sdspec score estimates campaign performance.</p>
        </section>
      </article>

      {_chart_figure(chart_id="cohort-score", title="Within-EE by evaluation cohort", subtitle="R2 is complete; the sdspec campaign backfill is partial", fallback_svg=cohort_svg, note="Percentages use available scored records inside each cohort. The raw acceptance metric remains hits divided by the fixed campaign population.", source_id="chart-source-cohort", datasets=r2_dataset + "; " + sd_dataset)}

      <article class="reading">
        <section class="narrative">
          <h2>The retrieval breaks at moderate and high AOD, not at a cloud threshold</h2>
          <p>R2 is {tooltip.wrap(r2_low_label, r2_dataset)} within EE below AOD 0.1, but falls to
          {tooltip.wrap(r2_mid_label, r2_dataset)} at 0.4-0.6 and
          {tooltip.wrap(r2_high_label, r2_dataset)} at 1.0-1.5. Cost is a materially stronger severity signal than cloud fraction: Pearson correlation with absolute error is
          {tooltip.wrap(f"{cost_corr:+.3f}", r2_dataset)} for final cost and only
          {tooltip.wrap(f"{cloud_corr:+.3f}", r2_dataset)} for scene cloud fraction.</p>
        </section>
      </article>

      {_chart_figure(chart_id="aod-bin-score", title="Within-EE by AERONET AOD bin", subtitle="Complete R2 compared with the available exact sdspec outputs", fallback_svg=bins_svg, note="The campaign deliberately contains many moderate and high-AOD cases. The sdspec series is partial but already reproduces the same failure shape.", source_id="chart-source-aod", datasets=r2_dataset + "; " + sd_dataset)}

      <article class="reading">
        <section class="narrative">
          <h2>The reported surface uncertainty is overconfident</h2>
          <p>The current ensemble spread measures variability among prior realizations, not error against the target scene. It therefore stays small when every historical realization is consistently wrong. This is why flat cost curves and band disagreement are common while the solver still assigns high confidence.</p>
        </section>
        <section class="card table-card">
          <div class="card-head"><h3>Visible surface-prior calibration</h3><p>Station-pixel residual at the cost node nearest AERONET truth AOD</p></div>
          <div class="table-scroll"><table><thead><tr><th>Band</th><th>Median |prior-ref|</th><th>Median stated sigma</th><th>Median standardized error</th><th>Within 1 sigma</th><th>Within 2 sigma</th></tr></thead><tbody>{prior_rows}</tbody></table></div>
        </section>
        <section class="card table-card">
          <div class="card-head"><h3>Paired prior-builder checks</h3><p>Old prior -> candidate median absolute BOA error on exactly the same matchups</p></div>
          <div class="table-scroll"><table><thead><tr><th>Candidate prior</th><th>Common</th><th>B02 median error</th><th>B02 cases</th><th>B04 median error</th><th>B04 cases</th></tr></thead><tbody>{prior_quality_rows}</tbody></table></div>
        </section>

        <section class="narrative">
          <h2>The current clean-day prior is often built from non-clean observations</h2>
          <p>Across {tooltip.wrap(str(sidecar_windows), "l1c_test/seasonal/*/*_meta.json")} campaign sidecar windows, the median within-window 60th-percentile MAIAC AOD is
          {tooltip.wrap(sidecar_q60_label, "l1c_test/seasonal/*/*_meta.json")}. Consequently,
          {tooltip.wrap(selected_above_label, "l1c_test/seasonal/*/*_meta.json")} of observations admitted by the relative low-AOD fraction exceed AOD 0.15. A strict threshold alone is not sufficient either: among sites with metadata,
          {tooltip.wrap(str(zero_clean_sites), "l1c_test/seasonal/*/*_meta.json")} have no historical scene at or below 0.15. The next builder therefore needs an absolute-clean preference with a minimum-sample fallback and an uncertainty penalty for residual aerosol contamination.</p>
          <p>The historical atmosphere correction also clamps
          {tooltip.wrap(wvp_outside_label, "l1c_test/seasonal/*/*_meta.json + custom_ac_phase2_sixs.py")} of saved Sen2Cor WVP values outside its
          {tooltip.wrap("0.5-4.5 cm", "custom_ac_phase2_sixs.py + bestpixel/atmosphere.py")} coefficient grid. Because the composite code uses scene-mean WVP, the simplest experiment is to compute one 6S coefficient set at each scene's actual WVP; this is both exact at that scalar and cheaper than evaluating the full TCWV node set per day.</p>
        </section>

        <section class="narrative">
          <h2>The final node is often chosen by aggregation and the AOD backstop</h2>
          <p>On the intentionally bad 88-case diagnostic set, the scene cost-curve minimum is within EE for
          {tooltip.wrap("26/88", local_dataset)}, but the final solver node is within EE for only
          {tooltip.wrap("6/88", local_dataset)}. On the exact sdspec outputs, the retrieved AOD lies within 0.05 of MAIAC for
          {tooltip.wrap(f"{sd_near_maiac_pct:.1f}%", sd_dataset)} of records with MAIAC coverage. This is not yet a genuinely surface-dominated solve.</p>
        </section>

        <section class="narrative">
          <h2>The hard residual requires both better surface information and different aerosol optics</h2>
          <p>Of the {tooltip.wrap("61", "campaign250_surface_unsolved61_mids.txt + " + FINDINGS.name)} matchups missed by every saved surface-only arm, most are under-retrievals. No single scalar gate separates them: high cost identifies severe smoke mismatch, while many moderate and dust failures have low cost and need better spectral constraints or prior calibration.</p>
        </section>
        <section class="card table-card">
          <div class="card-head"><h3>Unsolved failure anatomy</h3><p>Mutually exclusive diagnostic buckets over the surface-only unsolved set</p></div>
          <div class="table-scroll"><table><thead><tr><th>Failure group</th><th>Cases</th><th>Primary interpretation</th></tr></thead><tbody>{failure_rows}</tbody></table></div>
        </section>

        <section class="narrative" data-contract-section="scope-data-and-metric-definitions">
          <h2>Scope and metric</h2>
          <p><strong>Population:</strong> the fixed campaign-250 matchup list. <strong>Grain:</strong> one Sentinel-2 scene and AERONET matchup. <strong>Success:</strong> |retrieved AOD550 - AERONET AOD550| <= 0.05 + 0.15*AERONET. Raw scoring always divides by 250; failed or missing retrievals are misses. The exact sdspec series is reported on available outputs until all buildable cases finish and is never substituted for the raw denominator.</p>
        </section>

        <section class="narrative" data-contract-section="methodology">
          <h2>What was checked</h2>
          <p>The diagnosis reconciles complete R2 outputs, the exact acixThree-style 6S/sdspec reproduction, truth-AOD surface references, full-campaign local cost cubes, all saved surface-only arms, and the matched experiment table above. Driver checks cover AOD regime, cohort composition, prior bias and uncertainty, target and historical WVP, cost shape, B02/B04 disagreement, cloud fraction, extraction, aerosol family, and anchor iteration.</p>
          <p>Same-scene S2 L2A WVP was resolved for {tooltip.wrap("250/250", "s2_l2a_wvp_campaign250/*.json")}: 248 exact platform/time matches and two near-simultaneous S2A fallbacks for S2C records. The tested path holds raw WVP/1000 in centimetres through anchor correction and every 6S candidate; the undocumented 0.85 factor remains rejected because it has no independent closure or AOD benefit.</p>
          <p>The ACIX9 implementation audit found two concrete defects in the first smoke result: B06 was omitted from target loading, and nine-band composites were reshaped with a hard-coded seven-band count. The corrected run records non-null B02/B04 overrides from B06/B08/B11/B12 anchors. This fix makes the replication valid, but does not make it accurate.</p>
        </section>

        <section class="narrative experiment" data-contract-section="experimental-design">
          <h2>The staged tests do not support a full-campaign promotion</h2>
          <p>The development manifest contains all 89 R2 failures plus 54 AOD-stratified R2 hits. Percentages below use each arm's available records; the paired column is the promotion evidence because it compares exactly the same matchups. None of the controlled changes has positive paired lift large enough to approach the required score budget.</p>
        </section>
        <section class="card table-card experiment-table">
          <div class="card-head"><h3>Measured surface-prior and solver experiments</h3><p>6S RT, raw same-scene L2A WVP, no external AOD in the reported objective</p></div>
          <div class="table-scroll"><table><thead><tr><th>Controlled arm</th><th>Available score</th><th>Paired gains / losses</th><th>Decision</th></tr></thead><tbody>{experiment_table_rows}</tbody></table></div>
        </section>
        <section class="narrative">
          <p><strong>Water vapour consistency is retained for physical correctness, not score lift.</strong> Raw L2A WVP, the inherited 0.85 scaling, and the prior value produced no gains or losses on the 12 WVP-extreme cases. Exact historical-WVP 6S correction changes paired static B02/B04 median errors by less than 0.0002.</p>
          <p><strong>Solver and species alternatives were negative controls.</strong> Offline robust likelihood achieved {tooltip.wrap("112/250 = 44.8%", "robust_surface_likelihood_e3_20260709.json")}; fixed continental, biomass, desert, and maritime selection reached at best {tooltip.wrap("7/16", "canonical_aerosol_selector_e4_high12_20260709.json")}; anchor iteration had no material paired effect. The conservative high-spread band correction adds one hit on the 51-site surface-only sample, but the agreeing-band aggregation rule fails its independent 249-site check with {tooltip.wrap(aggregation_label, AGGREGATION_CONSENSUS.name)}.</p>
        </section>

        <section class="narrative" data-contract-section="limitations-uncertainty-and-robustness-checks">
          <h2>Limits and robustness checks</h2>
          <p>The exact sdspec campaign remains coverage-limited because only {tooltip.wrap(l1c_coverage_label, "l1c_seasonal_species_lut/*.npz")} original L1C priors are currently buildable and {tooltip.wrap(f"{sd_valid}/250", sd_dataset)} sdspec results existed at this snapshot. The new controlled arms are intentionally smaller stop/go screens, so their raw percentages are not campaign estimates. Their paired gains and losses, prior-quality gates, and independent 249-site aggregation check control the promotion decision.</p>
          <p>Prior-reference residuals use AERONET AOD to diagnose the surface model and must never become runtime features. The campaign has already been inspected repeatedly, so any future raw 250-site improvement would not be independent validation; site-grouped folds and an external surface-prior holdout remain required.</p>
          <p>The uncertainty calibration samples the cost-cube node nearest AERONET truth rather than interpolating BOA at the exact truth. It therefore includes small AOD-grid quantization as well as surface-prior and RT-model discrepancy. That is conservative for likelihood calibration, but it should not be interpreted as a pure surface-composite error measurement.</p>
          <p>Raw L2A WVP is retained because it removes an atmosphere inconsistency, but the extreme-case ablation shows no AOD lift. Prescribed biomass aerosol previously recovered only {tooltip.wrap("3/8", FINDINGS.name)} catastrophic smoke cases, while the four-family selector regressed controls. These are negative results, not evidence that another threshold on the same inputs will reach the target.</p>
        </section>

        <section class="narrative" data-contract-section="recommended-next-steps">
          <h2>Require new surface information before another AOD sweep</h2>
          <ol class="method-list">
            <li>Freeze the tested configurations as negative controls. Do not run E5 or tune more band/cost thresholds on campaign truth.</li>
            <li>Keep raw same-scene L2A WVP and the corrected ACIX9 loader as implementation fixes, but do not claim retrieval improvement from either.</li>
            <li>Develop the next prior against a site-grouped, leave-year-out BOA target before solving AOD. It must reduce B02 and B04 median error to at most 0.010 and produce calibrated one- and two-sigma coverage on held-out scenes.</li>
            <li>If B8A/B11/B12-conditioned history cannot meet that prior gate, add genuinely new surface information, such as a geometrically matched nearest clean observation or BRDF/land-cover constraint. More selection among the existing saved priors is bounded by the 75.6% oracle.</li>
            <li>Only after the prior passes should the robust covariance likelihood and penalized smoke/dust family be retested. Then run one frozen campaign-250 configuration and an external holdout.</li>
          </ol>
        </section>

        <section class="narrative" data-contract-section="further-questions">
          <h2>Decision-relevant questions left open</h2>
          <ul class="method-list">
            <li>Which observable marks a historical surface observation as spectrally representative of the target date without using an aerosol truth product?</li>
            <li>Can orbit, view geometry, phenology, snow state, and land-cover matching reduce the remaining target-scene B02/B04 prior error below 0.010?</li>
            <li>After that prior improvement, does full residual covariance make B01/B03 useful without repeating the known global-B01 regression?</li>
            <li>Is the 80% target identifiable from Sentinel-2 history alone, or does the hard smoke/dust residual require an independent aerosol or surface constraint?</li>
          </ul>
        </section>

        <section class="caveat">
          <strong>Interpretation.</strong> The current best complete method remains {tooltip.wrap(f"{100 * r2_hits / 250:.1f}%", r2_dataset)} raw. The target requires forty additional hits, while the oracle over all already-saved surface arms reaches only {tooltip.wrap(f"{surface_oracle_pct:.1f}%", "phaseD_results_campaign250_K2 through R2 plus fixed L2A ensembles")}. The tested changes do not close that information gap, so stopping before a full-campaign rerun is the statistically and computationally correct result.
        </section>
      </article>
    </main>
    """

    source = SHELL.read_text(encoding="utf-8")
    source = re.sub(
        r"\s*<!-- Deliver this HTML file itself\..*?-->\s*",
        "\n",
        source,
        flags=re.DOTALL,
    )
    source = source.replace("{{TITLE}}", title)
    source = source.replace("{{SOURCE_AND_DATE}}", f"Saved SIAC campaign artifacts | {generated}")
    source = source.replace("{{REPORT_AUDIENCE}}", "technical")
    source = re.sub(
        r'<main data-report-audience="[^"]+">.*?</main>',
        main,
        source,
        flags=re.DOTALL,
    )
    extra_css = """
    .method-list { margin: 10px 0 0; padding-left: 22px; color: var(--secondary); }
    .method-list li { margin: 0 0 12px; }
    .experiment { padding-top: 8px; border-top: 1px solid var(--border-strong); }
    td:nth-child(3) { white-space: normal; text-align: left; min-width: 310px; }
    .experiment-table td:nth-child(3) { min-width: 150px; }
    .experiment-table td:nth-child(4) { white-space: normal; text-align: left; min-width: 300px; }
    @media (max-width: 800px) { td:nth-child(3) { min-width: 260px; } }
    """
    source = source.replace("  </style>", extra_css + "  </style>")
    shell_path = OUT / "report-shell.html"
    shell_path.write_text(source, encoding="utf-8")

    payload = {
        "charts": [
            _payload_chart(
                chart_id="cohort-score",
                title="Within-EE by evaluation cohort",
                data=cohort_chart,
                category="cohort",
            ),
            _payload_chart(
                chart_id="aod-bin-score",
                title="Within-EE by AERONET AOD bin",
                data=bin_chart,
                category="aod_bin",
            ),
        ]
    }
    payload_path = OUT / "report-payload.json"
    payload_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    source_notes = {
        "audience": "technical",
        "delivery_mode": "html",
        "question": "Why does the surface-driven AOD retrieval fail on campaign-250, and did the staged 6S surface-prior experiments justify a full 250-site run?",
        "required_section_map": [
            "title",
            "technical-summary",
            "key-findings",
            "scope-data-and-metric-definitions",
            "methodology",
            "experimental-design",
            "limitations-uncertainty-and-robustness-checks",
            "recommended-next-steps",
            "further-questions",
        ],
        "sources": [
            str(R2_DIR),
            str(SDSPEC_DIR / SDSPEC_PATTERN),
            str(PQ_DIR),
            str(R2_LOCAL_DIR),
            str(R2_COST_DIR),
            str(FINDINGS),
            str(ROOT / "campaign250_mids.txt"),
            str(ROOT / "all62_mids.txt"),
            str(ROOT / "l1c_test" / "seasonal" / "*" / "*_meta.json"),
            str(EXPERIMENT_SUMMARY),
            str(PRIOR_QUALITY_FINAL),
            str(BAND_CONSENSUS_FINAL),
            str(AGGREGATION_CONSENSUS),
        ],
        "chart_map": [
            {
                "section": "cohort generalization",
                "question": "Does the selected all-62 result generalize?",
                "family": "comparison",
                "type": "grouped bar",
                "fields": ["cohort", "series", "within_ee"],
                "claim": "Both R2 and exact sdspec collapse on the campaign remainder.",
                "palette": {"R2": "#0169cc", "acixThree sdspec": "#e25507"},
            },
            {
                "section": "AOD regime",
                "question": "Where does performance fail along the AOD axis?",
                "family": "comparison",
                "type": "grouped bar",
                "fields": ["aod_bin", "series", "within_ee"],
                "claim": "Performance falls sharply from moderate AOD onward.",
                "palette": {"R2": "#0169cc", "acixThree sdspec": "#e25507"},
            },
        ],
        "omissions": {
            "full_campaign_E5": "No staged candidate passed its paired promotion gate, so a new 250-site retrieval was intentionally not run.",
            "staged_score_chart": "Experiment cohorts have different small denominators; an exact table with paired gains/losses is less misleading than a percentage-only chart.",
            "causal_claim": "All driver results are diagnostic associations or controlled ablations; no causal claim is made from correlation alone.",
        },
        "qa": {
            "structural": "Required technical sections, chart hosts, static fallbacks, source tooltip ids, placeholders, and external dependencies checked.",
            "static_visual": "Both same-data SVG fallbacks were rasterized and inspected at 960x420.",
            "live_visual": "Packaged Recharts runtime embedded; no browser executable was available in this environment for live-mount visual inspection.",
        },
    }
    (OUT / "source-notes.json").write_text(json.dumps(source_notes, indent=2), encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(EMBED),
            "--input",
            str(shell_path),
            "--payload",
            str(payload_path),
            "--output",
            str(OUT / "report.html"),
        ],
        check=True,
    )
    print(OUT / "report.html")


if __name__ == "__main__":
    build()

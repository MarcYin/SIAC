"""Build the paired B01+B02+B03+B04 versus B02+B03+B04 AOD report."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
COHORT = ROOT / "campaign250_lowcloud20_mids.txt"
BASELINE_DIR = ROOT / "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711"
FOUR_BAND_DIR = ROOT / "phaseD_results_lowcloud20_singleprior_b01b03_chi2_full_20260712"
DEFAULT_OUTPUT = ROOT / "reports/aod-b01-band-comparison-20260712"

SOURCE_ID = "b01-band-comparison"
TITLE = "Adding Sentinel-2 B01: paired gains and losses"


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        matchup_id = str(row.get("matchup_id") or path.stem)
        if matchup_id in records:
            raise ValueError(f"Duplicate matchup ID {matchup_id!r} in {directory}")
        records[matchup_id] = row
    return records


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _metric(record: dict[str, Any], key: str) -> float | None:
    return _finite(record.get("solver", {}).get(key))


def _nested(record: dict[str, Any], container: str, key: str) -> float | None:
    return _finite(record.get(container, {}).get(key))


def _truth_bin(value: float) -> str:
    if value < 0.1:
        return "<0.1"
    if value < 0.2:
        return "0.1-0.2"
    if value < 0.4:
        return "0.2-0.4"
    if value < 0.6:
        return "0.4-0.6"
    if value < 1.0:
        return "0.6-1.0"
    return ">=1.0"


def _transition(baseline_hit: bool, four_band_hit: bool) -> tuple[str, str]:
    if not baseline_hit and four_band_hit:
        return "gain", "Improved to within EE"
    if baseline_hit and not four_band_hit:
        return "loss", "Regressed outside EE"
    if baseline_hit:
        return "retained_hit", "Within EE in both"
    return "retained_miss", "Outside EE in both"


def _fmt_minima(record: dict[str, Any], bands: tuple[str, ...]) -> str:
    parts = []
    for band in bands:
        value = _metric(record, f"surface_band_{band}_argmin_aot")
        parts.append(f"{band}={value:.3f}" if value is not None else f"{band}=NA")
    return "; ".join(parts)


def _case_row(
    matchup_id: str,
    baseline: dict[str, Any],
    four_band: dict[str, Any],
) -> dict[str, Any]:
    if baseline.get("status") != "OK" or four_band.get("status") != "OK":
        raise ValueError(f"Both methods must be status OK for {matchup_id}")
    truth = float(baseline["truth"])
    if not math.isclose(truth, float(four_band["truth"]), abs_tol=1e-12):
        raise ValueError(f"Truth mismatch for {matchup_id}")
    baseline_aod = float(baseline["retrieved"])
    four_band_aod = float(four_band["retrieved"])
    threshold = 0.05 + 0.15 * truth
    baseline_error = baseline_aod - truth
    four_band_error = four_band_aod - truth
    baseline_abs_error = abs(baseline_error)
    four_band_abs_error = abs(four_band_error)
    baseline_hit = baseline_abs_error <= threshold + 1e-12
    four_band_hit = four_band_abs_error <= threshold + 1e-12
    if bool(baseline.get("within_ee")) != baseline_hit:
        raise ValueError(f"Baseline EE flag mismatch for {matchup_id}")
    if bool(four_band.get("within_ee")) != four_band_hit:
        raise ValueError(f"Four-band EE flag mismatch for {matchup_id}")
    transition, transition_label = _transition(baseline_hit, four_band_hit)

    baseline_pass1 = _nested(baseline, "anchor_iterate", "pass1_scene_mean")
    four_band_pass1 = _nested(four_band, "anchor_iterate", "pass1_scene_mean")
    b01_min = _metric(four_band, "surface_band_B01_argmin_aot")
    non_b01_minima = [
        value
        for band in ("B02", "B03", "B04")
        if (value := _metric(four_band, f"surface_band_{band}_argmin_aot")) is not None
    ]
    visible_median = statistics.median(non_b01_minima) if non_b01_minima else None
    b01_pull = (
        b01_min - visible_median if b01_min is not None and visible_median is not None else None
    )
    baseline_support = _metric(baseline, "surface_solved_pixel_fraction")
    four_band_support = _metric(four_band, "surface_solved_pixel_fraction")

    return {
        "matchup_id": matchup_id,
        "site": str(baseline.get("site") or matchup_id.split("__", 1)[0]),
        "product_id": str(baseline.get("product_id") or ""),
        "transition": transition,
        "transition_label": transition_label,
        "truth_bin": _truth_bin(truth),
        "regime": str(baseline.get("regime") or ""),
        "truth_aod": truth,
        "ee_threshold": threshold,
        "cloud_fraction": _finite(baseline.get("cloud_frac")),
        "baseline_aod": baseline_aod,
        "baseline_error": baseline_error,
        "baseline_abs_error": baseline_abs_error,
        "baseline_norm_error": baseline_abs_error / threshold,
        "baseline_within_ee": baseline_hit,
        "four_band_aod": four_band_aod,
        "four_band_error": four_band_error,
        "four_band_abs_error": four_band_abs_error,
        "four_band_norm_error": four_band_abs_error / threshold,
        "four_band_within_ee": four_band_hit,
        "aod_shift": four_band_aod - baseline_aod,
        "abs_error_change": four_band_abs_error - baseline_abs_error,
        "norm_error_change": (four_band_abs_error - baseline_abs_error) / threshold,
        "baseline_pass1_aod": baseline_pass1,
        "four_band_pass1_aod": four_band_pass1,
        "pass1_aod_shift": (
            four_band_pass1 - baseline_pass1
            if four_band_pass1 is not None and baseline_pass1 is not None
            else None
        ),
        "baseline_curve_min_aod": _metric(baseline, "surface_cost_curve_min_aot"),
        "four_band_curve_min_aod": _metric(four_band, "surface_cost_curve_min_aot"),
        "b01_argmin_aod": b01_min,
        "four_band_non_b01_median_argmin": visible_median,
        "b01_argmin_pull": b01_pull,
        "baseline_band_minima": _fmt_minima(baseline, ("B02", "B03", "B04")),
        "four_band_band_minima": _fmt_minima(four_band, ("B01", "B02", "B03", "B04")),
        "baseline_band_spread": _metric(baseline, "surface_band_argmin_spread"),
        "four_band_band_spread": _metric(four_band, "surface_band_argmin_spread"),
        "baseline_cost_per_band": _metric(baseline, "cost_final_per_band"),
        "four_band_cost_per_band": _metric(four_band, "cost_final_per_band"),
        "baseline_support": baseline_support,
        "four_band_support": four_band_support,
        "support_change": (
            four_band_support - baseline_support
            if four_band_support is not None and baseline_support is not None
            else None
        ),
        "four_band_b01_prior_boa": _nested(four_band, "prior_boa", "B01"),
        "four_band_b01_prior_unc": _nested(four_band, "prior_unc", "B01"),
    }


def _method_summary(rows: list[dict[str, Any]], prefix: str) -> dict[str, Any]:
    errors = [float(row[f"{prefix}_error"]) for row in rows]
    absolute = [abs(value) for value in errors]
    hits = sum(bool(row[f"{prefix}_within_ee"]) for row in rows)
    return {
        "hits": hits,
        "rate": hits / len(rows),
        "rmse": math.sqrt(statistics.fmean(value * value for value in errors)),
        "mae": statistics.fmean(absolute),
        "bias": statistics.fmean(errors),
        "median_abs_error": statistics.median(absolute),
    }


def _mcnemar_exact(gains: int, losses: int) -> float:
    discordant = gains + losses
    if discordant == 0:
        return 1.0
    lower = min(gains, losses)
    probability = (
        2.0 * sum(math.comb(discordant, index) for index in range(lower + 1)) / (2**discordant)
    )
    return min(1.0, probability)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty dataset to {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _source(generated_at: str) -> dict[str, Any]:
    return {
        "id": SOURCE_ID,
        "label": "Paired 152-case Sentinel-2 B01 band experiment",
        "path": "all-cases.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('all-cases.csv', header = true);",
            "description": (
                "Load the bounded paired rows recomputed from the terminal three-band and "
                "four-band experiment JSON records."
            ),
            "executed_at": generated_at,
            "tables_used": ["all-cases.csv"],
            "filters": [
                "Fixed cohort contains 152 unique campaign matchups with cloud fraction below 0.20",
                "Both experiment arms require terminal status OK and finite truth and retrieved AOD",
                "Baseline predicts and solves B02, B03, B04; candidate predicts and solves B01, B02, B03, B04",
                "Surface prior, ExtraTrees model, chi-square cost, iteration, masks, backstop, LUT RT, and scene-mean extraction are held fixed",
            ],
            "metric_definitions": [
                "Within expected error: abs(retrieved AOD - AERONET AOD) <= 0.05 + 0.15 * AERONET AOD",
                "Normalized EE error: absolute AOD error divided by the matchup-specific expected-error threshold",
                "Gain: baseline outside EE and four-band candidate within EE",
                "Loss: baseline within EE and four-band candidate outside EE",
                "Absolute-error change: four-band absolute error minus three-band absolute error; negative values improve",
                "B01 argmin pull: four-band B01-only AOD minimum minus the median B02/B03/B04 AOD minimum",
                "McNemar exact p-value: two-sided exact binomial test over paired gains and losses",
            ],
        },
    }


def _cards() -> list[dict[str, Any]]:
    return [
        {
            "id": "four-band-rate-card",
            "dataset": "headline",
            "sourceId": SOURCE_ID,
            "description": "Strict within-EE rate on the fixed 152-matchup denominator.",
            "metrics": [
                {"label": "Four-band within EE", "field": "four_band_rate", "format": "percent"},
                {"label": "Three-band control", "field": "baseline_rate", "format": "percent"},
                {
                    "label": "Rate change",
                    "field": "rate_change",
                    "format": "percent",
                    "signed": True,
                },
            ],
        },
        {
            "id": "paired-gains-card",
            "dataset": "headline",
            "sourceId": SOURCE_ID,
            "description": "Cases crossing into EE after B01 is added.",
            "metrics": [
                {"label": "Paired gains", "field": "gains", "format": "number"},
                {"label": "Paired losses", "field": "losses", "format": "number"},
                {
                    "label": "Net hits",
                    "field": "net_hits",
                    "format": "number",
                    "signed": True,
                },
            ],
        },
        {
            "id": "rmse-card",
            "dataset": "headline",
            "sourceId": SOURCE_ID,
            "description": "RMSE is sensitive to the largest high-AOD errors.",
            "metrics": [
                {"label": "Four-band RMSE", "field": "four_band_rmse", "format": "number"},
                {"label": "Three-band RMSE", "field": "baseline_rmse", "format": "number"},
                {
                    "label": "RMSE change",
                    "field": "rmse_change",
                    "format": "number",
                    "signed": True,
                },
            ],
        },
        {
            "id": "median-error-card",
            "dataset": "headline",
            "sourceId": SOURCE_ID,
            "description": "Median absolute error describes the typical matchup.",
            "metrics": [
                {
                    "label": "Four-band median |error|",
                    "field": "four_band_median_abs_error",
                    "format": "number",
                },
                {
                    "label": "Three-band median |error|",
                    "field": "baseline_median_abs_error",
                    "format": "number",
                },
                {
                    "label": "Median change",
                    "field": "median_abs_error_change",
                    "format": "number",
                    "signed": True,
                },
            ],
        },
    ]


def _charts() -> list[dict[str, Any]]:
    return [
        {
            "id": "method-score-chart",
            "title": "Strict within-EE rate by solve-band set",
            "subtitle": "Same 152 low-cloud matchups and retrieval recipe; only B01 inclusion changes",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "How does global B01 inclusion change the primary score?",
            "rationale": "Two fixed-denominator bars show the aggregate score without mixing in partial screens.",
            "comparisonContext": {
                "baseline": "B02+B03+B04 chi-square",
                "denominator": "152 low-cloud matchups",
                "grain": "solve-band set",
                "unit": "strict within-EE fraction",
            },
            "dataset": "method_scores",
            "sourceId": SOURCE_ID,
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Solve bands"},
                "y": {
                    "field": "strict_rate",
                    "type": "quantitative",
                    "label": "Strict within EE",
                },
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "rmse", "type": "quantitative"},
                    {"field": "mae", "type": "quantitative"},
                    {"field": "bias", "type": "quantitative"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "settings": {"showValues": True, "categoryLabelPolicy": "wrap"},
        },
        {
            "id": "paired-error-chart",
            "title": "EE-normalized error before and after adding B01",
            "subtitle": "Each point is one matchup; reference lines at 1.0 mark the expected-error boundary",
            "type": "scatter",
            "intent": "relationship",
            "question": "Which cases cross the EE boundary when B01 is added?",
            "rationale": "A paired-error scatter exposes gains, losses, retained hits, and retained misses on one common scale.",
            "comparisonContext": {
                "baseline": "Three-band normalized error",
                "denominator": "152 low-cloud matchups",
                "grain": "matchup",
                "normalization": "absolute error divided by matchup-specific EE threshold",
                "unit": "EE-normalized absolute error",
            },
            "dataset": "comparison_cases",
            "sourceId": SOURCE_ID,
            "encodings": {
                "x": {
                    "field": "baseline_norm_error",
                    "type": "quantitative",
                    "label": "Three-band |error| / EE",
                },
                "y": {
                    "field": "four_band_norm_error",
                    "type": "quantitative",
                    "label": "Four-band |error| / EE",
                },
                "color": {
                    "field": "transition_label",
                    "type": "nominal",
                    "label": "Transition",
                },
                "tooltip": [
                    {"field": "site", "type": "nominal"},
                    {"field": "truth_aod", "type": "quantitative"},
                    {"field": "baseline_aod", "type": "quantitative"},
                    {"field": "four_band_aod", "type": "quantitative"},
                    {"field": "cloud_fraction", "type": "quantitative", "format": "percent"},
                ],
            },
            "combinationRationale": "Transition color identifies the four EE quadrants while point tooltips retain matchup identity.",
            "layout": "full",
            "palette": {"kind": "categorical"},
            "legend": {
                "interactive": True,
                "position": "bottom",
                "sort": "labelAsc",
                "title": "EE transition",
            },
            "referenceLines": [
                {"axis": "x", "value": 1.0, "label": "Three-band EE boundary", "color": "neutral"},
                {"axis": "y", "value": 1.0, "label": "Four-band EE boundary", "color": "neutral"},
            ],
        },
        {
            "id": "transition-chart",
            "title": "Paired EE boundary crossings",
            "subtitle": "Only cases whose within-EE status changes; retained hits and misses are excluded",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "Do B01 gains outnumber B01 losses?",
            "rationale": "Two transition counts directly show the asymmetry behind the net hit change.",
            "comparisonContext": {
                "baseline": "Three-band within-EE status",
                "denominator": "32 discordant matchups",
                "grain": "transition direction",
                "unit": "matchups",
            },
            "dataset": "transition_changes",
            "sourceId": SOURCE_ID,
            "encodings": {
                "x": {"field": "transition", "type": "nominal", "label": "Transition"},
                "y": {"field": "count", "type": "quantitative", "label": "Matchups"},
            },
            "valueFormat": "number",
            "unit": "matchups",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "settings": {"showValues": True, "categoryLabelPolicy": "wrap"},
        },
        {
            "id": "truth-strata-chart",
            "title": "Net within-EE hit change by AERONET AOD",
            "subtitle": "Four-band hits minus three-band hits within each fixed truth-AOD bin",
            "type": "bar",
            "intent": "decomposition",
            "question": "Where in the AOD range do the gains and losses occur?",
            "rationale": "Ordered signed bars preserve the physical AOD range and expose segment-specific net movement.",
            "comparisonContext": {
                "baseline": "Three-band hits in each truth bin",
                "denominator": "Matchups in each AERONET AOD bin",
                "grain": "truth-AOD bin",
                "unit": "net within-EE hits",
            },
            "dataset": "truth_strata",
            "sourceId": SOURCE_ID,
            "encodings": {
                "x": {"field": "truth_bin", "type": "ordinal", "label": "AERONET AOD"},
                "y": {"field": "net_hits", "type": "quantitative", "label": "Net hits"},
                "tooltip": [
                    {"field": "count", "type": "quantitative"},
                    {"field": "gains", "type": "quantitative"},
                    {"field": "losses", "type": "quantitative"},
                    {"field": "baseline_hits", "type": "quantitative"},
                    {"field": "four_band_hits", "type": "quantitative"},
                ],
            },
            "valueFormat": "number",
            "unit": "net hits",
            "layout": "full",
            "palette": {"kind": "diverging", "midpoint": 0},
            "referenceLines": [
                {"axis": "y", "value": 0, "label": "No net change", "color": "neutral"}
            ],
            "settings": {"showValues": True},
        },
    ]


def _case_columns() -> list[dict[str, Any]]:
    return [
        {"field": "site", "label": "Site", "type": "text"},
        {"field": "truth_aod", "label": "AERONET", "format": "number"},
        {"field": "cloud_fraction", "label": "Cloud", "format": "percent"},
        {"field": "baseline_aod", "label": "3-band AOD", "format": "number"},
        {"field": "baseline_norm_error", "label": "3-band |err|/EE", "format": "number"},
        {"field": "four_band_aod", "label": "4-band AOD", "format": "number"},
        {"field": "four_band_norm_error", "label": "4-band |err|/EE", "format": "number"},
        {
            "field": "abs_error_change",
            "label": "Change in |error|",
            "format": "number",
            "signed": True,
            "movement": True,
        },
        {
            "field": "aod_shift",
            "label": "4-band - 3-band",
            "format": "number",
            "signed": True,
        },
        {"field": "baseline_pass1_aod", "label": "3-band pass 1", "format": "number"},
        {"field": "four_band_pass1_aod", "label": "4-band pass 1", "format": "number"},
        {"field": "b01_argmin_aod", "label": "B01 min", "format": "number"},
        {
            "field": "b01_argmin_pull",
            "label": "B01 min - B02/03/04 median",
            "format": "number",
            "signed": True,
        },
        {"field": "baseline_band_minima", "label": "3-band minima", "type": "text"},
        {"field": "four_band_band_minima", "label": "4-band minima", "type": "text"},
        {"field": "baseline_band_spread", "label": "3-band spread", "format": "number"},
        {"field": "four_band_band_spread", "label": "4-band spread", "format": "number"},
        {"field": "four_band_b01_prior_boa", "label": "B01 surface prior", "format": "number"},
        {"field": "four_band_b01_prior_unc", "label": "B01 prior unc", "format": "number"},
        {"field": "baseline_support", "label": "3-band support", "format": "percent"},
        {"field": "four_band_support", "label": "4-band support", "format": "percent"},
    ]


def _tables() -> list[dict[str, Any]]:
    return [
        {
            "id": "gain-table",
            "title": "Cases improved to within expected error",
            "subtitle": "All 10 paired gains, ordered by reduction in absolute AOD error",
            "dataset": "gain_cases",
            "sourceId": SOURCE_ID,
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "abs_error_change", "direction": "asc"},
            "columns": _case_columns(),
        },
        {
            "id": "loss-table",
            "title": "Cases moved outside expected error",
            "subtitle": "All 22 paired losses, ordered by increase in absolute AOD error",
            "dataset": "loss_cases",
            "sourceId": SOURCE_ID,
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "abs_error_change", "direction": "desc"},
            "columns": _case_columns(),
        },
    ]


def build_report(output: Path) -> dict[str, Any]:
    cohort = [
        line.strip() for line in COHORT.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if len(cohort) != 152 or len(set(cohort)) != 152:
        raise ValueError("Expected exactly 152 unique low-cloud matchup IDs")
    baseline_records = _load_records(BASELINE_DIR)
    four_band_records = _load_records(FOUR_BAND_DIR)
    missing_baseline = sorted(set(cohort) - set(baseline_records))
    missing_four_band = sorted(set(cohort) - set(four_band_records))
    if missing_baseline or missing_four_band:
        raise ValueError(
            f"Incomplete inputs: baseline missing={len(missing_baseline)}, "
            f"four-band missing={len(missing_four_band)}"
        )

    rows = [
        _case_row(matchup_id, baseline_records[matchup_id], four_band_records[matchup_id])
        for matchup_id in cohort
    ]
    baseline_summary = _method_summary(rows, "baseline")
    four_band_summary = _method_summary(rows, "four_band")
    gains = [row for row in rows if row["transition"] == "gain"]
    losses = [row for row in rows if row["transition"] == "loss"]
    retained_hits = [row for row in rows if row["transition"] == "retained_hit"]
    retained_misses = [row for row in rows if row["transition"] == "retained_miss"]
    gains.sort(key=lambda row: float(row["abs_error_change"]))
    losses.sort(key=lambda row: float(row["abs_error_change"]), reverse=True)
    exact_p = _mcnemar_exact(len(gains), len(losses))

    bin_order = ("<0.1", "0.1-0.2", "0.2-0.4", "0.4-0.6", "0.6-1.0", ">=1.0")
    truth_strata = []
    for truth_bin in bin_order:
        subset = [row for row in rows if row["truth_bin"] == truth_bin]
        baseline_hits = sum(bool(row["baseline_within_ee"]) for row in subset)
        four_band_hits = sum(bool(row["four_band_within_ee"]) for row in subset)
        truth_strata.append(
            {
                "truth_bin": truth_bin,
                "count": len(subset),
                "baseline_hits": baseline_hits,
                "four_band_hits": four_band_hits,
                "gains": sum(row["transition"] == "gain" for row in subset),
                "losses": sum(row["transition"] == "loss" for row in subset),
                "net_hits": four_band_hits - baseline_hits,
            }
        )

    generated_at = datetime.now(timezone.utc).isoformat()
    headline = {
        "cohort_count": len(rows),
        "baseline_hits": baseline_summary["hits"],
        "baseline_rate": baseline_summary["rate"],
        "baseline_rmse": baseline_summary["rmse"],
        "baseline_mae": baseline_summary["mae"],
        "baseline_bias": baseline_summary["bias"],
        "baseline_median_abs_error": baseline_summary["median_abs_error"],
        "four_band_hits": four_band_summary["hits"],
        "four_band_rate": four_band_summary["rate"],
        "four_band_rmse": four_band_summary["rmse"],
        "four_band_mae": four_band_summary["mae"],
        "four_band_bias": four_band_summary["bias"],
        "four_band_median_abs_error": four_band_summary["median_abs_error"],
        "rate_change": four_band_summary["rate"] - baseline_summary["rate"],
        "rmse_change": four_band_summary["rmse"] - baseline_summary["rmse"],
        "mae_change": four_band_summary["mae"] - baseline_summary["mae"],
        "median_abs_error_change": (
            four_band_summary["median_abs_error"] - baseline_summary["median_abs_error"]
        ),
        "gains": len(gains),
        "losses": len(losses),
        "net_hits": len(gains) - len(losses),
        "retained_hits": len(retained_hits),
        "retained_misses": len(retained_misses),
        "mcnemar_exact_p": exact_p,
    }
    method_scores = [
        {
            "method": "B02+B03+B04",
            "hits": baseline_summary["hits"],
            "strict_rate": baseline_summary["rate"],
            "rmse": baseline_summary["rmse"],
            "mae": baseline_summary["mae"],
            "bias": baseline_summary["bias"],
        },
        {
            "method": "B01+B02+B03+B04",
            "hits": four_band_summary["hits"],
            "strict_rate": four_band_summary["rate"],
            "rmse": four_band_summary["rmse"],
            "mae": four_band_summary["mae"],
            "bias": four_band_summary["bias"],
        },
    ]
    transition_changes = [
        {"transition": "Improved to within EE", "count": len(gains)},
        {"transition": "Regressed outside EE", "count": len(losses)},
    ]

    source = _source(generated_at)
    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
        {
            "id": "technical-summary",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Technical summary\n\n"
                "**Adding B01 globally reduces strict within-EE performance from "
                f"{baseline_summary['hits']}/152 = {100 * baseline_summary['rate']:.1f}% to "
                f"{four_band_summary['hits']}/152 = {100 * four_band_summary['rate']:.1f}%.** "
                f"The paired comparison contains {len(gains)} gains and {len(losses)} losses, "
                f"for a net change of {len(gains) - len(losses):+d} hits.\n\n"
                "**The aggregate error metrics move in different directions.** RMSE improves "
                f"from {baseline_summary['rmse']:.3f} to {four_band_summary['rmse']:.3f}, while "
                f"MAE worsens from {baseline_summary['mae']:.3f} to {four_band_summary['mae']:.3f} "
                f"and median absolute error rises from {baseline_summary['median_abs_error']:.3f} "
                f"to {four_band_summary['median_abs_error']:.3f}. This means a few large errors "
                "shrink while the typical matchup becomes less accurate.\n\n"
                f"**The two-sided exact McNemar p-value is {exact_p:.4f}.** This tests only the "
                "asymmetry of the 10/22 boundary crossings; it does not establish why B01 changes "
                "the retrieval. The tables below retain every gain and loss for inspection."
            ),
        },
        {
            "id": "headline-strip",
            "type": "metric-strip",
            "cardIds": [
                "four-band-rate-card",
                "paired-gains-card",
                "rmse-card",
                "median-error-card",
            ],
        },
        {
            "id": "aggregate-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Global B01 inclusion loses 12 net EE hits\n\n"
                "Both methods use the same 152 low-cloud matchups, S2 monthly SWIR/NIR-anchored "
                "ExtraTrees surface prior, chi-square objective, anchor iteration, atmospheric "
                "backstop, masking, LUT radiative transfer, and scene-mean extraction. The only "
                "planned change is predicting and solving B01 in addition to B02/B03/B04. The "
                "bars therefore use one fixed denominator and include every matchup."
            ),
        },
        {"id": "method-score-chart-block", "type": "chart", "chartId": "method-score-chart"},
        {
            "id": "paired-geometry-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## The paired error surface shows both recoveries and regressions\n\n"
                "Normalized error equals absolute AOD error divided by the matchup-specific EE "
                "threshold. Values at or below 1 are within EE. Points moving from right of the "
                "vertical boundary to below the horizontal boundary are gains; points moving from "
                "left of the vertical boundary to above the horizontal boundary are losses. The "
                "quadrants preserve all 152 outcomes, not only the 32 changed statuses."
            ),
        },
        {"id": "paired-error-chart-block", "type": "chart", "chartId": "paired-error-chart"},
        {
            "id": "transition-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Losses outnumber gains by more than two to one\n\n"
                f"B01 recovers {len(gains)} cases that the three-band solve misses, but moves "
                f"{len(losses)} previously successful cases outside EE. The remaining "
                f"{len(retained_hits)} cases stay within EE and {len(retained_misses)} stay outside. "
                "The chart isolates the discordant cases; the case tables provide the exact paired "
                "retrieval and spectral diagnostics behind each count."
            ),
        },
        {"id": "transition-chart-block", "type": "chart", "chartId": "transition-chart"},
        {
            "id": "truth-strata-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Net improvement is confined to the 0.4-0.6 AOD bin\n\n"
                "The four-band solve gains three net hits for AERONET AOD 0.4-0.6. Every other "
                "truth bin is flat or negative: -5 below 0.1, -4 at 0.1-0.2, -3 at 0.2-0.4, "
                "-1 at 0.6-1.0, and -2 at or above 1.0. These are descriptive cohort cuts, not "
                "an independently validated rule for switching B01 on or off."
            ),
        },
        {"id": "truth-strata-chart-block", "type": "chart", "chartId": "truth-strata-chart"},
        {
            "id": "gain-cases-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## All 10 improvements, with three-band and four-band evidence side by side\n\n"
                "Negative change in absolute error means the four-band result moved closer to "
                "AERONET. The pass-1 columns expose iteration movement; band minima show the "
                "surface-only AOD preference of each visible band; and B01 pull compares the B01 "
                "minimum with the median B02/B03/B04 minimum in the same four-band run."
            ),
        },
        {"id": "gain-table-block", "type": "table", "tableId": "gain-table"},
        {
            "id": "loss-cases-finding",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## All 22 losses, using the same diagnostics and ordering\n\n"
                "Positive change in absolute error means the four-band result moved farther from "
                "AERONET. The table is deliberately symmetric with the gains table so differences "
                "in B01 pull, iteration, spectral disagreement, support, cloud fraction, and prior "
                "uncertainty can be inspected without changing definitions between groups."
            ),
        },
        {"id": "loss-table-block", "type": "table", "tableId": "loss-table"},
        {
            "id": "scope-definitions",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Scope, data, and metric definitions\n\n"
                "The population is the fixed 152-matchup campaign subset with recorded cloud "
                "fraction below 20%. Both arms have 152 terminal `OK` outputs and use the same "
                "AERONET AOD reference. Expected error is `|retrieved - truth| <= 0.05 + "
                "0.15 * truth`. Strict rates divide by all 152 planned matchups. A gain or loss "
                "requires crossing this same threshold between the paired methods."
            ),
        },
        {
            "id": "methodology",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Paired experimental design\n\n"
                "The control predicts and solves B02/B03/B04. The candidate predicts and solves "
                "B01/B02/B03/B04. Both use an S2 monthly SWIR/NIR-anchored ExtraTrees prior with "
                "20 estimators, one anchor iteration, a visible reflectance uncertainty floor of "
                "0.006, calibrated atmospheric backstop, three/four-band chi-square cost, 60 m "
                "solver grid, and scene-mean AOD extraction. Metrics are recomputed from terminal "
                "JSON records rather than copied from prior summaries."
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "sourceId": SOURCE_ID,
            "body": (
                "## Limitations and robustness checks\n\n"
                "The comparison is paired and complete, but it uses one selected low-cloud "
                "campaign rather than an independent site/year holdout. RMSE and within-EE rate "
                "answer different questions: RMSE is dominated by large errors, while EE counts "
                "threshold crossings. Band argmin and B01-pull fields are diagnostics of the saved "
                "surface cost, not causal attribution. The exact McNemar result is borderline and "
                "does not adjust for earlier method exploration on this cohort."
            ),
        },
        {
            "id": "next-steps",
            "type": "markdown",
            "body": (
                "## Recommended next checks\n\n"
                "1. Compare gain and loss distributions for B01 prediction uncertainty, B01 pull, "
                "and pass-1 movement using a site-grouped development split.\n"
                "2. Test any proposed B01 inclusion rule on a separate site/year holdout before "
                "treating the 0.4-0.6 pattern as operational.\n"
                "3. Keep the global four-band arm as a negative control when evaluating selective "
                "or uncertainty-weighted B01 formulations."
            ),
        },
        {
            "id": "further-questions",
            "type": "markdown",
            "body": (
                "## Further questions\n\n"
                "- Are losses concentrated where the B01 surface prior has low stated uncertainty "
                "but disagrees with B02/B03/B04?\n"
                "- Does B01 help only when its cost minimum agrees with the other visible bands, "
                "or when it supplies an independent correction?\n"
                "- Does the same gain/loss asymmetry persist with tile-wide CCI aerosol species "
                "and on an independent low-cloud cohort?"
            ),
        },
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "description": (
                "Paired technical comparison of global Sentinel-2 B01 inclusion over the fixed "
                "152-case low-cloud AOD cohort."
            ),
            "generatedAt": generated_at,
            "blocks": blocks,
            "cards": _cards(),
            "charts": _charts(),
            "tables": _tables(),
            "sources": [source],
        },
        "snapshot": {
            "version": 1,
            "status": "ready",
            "generatedAt": generated_at,
            "datasets": {
                "headline": [headline],
                "method_scores": method_scores,
                "comparison_cases": rows,
                "transition_changes": transition_changes,
                "truth_strata": truth_strata,
                "gain_cases": gains,
                "loss_cases": losses,
            },
        },
        "sources": [source],
    }

    analysis = {
        "generated_at": generated_at,
        "cohort_count": len(rows),
        "input_directories": {
            "baseline": str(BASELINE_DIR),
            "four_band": str(FOUR_BAND_DIR),
        },
        "headline": headline,
        "truth_strata": truth_strata,
        "gain_matchup_ids": [row["matchup_id"] for row in gains],
        "loss_matchup_ids": [row["matchup_id"] for row in losses],
    }
    source_notes = (
        "# Source and chart notes\n\n"
        "Audience: technical. Delivery mode: portable HTML.\n\n"
        "Required structure mapping: title; technical summary; aggregate and paired visual "
        "findings; scope/definitions; paired design; limitations; next checks; further questions.\n\n"
        "Chart map:\n"
        "- Aggregate score: horizontal bar, solve-band set -> strict rate.\n"
        "- Paired geometry: scatter, three-band normalized error -> four-band normalized error, "
        "grouped by EE transition.\n"
        "- Boundary crossings: horizontal bar, transition -> case count.\n"
        "- Truth strata: diverging bar, AERONET AOD bin -> net hits.\n\n"
        "Omissions: no spatial imagery or raw cost curves were available for the full four-band "
        "run because cost-cube dumping was not enabled. The report retains scalar cost and band "
        "diagnostics from every terminal result instead.\n"
    )

    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "all-cases.csv", rows)
    _write_csv(output / "gain-cases.csv", gains)
    _write_csv(output / "loss-cases.csv", losses)
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    (output / "artifact.json").write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    (output / "source-notes.md").write_text(source_notes, encoding="utf-8")
    return analysis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    analysis = build_report(args.output)
    print(json.dumps(analysis["headline"], indent=2))
    print(args.output / "artifact.json")


if __name__ == "__main__":
    main()

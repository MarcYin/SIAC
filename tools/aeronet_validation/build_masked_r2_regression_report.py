"""Build the portable technical report for the masked-R2 score diagnosis."""

from __future__ import annotations

import csv
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
REPORT_DIR = REPO / "notes" / "reports" / "masked-r2-regression-20260711"
COMPARISON_PATH = REPORT_DIR / "comparison.json"
ROWS_PATH = REPORT_DIR / "comparison-rows.csv"


def _number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _pct(value: float) -> str:
    return f"{value:.1f}%"


def _pp(value: float) -> str:
    return f"{value:+.1f} pp"


def _source() -> dict[str, Any]:
    return {
        "id": "masked-r2-comparison",
        "label": "Masked R2 paired campaign comparison",
        "path": "notes/reports/masked-r2-regression-20260711/comparison.json",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": (
                "SELECT * FROM read_csv_auto("
                "'notes/reports/masked-r2-regression-20260711/comparison-rows.csv', "
                "header = true);"
            ),
            "description": "Load the paired matchup-level comparison rows used for the score decomposition and diagnostics.",
            "tables_used": ["notes/reports/masked-r2-regression-20260711/comparison-rows.csv"],
            "filters": [
                "campaign250_mids.txt defines the 250-matchup strict denominator",
                "status OK with finite truth and retrieved AOD defines a valid retrieval",
                "surface_e1_matched_mids.txt defines the 143-case mask-only control",
            ],
            "metric_definitions": [
                "Within EE: abs(retrieved AOD - AERONET AOD) <= 0.05 + 0.15 * AERONET AOD",
                "Strict score: within-EE hits divided by all 250 expected matchups; failed and no-valid records are misses",
                "Valid-only score: within-EE hits divided by finite status-OK retrievals",
            ],
        },
    }


def build() -> Path:
    comparison = json.loads(COMPARISON_PATH.read_text(encoding="utf-8"))
    with ROWS_PATH.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))

    generated_at = datetime.now(UTC).isoformat()
    old = comparison["old"]
    new = comparison["new"]
    decomposition = comparison["decomposition"]
    diagnostics = comparison["common_diagnostics"]
    break_diagnostics = comparison["break_diagnostics"]
    mask_control = comparison["mask_only_control"]
    abstention_diagnostics = comparison["new_abstention_diagnostics"]

    summary_rows = [
        {
            "new_strict": _pct(new["strict_pct"]),
            "old_strict": _pct(old["strict_pct"]),
            "strict_delta": _pp(new["strict_pct"] - old["strict_pct"]),
            "new_valid": _pct(new["valid_pct"]),
            "old_valid": _pct(old["valid_pct"]),
            "valid_delta": _pp(new["valid_pct"] - old["valid_pct"]),
            "new_coverage": _pct(new["coverage_pct"]),
            "old_coverage": _pct(old["coverage_pct"]),
            "coverage_delta": _pp(new["coverage_pct"] - old["coverage_pct"]),
            "no_valid": new["statuses"].get("NO_VALID_OBSERVATION", 0),
        }
    ]
    waterfall_rows = [
        {
            "component": "Historical R2 hits",
            "hit_change": old["hits"],
            "stage": "baseline",
            "common_fixes": decomposition["fixes"],
            "common_breaks": decomposition["breaks"],
        },
        {
            "component": "Common-valid changes",
            "hit_change": decomposition["common_net"],
            "stage": "retrieval",
            "common_fixes": decomposition["fixes"],
            "common_breaks": decomposition["breaks"],
        },
        {
            "component": "Coverage abstentions",
            "hit_change": decomposition["coverage_net"],
            "stage": "coverage",
            "common_fixes": decomposition["fixes"],
            "common_breaks": decomposition["breaks"],
        },
    ]

    common_rows: list[dict[str, Any]] = []
    for row in source_rows:
        if row["old_status"] != "OK" or row["new_status"] != "OK":
            continue
        common_rows.append(
            {
                "matchup_id": row["matchup_id"],
                "site": row["site"],
                "truth_aod": _number(row["truth"]),
                "outcome": row["outcome"],
                "aod_delta": _number(row["aod_delta"]),
                "prior_b02_delta": _number(row["prior_b02_delta"]),
                "prior_b04_delta": _number(row["prior_b04_delta"]),
                "cloud_fraction": _number(row["cloud_frac"]),
                "valid_support_fraction": _number(row["valid_support_frac"]),
            }
        )

    break_rows = []
    for row in comparison["breaks"]:
        break_rows.append(
            {
                "site": row["site"],
                "truth_aod": row["truth"],
                "old_aod": row["old_aod"],
                "new_aod": row["new_aod"],
                "aod_delta": row["aod_delta"],
                "prior_b02_delta": row["prior_b02_delta"],
                "cloud_fraction": row["cloud_frac"],
                "valid_support_fraction": row["valid_support_frac"],
            }
        )

    mask_rows = [
        {
            "mode": "Proper mask",
            "hits": mask_control["masked"]["hits"],
            "valid": mask_control["masked"]["valid"],
            "strict_pct": mask_control["masked"]["strict_pct"],
            "valid_pct": mask_control["masked"]["valid_pct"],
        },
        {
            "mode": "Cloud bypass",
            "hits": mask_control["cloud_bypass"]["hits"],
            "valid": mask_control["cloud_bypass"]["valid"],
            "strict_pct": mask_control["cloud_bypass"]["strict_pct"],
            "valid_pct": mask_control["cloud_bypass"]["valid_pct"],
        },
    ]
    abstention_classes = abstention_diagnostics["classes"]
    abstention_rows = [
        {
            "mechanism": "Fully cloud masked",
            "scenes": abstention_classes["fully_cloud_masked"],
        },
        {
            "mechanism": "Other zero support",
            "scenes": abstention_classes["zero_admissible_support_other"],
        },
        {
            "mechanism": "Pooling threshold",
            "scenes": abstention_classes["valid_support_but_pool_unsolved"],
        },
    ]

    source = _source()
    title = "Why masked R2 appeared to regress"
    manifest = {
        "version": 1,
        "surface": "report",
        "title": title,
        "description": "Technical decomposition of the campaign-250 score change and its mask-only control.",
        "generatedAt": generated_at,
        "sources": [source],
        "cards": [
            {
                "id": "strict-card",
                "dataset": "summary",
                "sourceId": source["id"],
                "description": "Raw within-EE score over all 250 expected matchups.",
                "metrics": [
                    {"label": "Masked R2 strict", "field": "new_strict"},
                    {"label": "Historical", "field": "old_strict"},
                    {"label": "Delta", "field": "strict_delta", "signed": True},
                ],
            },
            {
                "id": "valid-card",
                "dataset": "summary",
                "sourceId": source["id"],
                "description": "Within-EE rate among finite status-OK retrievals.",
                "metrics": [
                    {"label": "Masked R2 valid", "field": "new_valid"},
                    {"label": "Historical", "field": "old_valid"},
                    {"label": "Delta", "field": "valid_delta", "signed": True},
                ],
            },
            {
                "id": "coverage-card",
                "dataset": "summary",
                "sourceId": source["id"],
                "description": "Share of the 250 matchups with a finite status-OK retrieval.",
                "metrics": [
                    {"label": "Masked coverage", "field": "new_coverage"},
                    {"label": "Historical", "field": "old_coverage"},
                    {"label": "Delta", "field": "coverage_delta", "signed": True},
                ],
            },
        ],
        "charts": [
            {
                "id": "score-waterfall",
                "title": "Strict-score decomposition",
                "subtitle": "Campaign-250 hit count; 161 historical hits reconcile exactly to 148 masked-run hits",
                "type": "waterfall",
                "intent": "decomposition",
                "question": "Which components account for the 13-hit strict-score decline?",
                "rationale": "The components are additive and reconcile the historical and new hit counts exactly.",
                "comparisonContext": {
                    "baseline": "Historical R2",
                    "denominator": "250 expected matchups",
                    "grain": "matchup",
                    "unit": "within-EE hits",
                },
                "dataset": "decomposition",
                "sourceId": source["id"],
                "encodings": {
                    "x": {"field": "component", "type": "nominal", "label": "Component"},
                    "y": {"field": "hit_change", "type": "quantitative", "label": "Hit change"},
                },
                "valueFormat": "number",
                "unit": "hits",
                "layout": "full",
            },
            {
                "id": "prior-scatter",
                "title": "AOD change versus B02 surface-prior change",
                "subtitle": "232 common valid matchups; Pearson r = 0.59",
                "type": "scatter",
                "intent": "relationship",
                "question": "Do retrieval changes track surface-prior changes?",
                "rationale": "A matchup-level scatter exposes direction, concentration, and outliers without aggregating away the relationship.",
                "comparisonContext": {
                    "baseline": "Historical R2",
                    "denominator": "232 common valid matchups",
                    "grain": "matchup",
                    "unit": "reflectance and AOD delta",
                },
                "dataset": "common_matchups",
                "sourceId": source["id"],
                "encodings": {
                    "x": {
                        "field": "prior_b02_delta",
                        "type": "quantitative",
                        "label": "B02 prior delta",
                    },
                    "y": {"field": "aod_delta", "type": "quantitative", "label": "AOD delta"},
                    "label": {"field": "site", "type": "text"},
                    "tooltip": [
                        {"field": "site", "type": "text"},
                        {"field": "truth_aod", "type": "quantitative"},
                        {"field": "cloud_fraction", "type": "quantitative", "format": "percent"},
                        {
                            "field": "valid_support_fraction",
                            "type": "quantitative",
                            "format": "percent",
                        },
                    ],
                },
                "valueFormat": "number",
                "layout": "full",
                "referenceLines": [
                    {"axis": "x", "value": 0, "label": "No prior change", "color": "neutral"},
                    {"axis": "y", "value": 0, "label": "No AOD change", "color": "neutral"},
                ],
            },
            {
                "id": "mask-control",
                "title": "Proper mask and cloud-bypass control",
                "subtitle": "Failure-enriched 143-matchup development set; strict within-EE hits",
                "type": "bar",
                "intent": "comparison",
                "question": "Does cloud masking itself reduce strict within-EE hits?",
                "rationale": "Two directly comparable modes are best shown as a simple absolute bar comparison.",
                "comparisonContext": {
                    "baseline": "Cloud bypass",
                    "denominator": "143 expected matchups",
                    "grain": "mask mode",
                    "unit": "within-EE hits",
                },
                "dataset": "mask_control",
                "sourceId": source["id"],
                "encodings": {
                    "x": {"field": "mode", "type": "nominal", "label": "Mode"},
                    "y": {"field": "hits", "type": "quantitative", "label": "Within-EE hits"},
                    "tooltip": [
                        {"field": "valid", "type": "quantitative"},
                        {"field": "valid_pct", "type": "quantitative", "unit": "%"},
                    ],
                },
                "valueFormat": "number",
                "unit": "hits",
                "layout": "full",
            },
            {
                "id": "abstention-breakdown",
                "title": "Why 17 scenes returned no finite AOD",
                "subtitle": "Eight are fully cloud masked; six retain valid samples but fail the local pooling threshold",
                "type": "bar",
                "intent": "composition",
                "question": "Do all no-valid statuses mean the mask removed every observation?",
                "rationale": "The mutually exclusive mechanisms separate expected masking from a downstream solve failure.",
                "comparisonContext": {
                    "baseline": "Masked R2",
                    "denominator": "17 no-finite-AOD scenes",
                    "grain": "failure mechanism",
                    "unit": "scenes",
                },
                "dataset": "abstentions",
                "sourceId": source["id"],
                "encodings": {
                    "x": {"field": "mechanism", "type": "nominal", "label": "Mechanism"},
                    "y": {"field": "scenes", "type": "quantitative", "label": "Scenes"},
                },
                "valueFormat": "number",
                "unit": "scenes",
                "layout": "full",
            },
        ],
        "tables": [
            {
                "id": "break-table",
                "title": "Sixteen common-valid regressions",
                "subtitle": "Cases that moved from within EE historically to outside EE in the new run",
                "dataset": "breaks",
                "sourceId": source["id"],
                "layout": "full",
                "density": "dense",
                "defaultSort": {"field": "aod_delta", "direction": "asc"},
                "columns": [
                    {"field": "site", "label": "Site", "type": "text"},
                    {"field": "truth_aod", "label": "Truth AOD", "format": "number"},
                    {"field": "old_aod", "label": "Old AOD", "format": "number"},
                    {"field": "new_aod", "label": "New AOD", "format": "number"},
                    {
                        "field": "aod_delta",
                        "label": "AOD delta",
                        "format": "number",
                        "movement": True,
                    },
                    {
                        "field": "prior_b02_delta",
                        "label": "B02-prior delta",
                        "format": "number",
                        "movement": True,
                    },
                    {"field": "cloud_fraction", "label": "Cloud fraction", "format": "percent"},
                    {
                        "field": "valid_support_fraction",
                        "label": "Valid support",
                        "format": "percent",
                    },
                ],
            }
        ],
        "blocks": [
            {"id": "title", "type": "markdown", "body": f"# {title}"},
            {
                "id": "technical-summary",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Technical summary\n\n"
                    "**The 5.2-point headline decline is not a cloud-mask-only result.** The 64.4% historical and 59.2% new strict scores compare different retrieval systems: the new run also changed libRadtran LUT to native 6S, changed atmospheric-provider TCWV to same-scene L2A WVP, enforced finite full-band support, and stopped filling invalid output pixels.\n\n"
                    "Of the 13 lost hits, 6 are newly abstained scenes and 7 are the net of 9 fixes and 16 breaks among 232 common valid matchups. On that common set, the score changes from 155/232 (66.8%) to 148/232 (63.8%). In the controlled 143-case mask comparison, proper masking scores 51 hits versus 50 for cloud bypass."
                ),
            },
            {
                "id": "headline-metrics",
                "type": "metric-strip",
                "cardIds": ["strict-card", "valid-card", "coverage-card"],
            },
            {
                "id": "decomposition-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Coverage explains six lost hits; common-case changes explain seven\n\n"
                    "The strict denominator treats no-finite-AOD results as misses. Seventeen scenes now abstain, and six of those had happened to fall within EE in the historical all-pixel solve. The remaining seven-hit net loss comes from common valid cases. This distinction is why the strict score falls 5.2 points while valid-only accuracy falls only 0.9 points."
                ),
            },
            {"id": "decomposition-chart", "type": "chart", "chartId": "score-waterfall"},
            {
                "id": "abstention-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## No-valid status currently mixes masking with a pooling failure\n\n"
                    "Eight of the 17 scenes have 100% raw cloud coverage and correctly abstain. Three partially cloudy coastal or island scenes have zero admissible cloud/water/finite-prior support. The other six retain 34-169 valid samples and all 68 AOD cost nodes, but no local 20x20-cell window reaches the required 80 samples. Their AOD rasters therefore remain all-NaN. Those six should be diagnosed as `VALID_SUPPORT_BUT_UNSOLVED`, not `NO_VALID_OBSERVATION`; a sparse-support fallback must be tested rather than weakening the cloud mask."
                ),
            },
            {"id": "abstention-chart", "type": "chart", "chartId": "abstention-breakdown"},
            {
                "id": "mask-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## The controlled mask test slightly favors proper masking\n\n"
                    "Holding native 6S, L2A WVP, the surface prior, and solver settings fixed, cloud bypass produces 3 fixes but 7 breaks on common valid cases. It recovers 3 additional hits from scenes that proper masking abstains on, leaving proper masking ahead by one strict hit overall. Reverting the cloud mask would therefore reduce staged performance and restore low-confidence cloudy retrievals."
                ),
            },
            {"id": "mask-chart", "type": "chart", "chartId": "mask-control"},
            {
                "id": "prior-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Surface-prior movement, not cloud fraction, tracks the AOD movement\n\n"
                    f"Across common valid matchups, AOD change correlates with B02 and B04 prior change at r={diagnostics['corr_aod_delta_prior_b02_delta']:.2f} and r={diagnostics['corr_aod_delta_prior_b04_delta']:.2f}. The corresponding correlations with cloud fraction and valid-support fraction are {diagnostics['corr_aod_delta_cloud_frac']:.2f} and {diagnostics['corr_aod_delta_valid_support_frac']:.2f}. The 16 breaks have a median AOD shift of {break_diagnostics['median_aod_delta']:.3f} and median B02-prior shift of {break_diagnostics['median_prior_b02_delta']:.4f}. These are diagnostic associations; the submitted factorial is required to attribute the prior shift between RT backend and TCWV source."
                ),
            },
            {"id": "prior-chart", "type": "chart", "chartId": "prior-scatter"},
            {
                "id": "definitions",
                "type": "markdown",
                "body": (
                    "## Scope and metric definitions\n\n"
                    "The population is the fixed campaign of 250 Sentinel-2/AERONET matchups. A hit satisfies `|retrieved - truth| <= 0.05 + 0.15 * truth`. Strict accuracy divides hits by all 250 expected matchups. Valid-only accuracy divides by finite `status=OK` retrievals. The mask-only control uses the 143-case development set containing all 89 historical R2 failures plus 54 AOD-stratified controls."
                ),
            },
            {
                "id": "methodology",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Paired decomposition methodology\n\n"
                    "Records were joined by matchup identifier. Common valid records were classified as fixes, breaks, retained hits, or retained misses. Coverage loss counts historical hits whose new status is no-valid or failed. The two components reconcile exactly to the total hit delta. Prior and AOD deltas are new minus historical values sampled and aggregated by the campaign harness."
                ),
            },
            {"id": "breaks-table", "type": "table", "tableId": "break-table"},
            {
                "id": "limitations",
                "type": "markdown",
                "body": (
                    "## Limitations and robustness checks\n\n"
                    "The historical and new full runs are not a mask-only experiment. They differ in RT backend, TCWV source, invalid-support behavior, and some solver-grid extents after intervening code changes. Correlation with prior movement does not identify which change caused that movement. The 143-case mask control is failure-enriched and should be interpreted as an ablation, not a population estimate. Historical `cloud_frac` recorded a final mask rather than the raw M1 cloud mask, so cloud comparisons use the new run only."
                ),
            },
            {
                "id": "next-steps",
                "type": "markdown",
                "body": (
                    "## Recommended next steps\n\n"
                    "1. Keep proper masking as the default; do not restore all-pixel solving.\n"
                    "2. Use clean LUT factorial job `38944577` and node-isolated 6S job `38944641` to isolate RT, water-vapour, cloud, and water-mask effects on 66 affected and matched-control scenes.\n"
                    "3. Validate the sparse-valid-support fallback with saved-cube campaign replay job `38944665` while leaving fully cloudy scenes abstained.\n"
                    "4. Promote only a configuration whose paired gains survive the same strict coverage accounting.\n"
                    "5. Report strict score, valid-only score, and retrieval coverage together for every future experiment."
                ),
            },
            {
                "id": "further-questions",
                "type": "markdown",
                "body": (
                    "## Further questions\n\n"
                    "- Does native 6S or the L2A WVP lock cause the visible-prior shift?\n"
                    "- Why do several low-cloud scenes have limited full-band surface-prior support?\n"
                    "- Can a quality-aware fallback recover useful high-AOD scenes without admitting fully cloudy London-like observations?"
                ),
            },
        ],
    }

    artifact = {
        "surface": "report",
        "manifest": manifest,
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "summary": summary_rows,
                "decomposition": waterfall_rows,
                "common_matchups": common_rows,
                "mask_control": mask_rows,
                "abstentions": abstention_rows,
                "breaks": break_rows,
            },
        },
        "sources": [source],
    }
    artifact_path = REPORT_DIR / "artifact.json"
    artifact_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    source_notes = {
        "audience": "technical",
        "required_structure": [
            "title",
            "technical summary",
            "key findings with visual evidence",
            "scope, data, and metric definitions",
            "methodology",
            "limitations, uncertainty, and robustness checks",
            "recommended next steps",
            "further questions",
        ],
        "chart_map": [
            {
                "section": "Score decomposition",
                "question": "Which additive components explain the hit decline?",
                "family": "decomposition",
                "type": "waterfall",
                "fields": ["component", "hit_change"],
                "claim": "Seven hits are common-case changes and six are coverage loss.",
                "palette": "hard two-root cap with neutral baseline",
            },
            {
                "section": "Mask control",
                "question": "Does masking itself reduce hits?",
                "family": "comparison",
                "type": "bar",
                "fields": ["mode", "hits"],
                "claim": "Proper masking is ahead by one strict hit.",
                "palette": "single-root blue",
            },
            {
                "section": "Abstention mechanisms",
                "question": "Does every no-valid result mean zero usable input?",
                "family": "composition",
                "type": "bar",
                "fields": ["mechanism", "scenes"],
                "claim": "Six of 17 abstentions retain valid samples but fail pooling.",
                "palette": "single-root blue",
            },
            {
                "section": "Prior relationship",
                "question": "Do AOD changes track surface-prior changes?",
                "family": "relationship",
                "type": "scatter",
                "fields": ["prior_b02_delta", "aod_delta"],
                "claim": "AOD and B02-prior deltas correlate at r=0.59.",
                "palette": "single-root blue with neutral zero references",
            },
        ],
        "omissions": [],
        "factorial_jobs": ["38944577", "38944641"],
        "pooling_replay_job": "38944665",
    }
    (REPORT_DIR / "source-notes.json").write_text(
        json.dumps(source_notes, indent=2), encoding="utf-8"
    )
    return artifact_path


if __name__ == "__main__":
    print(build())

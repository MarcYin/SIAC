"""Build a rounded-regime multi-source AOD selector over saved outputs.

This is a lightweight saved-result generator: it reads already collected SIAC,
MAIAC, S2 L2A AOT, CAMS/MERRA, VIIRS, and aerosol-context JSON records, applies
a fixed shallow candidate-selection tree, and writes Phase-D-like validation
records. It does not run retrieval or contact remote data services.

The default v2 rule uses rounded, regime-scale thresholds instead of the exact
validation-fitted split values used by the first diagnostic tree.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("SIAC_REFACTOR_ROOT", "/gws/ssde/j25a/nceo_isp/public/siac_refactor"))
MATCHUPS = ROOT / "matchups" / "matchups.csv"
OUT = ROOT / os.environ.get(
    "MULTISOURCE_SELECTOR_OUTDIR",
    "phaseD_results_campaign250_multisource_tree_v2",
)

Q1_LOW_AOD_MAX = 0.40
R2_CLEAN_AOD_MAX = 0.20
MERRA_LOW_AOD_MAX = 0.15
VIIRS_DUST_LOW_AOD_MAX = 0.10
R2_BELOW_MODEL_MARGIN = -0.04
CAMS_SULPHATE_FRACTION_MAX = 0.20
R2_SPATIAL_STD_MAX = 0.30
VIIRS_DUST_TO_MODEL_MAX = 1.85

DIRS = {
    "R2": "phaseD_results_campaign250_R2_full",
    "Q1": "phaseD_results_campaign250_Q1_static_ext",
    "L2A_monthly": "phaseD_results_campaign250_l2a_monthly_median3_scene_mean",
    "L2A_pc": "phaseD_results_campaign250_l2a_pc_production_mean3_scene_mean",
    "S2": "s2_l2a_aot_campaign250",
    "CAMS": "gee_cams_aod_campaign250",
    "MERRA": "gee_merra_aod_campaign250",
    "GEE_mean": "gee_model_mean_aod_campaign250",
    "VIIRS0": "viirs_aod550_qc0_campaign250",
    "VIIRS2": "viirs_aod550_qc2_campaign250",
    "VDUST": "viirs_dust_qc2_campaign250",
    "VGEN": "viirs_generic_qc2_campaign250",
    "VURB": "viirs_urban_qc2_campaign250",
    "VSMOKE": "viirs_smoke_qc2_campaign250",
}

COMPACT = (
    "R2",
    "Q1",
    "L2A_monthly",
    "L2A_pc",
    "MAIAC",
    "S2",
    "GEE_mean",
    "MERRA",
    "VIIRS0",
    "VDUST",
    "VSMOKE",
)
EXT = (
    "MAIAC",
    "S2",
    "CAMS",
    "MERRA",
    "GEE_mean",
    "VIIRS0",
    "VIIRS2",
    "VDUST",
    "VGEN",
    "VURB",
    "VSMOKE",
)
NOMODEL = ("R2", "Q1", "L2A_monthly", "L2A_pc", "MAIAC", "S2")


def _safe_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _clip_aod(value: float) -> float:
    return min(5.0, max(0.0, value))


def _within_ee(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _linear_quantile(values: list[float], q: float) -> float | None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return None
    if len(clean) == 1:
        return clean[0]
    pos = q * (len(clean) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return clean[lo]
    frac = pos - lo
    return clean[lo] * (1.0 - frac) + clean[hi] * frac


def _leq(value: float | None, threshold: float) -> bool:
    return value is not None and math.isfinite(value) and value <= threshold


def _read_result_value(root: Path, dirname: str, matchup_id: str) -> tuple[float | None, dict[str, object]]:
    path = root / dirname / f"{matchup_id}.json"
    if not path.exists():
        return None, {}
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status", "OK") != "OK":
        return None, record
    return _safe_float(record.get("retrieved")), record


def _load_inputs(root: Path, matchup_id: str) -> tuple[dict[str, float], dict[str, object]]:
    values: dict[str, float] = {}
    records: dict[str, object] = {}
    for name, dirname in DIRS.items():
        value, record = _read_result_value(root, dirname, matchup_id)
        records[name] = record
        if value is not None:
            values[name] = value

    maiac_path = root / "maiac_qa" / f"{matchup_id}.json"
    if maiac_path.exists():
        maiac = json.loads(maiac_path.read_text(encoding="utf-8"))
        records["MAIAC"] = maiac
        value = _safe_float(maiac.get("aot"))
        if value is not None:
            values["MAIAC"] = value

    context_path = root / "gee_aerosol_context_campaign250" / f"{matchup_id}.json"
    context: dict[str, object] = {}
    if context_path.exists():
        context_record = json.loads(context_path.read_text(encoding="utf-8"))
        context = context_record.get("values", {}) if isinstance(context_record.get("values"), dict) else {}
    records["context"] = context
    return values, records


def _fallback(values: dict[str, float]) -> float | None:
    for key in ("GEE_mean", "R2", "MERRA", "CAMS", "S2"):
        if key in values:
            return values[key]
    return None


def _candidate(values: dict[str, float], name: str) -> float | None:
    fallback = _fallback(values)
    if name in values:
        return values[name]
    if name == "compact_q0.6":
        return _linear_quantile([values[k] for k in COMPACT if k in values], 0.6) or fallback
    if name == "compact_q0.85":
        return _linear_quantile([values[k] for k in COMPACT if k in values], 0.85) or fallback
    if name == "ext_q0.25":
        return _linear_quantile([values[k] for k in EXT if k in values], 0.25) or fallback
    if name == "nomodel_q1":
        available = [values[k] for k in NOMODEL if k in values]
        return max(available) if available else fallback
    if name == "models_q1":
        available = [values[k] for k in ("CAMS", "MERRA", "GEE_mean") if k in values]
        return max(available) if available else fallback
    return fallback


def _select(values: dict[str, float], records: dict[str, object]) -> tuple[str, float | None]:
    q1 = values.get("Q1")
    r2 = values.get("R2")
    gee = values.get("GEE_mean")
    merra = values.get("MERRA")
    vdust = values.get("VDUST")
    context = records.get("context", {})
    sulphate_frac = None
    if isinstance(context, dict):
        sulphate_frac = _safe_float(
            context.get("cams_sulphate_aerosol_optical_depth_at_550nm_surface_frac")
        )
    r2_std = None
    r2_record = records.get("R2")
    if isinstance(r2_record, dict):
        solver = r2_record.get("solver")
        if isinstance(solver, dict):
            r2_std = _safe_float(solver.get("aot_std"))

    if _leq(q1, Q1_LOW_AOD_MAX):
        if _leq(r2, R2_CLEAN_AOD_MAX):
            if _leq(merra, MERRA_LOW_AOD_MAX):
                if _leq(vdust, VIIRS_DUST_LOW_AOD_MAX):
                    value = _candidate(values, "compact_q0.6")
                    return "1.1*compact_q0.6", None if value is None else _clip_aod(1.1 * value)
                return "ext_q0.25", _candidate(values, "ext_q0.25")
            return "nomodel_q1", _candidate(values, "nomodel_q1")
        r2_minus_gee = None if r2 is None or gee is None else r2 - gee
        if _leq(r2_minus_gee, R2_BELOW_MODEL_MARGIN):
            value = _candidate(values, "GEE_mean")
            return "1.1*GEE_mean", None if value is None else _clip_aod(1.1 * value)
        if _leq(sulphate_frac, CAMS_SULPHATE_FRACTION_MAX):
            value = _candidate(values, "GEE_mean")
            return "1.5*GEE_mean", None if value is None else _clip_aod(1.5 * value)
        value = _candidate(values, "GEE_mean")
        return "1.1*GEE_mean", None if value is None else _clip_aod(1.1 * value)

    if _leq(r2_std, R2_SPATIAL_STD_MAX):
        vdust_over_gee = None if vdust is None or gee is None or abs(gee) < 1e-12 else vdust / gee
        if _leq(vdust_over_gee, VIIRS_DUST_TO_MODEL_MAX):
            return "compact_q0.85", _candidate(values, "compact_q0.85")
        return "R2", _candidate(values, "R2")

    value = _candidate(values, "MERRA")
    return "1.25*MERRA", None if value is None else _clip_aod(1.25 * value)


def build_one(root: Path, rows: dict[str, dict[str, str]], matchup_id: str) -> dict[str, object]:
    row = rows[matchup_id]
    values, records = _load_inputs(root, matchup_id)
    selected, retrieved = _select(values, records)
    truth = float(row["aeronet_aod550_mean"])
    base: dict[str, object] = {
        "matchup_id": matchup_id,
        "site": row.get("site") or matchup_id.split("__")[0],
        "source": "multisource_tree_v2",
        "truth": truth,
        "selector": "rounded_regime_candidate_tree_v2",
        "selected_candidate": selected,
        "available_sources": sorted(values),
        "note": "Rounded multi-source candidate selector; verify on an independent holdout before production use.",
    }
    if retrieved is None:
        base.update({"status": "FAILED", "error_type": "MissingAOD", "reason": "selector produced null AOD"})
        return base
    base.update(
        {
            "status": "OK",
            "retrieved": retrieved,
            "err": retrieved - truth,
            "within_ee": _within_ee(retrieved, truth),
        }
    )
    return base


def main(argv: list[str]) -> int:
    mids = argv or [line.strip() for line in (ROOT / "campaign250_mids.txt").read_text().splitlines() if line.strip()]
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    OUT.mkdir(parents=True, exist_ok=True)
    for matchup_id in mids:
        record = build_one(ROOT, rows, matchup_id)
        (OUT / f"{matchup_id}.json").write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {len(mids)} records to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

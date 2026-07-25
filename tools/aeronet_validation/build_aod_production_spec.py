"""Build the production reproduction specification for the validated AOD method."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import importlib.metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tarfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_OUTPUT = ROOT / "reports/aod-production-reproduction-spec-20260713"
WEB_ASSETS = Path(__file__).with_name("aod_production_spec")

MODEL_MANIFEST = ROOT / "reports/aod-et35-global-offset-model-manifest-20260713.json"
MODEL = ROOT / "models/aod_et35_global_offset_20260713.joblib"
SELECTION_REPORT = ROOT / "reports/aod-target-domain-recipe-selection-20260713.json"
OFFSET_REPORT = ROOT / "reports/aod-global-offset-validation-20260713.json"
SEED_REPORT = ROOT / "reports/aod-selected-et-seed-robustness-20260713.json"
ABLATION_REPORT = ROOT / "reports/aod-generic-feature-ablation-20260713.json"
FINAL_DASHBOARD = ROOT / "reports/aod-final-performance-dashboard-20260713"

TARGET_RESULTS = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
TARGET_MAIAC = ROOT / "maiac_qa_lowcloud20_native_adaptive"
TARGET_CAMS = ROOT / "gee_aerosol_context_campaign250"
TARGET_MIDS = ROOT / "campaign250_lowcloud20_mids.txt"
TARGET_COHORT_MANIFEST = ROOT / "campaign250_lowcloud20_mids.json"
MATCHUPS = ROOT / "matchups/matchups.csv"

TRAIN_RESULTS = ROOT / "phaseD_results_external_lowcloud588_native_maiac_b03_chi2_20260713"
TRAIN_MAIAC = ROOT / "maiac_qa_external_native_lowcloud588"
TRAIN_CAMS = ROOT / "gee_aerosol_context_train1000"
TRAIN_MIDS = ROOT / "external_train1000_lowcloud20_mids.txt"

LUT = ROOT.parent / "libradtran_continental_average_lut_1nm.zarr.zip"
LUT_PUBLIC_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)
REPORT_PUBLIC_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/siac_refactor/"
    "reports/aod-production-reproduction-spec-20260713/"
)

EXACT_ENV = {
    "MAIAC_AGG": "median",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "RAYON_NUM_THREADS": "1",
    "POLARS_MAX_THREADS": "1",
    "GDAL_NUM_THREADS": "1",
    "MALLOC_ARENA_MAX": "2",
    "PHASE_D_CONFIG_PRESET": "l2a_monthly_predictor",
    "PHASE_D_BESTPIXEL_ENDPOINT": "pc",
    "PHASE_D_BESTPIXEL_SOURCE": "l2a",
    "PHASE_D_USE_PIXI_BESTPIXEL": "1",
    "PHASE_D_PREDICT_VISIBLE_BANDS": "B02,B03,B04",
    "PHASE_D_SOLVE_BANDS": "B02,B03,B04",
    "PHASE_D_PREDICT_VISIBLE_FLOOR": "0.006",
    "PHASE_D_PREDICTOR_MODEL": "et20",
    "PHASE_D_ANCHOR_ITERATE": "1",
    "PHASE_D_ANCHOR_ITERATE_MODE": "anchor",
    "PHASE_D_TAU_PRIOR": "1",
    "PHASE_D_TAU_GATE": "8.0",
    "PHASE_D_BACKSTOP_CAL": "1",
    "PHASE_D_MAIAC_DIR": TARGET_MAIAC.name,
    "PHASE_D_ALLOW_CLOUD_RETRIEVAL": "0",
    "PHASE_D_IGNORE_CLOUD_WATER": "0",
    "PHASE_D_RETRIEVAL_SPATIAL_POOL": "scene_mean",
    "PHASE_D_SKIP_CORRECTION": "1",
    "PHASE_D_MAX_WORKERS": "1",
    "PHASE_D_SERIAL_POOLS": "1",
    "PHASE_D_COST_MODE": "chi2",
    "PHASE_D_LUT_PATH": str(LUT),
}

AOT_AXIS = [
    *[round(0.01 * value, 6) for value in range(1, 20)],
    *[round(0.2 + 0.025 * value, 6) for value in range(12)],
    *[round(0.5 + 0.05 * value, 6) for value in range(20)],
    *[round(1.5 + 0.1 * value, 6) for value in range(11)],
    *[round(2.75 + 0.25 * value, 6) for value in range(6)],
]

LUT_AXES = {
    "wavelength_nm": {"start": 340, "stop": 2600, "step": 1, "nodes": 2261},
    "sza_deg": [0, 5, 10, 15, 20, 25, 30, 35, 45, 55, 65, 75, 85],
    "vza_deg": [0, 5, 10, 15, 20, 25, 30],
    "raa_deg": [0, 20, 40, 60, 80, 100, 120, 140, 160, 180],
    "aot_550": [0.001, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5],
    "tcwv_mm": [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100],
    "ozone_du": [200, 250, 300, 350, 400, 450, 500, 550, 600],
    "altitude_km": [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9],
}

REPO_SOURCE_FILES = (
    "python/siac/adapters/satellite/sentinel2.py",
    "python/siac/adapters/bestpixel.py",
    "python/siac/adapters/bestpixel_surface_prior.py",
    "python/siac/adapters/atmo/mcd19_earthaccess.py",
    "python/siac/adapters/atmo/maiac_day_aod.py",
    "python/siac/algorithms/grid/assembler.py",
    "python/siac/algorithms/rt/lut/backend.py",
    "python/siac/algorithms/rt/lut/_spectral.py",
    "python/siac/algorithms/rt/lut/_spectral_math.py",
    "python/siac/algorithms/solver/surface_driven.py",
    "python/siac/algorithms/surface/seasonal_predictor.py",
    "python/siac/algorithms/aod_calibration.py",
    "src/siac_rs/src/surface_driven.rs",
    "tools/aeronet_validation/matchup.py",
    "tools/aeronet_validation/build_low_cloud_cohort.py",
    "tools/aeronet_validation/stage_maiac_prior.py",
    "tools/aeronet_validation/collect_gee_aerosol_context.py",
    "tools/aeronet_validation/aod_residual_calibration.py",
    "tools/aeronet_validation/select_generic_aod_calibrator.py",
    "tools/aeronet_validation/evaluate_global_aod_offset.py",
    "tools/aeronet_validation/evaluate_selected_aod_seed_robustness.py",
    "tools/aeronet_validation/export_generic_aod_calibrator.py",
    "tools/aeronet_validation/verify_aod_reproduction_package.py",
    "tools/aeronet_validation/campaign250_lowcloud20_single_prior_adaptive_submit.sbatch",
)

HARNESS_SOURCE_FILES = (
    ROOT / "phaseD_run_one_site.py",
    ROOT / "phaseD_array_runner.py",
    ROOT / "phaseD_prebuilt_maiac.py",
    ROOT / "phaseD_gee_sources.py",
)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, records: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")))
            stream.write("\n")


def _read_ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _records(directory: Path, ids: Sequence[str]) -> list[dict[str, Any]]:
    records = []
    for matchup_id in ids:
        path = directory / f"{matchup_id}.json"
        if not path.exists():
            raise FileNotFoundError(path)
        record = _load_json(path)
        if str(record.get("matchup_id") or matchup_id) != matchup_id:
            raise ValueError(f"matchup_id mismatch in {path}")
        records.append(record)
    return records


def _sanitize_context(record: Mapping[str, Any]) -> dict[str, Any]:
    # Inference reads only the operational fields. Remove evaluation labels from
    # the frozen context package so production integrations cannot consume them.
    return {key: value for key, value in record.items() if key not in {"truth", "regime"}}


def _finite(values: Iterable[Any]) -> list[float]:
    output: list[float] = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            output.append(number)
    return output


def _quantiles(values: Sequence[float]) -> dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        return {}

    def pick(fraction: float) -> float:
        index = fraction * (len(ordered) - 1)
        low = int(math.floor(index))
        high = int(math.ceil(index))
        if low == high:
            return ordered[low]
        return ordered[low] + (ordered[high] - ordered[low]) * (index - low)

    return {
        "min": ordered[0],
        "p25": pick(0.25),
        "median": pick(0.5),
        "p75": pick(0.75),
        "max": ordered[-1],
    }


def _target_summary(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    statuses: dict[str, int] = {}
    for record in records:
        status = str(record.get("status"))
        statuses[status] = statuses.get(status, 0) + 1
    return {
        "count": len(records),
        "statuses": statuses,
        "physical_hits": sum(bool(record.get("within_ee")) for record in records),
        "tau_gate_fired": sum(
            bool((record.get("solver") or {}).get("surface_tau_gate_fired")) for record in records
        ),
        "mask_fallback": [
            str(record["matchup_id"])
            for record in records
            if bool((record.get("solver") or {}).get("surface_cloud_mask_bypassed"))
            or bool((record.get("solver") or {}).get("surface_water_mask_bypassed"))
        ],
        "anchor_iteration_count": sum(isinstance(record.get("anchor_iterate"), Mapping) for record in records),
        "runtime_seconds": _quantiles(_finite(record.get("runtime_s") for record in records)),
        "valid_aot_fraction": _quantiles(
            _finite(record.get("valid_aot_fraction") for record in records)
        ),
    }


def _maiac_summary(records: Sequence[Mapping[str, Any]], matchup_rows: Mapping[str, Mapping[str, str]]) -> dict[str, Any]:
    valid = [record for record in records if record.get("status") == "OK"]
    offsets = []
    for record in valid:
        when = str(matchup_rows[str(record["matchup_id"])]["sensing_time_utc"])[:10]
        day = str(record.get("day") or when)
        offsets.append(abs((datetime.fromisoformat(day) - datetime.fromisoformat(when)).days))
    return {
        "count": len(records),
        "ok": len(valid),
        "no_valid": len(records) - len(valid),
        "same_day": sum(value == 0 for value in offsets),
        "one_day_offset": sum(value == 1 for value in offsets),
        "half_degree_counts": {
            str(value): sum(float(record.get("half_deg", 0)) == value for record in records)
            for value in (0.12, 0.25)
        },
        "n_good": _quantiles(_finite(record.get("n_good") for record in valid)),
    }


def _load_matchups() -> tuple[list[str], dict[str, dict[str, str]]]:
    with MATCHUPS.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return fieldnames, {str(row["matchup_id"]): row for row in rows}


def _write_matchup_subset(
    path: Path,
    fieldnames: Sequence[str],
    rows: Mapping[str, Mapping[str, str]],
    ids: Sequence[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for matchup_id in ids:
            writer.writerow(rows[matchup_id])


@contextmanager
def _temporary_environment(values: Mapping[str, str]):
    previous = {key: os.environ.get(key) for key in values}
    os.environ.update(values)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _resolved_config() -> dict[str, Any]:
    for path in (str(REPO_ROOT / "python"), str(ROOT)):
        if path not in sys.path:
            sys.path.insert(0, path)
    with _temporary_environment(EXACT_ENV):
        from phaseD_run_one_site import build_config

        config = build_config(5)
    if hasattr(config, "model_dump"):
        return dict(config.model_dump(mode="json"))
    raise TypeError("Resolved SIAC config does not expose model_dump")


def _git(command: Sequence[str]) -> str:
    return subprocess.run(
        ["git", *command],
        cwd=REPO_ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()


def _environment_payload() -> dict[str, Any]:
    packages = {}
    for name in (
        "numpy",
        "scipy",
        "xarray",
        "rioxarray",
        "rasterio",
        "scikit-learn",
        "joblib",
        "pydantic",
        "bestpixel",
    ):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    status = _git(["status", "--short"])
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "git_commit": _git(["rev-parse", "HEAD"]),
        "git_worktree_dirty": bool(status),
        "git_status_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "pixi_lock_sha256": _sha256(REPO_ROOT / "pixi.lock"),
        "packages": packages,
        "determinism_environment": EXACT_ENV,
        "note": (
            "The validated source tree was dirty. The critical-source archive and per-file "
            "hashes, not the commit alone, define the tested implementation."
        ),
    }


def _feature_definition(name: str) -> tuple[str, str, str]:
    if name in {"atmo_aot_mean", "atmo_aot_median"}:
        return (
            "MAIAC",
            "AOT550",
            "Mean or median of the staged tile-constant MAIAC AOT field. The two are equal here.",
        )
    if name in {"atmo_aot_unc_mean", "atmo_aot_unc_median"}:
        return (
            "MAIAC",
            "AOT550",
            "Staged calibrated uncertainty: clip(max(native_unc, 0.04) + 0.18*max(AOT-0.3, 0) + sparse, 0.04, 0.6), sparse=0.05 when n_good<8.",
        )
    if name == "consistency_atmo_minus_siac":
        return ("Contrast", "AOT550", "staged_MAIAC_AOT - physical_scene_mean_AOT")
    if name == "consistency_atmo_to_siac":
        return ("Contrast", "ratio", "staged_MAIAC_AOT / physical_scene_mean_AOT; missing when denominator magnitude < 1e-8")
    if name == "consistency_cams_minus_siac":
        return ("Contrast", "AOT550", "CAMS_total_AOT550 - physical_scene_mean_AOT")
    if name == "consistency_cams_to_siac":
        return ("Contrast", "ratio", "CAMS_total_AOT550 / physical_scene_mean_AOT; missing when denominator magnitude < 1e-8")
    if name.startswith("context_cams_"):
        unit = "fraction" if name.endswith("_frac") else "AOT"
        return ("CAMS", unit, "CAMS NRT field averaged over the scene-centred 0.12 x 0.12 degree box and +/-3 hour window.")
    if name == "solver_cost_final":
        return ("Solver", "chi-square", "Median finite observation-only pooled minimum cost (jmin) from the returned static or tau solve.")
    if name == "solver_cost_final_per_band":
        return ("Solver", "chi-square/band", "solver_cost_final divided by the three solve bands.")
    if name.endswith("_argmin_cost"):
        return ("Solver", "chi-square", "Minimum of the scene-median single-band cost curve over the 68 AOT nodes.")
    if name.endswith("_cost_final_node"):
        return ("Solver", "chi-square", "Scene-median single-band cost evaluated at the final scene AOT node.")
    if name.endswith("_residual_final_node"):
        return ("Solver", "reflectance", "Scene-median absolute BOA-minus-prior residual at the final scene AOT node.")
    if name == "solver_surface_band_argmin_spread":
        return ("Solver", "AOT550", "Maximum minus minimum of the B02/B03/B04 single-band argmin AOT values.")
    if name == "solver_surface_cost_curve_min_cost":
        return ("Solver", "chi-square", "Minimum of the scene-median combined surface cost curve over AOT nodes.")
    if name == "solver_surface_cost_curve_second_delta":
        return ("Solver", "chi-square", "Second-lowest minus lowest scene-median combined surface cost.")
    if name == "solver_surface_static_cost_per_band":
        return ("Solver", "chi-square/band", "Static first-solve final_cost divided by three; always present when the tau gate is configured.")
    if name == "solver_surface_tau_cost_per_band":
        return ("Solver", "chi-square/band", "Tau-dependent rerun final_cost divided by three; missing and median-imputed when the gate does not fire.")
    if name in {"solver_surface_tau_gate_fired", "solver_surface_tau_prior_enabled"}:
        return ("Solver", "boolean", "Numeric 0/1 solver diagnostic for the cost-gated tau-dependent surface-prior path.")
    if name == "time_utc_sin":
        return ("Time", "unitless", "sin(2*pi*(UTC hour + minute/60 + second/3600)/24).")
    return ("Other", "varies", "Operational feature emitted by the frozen extractor.")


def _source_map() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for relative in REPO_SOURCE_FILES:
        path = REPO_ROOT / relative
        if not path.exists():
            raise FileNotFoundError(path)
        records.append(
            {
                "component": relative,
                "origin": "repository",
                "path": str(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    for path in HARNESS_SOURCE_FILES:
        if not path.exists():
            raise FileNotFoundError(path)
        records.append(
            {
                "component": path.name,
                "origin": "campaign harness",
                "path": str(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    return records


def _write_source_archive(path: Path, source_map: Sequence[Mapping[str, Any]]) -> None:
    with tarfile.open(path, "w:gz") as archive:
        for record in source_map:
            source = Path(str(record["path"]))
            prefix = "repo" if record["origin"] == "repository" else "harness"
            arcname = f"{prefix}/{record['component']}"
            archive.add(source, arcname=arcname, recursive=False)


def _copy(path: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _recipe_toml() -> str:
    env_lines = "\n".join(f'{key} = "{value}"' for key, value in EXACT_ENV.items())
    return f'''# Normalized reproduction recipe. resolved-config.json is the authoritative SIAC config.
[benchmark]
cohort = "campaign250_lowcloud20_mids.txt"
cloud_definition = "OmniCloudMask cloud_frac < 0.20"
expected_error = "abs(retrieved - aeronet) <= 0.05 + 0.15*aeronet"
retrieval_extraction = "scene_mean"

[surface]
type = "single_s2_l2a_bestpixel"
lookback_years = 5
seasonal_window_months = 1
top_k = 15
low_aod_day_fraction = 0.6
robust_clip_mad = 1.5
anchors = ["B8A", "B11", "B12"]
predicted_bands = ["B02", "B03", "B04"]
anchor_regressor = "ExtraTreesRegressor(n_estimators=20,min_samples_leaf=5,random_state=0,n_jobs=1)"
anchor_fixed_point_passes = 1

[solver]
bands = ["B02", "B03", "B04"]
grid_m = 60
cost = "chi2"
aot_axis = "acixthree_68"
tcwv_cm = 2.0
pool_radius_m = 600
pool_window_cells = 20
pool_min_count = 80
tau_gate_cost_per_band = 8.0
aerosol_species = "none"

[calibration]
artifact = "aod-calibrator.joblib"
feature_schema_count = 137
selected_feature_count = 35
trees = 1500
criterion = "absolute_error"
random_state = 20260713
min_samples_leaf = 4
max_features = 0.5
global_log_offset = -0.0125

[environment]
{env_lines}
'''


def _pseudocode() -> str:
    return '''INPUT: Sentinel-2 L1C SAFE, AOI, sensing time, frozen model, data-service credentials
1. Read MSI bands and apply (DN + RADIO_ADD_OFFSET) / QUANTIFICATION_VALUE; clip [0, 1.5].
2. Build OmniCloudMask at 10 m and buffered land/water exclusion; parse mean scene geometry.
3. Stage QA-best MCD19A2 AOT near the sensing day; use it as the tile-constant AOT backstop.
4. For each of five complete prior years and scene month +/-1:
   a. get QA-MAIAC AOT by day; retain the lowest ceil(0.6*N) days in each year-month;
   b. build one top-15 Planetary Computer Sentinel-2 L2A bestpixel composite.
5. Reproject the up-to-15 seven-band composites to the observation-aligned grid.
6. Robustly reduce them to temporal median and RMSE uncertainty after 1.5*MAD clipping.
7. Correct scene B8A/B11/B12 TOA to BOA at staged MAIAC AOT through the same RT LUT.
8. Train one 20-tree ExtraTrees model per historical realization from historical anchors plus
   four-band localizer to historical B02/B03/B04; predict each visible band and median ensemble.
9. Add per-band intercept+slope*AOT debias, replace valid predictions in (0.001, 0.6), floor sigma 0.006.
10. Area-resample native B02/B03/B04 TOA and the surface prior to the 60 m retrieval grid.
11. For each of 68 AOT nodes, hold TCWV=2 cm, O3=0.30 atm-cm, elevation=0, use scene-mean
    geometry, linearly interpolate the continental-average libRadtran LUT, convolve with S2 RSRF,
    correct TOA to BOA, and compute the three-band chi-square surface cost.
12. For every AOT node, take the centred 20x20 rolling median (>=80 finite cells). Add the calibrated
    Gaussian MAIAC backstop and choose the AOT node with minimum total cost per valid pixel.
13. If static median observation cost / 3 > 8, rebuild the visible prior at every candidate AOT and
    rerun the same solve using the tau-dependent prior.
14. Compute pass-1 scene-mean AOT. Rebuild only the anchor-conditioned surface prior at this AOT,
    repeat steps 10-13 once against the original MAIAC backstop, and retain pass 2.
15. Emit physical AOT as mean of all finite AOI grid cells and collect frozen operational diagnostics.
16. Query CAMS NRT +/-3 h over a 0.12 x 0.12 degree box; extract the 137-feature frozen schema.
17. Inference pipeline: median imputation -> standardization -> frozen F-test 35-feature selection ->
    1500-tree ExtraTrees L1 model predicting log((AERONET+1/3)/(physical+1/3)).
18. raw=clip((physical+1/3)*exp(model_output)-1/3,0,4);
    final=clip((raw+1/3)*exp(-0.0125)-1/3,0,4).
OUTPUT: physical AOT map, uncertainty/cost diagnostics, scene physical AOT, calibrated scene AOT.
'''


def _h(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _fmt(value: Any, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "n/a"
    return f"{number:.{digits}f}"


def _table(headers: Sequence[str], rows: Iterable[Sequence[Any]], *, css_class: str = "") -> str:
    head = "".join(f"<th scope=\"col\">{_h(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>" for row in rows
    )
    return f'<div class="table-wrap"><table class="{_h(css_class)}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _render_html(
    *,
    model_manifest: Mapping[str, Any],
    selection: Mapping[str, Any],
    offset: Mapping[str, Any],
    target_summary: Mapping[str, Any],
    maiac_summary: Mapping[str, Any],
    source_map: Sequence[Mapping[str, Any]],
    artifacts: Sequence[Mapping[str, Any]],
    generated_utc: str,
) -> str:
    selected = list(model_manifest["features"]["selected"])
    importance = {
        row["feature"]: float(row["importance"])
        for row in model_manifest["features"]["importance"]
    }
    feature_rows = []
    for index, name in enumerate(selected, start=1):
        group, unit, definition = _feature_definition(name)
        feature_rows.append(
            f'<tr data-feature-row data-search="{_h((name + " " + group + " " + definition).lower())}" data-group="{_h(group)}">'
            f"<td>{index}</td><td><code>{_h(name)}</code></td><td>{_h(group)}</td>"
            f"<td>{_h(unit)}</td><td>{_h(definition)}</td><td>{importance.get(name, 0):.4f}</td></tr>"
        )
    source_rows = "".join(
        f'<tr data-source-row data-search="{_h((str(row["component"]) + " " + str(row["origin"])).lower())}">'
        f'<td><code>{_h(row["component"])}</code></td><td>{_h(row["origin"])}</td>'
        f'<td class="hash"><code>{_h(row["sha256"])}</code></td><td>{int(row["bytes"]):,}</td></tr>'
        for row in source_map
    )
    artifact_rows = "".join(
        f'<tr><td><a href="{_h(row["href"])}">{_h(row["name"])}</a></td>'
        f'<td>{_h(row["purpose"])}</td><td>{int(row["bytes"]):,}</td>'
        f'<td class="hash"><code>{_h(row["sha256"])}</code></td></tr>'
        for row in artifacts
    )

    target_metrics = model_manifest["validation_reproduction"]["target"]
    confirmation = model_manifest["validation_reproduction"]["target_confirmation"]
    external = model_manifest["validation_reproduction"]["external_holdout"]
    unadjusted_block = offset["target_unadjusted"]
    unadjusted = unadjusted_block["candidate"]
    physical_hits = int(target_summary["physical_hits"])
    physical_pct = 100 * physical_hits / int(target_summary["count"])
    mask_cases = target_summary["mask_fallback"]
    mask_case_list = ", ".join(f"<code>{_h(value)}</code>" for value in mask_cases)
    recipe = model_manifest["recipe"]
    policy = selection["selection_policy"]

    input_rows = [
        ("Sentinel-2 L1C", "GCS public bucket", "One fixed product ID and AOI; B02/B03/B04 solve, B8A/B11/B12 anchor", "Product IDs in matchups subsets"),
        ("Sentinel-2 L2A history", "Planetary Computer / bestpixel", "Five complete prior years; scene month +/-1; top_k=15; catalog cloud <=90%", "Source service/cache; full pixel cubes not embedded"),
        ("MAIAC MCD19A2", "Earthdata native tiles + GEE day gate", "Best-quality QA; nearest day within +/-2; tile-fixed backstop; clean-day surface gate", "All scalar backstops frozen as JSONL"),
        ("CAMS NRT", "ECMWF/CAMS/NRT via Earth Engine", "+/-3 hours; mean over 0.12 x 0.12 degree box; 40 km reduction", "All selected target/training contexts frozen as JSONL"),
        ("AERONET V3", "matchups.csv", "+/-30 minutes around overpass; arithmetic mean AOD550; minimum one observation", "Target and training matchup subsets"),
        ("RT LUT", "libRadtran continental average", "1 nm spectral LUT; Sentinel-2 RSRF convolution; linear state interpolation", "External 198 GB public artifact; identity recorded"),
        ("Calibrator", "Frozen joblib artifact", "137-schema operational input; 35 selected features; no site/geography/case routing", "Model embedded in downloads"),
    ]

    validation_rows = [
        ("Physical retrieval", "152", f"{physical_hits}", f"{physical_pct:.2f}%", "Before supervised calibration"),
        ("Selected model, no offset", str(unadjusted["count"]), str(unadjusted["hits"]), f'{unadjusted["within_ee_percent"]:.2f}%', "Frozen recipe, descriptive all-target result"),
        ("Final exported artifact", str(target_metrics["count"]), str(target_metrics["hits"]), f'{target_metrics["within_ee_percent"]:.2f}%', "Primary 152-case result"),
        ("Confirmation fold 4", str(confirmation["count"]), str(confirmation["hits"]), f'{confirmation["within_ee_percent"]:.2f}%', "Labels not used for recipe or offset selection"),
        ("External holdout", str(external["count"]), str(external["hits"]), f'{external["within_ee_percent"]:.2f}%', "Distribution-shift warning"),
    ]

    axis_rows = [
        ("Wavelength", "340-2600 nm", "1 nm", "2261"),
        ("SZA", "0-85 deg", "13 explicit nodes", "13"),
        ("VZA", "0-30 deg", "5 deg", "7"),
        ("RAA", "0-180 deg", "20 deg", "10"),
        ("AOT550", "0.001-5", "16 explicit nodes", "16"),
        ("TCWV", "0-100 mm", "12 explicit nodes", "12"),
        ("Ozone", "200-600 DU", "50 DU", "9"),
        ("Altitude", "0-9 km", "15 explicit nodes", "15"),
    ]

    return f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Production reproduction specification for the validated SIAC low-cloud Sentinel-2 AOD algorithm.">
  <title>SIAC AOD production reproduction specification</title>
  <link rel="stylesheet" href="assets/app.css">
  <script defer src="assets/app.js"></script>
</head>
<body>
  <header class="topbar">
    <a class="brand" href="#top" aria-label="Go to top"><span class="brand-mark">S</span><span>SIAC</span></a>
    <span class="topbar-context">Production reproduction specification</span>
    <a class="topbar-link" href="../aod-final-performance-dashboard-20260713/">Performance dashboard</a>
  </header>

  <div class="page-shell" id="top">
    <aside class="toc" aria-label="Document navigation">
      <p class="toc-title">Contents</p>
      <nav>
        <a href="#overview">System overview</a>
        <a href="#contracts">Input contracts</a>
        <a href="#benchmark">Benchmark</a>
        <a href="#ingestion">L1C ingestion</a>
        <a href="#atmosphere">Atmospheric prior</a>
        <a href="#surface">Surface prior</a>
        <a href="#rt">RT model</a>
        <a href="#solver">AOD solver</a>
        <a href="#iteration">Iteration and result</a>
        <a href="#calibration">Calibration</a>
        <a href="#features">Selected features</a>
        <a href="#validation">Validation</a>
        <a href="#production">Production controls</a>
        <a href="#artifacts">Reproduction files</a>
        <a href="#source-map">Source map</a>
        <a href="#limits">Limits</a>
        <a href="#checklist">Acceptance checklist</a>
      </nav>
    </aside>

    <main>
      <section class="title-band" aria-labelledby="page-title">
        <div class="title-copy">
          <p class="eyebrow">Validated method / frozen 13 July 2026</p>
          <h1 id="page-title">Low-cloud Sentinel-2 AOD retrieval</h1>
          <p class="lede">Complete production contract for data acquisition, radiometric preparation, the single S2 L2A surface prior, MAIAC backstop, libRadtran LUT, physical inversion, one fixed-point rerun, CAMS diagnostics, and the frozen generic AOD calibrator.</p>
        </div>
        <div class="metric-strip" aria-label="Primary validation metrics">
          <div><strong>88.16%</strong><span>within expected error</span></div>
          <div><strong>134 / 152</strong><span>final target hits</span></div>
          <div><strong>111 / 152</strong><span>physical hits</span></div>
          <div><strong>152 / 152</strong><span>terminal physical jobs</span></div>
        </div>
        <div class="claim-boundary"><strong>Claim boundary.</strong> The 88.16% result applies to the fixed 152-image cohort with raw OmniCloudMask cloud fraction strictly below 20%. It is not a global accuracy claim. External holdout performance is {external['within_ee_percent']:.2f}%.</div>
      </section>

      <section id="overview" class="doc-section">
        <div class="section-heading"><p>01</p><div><h2>System overview</h2><p>One uniform physical retrieval followed by one uniform learned calibration. There is no per-case source routing.</p></div></div>
        <div class="flow" aria-label="End to end algorithm flow">
          <div><span>1</span><strong>S2 L1C</strong><small>calibrate, geometry, masks</small></div><i></i>
          <div><span>2</span><strong>S2 L2A history</strong><small>five-year clean-day dictionary</small></div><i></i>
          <div><span>3</span><strong>Surface prior</strong><small>NIR/SWIR anchored visible prediction</small></div><i></i>
          <div><span>4</span><strong>RT sweep</strong><small>68 AOT nodes, 3 visible bands</small></div><i></i>
          <div><span>5</span><strong>Physical AOT</strong><small>pool, backstop, one rerun</small></div><i></i>
          <div><span>6</span><strong>Calibration</strong><small>CAMS + diagnostics, one frozen model</small></div>
        </div>
        <div class="two-col">
          <div>
            <h3>Operational output</h3>
            <p>The physical stage returns a 60 m AOT field, uncertainty, observation-only cost, masks, and solver diagnostics. The validated scalar is the arithmetic mean of every finite AOT cell in the requested AOI after the second anchor pass. The calibrator then returns one scalar final AOD in [0, 4].</p>
          </div>
          <div>
            <h3>Reproduction levels</h3>
            <ol class="compact-list">
              <li><strong>Exact scalar replay:</strong> frozen 152 result records, contexts, feature schema, model, and expected predictions are embedded.</li>
              <li><strong>Deterministic physical replay:</strong> requires the same S2 pixels, bestpixel composites, staged MAIAC scalars, LUT, code hashes, and serialized warps.</li>
              <li><strong>Source refresh:</strong> reruns catalog queries and is scientifically comparable, but mutable upstream services can change bytes.</li>
            </ol>
          </div>
        </div>
        <div class="callout warning"><strong>Implementation identity.</strong> Git commit <code>{_h(_git(['rev-parse', 'HEAD']))}</code> is not sufficient because the validated tree was dirty. <a href="downloads/critical-source-snapshot.tar.gz">The critical-source snapshot</a> and <a href="downloads/source-map.json">per-file SHA-256 map</a> define the tested implementation.</div>
      </section>

      <section id="contracts" class="doc-section">
        <div class="section-heading"><p>02</p><div><h2>Input contracts</h2><p>Required data, exact selection rules, and what is frozen in this package.</p></div></div>
        {_table(("Input", "Source", "Selection and units", "Frozen handoff"), input_rows)}
        <h3>AOI and units</h3>
        <ul class="fact-grid">
          <li><strong>Target AOI</strong><span>station longitude +/-0.045 deg; latitude +/-0.0351 deg</span></li>
          <li><strong>AOD wavelength</strong><span>550 nm, unitless</span></li>
          <li><strong>Reflectance</strong><span>unitless TOA or BOA</span></li>
          <li><strong>TCWV</strong><span>solver input 2.0 cm; LUT coordinate converted to 20 mm</span></li>
          <li><strong>Ozone</strong><span>0.30 atm-cm; converted to 300 DU</span></li>
          <li><strong>Elevation</strong><span>0 km in the selected benchmark; DEM path disabled</span></li>
        </ul>
      </section>

      <section id="benchmark" class="doc-section">
        <div class="section-heading"><p>03</p><div><h2>Benchmark construction</h2><p>The 152 denominator is a fixed subset of a prior 250-case campaign, not 528 and not a fresh cloud-catalog query.</p></div></div>
        <div class="three-col metrics-detail">
          <div><span>Cohort</span><strong>152 / 250</strong><p>Strictly <code>cloud_frac &lt; 0.20</code> from raw OmniCloudMask output in the fixed campaign run. One of 250 was unclassified.</p></div>
          <div><span>AERONET truth</span><strong>&plusmn;30 min</strong><p>Arithmetic mean of available AOD550 observations around the S2 sensing time; at least one observation.</p></div>
          <div><span>Expected error</span><strong>0.05 + 0.15A</strong><p>A retrieval is a hit when <code>|AOD_retrieved - AOD_AERONET| &le; 0.05 + 0.15*AOD_AERONET</code>.</p></div>
        </div>
        <h3>Product deduplication</h3>
        <p>The campaign catalog searched Sentinel-2 L1C over 2024. Products at the same sensing minute were deduplicated by choosing the highest processing baseline and then the MGRS tile. The 152-product IDs are frozen in <a href="downloads/target-matchups.csv">target-matchups.csv</a>.</p>
        <h3>Label-use audit</h3>
        <div class="audit-grid">
          <div><strong>External labels</strong><span>Used to fit the calibrator on 588 records from 234 sites.</span></div>
          <div><strong>Target folds 0-3</strong><span>Used to choose the recipe and -0.0125 global log offset.</span></div>
          <div><strong>Target fold 4</strong><span>26 records; labels held back from recipe and offset selection.</span></div>
          <div><strong>Inference</strong><span>No AERONET target, site, latitude/longitude, geography, or case route is an input.</span></div>
        </div>
        <p class="footnote">Site fold is deterministic: <code>int(sha256(site_utf8).hexdigest()[:8], 16) % 5</code>. Development folds are 0-3; confirmation fold is 4.</p>
      </section>

      <section id="ingestion" class="doc-section">
        <div class="section-heading"><p>04</p><div><h2>Sentinel-2 L1C ingestion</h2><p>Radiometry, angle handling, masks, and grid alignment.</p></div></div>
        <div class="formula"><span>TOA calibration</span><code>rho_TOA = clip((DN + RADIO_ADD_OFFSET_band) / QUANTIFICATION_VALUE, 0, 1.5)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <ul class="method-list">
          <li><strong>Metadata.</strong> Read <code>QUANTIFICATION_VALUE</code> and per-band <code>RADIO_ADD_OFFSET</code> from <code>MTD_MSIL1C.xml</code>. If quantification is absent, the adapter warns and uses 10000. A production run should fail instead if exact reproducibility is required.</li>
          <li><strong>Bands.</strong> Native B02, B03, B04 are loaded for the solve. B8A, B11, and B12 are additionally loaded for scene anchoring. The adapter can align mixed-resolution bands bilinearly for its observation bundle, but grid assembly deliberately reloads each native band and area-resamples it once to the 60 m solver grid.</li>
          <li><strong>Angles.</strong> Parse native 23 x 23, approximately 5 km sun-angle grids and per-band/per-detector view-angle grids. Azimuths are unwrapped before interpolation; all runtime angles are radians. The physical sweep replaces each field with its finite AOI mean, then derives relative azimuth and normalizes it to 0-180 degrees.</li>
          <li><strong>Cloud.</strong> OmniCloudMask runs on CPU at 10 m. Class output is converted to a boolean mask where true means excluded. Downsampling is conservative: a target cell is excluded if its contributing source footprint contains cloud.</li>
          <li><strong>Water.</strong> The landWater2020 VRT is subset to the AOI and dilated by 32 native mask pixels before exclusion on the retrieval grid.</li>
          <li><strong>Resampling.</strong> Continuous TOA, surface reflectance, uncertainty, and localizer planes use area averaging onto the 60 m retrieval grid. Atmospheric state uses bilinear interpolation. All resampling workers are fixed to one.</li>
        </ul>
        <div class="callout neutral"><strong>PSF status.</strong> The resolved preset contains historical PSF parameters, but <code>build_toa_psf_config</code> enables observation-side PSF only for <code>kernel_model</code> and <code>whittaker</code>. The selected <code>bestpixel</code> method does not execute PSF convolution or shift fitting.</div>
      </section>

      <section id="atmosphere" class="doc-section">
        <div class="section-heading"><p>05</p><div><h2>MAIAC atmospheric prior</h2><p>One tile-fixed aerosol centre regularizes the physical solve and initializes the surface anchor.</p></div></div>
        <h3>Native staging</h3>
        <ol class="numbered-method">
          <li><span>1</span><div>Download MCD19A2 granules from the sensing day +/-2 days with <code>best_quality_qa=true</code>, at most 40 granules.</div></li>
          <li><span>2</span><div>On each native product grid retain pixels with cloud bits 0-2 equal to 1, adjacency bits 5-7 equal to 0, and AOD-QA bits 8-11 equal to 0. Apply product scale and fill masking before aggregation.</div></li>
          <li><span>3</span><div>Select native pixel centres inside a station-centred +/-0.12 degree box. Three sparse cases use +/-0.25 degrees. Group pixels by granule day.</div></li>
          <li><span>4</span><div>Select the nearest valid day; ties prefer the day with more pixels. Aggregate AOT and native uncertainty by median.</div></li>
          <li><span>5</span><div>Construct a tile-constant 60 m atmospheric state with AOT from MAIAC, TCWV=2.0 cm, O3=0.30 atm-cm, elevation=0, TCWV uncertainty=0.3, and ozone uncertainty=0.03.</div></li>
        </ol>
        <div class="formula"><span>Staged MAIAC uncertainty retained for calibration features</span><code>sigma_stage = clip(max(sigma_native, 0.04) + 0.18*max(AOT-0.3, 0) + (0.05 if n_good&lt;8 else 0), 0.04, 0.6)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <div class="formula"><span>Physical-solver backstop width</span><code>loose=max(0.5m,0.02); mid=max(0.07,0.5m/(1+exp(-(m-0.5)/0.15))); sigma=loose if m&lt;0.15 else mid</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <p>The solver ignores <code>sigma_stage</code> when calibrated backstop mode is enabled and recomputes the second formula from centre <code>m</code>. The calibrator feature extractor, however, receives the staged sidecar and therefore uses <code>sigma_stage</code>. Keeping these two uncertainty contracts separate is necessary for exact replay.</p>
        <ul class="fact-grid">
          <li><strong>Staged OK</strong><span>{maiac_summary['ok']} / {maiac_summary['count']}</span></li>
          <li><strong>Same day</strong><span>{maiac_summary['same_day']}</span></li>
          <li><strong>Previous/next day</strong><span>{maiac_summary['one_day_offset']}</span></li>
          <li><strong>Native good pixels</strong><span>min {_fmt(maiac_summary['n_good']['min'], 0)}, median {_fmt(maiac_summary['n_good']['median'], 0)}, max {_fmt(maiac_summary['n_good']['max'], 0)}</span></li>
        </ul>
        <p>Four staged records have no valid native AOD and defer to the live QA-filtered GEE provider: Barrow, MCO-Hanimaadhoo, Tai_Ping, and Yellowknife. The frozen result receipts preserve their exact backstop values.</p>
      </section>

      <section id="surface" class="doc-section">
        <div class="section-heading"><p>06</p><div><h2>Single S2 L2A surface prior</h2><p>Every scene uses the same source type and construction. MODIS supplies aerosol gating, not surface reflectance.</p></div></div>
        <div class="callout success"><strong>Source contract.</strong> Surface reflectance comes only from Sentinel-2 L2A Planetary Computer bestpixel composites. There is no MODIS surface branch, source selector, per-case route, or aerosol-species selector in the final method.</div>
        <h3>Historical dictionary</h3>
        <ol class="numbered-method">
          <li><span>1</span><div>Use the five complete years before the observation year. For each year query the scene month and its two adjacent months with calendar wrapping, producing at most 15 year-month realizations.</div></li>
          <li><span>2</span><div>For every year-month obtain QA-MAIAC AOT by day over the AOI. With no absolute AOD ceiling configured, sort scored days and retain the lowest <code>ceil(0.6*N)</code>. Pass the retained map to bestpixel with <code>reject_unknown=true</code>. If the day gate fails entirely, retain all candidate days and record the degraded provenance.</div></li>
          <li><span>3</span><div>Build one top-15 monthly composite from Planetary Computer S2 L2A with catalog scene cloud <=90%, UTM output, uncertainty emission, and canonical bands coastal, blue, green, red, nir, swir16, swir22.</div></li>
          <li><span>4</span><div>Reproject every realization to the observation-aligned target grid. S2 canonical identities map directly to B01, B02, B03, B04, B8A, B11, B12; the configured KNN mapper is not invoked for this same-sensor path.</div></li>
          <li><span>5</span><div>For each band and pixel compute the initial temporal median, MAD x 1.4826, discard realization values beyond 1.5 MAD (MAD=0 values remain), then recompute the median.</div></li>
        </ol>
        <div class="formula"><span>Base prior uncertainty</span><code>sigma_BOA = max(sqrt(mean((rho_i - median(rho_i))^2)), 0.006)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <p>With one realization the first composite uncertainty replaces the temporal spread. Uncertainty is NaN wherever the median BOA is NaN.</p>

        <h3>NIR/SWIR anchored visible prediction</h3>
        <p>The base temporal median remains the fallback. B02, B03, and B04 are replaced only where the generic scene-conditioned predictor returns a finite value in (0.001, 0.6).</p>
        <ol class="numbered-method">
          <li><span>A</span><div>Area-resample target-scene B8A/B11/B12 TOA to the prior grid and retain pixels where all anchors are finite and in (0, 1.2).</div></li>
          <li><span>B</span><div>Correct these anchors to BOA at the median tile MAIAC AOT through the exact same RT backend, using AOI-mean geometry, TCWV=2 cm, O3=0.30 atm-cm, and elevation=0.</div></li>
          <li><span>C</span><div>Build a four-plane localizer as the realization mean of historical B01/B02/B03/B04 at each pixel.</div></li>
          <li><span>D</span><div>For each usable historical realization with at least 200 finite training pixels, fit <code>ExtraTreesRegressor(n_estimators=20, min_samples_leaf=5, random_state=0, n_jobs=1)</code>. Inputs are that realization's B8A/B11/B12 plus the four localizer planes; outputs are its B02/B03/B04.</div></li>
          <li><span>E</span><div>Predict target visible BOA from target corrected anchors plus the same localizer. Aggregate predictions by per-pixel median across realizations. Prediction uncertainty is MAD x 1.4826, floored at 0.006.</div></li>
          <li><span>F</span><div>Add the fixed per-band debias <code>intercept + slope*AOT_anchor</code>: B02=(-0.0003,0.0243), B03=(-0.0006,0.0235), B04=(-0.0011,0.0223).</div></li>
        </ol>
        <div class="callout warning"><strong>Harness detail promoted to contract.</strong> The package function nominally imports <code>ExtraTreeRegressor</code>. The validated runner patches that symbol to the 20-tree forest above. A production implementation must configure the forest explicitly; relying on runtime monkeypatching is not an acceptable deployment design.</div>
      </section>

      <section id="rt" class="doc-section">
        <div class="section-heading"><p>07</p><div><h2>Radiative-transfer model</h2><p>Continental-average libRadtran LUT with Sentinel-2 spectral integration.</p></div></div>
        <div class="two-col">
          <div>
            <h3>LUT contents</h3>
            <p>The external Zarr v3 archive stores <code>TOA_rho1</code>, <code>TOA_rho2</code>, <code>Eg_rho1</code>, and <code>Eg_rho2</code> at Lambertian reference reflectances rho1=0.15 and rho2=0.5. Aerosol is one continental-average libRadtran profile for the entire tile. Final config <code>surface_driven_aerosol_species=none</code> disables CCI or per-pixel species selection.</p>
          </div>
          <div>
            <h3>Interpolation</h3>
            <p>AOI-mean SZA/VZA/RAA are snapped to their nearest geometry axes. The scene subset is then linearly interpolated in AOT and TCWV for each candidate. O3=300 DU and altitude=0 km are fixed. Dense 1 nm coefficients are convolved against the tabulated Sentinel-2 RSRF using normalized trapezoidal response weighting.</p>
          </div>
        </div>
        {_table(("Axis", "Range", "Spacing", "Nodes"), axis_rows)}
        <div class="formula-stack">
          <div class="formula"><span>Derived coefficients</span><code>s=(Eg2-Eg1)/(rho2*Eg2-rho1*Eg1); path=(TOA2*rho1*Eg1-TOA1*rho2*Eg2)/(rho1*Eg1-rho2*Eg2); t_up=(TOA2-TOA1)/(rho2*Eg2-rho1*Eg1)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
          <div class="formula"><span>Coefficient conversion</span><code>Eg0=Eg1*(1-rho1*s); t_total=max(Eg0*t_up,1e-10); xap=1/t_total; xbp=path/t_total; xcp=s</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
          <div class="formula"><span>TOA to BOA</span><code>y=xap*rho_TOA-xbp; rho_BOA=y/(1+xcp*y)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        </div>
        <p>The LUT archive is {LUT.stat().st_size:,} bytes. Its whole-file SHA-256 was not computed because the artifact is approximately 198 GB. Identity is pinned by URL, byte size, modification timestamp, Zarr root metadata SHA-256 <code>8c3344d510ea5368074ed60d6f533be88b0343188694626a02da047d07404d61</code>, and ZIP central-listing SHA-256 <code>5f4fd7ac3be9e46513fa738997395cba9bd65e9d3cfd3c9fb7e49c4876479fee</code>. A production release should add a complete archive digest or content-addressed object version.</p>
      </section>

      <section id="solver" class="doc-section">
        <div class="section-heading"><p>08</p><div><h2>Physical AOD solver</h2><p>Three-band chi-square inversion, spatial evidence pooling, and a MAIAC Gaussian backstop.</p></div></div>
        <h3>Candidate grid and observation cost</h3>
        <p>The AOT axis has 68 float32 nodes: 0.01 steps below 0.20, 0.025 steps from 0.20 to below 0.50, 0.05 steps from 0.50 to below 1.50, 0.10 steps from 1.50 through 2.50, then 0.25 steps from 2.75 through 4.00. <a href="downloads/aot-axis.csv">The exact decimal sequence is downloadable.</a></p>
        <div class="formula"><span>Per-pixel, per-node surface cost</span><code>J_obs(k,p) = sum_b ((rho_BOA(k,p,b) - rho_prior(p,b)) / max(sigma_prior(p,b), 0.006))^2, b in B02/B03/B04</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <p>At least one finite band is accepted by the generic cost builder; the validated records provide all three solve bands. Costs are float64 and do not include a 0.5 factor.</p>

        <h3>Pooling and posterior node</h3>
        <div class="formula"><span>Spatially pooled total cost</span><code>J_total(k,p) = rolling_median_20x20(J_obs(k), min_count=80) + (tau_k - m_MAIAC)^2 / sigma_backstop^2</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <ul class="method-list">
          <li><strong>Validity.</strong> Start from finite TOA/surface support and exclude cloud and buffered water. A pixel is solved only if the pooled cost is finite for every one of the 68 AOT nodes.</li>
          <li><strong>Rolling window.</strong> The centred window edge is <code>round(2*600/60)=20</code> cells, clamped at scene edges. Even-count medians average the two central values. At least <code>max(1,20,20*20/5)=80</code> finite samples are required per node.</li>
          <li><strong>Argmin.</strong> Choose the first minimum of total cost in ascending AOT order. Report AOT uncertainty as half the range of nodes satisfying <code>J_total &lt;= J_min + 0.5</code>, floored at 0.02.</li>
          <li><strong>Cost output.</strong> <code>jmin</code> and <code>final_cost</code> are observation-only pooled cost, not total cost including MAIAC. The scene final cost is the median finite <code>jmin</code>.</li>
          <li><strong>No smoothness solve.</strong> The selected path has one node sweep, no iterative optimizer, no multigrid regularization despite inherited config fields, no backstop escape, and no aerosol-species branch.</li>
        </ul>

        <h3>Cost-gated tau-dependent prior</h3>
        <p>Run the static anchored prior first. If <code>final_cost / 3 &gt; 8</code> and the forest payload exists, recompute target B02/B03/B04 prior values at every candidate AOT by correcting B8A/B11/B12 at that node and evaluating the same forests, localizer, and debias. Then rerun the entire inversion and return the tau result. The algorithm does not compare static and tau outcomes after firing. The gate fired in {target_summary['tau_gate_fired']} of 152 cases.</p>

        <h3>Recorded no-valid recovery</h3>
        <p>Strict cloud and water masks are the default. If the job returns <code>NO_VALID_OBSERVATION</code>, the array wrapper archives that receipt and retries once with both cloud and water exclusions bypassed. This deterministic recovery occurred for {len(mask_cases)} cases: {mask_case_list}. It is a fallback state, not a hidden change to the normal objective.</p>
      </section>

      <section id="iteration" class="doc-section">
        <div class="section-heading"><p>09</p><div><h2>One fixed-point rerun and scalar result</h2><p>Iteration updates the surface anchor only; it does not replace the atmospheric backstop.</p></div></div>
        <ol class="numbered-method">
          <li><span>1</span><div>Build the surface prior with scene anchors corrected at staged MAIAC AOT and complete the static/tau-gated physical solve.</div></li>
          <li><span>2</span><div>Take the arithmetic mean of all finite pass-1 AOT grid cells. If it is finite and positive, construct a constant anchor atmosphere at that AOT with TCWV=2, O3=0.30, elevation=0.</div></li>
          <li><span>3</span><div>Rebuild the same five-year bestpixel surface prior and visible prediction, but correct B8A/B11/B12 at the pass-1 mean. Do not change the MAIAC solver backstop.</div></li>
          <li><span>4</span><div>Repeat the full physical solve once and return pass 2. No convergence tolerance or additional iteration is used. All {target_summary['anchor_iteration_count']} target records contain both pass receipts.</div></li>
        </ol>
        <div class="formula"><span>Validated physical scalar</span><code>AOD_physical = mean(AOT_grid[p] for every finite AOI cell p after pass 2)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        <p>Station-nearest AOT and +/-1.5 km window median are recorded as diagnostics but are not the validated model input. BOA correction is disabled (<code>skip_correction=true</code>) for this benchmark; production BOA generation is a separate downstream mode and is not part of the 88.16% result.</p>
      </section>

      <section id="calibration" class="doc-section">
        <div class="section-heading"><p>10</p><div><h2>Generic post-retrieval calibration</h2><p>One frozen model applies uniformly to every retrieval after physical inversion.</p></div></div>
        <h3>CAMS operational context</h3>
        <p>Query <code>ECMWF/CAMS/NRT</code> from sensing time -3 hours to +3 hours, average images, then reduce by mean over a station-centred +/-0.06 degree box at 40 km. Required fields are total AOD at 469/550/670/865 nm and 550 nm black-carbon, dust, organic-matter components plus dust fraction. MERRA fields were collected during diagnosis but are explicitly excluded from the final feature extractor.</p>

        <h3>Training data and target</h3>
        <ul class="method-list">
          <li><strong>Training set.</strong> 588 complete external low-cloud physical retrievals from 234 sites. The external list uses catalog scene cloud cover below 20%. The ordered ID digest is <code>{_h(model_manifest['training']['matchup_id_sha256'])}</code>.</li>
          <li><strong>Feature schema.</strong> Build the sorted union of all operational feature names across training records. Missing values remain NaN. The frozen artifact stores all 137 names in order.</li>
          <li><strong>Leakage denylist.</strong> Feature names containing AERONET, truth, within_ee, ee_threshold, error, err, regime, runtime, or total_s fail extraction. Geography is disabled.</li>
          <li><strong>Regression target.</strong> <code>y = log((AOD_AERONET + 1/3) / (AOD_physical + 1/3))</code>.</li>
        </ul>

        <h3>Frozen estimator</h3>
        <div class="pipeline" aria-label="Calibration pipeline">
          <div><strong>137 inputs</strong><span>frozen schema order</span></div><b></b>
          <div><strong>Median imputer</strong><span>fit on 588 training rows</span></div><b></b>
          <div><strong>StandardScaler</strong><span>training mean and variance</span></div><b></b>
          <div><strong>SelectKBest</strong><span>f_regression, k=35</span></div><b></b>
          <div><strong>ExtraTrees</strong><span>1500 trees, L1 splits</span></div>
        </div>
        <ul class="fact-grid">
          <li><strong>Recipe</strong><span><code>{_h(recipe['name'])}</code></span></li>
          <li><strong>Trees</strong><span>{model_manifest['tree_count']}, random_state={model_manifest['random_state']}</span></li>
          <li><strong>Criterion</strong><span>absolute_error</span></li>
          <li><strong>Tree controls</strong><span>min_samples_leaf=4, max_features=0.5, n_jobs=4</span></li>
          <li><strong>Model versions</strong><span>NumPy {model_manifest['software']['numpy']}; scikit-learn {model_manifest['software']['scikit_learn']}; joblib {model_manifest['software']['joblib']}</span></li>
          <li><strong>Model SHA-256</strong><span class="hash-inline"><code>{_h(model_manifest['model_sha256'])}</code></span></li>
        </ul>
        <div class="formula-stack">
          <div class="formula"><span>Model reconstruction</span><code>AOD_raw = clip((AOD_physical + 1/3)*exp(y_hat) - 1/3, 0, 4)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
          <div class="formula"><span>Global correction</span><code>AOD_final = clip((AOD_raw + 1/3)*exp(-0.0125) - 1/3, 0, 4)</code><button type="button" data-copy aria-label="Copy formula">Copy</button></div>
        </div>
        <p>The recipe was ranked on target development folds by within-EE hits, then MAE and RMSE. The global offset was selected from -0.05 to +0.05 in 0.0025 steps using folds 0-3, ranked by hits, MAE, RMSE, then smallest absolute offset. These target labels are not used at inference, but they are part of the reported target-domain model-selection procedure.</p>
      </section>

      <section id="features" class="doc-section">
        <div class="section-heading"><p>11</p><div><h2>Selected 35 features</h2><p>The fitted artifact still expects the full 137-name schema; its frozen selector forwards these 35 columns.</p></div></div>
        <div class="filterbar">
          <label><span>Search features</span><input id="feature-search" type="search" placeholder="name or definition" autocomplete="off"></label>
          <label><span>Group</span><select id="feature-group"><option value="all">All groups</option><option>MAIAC</option><option>Contrast</option><option>CAMS</option><option>Solver</option><option>Time</option></select></label>
          <p id="feature-count">35 of 35 features</p>
        </div>
        <div class="table-wrap"><table class="feature-table"><thead><tr><th>#</th><th>Frozen name</th><th>Group</th><th>Unit</th><th>Exact operational definition</th><th>Importance</th></tr></thead><tbody>{''.join(feature_rows)}</tbody></table></div>
        <p class="footnote">Importance is the fitted ExtraTrees impurity importance and is descriptive, not a causal attribution. <a href="downloads/feature-schema.json">Download all 137 names</a> or <a href="downloads/selected-features.csv">the selected table as CSV</a>.</p>
      </section>

      <section id="validation" class="doc-section">
        <div class="section-heading"><p>12</p><div><h2>Validation result and robustness</h2><p>Measured performance, model-selection separation, and known generalization limits.</p></div></div>
        {_table(("Evaluation", "N", "Hits", "Within EE", "Interpretation"), validation_rows)}
        <div class="three-col metrics-detail">
          <div><span>Final RMSE</span><strong>{target_metrics['rmse']:.3f}</strong><p>Target 152</p></div>
          <div><span>Final MAE</span><strong>{target_metrics['mae']:.3f}</strong><p>Target 152</p></div>
          <div><span>Final bias</span><strong>{target_metrics['bias']:.3f}</strong><p>Target 152</p></div>
        </div>
        <ul class="method-list">
          <li><strong>Threshold.</strong> At least 133/152 hits are required to exceed 87%. The exported artifact has 134/152.</li>
          <li><strong>Cross-host artifact replay.</strong> Maximum absolute prediction difference from the selection receipts is {_fmt(model_manifest['validation_reproduction']['max_abs_prediction_delta'], 6)}, below the 0.005 tolerance, with zero within-EE classification disagreements.</li>
          <li><strong>Seed/tree robustness.</strong> Every ten tested full-schema seed/tree variant reached at least 133 hits. The frozen production variant is seed 20260713 with 1500 trees.</li>
          <li><strong>Transitions.</strong> Relative to physical retrieval there are 28 gains, 5 losses, 13 remaining misses, and 106 retained hits.</li>
          <li><strong>Unseen sites.</strong> Target records from sites absent in external training reach 46/53 = 86.79%, one hit below the target. External holdout reaches 83/123 = 67.48%.</li>
        </ul>
        <div class="callout warning"><strong>Selection caveat.</strong> The 88.16% target result includes recipe and global-offset choices made with labels from target folds 0-3. Fold 4 is held back from those choices, but the full 152 is not an untouched external test. Use the external and unseen-site results when assessing deployment risk.</div>
        <a class="command-link" href="../aod-final-performance-dashboard-20260713/"><span>Open detailed performance evidence</span><strong>Case explorer, spatial plots, spectra, costs, gains, losses, and experiment ledger</strong></a>
      </section>

      <section id="production" class="doc-section">
        <div class="section-heading"><p>13</p><div><h2>Production controls</h2><p>Determinism, orchestration, retries, monitoring, and failure semantics.</p></div></div>
        <h3>Execution contract</h3>
        <ul class="method-list">
          <li><strong>Idempotency.</strong> One JSON receipt per matchup ID. Existing <code>OK</code>, <code>NO_VALID_OBSERVATION</code>, or <code>DATA_UNAVAILABLE</code> receipts are terminal unless force-rerun is explicitly enabled.</li>
          <li><strong>Process retries.</strong> The array runner attempts a site up to three times. Resource/thread failures switch to single-thread mode. The SLURM wrapper also permits three process invocations with staggered delays.</li>
          <li><strong>Determinism.</strong> Set OMP, OpenBLAS, MKL, NumExpr, Rayon, Polars, and GDAL threads to one; set <code>PHASE_D_MAX_WORKERS=1</code> and <code>PHASE_D_SERIAL_POOLS=1</code>. The latter patches every <code>ThreadPoolExecutor</code> to one worker because concurrent GDAL warps moved flat-cost argmins.</li>
          <li><strong>Resource envelope.</strong> Validated SLURM allocation: 4 CPUs, 28 GB memory, 8-hour limit, at most eight concurrent array tasks. Observed runtime was {_fmt(target_summary['runtime_seconds']['min'], 1)}-{_fmt(target_summary['runtime_seconds']['max'], 1)} seconds, median {_fmt(target_summary['runtime_seconds']['median'], 1)} seconds per target scene.</li>
          <li><strong>Failure records.</strong> Every exception becomes a structured receipt with status, reason, type, traceback tail, and runtime. Missing output is not counted as success.</li>
          <li><strong>Credentials.</strong> GEE service-account and Earthdata credentials are runtime secrets. This public package records required variable names and source contracts but never embeds a private key.</li>
        </ul>
        <h3>Recommended service boundaries</h3>
        <div class="service-grid">
          <div><strong>Acquisition</strong><span>Resolve immutable S2 product ID; materialize SAFE and L2A source IDs; checksum cache.</span></div>
          <div><strong>Prior builder</strong><span>Materialize MAIAC day gate and S2 L2A composite stack with source receipts.</span></div>
          <div><strong>Physical worker</strong><span>Load frozen LUT, build one prior, run two passes, emit map plus diagnostics.</span></div>
          <div><strong>Context worker</strong><span>Fetch CAMS fields and enforce temporal/spatial completeness independently.</span></div>
          <div><strong>Calibrator</strong><span>Load model by SHA-256, validate 137 schema, reject nonfinite physical AOD.</span></div>
          <div><strong>Validator</strong><span>Reconcile counts, statuses, source IDs, config hash, model hash, and golden predictions.</span></div>
        </div>
        <details>
          <summary>Exact run environment</summary>
          <pre><code>{_h(chr(10).join(f'export {key}={value}' for key, value in EXACT_ENV.items()))}</code></pre>
        </details>
        <details>
          <summary>Reference pseudocode</summary>
          <pre><code>{_h(_pseudocode())}</code></pre>
        </details>
      </section>

      <section id="artifacts" class="doc-section">
        <div class="section-heading"><p>14</p><div><h2>Reproduction files</h2><p>Frozen artifacts and machine-readable contracts generated with this page.</p></div></div>
        <div class="table-wrap"><table><thead><tr><th>File</th><th>Purpose</th><th>Bytes</th><th>SHA-256</th></tr></thead><tbody>{artifact_rows}</tbody></table></div>
        <p>Use <a href="downloads/reproduction-manifest.json">reproduction-manifest.json</a> as the package inventory and <a href="build-receipt.json">build-receipt.json</a> as the outer receipt. The 198 GB LUT is linked rather than copied: <a href="{_h(LUT_PUBLIC_URL)}">public LUT archive</a>.</p>
      </section>

      <section id="source-map" class="doc-section">
        <div class="section-heading"><p>15</p><div><h2>Critical source map</h2><p>Exact tested bytes for implementation modules and campaign orchestration.</p></div></div>
        <div class="filterbar source-filter"><label><span>Search source</span><input id="source-search" type="search" placeholder="component or origin" autocomplete="off"></label><p id="source-count">{len(source_map)} of {len(source_map)} files</p></div>
        <div class="table-wrap"><table><thead><tr><th>Component</th><th>Origin</th><th>SHA-256</th><th>Bytes</th></tr></thead><tbody>{source_rows}</tbody></table></div>
      </section>

      <section id="limits" class="doc-section">
        <div class="section-heading"><p>16</p><div><h2>Known limits and unresolved production risks</h2><p>These are part of the handoff and should not be hidden by the headline metric.</p></div></div>
        <ul class="risk-list">
          <li><span>High</span><div><strong>External generalization is substantially lower.</strong><p>83/123 external holdout and 46/53 unseen target sites do not independently demonstrate 87% deployment accuracy.</p></div></li>
          <li><span>High</span><div><strong>Mutable upstream services prevent byte-identical refreshes.</strong><p>Planetary Computer STAC, bestpixel implementation, Earthdata granules, GEE reductions, and CAMS NRT collections can change. Preserve source IDs and materialized arrays for regulated replay.</p></div></li>
          <li><span>Medium</span><div><strong>The LUT lacks a full content digest in this campaign.</strong><p>Metadata and ZIP inventory identities are recorded, but release engineering should compute and publish a complete archive SHA-256 once.</p></div></li>
          <li><span>Medium</span><div><strong>Three cases require mask bypass.</strong><p>The recovery prevents job failure but permits cloud/water-contaminated evidence. Surface this flag in every downstream quality record.</p></div></li>
          <li><span>Medium</span><div><strong>Target-domain label selection is real.</strong><p>The global -0.0125 offset and recipe were selected on target development labels. Preserve fold 4 and acquire a new prospective test before claiming production accuracy.</p></div></li>
          <li><span>Low</span><div><strong>Forest anchor predictor is currently injected by a shim.</strong><p>Promote it to a typed config field and unit-test serialization before deployment.</p></div></li>
        </ul>
      </section>

      <section id="checklist" class="doc-section final-section">
        <div class="section-heading"><p>17</p><div><h2>Acceptance checklist</h2><p>A production reimplementation is equivalent only when every mandatory check passes.</p></div></div>
        <div class="checklist" data-checklist>
          <label><input type="checkbox"><span>Model SHA-256 equals <code>{_h(model_manifest['model_sha256'])}</code>.</span></label>
          <label><input type="checkbox"><span>Resolved config semantically matches <a href="downloads/resolved-config.json">resolved-config.json</a>; all inert fields are identified.</span></label>
          <label><input type="checkbox"><span>The 68-node float32 AOT axis matches <a href="downloads/aot-axis.csv">aot-axis.csv</a> exactly.</span></label>
          <label><input type="checkbox"><span>S2 product IDs, target cohort digest, MAIAC/CAMS context digests, and source IDs reconcile before execution.</span></label>
          <label><input type="checkbox"><span>All target physical receipts are terminal: 152 OK, zero missing, zero FAILED.</span></label>
          <label><input type="checkbox"><span>Physical replay returns 111/152 within EE and all 152 records contain anchor pass 1 and pass 2.</span></label>
          <label><input type="checkbox"><span>Artifact inference uses the 137-name schema and reproduces reference predictions within absolute tolerance 0.005.</span></label>
          <label><input type="checkbox"><span>Within-EE classifications agree for all 152 cases and final count is 134/152.</span></label>
          <label><input type="checkbox"><span>Mask bypass, gate failure, live fallback, missing CAMS, and model-schema errors are explicit quality states.</span></label>
          <label><input type="checkbox"><span>A prospective site-held-out validation is completed before extending the 88.16% claim beyond this cohort.</span></label>
        </div>
        <div class="generated">Generated {_h(generated_utc)}. Page URL: <a href="{_h(REPORT_PUBLIC_URL)}">{_h(REPORT_PUBLIC_URL)}</a></div>
      </section>
    </main>
  </div>
</body>
</html>
'''


def build(output: Path) -> dict[str, Any]:
    required = (
        MODEL_MANIFEST,
        MODEL,
        SELECTION_REPORT,
        OFFSET_REPORT,
        SEED_REPORT,
        ABLATION_REPORT,
        TARGET_MIDS,
        TARGET_COHORT_MANIFEST,
        MATCHUPS,
        TRAIN_MIDS,
        LUT,
        REPO_ROOT / "pixi.lock",
        REPO_ROOT / "pixi.toml",
    )
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)
    if len(AOT_AXIS) != 68 or AOT_AXIS[0] != 0.01 or AOT_AXIS[-1] != 4.0:
        raise ValueError("AOT axis contract is invalid")

    output.mkdir(parents=True, exist_ok=True)
    downloads = output / "downloads"
    assets = output / "assets"
    downloads.mkdir(exist_ok=True)
    assets.mkdir(exist_ok=True)

    target_ids = _read_ids(TARGET_MIDS)
    train_ids = _read_ids(TRAIN_MIDS)
    if len(target_ids) != 152 or len(set(target_ids)) != 152:
        raise ValueError(f"Expected 152 unique target IDs, found {len(target_ids)}")
    if len(train_ids) != 588 or len(set(train_ids)) != 588:
        raise ValueError(f"Expected 588 unique training IDs, found {len(train_ids)}")

    target_results = _records(TARGET_RESULTS, target_ids)
    target_maiac = _records(TARGET_MAIAC, target_ids)
    target_cams = _records(TARGET_CAMS, target_ids)
    train_results = _records(TRAIN_RESULTS, train_ids)
    train_maiac = _records(TRAIN_MAIAC, train_ids)
    train_cams = _records(TRAIN_CAMS, train_ids)
    if any(record.get("status") != "OK" for record in target_results):
        raise ValueError("Target physical result package is not 152/152 OK")
    if any(record.get("status") != "OK" for record in train_results):
        raise ValueError("Training physical result package is not 588/588 OK")
    if any(record.get("status") != "OK" for record in train_cams):
        raise ValueError("Training CAMS package is not 588/588 OK")

    model_manifest = _load_json(MODEL_MANIFEST)
    selection = _load_json(SELECTION_REPORT)
    offset = _load_json(OFFSET_REPORT)
    if _sha256(MODEL) != model_manifest["model_sha256"]:
        raise ValueError("Model digest does not match its manifest")
    if model_manifest["validation_reproduction"]["target"]["hits"] != 134:
        raise ValueError("Frozen model manifest no longer reports 134 target hits")

    training_digest = hashlib.sha256("\n".join(train_ids).encode("utf-8")).hexdigest()
    if training_digest != model_manifest["training"]["matchup_id_sha256"]:
        raise ValueError("Training ID order does not match the exported artifact")

    fieldnames, matchup_rows = _load_matchups()
    missing_matchups = [mid for mid in [*target_ids, *train_ids] if mid not in matchup_rows]
    if missing_matchups:
        raise ValueError(f"Missing {len(missing_matchups)} matchup rows")

    target_summary = _target_summary(target_results)
    if target_summary["physical_hits"] != 111:
        raise ValueError(f"Expected 111 physical hits, found {target_summary['physical_hits']}")
    maiac_summary = _maiac_summary(target_maiac, matchup_rows)

    _copy(MODEL, downloads / "aod-calibrator.joblib")
    _copy(MODEL_MANIFEST, downloads / "model-manifest.json")
    _copy(SELECTION_REPORT, downloads / "recipe-selection.json")
    _copy(OFFSET_REPORT, downloads / "global-offset-validation.json")
    _copy(SEED_REPORT, downloads / "seed-robustness.json")
    _copy(ABLATION_REPORT, downloads / "feature-ablation.json")
    _copy(TARGET_MIDS, downloads / "campaign250-lowcloud20-mids.txt")
    _copy(TARGET_COHORT_MANIFEST, downloads / "cohort-manifest.json")
    _copy(TRAIN_MIDS, downloads / "external-training-ids.txt")
    _copy(REPO_ROOT / "pixi.lock", downloads / "pixi.lock")
    _copy(REPO_ROOT / "pixi.toml", downloads / "pixi.toml")
    _copy(
        REPO_ROOT / "tools/aeronet_validation/verify_aod_reproduction_package.py",
        downloads / "verify-reproduction.py",
    )
    _copy(FINAL_DASHBOARD / "data/all-cases.csv", downloads / "reference-predictions.csv")

    _write_jsonl(downloads / "target-physical-results.jsonl", target_results)
    _write_jsonl(downloads / "target-maiac-inputs.jsonl", target_maiac)
    _write_jsonl(downloads / "target-cams-context.jsonl", map(_sanitize_context, target_cams))
    _write_jsonl(downloads / "training-physical-results.jsonl", train_results)
    _write_jsonl(downloads / "training-maiac-inputs.jsonl", train_maiac)
    _write_jsonl(downloads / "training-cams-context.jsonl", map(_sanitize_context, train_cams))
    _write_matchup_subset(downloads / "target-matchups.csv", fieldnames, matchup_rows, target_ids)
    _write_matchup_subset(downloads / "training-matchups.csv", fieldnames, matchup_rows, train_ids)

    _write_json(downloads / "resolved-config.json", _resolved_config())
    _write_json(downloads / "environment.json", _environment_payload())
    _write_json(downloads / "lut-axes.json", LUT_AXES)
    (downloads / "production-recipe.toml").write_text(_recipe_toml(), encoding="utf-8")
    (downloads / "algorithm-pseudocode.txt").write_text(_pseudocode(), encoding="utf-8")

    _write_csv(
        downloads / "aot-axis.csv",
        ("index", "aot_550"),
        [{"index": index, "aot_550": f"{value:.6f}"} for index, value in enumerate(AOT_AXIS)],
    )
    schema = list(model_manifest["features"]["schema"])
    selected = list(model_manifest["features"]["selected"])
    _write_json(downloads / "feature-schema.json", {"count": len(schema), "features": schema})
    selected_rows = []
    importance = {
        row["feature"]: row["importance"] for row in model_manifest["features"]["importance"]
    }
    for index, name in enumerate(selected, start=1):
        group, unit, definition = _feature_definition(name)
        selected_rows.append(
            {
                "index": index,
                "name": name,
                "group": group,
                "unit": unit,
                "definition": definition,
                "importance": importance[name],
            }
        )
    _write_csv(
        downloads / "selected-features.csv",
        ("index", "name", "group", "unit", "definition", "importance"),
        selected_rows,
    )

    input_contract = {
        "schema_version": "siac-aod-production-input-contract-v1",
        "target_cohort": {
            "count": 152,
            "cloud_field": "raw OmniCloudMask cloud_frac",
            "predicate": "cloud_frac < 0.20",
            "mids_sha256": _sha256(TARGET_MIDS),
        },
        "training": {
            "count": 588,
            "site_count": 234,
            "ordered_ids_sha256": training_digest,
            "catalog_cloud_predicate": "scene_cloud_cover < 20 percent",
        },
        "sentinel2": {
            "processing_level": "L1C",
            "solve_bands": ["B02", "B03", "B04"],
            "anchor_bands": ["B8A", "B11", "B12"],
            "reflectance_formula": "clip((DN + RADIO_ADD_OFFSET) / QUANTIFICATION_VALUE, 0, 1.5)",
        },
        "surface": {
            "source_type": "Sentinel-2 L2A bestpixel only",
            "endpoint": "Planetary Computer",
            "lookback_years": 5,
            "months_relative": [-1, 0, 1],
            "low_aod_day_fraction": 0.6,
        },
        "atmosphere": {
            "aot": "MCD19A2 QA-best native median nearest day +/-2",
            "tcwv_cm": 2.0,
            "ozone_atm_cm": 0.3,
            "elevation_km": 0.0,
        },
        "rt": {
            "artifact": str(LUT),
            "public_url": LUT_PUBLIC_URL,
            "byte_size": LUT.stat().st_size,
            "zarr_root_metadata_sha256": "8c3344d510ea5368074ed60d6f533be88b0343188694626a02da047d07404d61",
            "zip_listing_sha256": "5f4fd7ac3be9e46513fa738997395cba9bd65e9d3cfd3c9fb7e49c4876479fee",
        },
        "calibration": {
            "model_sha256": model_manifest["model_sha256"],
            "feature_schema_count": 137,
            "selected_feature_count": 35,
        },
    }
    _write_json(downloads / "input-contract.json", input_contract)

    source_map = _source_map()
    _write_json(downloads / "source-map.json", source_map)
    _write_source_archive(downloads / "critical-source-snapshot.tar.gz", source_map)

    shutil.copy2(WEB_ASSETS / "app.css", assets / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", assets / "app.js")

    purpose = {
        "aod-calibrator.joblib": "Frozen deployable GenericAODCalibrator artifact.",
        "model-manifest.json": "Model identity, schema, selected features, and replay metrics.",
        "resolved-config.json": "Full SIAC configuration after preset and environment overrides.",
        "production-recipe.toml": "Compact normalized operational recipe and environment.",
        "input-contract.json": "Machine-readable source, unit, and cohort contract.",
        "algorithm-pseudocode.txt": "Implementation-neutral end-to-end pseudocode.",
        "feature-schema.json": "All 137 frozen model input names in order.",
        "selected-features.csv": "Selected 35 features with definitions and importance.",
        "aot-axis.csv": "Exact 68-node physical AOT search grid.",
        "lut-axes.json": "Radiative-transfer LUT axes and node values.",
        "reference-predictions.csv": "152 physical, unadjusted, and final scalar predictions.",
        "target-physical-results.jsonl": "Full 152 physical receipts and solver diagnostics.",
        "target-maiac-inputs.jsonl": "Frozen 152 MAIAC staging records.",
        "target-cams-context.jsonl": "Frozen 152 operational CAMS contexts without labels.",
        "target-matchups.csv": "Fixed target product, time, station, and truth rows.",
        "campaign250-lowcloud20-mids.txt": "Ordered 152-case target denominator.",
        "cohort-manifest.json": "Cloud cohort provenance and digest.",
        "training-physical-results.jsonl": "Full 588 external physical training receipts.",
        "training-maiac-inputs.jsonl": "Frozen 588 external MAIAC staging records.",
        "training-cams-context.jsonl": "Frozen 588 external CAMS contexts without labels.",
        "training-matchups.csv": "Fixed 588 external product, time, and target rows.",
        "external-training-ids.txt": "Ordered training IDs matching the model digest.",
        "recipe-selection.json": "Target-development recipe selection audit.",
        "global-offset-validation.json": "Offset selection, nested replay, and case receipts.",
        "seed-robustness.json": "Seed and tree-count stability evaluation.",
        "feature-ablation.json": "Operational feature family ablation.",
        "critical-source-snapshot.tar.gz": "Exact algorithm and harness source bytes.",
        "source-map.json": "Per-file source origin and SHA-256 identity.",
        "environment.json": "Python, platform, package, git, and deterministic runtime identity.",
        "pixi.lock": "Resolved dependency lock.",
        "pixi.toml": "Project environment declaration.",
        "verify-reproduction.py": "Offline checksum, feature, model, and 134/152 replay verifier.",
    }
    artifacts = []
    for path in sorted(downloads.iterdir()):
        if path.name == "reproduction-manifest.json":
            continue
        artifacts.append(
            {
                "name": path.name,
                "href": f"downloads/{path.name}",
                "purpose": purpose.get(path.name, "Frozen reproduction artifact."),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )

    generated_utc = datetime.now(timezone.utc).isoformat()
    index = _render_html(
        model_manifest=model_manifest,
        selection=selection,
        offset=offset,
        target_summary=target_summary,
        maiac_summary=maiac_summary,
        source_map=source_map,
        artifacts=artifacts,
        generated_utc=generated_utc,
    )
    (output / "index.html").write_text(index, encoding="utf-8")

    inventory = {}
    for path in sorted(output.rglob("*")):
        if not path.is_file() or path.name in {"reproduction-manifest.json", "build-receipt.json"}:
            continue
        relative = str(path.relative_to(output))
        inventory[relative] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
    reproduction_manifest = {
        "schema_version": "siac-aod-production-reproduction-package-v1",
        "generated_utc": generated_utc,
        "public_url": REPORT_PUBLIC_URL,
        "target_result": model_manifest["validation_reproduction"]["target"],
        "physical_result": {
            "count": 152,
            "hits": 111,
            "within_ee_percent": 100 * 111 / 152,
        },
        "model_sha256": model_manifest["model_sha256"],
        "files": inventory,
        "external_artifacts": [input_contract["rt"]],
    }
    _write_json(downloads / "reproduction-manifest.json", reproduction_manifest)
    receipt = {
        "schema_version": "siac-aod-production-spec-build-v1",
        "generated_utc": generated_utc,
        "public_url": REPORT_PUBLIC_URL,
        "index_sha256": _sha256(output / "index.html"),
        "manifest_sha256": _sha256(downloads / "reproduction-manifest.json"),
        "model_sha256": model_manifest["model_sha256"],
        "target_records": len(target_results),
        "training_records": len(train_results),
        "critical_source_files": len(source_map),
    }
    _write_json(output / "build-receipt.json", receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    receipt = build(args.output)
    print(json.dumps(receipt, indent=2))
    print(args.output)
    print(REPORT_PUBLIC_URL)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/bin/bash
#SBATCH --job-name=c250maskmodis
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-264%2
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250maskmodis_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250maskmodis_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="${PHASE_D_MIDS_PATH:-$ROOT/masked_r2_regression_diagnostic_mids.txt}"
N_SITES="${PHASE_D_N_SITES:-$(sed '/^[[:space:]]*$/d' "$MIDS" | wc -l)}"
OPTIONS=(
  modis_single_day
  modis_timeseries_smooth
  modis_monthly_swir_nir_anchored
  s2_swir_nir_anchored
)

if [[ -n "${PHASE_D_MODIS_OPTION:-}" ]]; then
  SITE_INDEX=$SLURM_ARRAY_TASK_ID
  OPTION=$PHASE_D_MODIS_OPTION
else
  OPTION_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) / N_SITES))
  SITE_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) % N_SITES + 1))
  OPTION=${OPTIONS[$OPTION_INDEX]}
fi
MID=$(sed -n "${SITE_INDEX}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for option=$OPTION site index=$SITE_INDEX" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
export POLARS_MAX_THREADS=1
export GDAL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
export MAIAC_AGG=median
export TMPDIR="$ROOT/spx_tmp_masked_modis_prior_followups/${OPTION}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_campaign250_masked_options"

export PHASE_D_CONFIG_PRESET=example
export PHASE_D_BRDF_KIND=mcd43_gee
export PHASE_D_BRDF_CACHE_DIR="$ROOT/phaseD_cache/mcd43_gee_masked_followups"
export PHASE_D_BRDF_DATA_PATH=/home/users/marcyin/SIAC/tools/edown_runtime
export PHASE_D_BRDF_TEMPORAL_WINDOW=16
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_COST_MODE=auto2
export PHASE_D_SOLVE_BANDS=B02,B04
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MAIAC_DIR=maiac_qa
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_TARGET_TCWV_MODE=prior
export PHASE_D_TARGET_TCWV_FIXED=2.0
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_RT_BACKEND=lut
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_MIDS_FILE="$MIDS"
export PHASE_D_INDEX=$SITE_INDEX
unset PHASE_D_TAU_GATE PHASE_D_TAU_PRIOR PHASE_D_ANCHOR_ITERATE
unset PHASE_D_MONTHLY_COMPOSITES_KIND
unset PHASE_D_MONTHLY_COMPOSITE_UNCERTAINTY_SCALE

"$PHASE_D_BRDF_DATA_PATH" --help >/dev/null

case "$OPTION" in
  modis_single_day)
    export PHASE_D_SURFACE_PRIOR_METHOD=kernel_model
    ;;
  modis_timeseries_smooth)
    export PHASE_D_SURFACE_PRIOR_METHOD=whittaker
    ;;
  modis_monthly_swir_nir_anchored)
    export PHASE_D_SURFACE_PRIOR_METHOD=monthly_database
    export PHASE_D_MONTHLY_COMPOSITES_KIND=generated_brdf
    ;;
  modis_monthly_swir_nir_anchored_loose_backstop)
    export PHASE_D_SURFACE_PRIOR_METHOD=monthly_database
    export PHASE_D_MONTHLY_COMPOSITES_KIND=generated_brdf
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  modis_monthly_swir_nir_anchored_spread_only)
    export PHASE_D_SURFACE_PRIOR_METHOD=monthly_database
    export PHASE_D_MONTHLY_COMPOSITES_KIND=generated_brdf
    export PHASE_D_MONTHLY_COMPOSITE_UNCERTAINTY_SCALE=0
    ;;
  s2_swir_nir_anchored|swir_nir_anchored)
    export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
    export PHASE_D_BESTPIXEL_ENDPOINT=pc
    export PHASE_D_BESTPIXEL_SOURCE=l2a
    export PHASE_D_USE_PIXI_BESTPIXEL=1
    export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
    export PHASE_D_PREDICTOR_MODEL=et20
    export PHASE_D_ANCHOR_ITERATE=1
    ;;
  *)
    echo "Unknown option: $OPTION" >&2
    exit 2
    ;;
esac

RESULTS_SUFFIX="${PHASE_D_RESULTS_SUFFIX:-20260711}"
export PHASE_D_RESULTS_DIR="phaseD_results_masked_modis_prior_${OPTION}_${RESULTS_SUFFIX}"
mkdir -p "$TMPDIR" "$PHASE_D_BRDF_CACHE_DIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"

echo "TASK $SLURM_ARRAY_TASK_ID -> option=$OPTION site_index=$SITE_INDEX mid=$MID"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
attempt=1
max_attempts=3
while true; do
  status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || status=$?
  if (( status == 0 )) && [[ -f "$RESULT_JSON" ]]; then
    if "$PY" - "$RESULT_JSON" <<'PY'
import json
import sys

path = sys.argv[1]
try:
    record = json.loads(open(path, encoding="utf-8").read())
except Exception:
    raise SystemExit(1)
status = str(record.get("status", "")).upper()
reason = " ".join(
    str(record.get(key, "")) for key in ("reason", "error_type", "traceback_tail")
)
transient = (
    status == "FAILED"
    and (
        "can't start new thread" in reason
        or "Resource temporarily unavailable" in reason
        or "ThreadPoolBuildError" in reason
    )
)
raise SystemExit(0 if transient else 1)
PY
    then
      status=75
      echo "Transient failed record detected in $RESULT_JSON" >&2
    else
      exit 0
    fi
  elif (( status == 0 )); then
    exit 0
  fi
  if (( attempt >= max_attempts )); then
    echo "Task failed after $attempt attempts: option=$OPTION site_index=$SITE_INDEX mid=$MID status=$status" >&2
    exit "$status"
  fi
  sleep_for=$((60 + (SLURM_ARRAY_TASK_ID % 60) + (attempt - 1) * 120))
  echo "Attempt $attempt failed with status=$status; retrying in ${sleep_for}s" >&2
  sleep "$sleep_for"
  attempt=$((attempt + 1))
done

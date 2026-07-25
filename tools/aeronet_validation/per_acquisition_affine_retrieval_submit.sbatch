#!/bin/bash
#SBATCH --job-name=l2apairret
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-36%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairret_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairret_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${PAIR_AFFINE_CASES:-$ROOT/reports/aod-medium-physics-20260713/downloads/cases.csv}
CASE_INDEX=${PAIR_AFFINE_CASE_INDEX:-$SLURM_ARRAY_TASK_ID}
MID=$(sed -n "$((CASE_INDEX + 1))p" "$CASES" | cut -d, -f1)
HISTORY_ROOT=${PAIR_AFFINE_HISTORY_ROOT:-analysis/l2a_l1c_per_acquisition_affine_20260715/daily_histories}
HISTORY_TAG=${PAIR_AFFINE_HISTORY_TAG:?set PAIR_AFFINE_HISTORY_TAG to a mapped-history directory name}
RUN_TAG=${PAIR_AFFINE_RUN_TAG:-${HISTORY_TAG}_physical_static}
RESULTS_DIR=${PAIR_AFFINE_RESULTS_DIR:-phaseD_results_lowcloud20_pair_affine_${RUN_TAG}_mediumdev_20260715}
if [[ -z "${MID:-}" || "$MID" == "matchup_id" ]]; then
  echo "No matchup for task $CASE_INDEX in $CASES" >&2
  exit 2
fi
HISTORY_FILE="$ROOT/$HISTORY_ROOT/$HISTORY_TAG/${MID}.npz"
PRIOR_QUALITY_FILE="$ROOT/prior_quality/${MID}.json"
if [[ ! -f "$HISTORY_FILE" || ! -f "$PRIOR_QUALITY_FILE" ]]; then
  echo "Required seasonal inputs are missing: history=$HISTORY_FILE quality=$PRIOR_QUALITY_FILE" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT=python-gee@gee-marc.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/home/users/marcyin/gee-service-account.json
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
export MALLOC_CONF=background_thread:false,narenas:1
export MAIAC_AGG=median

# The mapped history remains the single L2A prior source.  Every target-scene
# RT input is shared by the anchor correction and the AOD solve; only AOD is
# swept by the solver.
export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_SEASONAL_PREDICTOR=1
export PHASE_D_SEASONAL_DIR="$HISTORY_ROOT/$HISTORY_TAG"
export PHASE_D_SEASONAL_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_BANDS=B02,B03,B04
export PHASE_D_SOLVE_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICT_VISIBLE_ANCHOR_SOURCE=atmo_prior
export PHASE_D_SEASONAL_ANCHOR_AOD_SOURCE=atmo_prior
export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
export PHASE_D_SEASONAL_DEBIAS=0
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_BANDS=B8A,B11,B12
export PHASE_D_ANCHOR_ITERATE=0
export PHASE_D_TAU_PRIOR=0
unset PHASE_D_TAU_GATE
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_adaptive
export PHASE_D_TARGET_TCWV_MODE=l2a
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_REFERENCE_TCWV=prior
export PHASE_D_DEM=glo30
export PHASE_D_DEM_PRE_M3=1
export PHASE_D_SCENE_MEAN_GEOMETRY=0
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_COST_MODE=chi2
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export PHASE_D_MATCHUP_ID="$MID"
export PHASE_D_RESULTS_DIR="$RESULTS_DIR"

export TMPDIR="$ROOT/spx_tmp_lowcloud20_pair_affine/$RUN_TAG/${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}"
mkdir -p "$TMPDIR" "$ROOT/$RESULTS_DIR"
RESULT_JSON="$ROOT/$RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
if [[ -f "$RESULT_JSON" && "${PAIR_AFFINE_FORCE_RERUN:-0}" != "1" ]]; then
  status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$RESULT_JSON" 2>/dev/null || echo malformed)
  if [[ "$status" == "ok" || "$status" == "data_unavailable" ]]; then
    echo "Existing terminal result: history=$HISTORY_TAG mid=$MID status=$status"
    exit 0
  fi
fi

for attempt in 1 2 3; do
  process_status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || process_status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$RESULT_JSON" 2>/dev/null || echo malformed)
    if [[ "$status" == "no_valid_observation" && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_pair_affine_no_valid/$RUN_TAG"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}.json"
      export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
      export PHASE_D_IGNORE_CLOUD_WATER=1
      status=mask_fallback
      continue
    fi
    if [[ "$status" == "ok" || "$status" == "data_unavailable" ]]; then
      exit 0
    fi
  else
    status=missing
  fi
  if (( attempt < 3 )); then
    sleep_for=$((60 + CASE_INDEX % 60 + (attempt - 1) * 120))
    echo "Attempt $attempt ended process=$process_status result=$status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done
echo "Retrieval failed after three attempts: history=$HISTORY_TAG mid=$MID result=$status" >&2
exit 1

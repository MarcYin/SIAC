#!/bin/bash
#SBATCH --job-name=l2aharmret
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=28G
#SBATCH --array=1-36%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aharmret_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aharmret_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=$ROOT/reports/aod-medium-physics-20260713/downloads/cases.csv
CASE_INDEX=${HARM_CASE_INDEX:-$SLURM_ARRAY_TASK_ID}
MID=$(sed -n "$((CASE_INDEX + 1))p" "$CASES" | cut -d, -f1)
HISTORY_TAG=${HARM_HISTORY_TAG:?set HARM_HISTORY_TAG to a daily-history directory name}
HISTORY_ROOT=${HARM_HISTORY_ROOT:-analysis/l2a_l1c_harmonizer_mediumdev_20260713/daily_histories}
DEBIAS=${HARM_DEBIAS:-off}
DIRECT_HISTORY=${HARM_DIRECT_HISTORY:-0}
TAU_MODE=${HARM_TAU_MODE:-off}
BACKSTOP_CAL=${HARM_BACKSTOP_CAL:-1}
VISIBLE_FLOOR=${HARM_VISIBLE_FLOOR:-0.006}
DUMP_COST=${HARM_DUMP_COST:-0}
DEM=${HARM_DEM:-none}
DEM_PRE_M3=${HARM_DEM_PRE_M3:-0}
RUN_TAG=${HARM_RUN_TAG:-${HISTORY_TAG}_${DEBIAS}_${TAU_MODE}}
if [[ -z "${MID:-}" ]]; then
  echo "No matchup for task $CASE_INDEX" >&2
  exit 2
fi
if [[ "$DEBIAS" != "on" && "$DEBIAS" != "off" ]]; then
  echo "HARM_DEBIAS must be 'on' or 'off', got '$DEBIAS'" >&2
  exit 2
fi
if [[ "$DIRECT_HISTORY" != "0" && "$DIRECT_HISTORY" != "1" ]]; then
  echo "HARM_DIRECT_HISTORY must be 0 or 1, got '$DIRECT_HISTORY'" >&2
  exit 2
fi
if [[ "$TAU_MODE" != "off" && "$TAU_MODE" != "gated8" && "$TAU_MODE" != "always" ]]; then
  echo "HARM_TAU_MODE must be off, gated8, or always; got '$TAU_MODE'" >&2
  exit 2
fi
if [[ "$BACKSTOP_CAL" != "0" && "$BACKSTOP_CAL" != "1" ]]; then
  echo "HARM_BACKSTOP_CAL must be 0 or 1, got '$BACKSTOP_CAL'" >&2
  exit 2
fi
if ! awk -v value="$VISIBLE_FLOOR" 'BEGIN { exit !(value + 0 > 0 && value + 0 < 0.1) }'; then
  echo "HARM_VISIBLE_FLOOR must be numeric and in (0, 0.1), got '$VISIBLE_FLOOR'" >&2
  exit 2
fi
if [[ "$DUMP_COST" != "0" && "$DUMP_COST" != "1" ]]; then
  echo "HARM_DUMP_COST must be 0 or 1, got '$DUMP_COST'" >&2
  exit 2
fi
if [[ "$DEM" != "none" && "$DEM" != "glo30" ]]; then
  echo "HARM_DEM must be 'none' or 'glo30', got '$DEM'" >&2
  exit 2
fi
if [[ "$DEM_PRE_M3" != "0" && "$DEM_PRE_M3" != "1" ]]; then
  echo "HARM_DEM_PRE_M3 must be 0 or 1, got '$DEM_PRE_M3'" >&2
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

export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_SEASONAL_PREDICTOR=1
export PHASE_D_SEASONAL_DIR="$HISTORY_ROOT/$HISTORY_TAG"
export PHASE_D_SEASONAL_BANDS=B02,B03,B04
export PHASE_D_SEASONAL_ANCHOR_AOD_SOURCE=maiac
export PHASE_D_SEASONAL_DEBIAS=$([[ "$DEBIAS" == "on" ]] && echo 1 || echo 0)
export PHASE_D_PREDICT_VISIBLE_BANDS=B02,B03,B04
export PHASE_D_SOLVE_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_FLOOR=$VISIBLE_FLOOR
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_ITERATE=0
export PHASE_D_TAU_PRIOR=0
unset PHASE_D_TAU_GATE
export PHASE_D_SEASONAL_TAU_PREDICTOR=0
if [[ "$TAU_MODE" == "gated8" ]]; then
  export PHASE_D_TAU_PRIOR=1
  export PHASE_D_TAU_GATE=8.0
  export PHASE_D_SEASONAL_TAU_PREDICTOR=1
elif [[ "$TAU_MODE" == "always" ]]; then
  export PHASE_D_TAU_PRIOR=1
  export PHASE_D_SEASONAL_TAU_PREDICTOR=1
fi
export PHASE_D_BACKSTOP_CAL=$BACKSTOP_CAL
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_adaptive
export PHASE_D_TARGET_TCWV_MODE=prior
export PHASE_D_REFERENCE_TCWV=2.0
export PHASE_D_DEM=$DEM
export PHASE_D_DEM_PRE_M3=$DEM_PRE_M3
export PHASE_D_PREDICT_VISIBLE_DEBIAS=$DEBIAS
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_COST_MODE=chi2
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export PHASE_D_MATCHUP_ID=$MID
export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_l2aharm_${RUN_TAG}_mediumdev_20260713"
export PHASE_D_DUMP_COST_CUBE=$DUMP_COST
export PHASE_D_COST_CUBE_DIR="phaseD_cost_cubes_lowcloud20_l2aharm_${RUN_TAG}_mediumdev_20260713"

RUNNER="$ROOT/phaseD_array_runner.py"
if [[ "$DIRECT_HISTORY" == "1" ]]; then
  export PHASE_D_PREBUILT_DIR="$HISTORY_ROOT/$HISTORY_TAG"
  export PHASE_D_PREBUILT_BANDS=B01,B02,B03,B04
  RUNNER=/home/users/marcyin/SIAC/tools/aeronet_validation/run_harmonized_retrieval.py
fi

export TMPDIR="$ROOT/spx_tmp_lowcloud20_l2aharm/$RUN_TAG/${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}"
mkdir -p "$TMPDIR"
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
if [[ -f "$RESULT_JSON" && "${PHASE_D_FORCE_RERUN:-0}" != "1" ]]; then
  status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$RESULT_JSON" 2>/dev/null || echo malformed)
  if [[ "$status" == "ok" || "$status" == "data_unavailable" ]]; then
    echo "Existing terminal result: history=$HISTORY_TAG mid=$MID status=$status"
    exit 0
  fi
fi

for attempt in 1 2 3; do
  process_status=0
  "$PY" "$RUNNER" || process_status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$RESULT_JSON" 2>/dev/null || echo malformed)
    if [[ "$status" == "no_valid_observation" && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_l2aharm_no_valid/$RUN_TAG"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
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

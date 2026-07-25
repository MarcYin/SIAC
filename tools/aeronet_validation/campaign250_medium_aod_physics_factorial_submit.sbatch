#!/bin/bash
#SBATCH --job-name=medaodphys1
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=28G
#SBATCH --array=1-608%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodphys1_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodphys1_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS_FILE=campaign250_lowcloud20_mids.txt
MIDS="$ROOT/$MIDS_FILE"
N_CASES=152
VARIANT_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / N_CASES ))
CASE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % N_CASES + 1 ))
MID=$(sed -n "${CASE_INDEX}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for case ${CASE_INDEX}" >&2
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

export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_PREDICT_VISIBLE_BANDS=B02,B03,B04
export PHASE_D_SOLVE_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_ITERATE=1
export PHASE_D_ANCHOR_ITERATE_MODE=anchor
export PHASE_D_TAU_PRIOR=1
export PHASE_D_TAU_GATE=8.0
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_adaptive
export PHASE_D_TARGET_TCWV_MODE=prior
export PHASE_D_REFERENCE_TCWV=2.0
export PHASE_D_DEM=none
export PHASE_D_PREDICT_VISIBLE_DEBIAS=on
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_COST_MODE=chi2
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export PHASE_D_MIDS_FILE="$MIDS_FILE"
export PHASE_D_INDEX=$CASE_INDEX

case "$VARIANT_INDEX" in
  0)
    VARIANT=dem
    export PHASE_D_DEM=glo30
    ;;
  1)
    VARIANT=dem_l2awvp
    export PHASE_D_DEM=glo30
    export PHASE_D_TARGET_TCWV_MODE=l2a
    export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
    export PHASE_D_REFERENCE_TCWV=prior
    ;;
  2)
    VARIANT=product_unc
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  3)
    VARIANT=no_surface_debias
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    ;;
  *)
    echo "Unknown variant index ${VARIANT_INDEX}" >&2
    exit 2
    ;;
esac

export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_mediumphysics_${VARIANT}_20260713"
export TMPDIR="$ROOT/spx_tmp_lowcloud20_mediumphysics/${VARIANT}/${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}"
mkdir -p "$TMPDIR"
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

if [[ -f "$RESULT_JSON" && "${PHASE_D_FORCE_RERUN:-0}" != "1" ]]; then
  existing_status=$(
    "$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
      "$RESULT_JSON" 2>/dev/null || echo malformed
  )
  if [[ "$existing_status" == "ok" || "$existing_status" == "data_unavailable" ]]; then
    echo "Existing terminal result: variant=$VARIANT mid=$MID status=$existing_status"
    exit 0
  fi
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> variant=$VARIANT case=$CASE_INDEX mid=$MID"
for attempt in 1 2 3; do
  process_status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || process_status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    result_status=$(
      "$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
        "$RESULT_JSON" 2>/dev/null || echo malformed
    )
    if [[ "$result_status" == "no_valid_observation" \
      && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_mediumphysics_no_valid/$VARIANT"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" \
        "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
      export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
      export PHASE_D_IGNORE_CLOUD_WATER=1
      result_status=mask_fallback
      echo "No valid masked observation; retrying with recorded mask fallback" >&2
      continue
    fi
    if [[ "$result_status" == "ok" || "$result_status" == "data_unavailable" ]]; then
      exit 0
    fi
  else
    result_status=missing
  fi
  if (( attempt < 3 )); then
    sleep_for=$((60 + CASE_INDEX % 60 + (attempt - 1) * 120))
    echo "Attempt $attempt ended with process=$process_status result=$result_status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "Retrieval failed after three attempts: variant=$VARIANT mid=$MID result=$result_status" >&2
exit 1

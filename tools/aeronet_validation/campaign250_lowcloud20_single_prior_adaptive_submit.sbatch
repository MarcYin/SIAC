#!/bin/bash
#SBATCH --job-name=lc20s2adapt
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=28G
#SBATCH --array=1-152%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20s2adapt_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20s2adapt_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS_FILE="${SINGLE_PRIOR_MIDS_FILE:-campaign250_lowcloud20_mids.txt}"
MIDS="$ROOT/$MIDS_FILE"
OPTION="${SINGLE_PRIOR_OPTION:-b03_loo_s10}"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for task ${SLURM_ARRAY_TASK_ID}" >&2
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
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export TMPDIR="$ROOT/spx_tmp_lowcloud20_single_prior/${OPTION}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_PREDICT_VISIBLE_BANDS=B02,B03,B04
export PHASE_D_SOLVE_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_ITERATE=1
export PHASE_D_TAU_PRIOR=1
export PHASE_D_TAU_GATE=8.0
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_MAIAC_DIR="${PHASE_D_MAIAC_DIR:-maiac_qa}"
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
if [[ "${PHASE_D_FORCE_MASK_FALLBACK:-0}" == "1" ]]; then
  export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
  export PHASE_D_IGNORE_CLOUD_WATER=1
fi
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_MIDS_FILE="$MIDS_FILE"
export PHASE_D_INDEX=$SLURM_ARRAY_TASK_ID

case "$OPTION" in
  b03_chi2)
    export PHASE_D_COST_MODE=chi2
    ;;
  b03_shape)
    export PHASE_D_COST_MODE=shape
    ;;
  b03_auto2)
    export PHASE_D_COST_MODE=auto2
    ;;
  b03_profile_s05)
    export PHASE_D_COST_MODE=profile_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.05
    ;;
  b03_profile_s10)
    export PHASE_D_COST_MODE=profile_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.10
    ;;
  b03_profile_s20)
    export PHASE_D_COST_MODE=profile_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.20
    ;;
  b03_loo_s05)
    export PHASE_D_COST_MODE=loo_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.05
    ;;
  b03_loo_s10)
    export PHASE_D_COST_MODE=loo_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.10
    ;;
  b03_loo_s20)
    export PHASE_D_COST_MODE=loo_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.20
    ;;
  b03_trimmed)
    export PHASE_D_COST_MODE=trimmed_chi2
    ;;
  b03_trimmed_loose)
    export PHASE_D_COST_MODE=trimmed_chi2
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  b03_trimmed_bs3)
    export PHASE_D_COST_MODE=trimmed_chi2
    export PHASE_D_BACKSTOP_CAL=0
    export PHASE_D_BACKSTOP_UNC_SCALE=3.0
    ;;
  b03_trimmed_bs10)
    export PHASE_D_COST_MODE=trimmed_chi2
    export PHASE_D_BACKSTOP_CAL=0
    export PHASE_D_BACKSTOP_UNC_SCALE=10.0
    ;;
  b01b03_chi2)
    export PHASE_D_PREDICT_VISIBLE_BANDS=B01,B02,B03,B04
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=chi2
    ;;
  b01b03_profile_s10)
    export PHASE_D_PREDICT_VISIBLE_BANDS=B01,B02,B03,B04
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=profile_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.10
    ;;
  b01b03_loo_s10)
    export PHASE_D_PREDICT_VISIBLE_BANDS=B01,B02,B03,B04
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=loo_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.10
    ;;
  b01b03_trimmed)
    export PHASE_D_PREDICT_VISIBLE_BANDS=B01,B02,B03,B04
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=trimmed_chi2
    ;;
  b03b08_chi2)
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08
    export PHASE_D_COST_MODE=chi2
    ;;
  b03b08_trimmed)
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08
    export PHASE_D_COST_MODE=trimmed_chi2
    ;;
  b03anchors_chi2)
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08,B11,B12
    export PHASE_D_COST_MODE=chi2
    ;;
  b03anchors_profile_s10)
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08,B11,B12
    export PHASE_D_COST_MODE=profile_scale
    export PHASE_D_PROFILE_SCALE_SIGMA=0.10
    ;;
  b01monthly_b03_chi2)
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=chi2
    ;;
  b03_anchor_weighted)
    export PHASE_D_COST_MODE=chi2
    export PHASE_D_PREDICT_VISIBLE_AGGREGATION=anchor_weighted
    export PHASE_D_PREDICT_VISIBLE_ANCHOR_MATCH_SCALE=0.05
    export PHASE_D_TAU_PRIOR=0
    unset PHASE_D_TAU_GATE
    ;;
  *)
    echo "Unknown SINGLE_PRIOR_OPTION=$OPTION" >&2
    exit 2
    ;;
esac

export PHASE_D_RESULTS_DIR="${SINGLE_PRIOR_RESULTS_DIR:-phaseD_results_lowcloud20_singleprior_${OPTION}_20260711}"
if [[ -n "${SINGLE_PRIOR_COST_CUBE_DIR:-}" ]]; then
  export PHASE_D_DUMP_COST_CUBE=1
  export PHASE_D_COST_CUBE_DIR="$SINGLE_PRIOR_COST_CUBE_DIR"
fi
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

if [[ -f "$RESULT_JSON" && "${PHASE_D_FORCE_RERUN:-0}" != "1" ]]; then
  existing_status=$(
    "$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
      "$RESULT_JSON" 2>/dev/null || echo malformed
  )
  if [[ "$existing_status" == "ok" || "$existing_status" == "no_valid_observation" \
    || "$existing_status" == "data_unavailable" ]]; then
    echo "Existing terminal result: option=$OPTION mid=$MID status=$existing_status"
    exit 0
  fi
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> option=$OPTION mid=$MID"
for attempt in 1 2 3; do
  status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    result_status=$(
      "$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
        "$RESULT_JSON" 2>/dev/null || echo malformed
    )
    if [[ "$result_status" == "no_valid_observation" \
      && "${PHASE_D_AUTO_MASK_FALLBACK:-1}" == "1" \
      && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_singleprior_no_valid"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" \
        "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
      export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
      export PHASE_D_IGNORE_CLOUD_WATER=1
      result_status=mask_fallback
      echo "No valid masked observation; retrying with recorded mask fallback" >&2
      continue
    fi
    if [[ "$result_status" == "ok" || "$result_status" == "no_valid_observation" \
      || "$result_status" == "data_unavailable" ]]; then
      exit 0
    fi
  else
    result_status=missing
  fi
  if (( attempt < 3 )); then
    sleep_for=$((60 + SLURM_ARRAY_TASK_ID % 60 + (attempt - 1) * 120))
    echo "Attempt $attempt ended with process=$status result=$result_status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "Retrieval failed after three attempts: option=$OPTION mid=$MID result=$result_status" >&2
exit 1

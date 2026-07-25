#!/bin/bash
#SBATCH --job-name=xrt_v2_ret
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=04:00:00
# Retrieval pools and numerical libraries are explicitly single-threaded below.
#SBATCH --cpus-per-task=1
# Prior physical retrievals peaked at 2.1 GB RSS; nine variants run sequentially.
#SBATCH --mem=4G
#SBATCH --array=1-152%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_v2_ret_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_v2_ret_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${CROSS_RT_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
MODEL_ROOT=${CROSS_RT_OUTPUT:-$ROOT/analysis/cross_rt_harmonizer_lowcloud152_20260716}
HISTORY_ROOT=${CROSS_RT_HISTORY_ROOT:-analysis/cross_rt_nonlinear_v2_lowcloud152_20260717/daily_histories}
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES" | cut -d, -f1)
if [[ -z "${MID:-}" || "$MID" == "matchup_id" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID in $CASES" >&2
  exit 2
fi

HISTORY_TAGS=(
  identity_daily
  hgb_target_band_cams_o3_blue_cap0p030
  hgb_target_band_cams_o3_solver_cap0p030
)

export GEE_SERVICE_ACCOUNT=python-gee@gee-marc.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/home/users/marcyin/gee-service-account.json
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
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
# Solver ozone follows the real CAMS column — matches the v2 target frame.
export PHASE_D_TCO3_MODE=cams

# Keep the solver and atmospheric state identical to the current physical run.
# Only the historical L2A harmonisation tag changes inside the loop below.
export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_SEASONAL_PREDICTOR=1
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
export PHASE_D_TAU_PRIOR=1
unset PHASE_D_TAU_GATE
export PHASE_D_SEASONAL_TAU_PREDICTOR=1
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=0
unset PHASE_D_BACKSTOP_CONFLICT_Z
unset PHASE_D_ATMO_PROVIDER
export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_adaptive
export PHASE_D_TARGET_TCWV_MODE=l2a_spatial
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_REFERENCE_TCWV=prior
export PHASE_D_DEM=glo30
export PHASE_D_DEM_PRE_M3=1
export PHASE_D_SCENE_MEAN_GEOMETRY=1
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_COST_MODE=chi2
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export PHASE_D_MATCHUP_ID=$MID
export PHASE_D_DUMP_COST_CUBE=0

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
RUNNER=/home/users/marcyin/SIAC/tools/aeronet_validation/run_harmonized_retrieval.py
failed=0
case_mask_fallback=0
for history_tag in "${HISTORY_TAGS[@]}"; do
  history_file="$ROOT/$HISTORY_ROOT/$history_tag/${MID}.npz"
  if [[ ! -f "$history_file" ]]; then
    echo "Missing harmonized history: $history_file" >&2
    failed=1
    continue
  fi
  export PHASE_D_SEASONAL_DIR="$HISTORY_ROOT/$history_tag"
  export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_crossrt_${history_tag}_physical_o3v2_20260717"
  # Observation masks do not depend on the surface-history variant. Once one
  # variant proves the strict mask has no support, start all later variants in
  # the same recovered mode instead of repeating a known failed first attempt.
  export PHASE_D_ALLOW_CLOUD_RETRIEVAL=$case_mask_fallback
  export PHASE_D_IGNORE_CLOUD_WATER=$case_mask_fallback
  export TMPDIR="$ROOT/spx_tmp_lowcloud20_crossrt_o3v2/${history_tag}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  mkdir -p "$TMPDIR" "$ROOT/$PHASE_D_RESULTS_DIR"
  result_json="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
  if [[ -f "$result_json" && "${CROSS_RT_FORCE_RERUN:-0}" != "1" ]]; then
    status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$result_json" 2>/dev/null || echo malformed)
    if [[ "$status" == "ok" || "$status" == "data_unavailable" ]]; then
      mask_bypassed=$("$PY" -c 'import json,sys; d=json.load(open(sys.argv[1])); print(int(bool((d.get("solver") or {}).get("surface_cloud_mask_bypassed"))))' "$result_json" 2>/dev/null || echo 0)
      if [[ "$mask_bypassed" == "1" ]]; then
        case_mask_fallback=1
      fi
      echo "Existing terminal result: tag=$history_tag mid=$MID status=$status"
      continue
    fi
  fi

  terminal=0
  # Four attempts guarantee at least one run after a third-attempt mask fallback.
  for attempt in 1 2 3 4; do
    process_status=0
    "$PY" "$RUNNER" || process_status=$?
    status=missing
    if [[ -f "$result_json" ]]; then
      status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' "$result_json" 2>/dev/null || echo malformed)
      if [[ "$status" == "ok" || "$status" == "data_unavailable" ]]; then
        terminal=1
        break
      fi
      if [[ "$status" == "no_valid_observation" && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
        diagnostic="$ROOT/diagnostics/lowcloud20_crossrt_nl_no_valid/$history_tag"
        mkdir -p "$diagnostic"
        mv "$result_json" "$diagnostic/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
        export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
        export PHASE_D_IGNORE_CLOUD_WATER=1
        case_mask_fallback=1
      fi
    fi
    if (( attempt < 4 )); then
      sleep $((45 + SLURM_ARRAY_TASK_ID % 30 + (attempt - 1) * 60))
    fi
    echo "Retry tag=$history_tag mid=$MID attempt=$attempt process=$process_status status=$status" >&2
  done
  if (( terminal == 0 )); then
    echo "Cross-RT retrieval failed: tag=$history_tag mid=$MID" >&2
    failed=1
  fi
done
exit "$failed"

#!/bin/bash
#SBATCH --job-name=medaodcombo
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=28G
#SBATCH --array=1-152%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodcombo_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodcombo_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS_FILE="${MEDIUM_PHYSICS_MIDS_FILE:-campaign250_lowcloud20_mids.txt}"
MIDS="$ROOT/$MIDS_FILE"
OPTION="${MEDIUM_PHYSICS_OPTION:-dem_no_debias}"
CASE_INDEX=$SLURM_ARRAY_TASK_ID
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
export MALLOC_CONF=background_thread:false,narenas:1
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
export PHASE_D_DEM_PRE_M3=0
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

case "$OPTION" in
  baseline_costcube)
    ;;
  anchor_converged)
    export PHASE_D_ANCHOR_MAX_PASSES=5
    export PHASE_D_ANCHOR_CONVERGENCE_TOL=0.005
    ;;
  additive_offset)
    export PHASE_D_COST_MODE=additive_offset
    ;;
  b01_additive_offset)
    export PHASE_D_PREDICT_VISIBLE_BANDS=B01,B02,B03,B04
    export PHASE_D_SOLVE_BANDS=B01,B02,B03,B04
    export PHASE_D_COST_MODE=additive_offset
    ;;
  dem_no_debias)
    export PHASE_D_DEM=glo30
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    ;;
  dem_l2awvp_no_debias)
    export PHASE_D_DEM=glo30
    export PHASE_D_TARGET_TCWV_MODE=l2a
    export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
    export PHASE_D_REFERENCE_TCWV=prior
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    ;;
  no_debias_product_unc)
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  dem_no_debias_product_unc)
    export PHASE_D_DEM=glo30
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  dem_product_unc)
    export PHASE_D_DEM=glo30
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  dem_pre_m3)
    export PHASE_D_DEM=glo30
    export PHASE_D_DEM_PRE_M3=1
    ;;
  dem_pre_m3_product_unc)
    export PHASE_D_DEM=glo30
    export PHASE_D_DEM_PRE_M3=1
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  tau_always)
    unset PHASE_D_TAU_GATE
    ;;
  tau_product_unc)
    unset PHASE_D_TAU_GATE
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  dem_pre_m3_tau_always)
    export PHASE_D_DEM=glo30
    export PHASE_D_DEM_PRE_M3=1
    unset PHASE_D_TAU_GATE
    ;;
  dem_pre_m3_tau_product_unc)
    export PHASE_D_DEM=glo30
    export PHASE_D_DEM_PRE_M3=1
    unset PHASE_D_TAU_GATE
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  dem_l2awvp_no_debias_product_unc)
    export PHASE_D_DEM=glo30
    export PHASE_D_TARGET_TCWV_MODE=l2a
    export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
    export PHASE_D_REFERENCE_TCWV=prior
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  swir_anchor)
    export PHASE_D_ANCHOR_BANDS=B11,B12
    ;;
  swir_anchor_no_debias)
    export PHASE_D_ANCHOR_BANDS=B11,B12
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    ;;
  maiac_orbit_matched)
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_orbitmatched_mediumdev_20260713
    ;;
  maiac_orbit_matched_product_unc)
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_orbitmatched_mediumdev_20260713
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  maiac_footprint0045)
    # Match the native-grid MAIAC aggregation footprint to the S2 retrieval AOI.
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_footprint0045_mediumdev_20260713
    ;;
  maiac_footprint0045_product_unc)
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_footprint0045_mediumdev_20260713
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  tau_no_debias_no_iter)
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_ANCHOR_ITERATE=0
    unset PHASE_D_TAU_GATE
    ;;
  tau_no_debias_no_iter_product_unc)
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_ANCHOR_ITERATE=0
    unset PHASE_D_TAU_GATE
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  maiac_orbit_tau_no_debias_no_iter)
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_orbitmatched_mediumdev_20260713
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_ANCHOR_ITERATE=0
    unset PHASE_D_TAU_GATE
    ;;
  maiac_orbit_tau_no_debias_no_iter_product_unc)
    export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_orbitmatched_mediumdev_20260713
    export PHASE_D_PREDICT_VISIBLE_DEBIAS=off
    export PHASE_D_ANCHOR_ITERATE=0
    unset PHASE_D_TAU_GATE
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  maiac_native_spatial)
    unset PHASE_D_MAIAC_DIR
    export PHASE_D_ATMO_PROVIDER=mcd19_native
    ;;
  maiac_native_spatial_product_unc)
    unset PHASE_D_MAIAC_DIR
    export PHASE_D_ATMO_PROVIDER=mcd19_native
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  lut_linear_geometry)
    export SIAC_LUT_SCENE_GEOMETRY_MODE=linear
    ;;
  lut_linear_geometry_product_unc)
    export SIAC_LUT_SCENE_GEOMETRY_MODE=linear
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  tau_swir_anchor)
    export PHASE_D_ANCHOR_BANDS=B11,B12
    unset PHASE_D_TAU_GATE
    ;;
  tau_swir_anchor_no_iter)
    export PHASE_D_ANCHOR_BANDS=B11,B12
    unset PHASE_D_TAU_GATE
    export PHASE_D_ANCHOR_ITERATE=0
    ;;
  tau_swir_anchor_product_unc)
    export PHASE_D_ANCHOR_BANDS=B11,B12
    unset PHASE_D_TAU_GATE
    export PHASE_D_BACKSTOP_CAL=0
    ;;
  water_buffer0)
    export PHASE_D_WATER_BUFFER=0
    ;;
  water_buffer8)
    export PHASE_D_WATER_BUFFER=8
    ;;
  water_buffer16)
    export PHASE_D_WATER_BUFFER=16
    ;;
  allow_cloud_retrieval)
    # Diagnostic for aerosol that OmniCloudMask may label as cloud. Water and
    # no-data remain masked; only the cloud exclusion is relaxed.
    export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
    ;;
  nir_cost)
    # Keep one historical S2 prior, but score the directly corresponding B08
    # NIR band as well as using B8A to predict the visible surface.
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08
    ;;
  nir_swir_cost)
    # Directly score all three surface anchors in the same prior/RT solve.
    export PHASE_D_SOLVE_BANDS=B02,B03,B04,B08,B11,B12
    ;;
  ridge_predictor)
    # Lower-sensitivity polynomial ridge surface mapping, retaining the same
    # historical S2 inputs, visible bands, RT, and atmospheric constraint.
    export PHASE_D_PREDICTOR_MODEL=ridge
    ;;
  clean_day_aod015)
    # Build the same S2 monthly surface prior only from genuinely clean
    # historical days, rather than a relative within-month AOD quantile.
    export PHASE_D_BESTPIXEL_AOD_MAX=0.15
    ;;
  *)
    echo "Unknown MEDIUM_PHYSICS_OPTION=$OPTION" >&2
    exit 2
    ;;
esac

export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_mediumphysics_${OPTION}_mediumdev_20260713"
export TMPDIR="$ROOT/spx_tmp_lowcloud20_mediumphysics/$OPTION/${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}"
mkdir -p "$TMPDIR"
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

if [[ -f "$RESULT_JSON" && "${PHASE_D_FORCE_RERUN:-0}" != "1" ]]; then
  existing_status=$(
    "$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
      "$RESULT_JSON" 2>/dev/null || echo malformed
  )
  if [[ "$existing_status" == "ok" || "$existing_status" == "data_unavailable" ]]; then
    echo "Existing terminal result: option=$OPTION mid=$MID status=$existing_status"
    exit 0
  fi
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> option=$OPTION mid=$MID"
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
      fallback_dir="$ROOT/diagnostics/lowcloud20_mediumphysics_no_valid/$OPTION"
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

echo "Retrieval failed after three attempts: option=$OPTION mid=$MID result=$result_status" >&2
exit 1

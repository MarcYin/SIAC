#!/bin/bash
#SBATCH --job-name=lc20physanchor
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
# host1065 exhausted its process/thread limit; host1262 and host1263 have
# cancelled otherwise idle jobs with no process output or resource use.
#SBATCH --exclude=host1042,host1057,host1065,host1153,host1201,host1205,host1206,host1243,host1261,host1262,host1263
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-152%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20physanchor_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20physanchor_%A_%a.out

set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
TAG="${PHYSICAL_ANCHOR_TAG:-20260715}"
MIDS_FILE="${PHYSICAL_ANCHOR_MIDS_FILE:-campaign250_lowcloud20_mids.txt}"
if [[ "$MIDS_FILE" != /* ]]; then
  MIDS_FILE="$ROOT/$MIDS_FILE"
fi
CASE_INDEX=$SLURM_ARRAY_TASK_ID
MID=$(sed -n "${CASE_INDEX}p" "$MIDS_FILE")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for task ${CASE_INDEX} in $MIDS_FILE" >&2
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
export TOKIO_WORKER_THREADS=1
export POLARS_MAX_THREADS=1
export GDAL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
export MALLOC_CONF=background_thread:false,narenas:1
export MAIAC_AGG=median

# One physical atmospheric state is shared by the anchor and the solve:
# MAIAC AOD; same-scene Sen2Cor/ACIX-compatible L2A WVP; GLO-30 elevation;
# one AOI-mean Sentinel-2 acquisition geometry. AOD is the only swept RT
# dimension; water vapour and elevation retain their spatial fields.
export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_PREDICT_VISIBLE_BANDS=B02,B03,B04
export PHASE_D_SOLVE_BANDS=B02,B03,B04
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICT_VISIBLE_ANCHOR_SOURCE=atmo_prior
# The runner's legacy seasonal switch is set explicitly as well.  The anchor
# now receives the complete atmospheric-prior AOD field, rather than a CAMS
# scalar fallback, so this is intentionally not interchangeable with CAMS.
export PHASE_D_SEASONAL_ANCHOR_AOD_SOURCE=atmo_prior
export PHASE_D_PREDICT_VISIBLE_DEBIAS="${PHYSICAL_ANCHOR_DEBIAS:-off}"
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_BANDS=B8A,B11,B12
export PHASE_D_ANCHOR_ITERATE="${PHYSICAL_ANCHOR_ANCHOR_ITERATE:-0}"
export PHASE_D_ANCHOR_ITERATE_MODE="${PHYSICAL_ANCHOR_ANCHOR_ITERATE_MODE:-anchor}"
export PHASE_D_ANCHOR_MAX_PASSES="${PHYSICAL_ANCHOR_ANCHOR_MAX_PASSES:-2}"
export PHASE_D_ANCHOR_CONVERGENCE_TOL="${PHYSICAL_ANCHOR_ANCHOR_CONVERGENCE_TOL:-0.005}"
case "${PHYSICAL_ANCHOR_TAU_MODE:-always}" in
  always)
    export PHASE_D_TAU_PRIOR=1
    unset PHASE_D_TAU_GATE
    ;;
  gated)
    export PHASE_D_TAU_PRIOR=1
    export PHASE_D_TAU_GATE="${PHYSICAL_ANCHOR_TAU_GATE:-8.0}"
    ;;
  static)
    export PHASE_D_TAU_PRIOR=0
    unset PHASE_D_TAU_GATE
    ;;
  *)
    echo "Unknown PHYSICAL_ANCHOR_TAU_MODE=${PHYSICAL_ANCHOR_TAU_MODE}" >&2
    exit 2
    ;;
esac
export PHASE_D_BACKSTOP_CAL="${PHYSICAL_ANCHOR_BACKSTOP_CAL:-1}"
export PHASE_D_BACKSTOP_ESCAPE=0
if [[ -n "${PHYSICAL_ANCHOR_CONFLICT_Z:-}" ]]; then
  export PHASE_D_BACKSTOP_CONFLICT_Z="$PHYSICAL_ANCHOR_CONFLICT_Z"
else
  unset PHASE_D_BACKSTOP_CONFLICT_Z
fi
if [[ -n "${PHYSICAL_ANCHOR_ATMO_PROVIDER:-}" ]]; then
  export PHASE_D_ATMO_PROVIDER="$PHYSICAL_ANCHOR_ATMO_PROVIDER"
else
  unset PHASE_D_ATMO_PROVIDER
fi
export PHASE_D_MAIAC_DIR=maiac_qa_lowcloud20_native_adaptive
# ``l2a_spatial`` reads the matching Sen2Cor WVP raster onto the atmospheric
# grid; legacy ``l2a`` retains the collected station-buffer scalar.
export PHASE_D_TARGET_TCWV_MODE="${PHYSICAL_ANCHOR_TCWV_MODE:-l2a_spatial}"
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_REFERENCE_TCWV=prior
export PHASE_D_DEM=glo30
export PHASE_D_DEM_PRE_M3=1
export PHASE_D_SCENE_MEAN_GEOMETRY="${PHYSICAL_ANCHOR_SCENE_MEAN_GEOMETRY:-1}"
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

export PHASE_D_RESULTS_DIR="${PHYSICAL_ANCHOR_RESULTS_DIR:-phaseD_results_lowcloud20_physical_anchor_maiac_l2awvp_dem_scenegeom_${TAG}}"
export TMPDIR="$ROOT/spx_tmp_lowcloud20_physical_anchor/${TAG}/${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}"
mkdir -p "$TMPDIR" "$ROOT/$PHASE_D_RESULTS_DIR"
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

if [[ -f "$RESULT_JSON" && "${PHASE_D_FORCE_RERUN:-0}" != "1" ]]; then
  existing_status=$("$PY" -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["status"]).lower())' \
    "$RESULT_JSON" 2>/dev/null || echo malformed)
  if [[ "$existing_status" == "ok" || "$existing_status" == "no_valid_observation" \
    || "$existing_status" == "data_unavailable" ]]; then
    echo "Existing terminal result: mid=$MID status=$existing_status"
    exit 0
  fi
fi

echo "TASK $CASE_INDEX -> physical-anchor mid=$MID"
result_status=missing
for attempt in 1 2 3; do
  process_status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || process_status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    mapfile -t result_fields < <(
      "$PY" -c 'import json,sys; row=json.load(open(sys.argv[1])); print(str(row["status"]).lower()); print(str(row.get("reason", "")))' \
        "$RESULT_JSON" 2>/dev/null || printf 'malformed\n\n'
    )
    result_status="${result_fields[0]:-malformed}"
    result_reason="${result_fields[1]:-}"
    if [[ "$result_status" == "no_valid_observation" \
      && "${result_reason,,}" == *"bestpixel surface prior produced no usable composites"* ]]; then
      refresh_dir="$ROOT/diagnostics/lowcloud20_physical_anchor_surface_refresh"
      mkdir -p "$refresh_dir"
      mv "$RESULT_JSON" "$refresh_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}.attempt${attempt}.json"
      result_status=surface_refresh
      echo "Surface composites unavailable; retrying in a fresh process" >&2
      continue
    fi
    if [[ "$result_status" == "no_valid_observation" \
      && "${PHASE_D_AUTO_MASK_FALLBACK:-1}" == "1" \
      && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_physical_anchor_no_valid"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${CASE_INDEX}.json"
      export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
      export PHASE_D_IGNORE_CLOUD_WATER=1
      result_status=mask_fallback
      echo "No valid masked observation; retrying with the recorded cloud-mask fallback" >&2
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
    sleep_for=$((60 + CASE_INDEX % 60 + (attempt - 1) * 120))
    echo "Attempt $attempt ended with process=$process_status result=$result_status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "Retrieval failed after three attempts: mid=$MID result=$result_status" >&2
exit 1

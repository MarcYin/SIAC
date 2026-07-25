#!/bin/bash
#SBATCH --job-name=lc20species
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1092,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-152%4
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20species_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20species_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
OPTION="${AEROSOL_SPECIES_OPTION:-sixs_cci3}"
TAG="${AEROSOL_SPECIES_TAG:-20260712}"
MIDS_FILE="${AEROSOL_SPECIES_MIDS_FILE:-$ROOT/campaign250_lowcloud20_mids.txt}"
if [[ "$MIDS_FILE" != /* ]]; then
  MIDS_FILE="$ROOT/$MIDS_FILE"
fi
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS_FILE")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for task ${SLURM_ARRAY_TASK_ID} in $MIDS_FILE" >&2
  exit 2
fi

case "$OPTION" in
  lut_continental|sixs_continental|sixs_cci3|sixs_cci_exact|sixs_canonical|libradtran_continental)
    ;;
  *)
    echo "Unknown AEROSOL_SPECIES_OPTION=$OPTION" >&2
    exit 2
    ;;
esac

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
export TMPDIR="$ROOT/spx_tmp_lowcloud20_aerosol_species/${OPTION}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_lowcloud20_aerosol_species"

# Keep the measured three-band recipe fixed across all RT/species arms.
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
export PHASE_D_BACKSTOP_CAL="${AEROSOL_SPECIES_BACKSTOP_CAL:-1}"
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_COST_MODE=chi2
export PHASE_D_MAIAC_DIR="${AEROSOL_SPECIES_MAIAC_DIR:-maiac_qa}"
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
if [[ "${PHASE_D_FORCE_MASK_FALLBACK:-0}" == "1" ]]; then
  export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
  export PHASE_D_IGNORE_CLOUD_WATER=1
fi
export PHASE_D_TARGET_TCWV_MODE=prior
export PHASE_D_REFERENCE_TCWV="${AEROSOL_SPECIES_REFERENCE_TCWV:-2.0}"
export PHASE_D_DEM="${AEROSOL_SPECIES_DEM:-none}"
export PHASE_D_PREDICT_VISIBLE_DEBIAS="${AEROSOL_SPECIES_DEBIAS:-on}"
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_MIDS_FILE="$MIDS_FILE"
export PHASE_D_INDEX=$SLURM_ARRAY_TASK_ID

case "$OPTION" in
  lut_continental)
    export PHASE_D_RT_BACKEND=lut
    export PHASE_D_AEROSOL_PROFILE=continental_average
    ;;
  sixs_continental)
    export PHASE_D_RT_BACKEND=sixs
    export PHASE_D_AEROSOL_PROFILE=continental
    ;;
  sixs_cci3)
    export PHASE_D_RT_BACKEND=sixs
    export PHASE_D_AEROSOL_PROFILE=continental
    export PHASE_D_SURFDRV_SPECIES=cci_climatology
    export PHASE_D_SURFDRV_SPECIES_N=3
    ;;
  sixs_cci_exact)
    export PHASE_D_RT_BACKEND=sixs
    export PHASE_D_AEROSOL_PROFILE=continental
    export PHASE_D_SURFDRV_SPECIES=cci_climatology_exact
    ;;
  sixs_canonical)
    export PHASE_D_RT_BACKEND=sixs
    export PHASE_D_AEROSOL_PROFILE=continental
    export PHASE_D_SURFDRV_SPECIES=canonical_6s
    ;;
  libradtran_continental)
    export PHASE_D_RT_BACKEND=libradtran
    export PHASE_D_AEROSOL_PROFILE=continental_average
    # Direct libRadtran internally bounds its worker count by memory. The
    # harness-wide ThreadPoolExecutor patch would otherwise serialize uvspec.
    export PHASE_D_SERIAL_POOLS=0
    ;;
esac

if [[ "$PHASE_D_RT_BACKEND" == "sixs" ]]; then
  export PHASE_D_SIXS_MODULE_PATH="$ROOT/../siac_env/siac_cache/rt6s/release/_siac_rt6s_native.cpython-311-x86_64-linux-gnu.so"
  export PHASE_D_SIXS_THREADS="${AEROSOL_SPECIES_SIXS_THREADS:-1}"
  export PHASE_D_SIXS_WORKER_LIBRARIES="${AEROSOL_SPECIES_SIXS_WORKER_LIBRARIES:-1}"
  export PHASE_D_SIXS_RUN_CACHE_DIR="$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
fi

export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_aerosol_species_${OPTION}_${TAG}"
COST_DIR="$ROOT/phaseD_cost_cubes_lowcloud20_aerosol_species_${OPTION}_${TAG}"
RESULT_JSON="$ROOT/$PHASE_D_RESULTS_DIR/${MID}.json"
mkdir -p "$COST_DIR" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
if [[ "${AEROSOL_SPECIES_DUMP_COST_CUBE:-1}" == "1" ]]; then
  export SIAC_DUMP_COST_CUBE="$COST_DIR/${MID}.npz"
fi

echo "TASK $SLURM_ARRAY_TASK_ID -> option=$OPTION mid=$MID"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
result_status=missing
result_reason=
for attempt in 1 2 3; do
  status=0
  "$PY" "$ROOT/phaseD_array_runner.py" || status=$?
  if [[ -f "$RESULT_JSON" ]]; then
    mapfile -t result_fields < <(
      "$PY" -c \
        'import json,sys; row=json.load(open(sys.argv[1])); print(str(row["status"]).lower()); print(str(row.get("reason", "")))' \
        "$RESULT_JSON" 2>/dev/null || printf 'malformed\n\n'
    )
    result_status="${result_fields[0]:-malformed}"
    result_reason="${result_fields[1]:-}"
    if [[ "$result_status" == "no_valid_observation" \
      && "${result_reason,,}" == *"bestpixel surface prior produced no usable composites"* ]]; then
      refresh_dir="$ROOT/diagnostics/lowcloud20_aerosol_species_surface_refresh/$OPTION"
      mkdir -p "$refresh_dir"
      mv "$RESULT_JSON" \
        "$refresh_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.attempt${attempt}.json"
      result_status=surface_refresh
      echo "Surface composites unavailable; retrying in a fresh process with unchanged masks" >&2
      continue
    fi
    if [[ "$result_status" == "no_valid_observation" \
      && "${PHASE_D_AUTO_MASK_FALLBACK:-1}" == "1" \
      && "$PHASE_D_IGNORE_CLOUD_WATER" != "1" ]]; then
      fallback_dir="$ROOT/diagnostics/lowcloud20_aerosol_species_no_valid/$OPTION"
      mkdir -p "$fallback_dir"
      mv "$RESULT_JSON" \
        "$fallback_dir/${MID}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
      export PHASE_D_ALLOW_CLOUD_RETRIEVAL=1
      export PHASE_D_IGNORE_CLOUD_WATER=1
      result_status=mask_fallback
      echo "No valid masked observation; retrying with recorded mask fallback" >&2
      continue
    fi
    if (( status == 0 )) && [[ "$result_status" == "ok" ]]; then
      exit 0
    fi
  fi
  if (( attempt < 3 )); then
    sleep_for=$((20 + attempt * 20))
    echo "Attempt $attempt ended with process=$status result=$result_status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "Retrieval failed after three attempts: option=$OPTION mid=$MID result=$result_status" >&2
exit 1

#!/bin/bash
#SBATCH --job-name=xrt_s6path
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-2
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_s6path_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_s6path_%A_%a.out
set -euo pipefail
umask 0022

# Species-pathology diagnostics: re-run U_of_Wisconsin_SSEC (clean blowup
# 0.057->0.274 under cci_climatology_exact) and ICIMOD (overshoot 1.22->2.07)
# in both species arms with cost-cube + calibration dumps enabled, so the
# per-band cost curves and predicted prior can be compared offline.
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS=$ROOT/species_pathology_mids.txt
HISTORY_ROOT=analysis/cross_rt_nonlinear_v3_deterrain_lowcloud152_20260718/daily_histories
HISTORY_TAG=hgb_target_band_cams_o3_blue_cap0p030
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")
if [[ -z "${MID:-}" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

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
unset PHASE_D_BACKSTOP_UNC_SCALE
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
export PHASE_D_MATCHUP_ID=$MID
export PHASE_D_SEASONAL_DIR="$HISTORY_ROOT/$HISTORY_TAG"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
RUNNER=/home/users/marcyin/SIAC/tools/aeronet_validation/run_harmonized_retrieval.py

# 6S RT backend, shared run cache with the full A/B (mostly warm).
export PHASE_D_RT_BACKEND=sixs
export PHASE_D_AEROSOL_PROFILE=continental
export PHASE_D_SIXS_MODULE_PATH=/gws/ssde/j25a/nceo_isp/public/siac_env/siac_cache/rt6s/release/_siac_rt6s_native.cpython-311-x86_64-linux-gnu.so
export PHASE_D_SIXS_THREADS=4
export PHASE_D_SIXS_WORKER_LIBRARIES=1
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_species_probe"
export PHASE_D_SIXS_RUN_CACHE_DIR="$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
mkdir -p "$PHASE_D_SIXS_RUN_CACHE_DIR"

# Diagnostics ON for this job.
export PHASE_D_DUMP_COST_CUBE=1

failed=0
for arm in sixs_cont sixs_cci; do
  unset PHASE_D_SURFDRV_SPECIES
  if [[ "$arm" == "sixs_cci" ]]; then
    export PHASE_D_SURFDRV_SPECIES=cci_climatology_exact
  fi
  export PHASE_D_CALIB_DUMP_DIR="$ROOT/calib_dumps_pathology_${arm}_20260719"
  export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_crossrt_${arm}_pathology_20260719"
  export TMPDIR="$ROOT/spx_tmp_species_pathology/${arm}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  mkdir -p "$TMPDIR" "$ROOT/$PHASE_D_RESULTS_DIR" "$PHASE_D_CALIB_DUMP_DIR"
  cube_file="$ROOT/phaseD_cost_cubes/$PHASE_D_RESULTS_DIR/${MID}.npz"
  if [[ -f "$cube_file" ]]; then
    echo "Cube exists: arm=$arm mid=$MID"
    continue
  fi
  ok=0
  for attempt in 1 2 3; do
    if "$PY" "$RUNNER"; then
      ok=1
      break
    fi
    sleep $((60 * attempt))
  done
  if (( ok == 0 )); then
    echo "Pathology dump failed: arm=$arm mid=$MID" >&2
    failed=1
  fi
done
exit "$failed"

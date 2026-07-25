#!/bin/bash
#SBATCH --job-name=xrt_calib152
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --array=1-152%10
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_calib152_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_calib152_%A_%a.out
set -euo pipefail
umask 0022

# Surface-prior evaluation dumps: re-run the 71.1% crossrt physical_dt
# baseline for 20 representative cases with PHASE_D_CALIB_DUMP_DIR so the
# per-pixel predicted prior, scene TOA, and RT coefficient curves are saved.
# BOA at AERONET-truth AOD is then reconstructable offline, giving a
# truth-referenced visible-band prior accuracy assessment.
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv
HISTORY_ROOT=analysis/cross_rt_nonlinear_v3_deterrain_lowcloud152_20260718/daily_histories
HISTORY_TAG=hgb_target_band_cams_o3_blue_cap0p030
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES" | cut -d, -f1)
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
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export PHASE_D_MATCHUP_ID=$MID
export PHASE_D_DUMP_COST_CUBE=0
export PHASE_D_SEASONAL_DIR="$HISTORY_ROOT/$HISTORY_TAG"
export PHASE_D_CALIB_DUMP_DIR="$ROOT/calib_dumps_crossrt_dt_20260719"
export PHASE_D_RESULTS_DIR="phaseD_results_lowcloud20_crossrt_calibdump_20260719"
export TMPDIR="$ROOT/spx_tmp_lowcloud20_calibdump/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR" "$ROOT/$PHASE_D_RESULTS_DIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
RUNNER=/home/users/marcyin/SIAC/tools/aeronet_validation/run_harmonized_retrieval.py
dump_file="$ROOT/calib_dumps_crossrt_dt_20260719/${MID}.npz"
if [[ -f "$dump_file" ]]; then
  echo "Calibration dump already present: $dump_file"
  exit 0
fi
for attempt in 1 2 3; do
  if "$PY" "$RUNNER" && [[ -f "$dump_file" ]]; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((45 * attempt + SLURM_ARRAY_TASK_ID % 30))
  fi
done
echo "Calibration dump failed after three attempts: $MID" >&2
exit 1

#!/bin/bash
#SBATCH --job-name=c250maskr2
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-250%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250maskr2_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250maskr2_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MAIAC_AGG=median
export TMPDIR="$ROOT/spx_tmp_campaign250_masked_r2/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_campaign250_masked_options"

# Masked 6S + same-scene L2A-WVP campaign derived from the R2 solver recipe.
# This is not a mask-only reproduction of historical R2, which used the LUT
# backend and atmospheric-prior TCWV; use the paired factorial for attribution.
export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_ITERATE=1
export PHASE_D_TAU_PRIOR=1
export PHASE_D_TAU_GATE=8.0
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=0
export PHASE_D_COST_MODE=auto2
export PHASE_D_SOLVE_BANDS=B02,B04
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MAIAC_DIR=maiac_qa
export PHASE_D_ALLOW_CLOUD_RETRIEVAL=0
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_TARGET_TCWV_MODE=l2a
export PHASE_D_TARGET_TCWV_FIXED=2.0
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_RT_BACKEND=sixs
export PHASE_D_AEROSOL_PROFILE=continental
export PHASE_D_SIXS_MODULE_PATH="$ROOT/../siac_env/siac_cache/rt6s/release/_siac_rt6s_native.cpython-311-x86_64-linux-gnu.so"
export PHASE_D_SKIP_CORRECTION=1
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_MIDS_FILE=campaign250_mids.txt
export PHASE_D_RESULTS_DIR=phaseD_results_campaign250_masked_r2_l2awvp_6s_20260710

COST_DIR="$ROOT/phaseD_cost_cubes_campaign250_masked_r2_l2awvp_6s_20260710"
mkdir -p "$COST_DIR" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
export SIAC_DUMP_COST_CUBE="$COST_DIR/${MID}.npz"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/phaseD_array_runner.py"

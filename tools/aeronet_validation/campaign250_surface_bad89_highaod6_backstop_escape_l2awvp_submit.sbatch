#!/bin/bash
#SBATCH --job-name=c250highaod6
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-6%2
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6_l2awvp_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6_l2awvp_%A_%a.out
set -uo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_surface_bad89_highaod6_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8 MAIAC_AGG=median
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_surface_bad89_highaod6_l2awvp"
export TMPDIR="$ROOT/spx_tmp_surface_bad89_highaod6_l2awvp"

export PHASE_D_CONFIG_PRESET=l2a_monthly_predictor
export PHASE_D_BESTPIXEL_ENDPOINT=pc
export PHASE_D_BESTPIXEL_SOURCE=l2a
export PHASE_D_USE_PIXI_BESTPIXEL=1
export PHASE_D_SOLVE_BANDS=B02,B04
export PHASE_D_MAIAC_DIR=maiac_qa
export PHASE_D_ALLOW_CLOUD_RETRIEVAL="${PHASE_D_ALLOW_CLOUD_RETRIEVAL:-0}"
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_PREDICT_VISIBLE_FLOOR=0.006
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_ANCHOR_ITERATE=1
export PHASE_D_TAU_PRIOR=1
export PHASE_D_TAU_GATE=8.0
export PHASE_D_BACKSTOP_CAL=1
export PHASE_D_BACKSTOP_ESCAPE=1
export PHASE_D_BACKSTOP_ESCAPE_LOW_AOT=0.08
export PHASE_D_BACKSTOP_ESCAPE_COST_THRESHOLD=20.0
export PHASE_D_BACKSTOP_ESCAPE_DELTA=0.05
export PHASE_D_BACKSTOP_ESCAPE_BAND_SPREAD=0.25
export PHASE_D_BACKSTOP_ESCAPE_MIN_JUMP=0.35
export PHASE_D_BACKSTOP_ESCAPE_COST_RATIO=1.8
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MIDS_FILE=campaign250_surface_bad89_highaod6_mids.txt
export PHASE_D_RESULTS_DIR=phaseD_results_campaign250_surface_bad89_highaod6_backstop_escape_l2awvp_20260710
export PHASE_D_TARGET_TCWV_MODE=l2a
export PHASE_D_TARGET_TCWV_FIXED=2.0
export PHASE_D_TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PHASE_D_RT_BACKEND=sixs
# 6S-native aerosol family; usstd/us60 is for atmosphere profile, not RT aerosol.
export PHASE_D_AEROSOL_PROFILE=continental

COST_DIR="$ROOT/phaseD_cost_cubes_campaign250_surface_bad89_highaod6_backstop_escape_20260710_l2awvp"
mkdir -p "$COST_DIR" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT" "$SIAC_RT_RUN_CACHE_ROOT/rt6s" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache" 2>/dev/null || true
chmod u+rwx,g+rx,o+rx "$TMPDIR" 2>/dev/null || true
export SIAC_DUMP_COST_CUBE="$COST_DIR/${MID}.npz"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/phaseD_array_runner.py"

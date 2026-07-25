#!/bin/bash
#SBATCH --job-name=c250fulllocal
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-250%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250fulllocal_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250fulllocal_%A_%a.out
set -uo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8 MAIAC_AGG=median

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
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MIDS_FILE=campaign250_mids.txt
export PHASE_D_RESULTS_DIR=phaseD_results_campaign250_R2_full_localdiag_20260705

COST_DIR="$ROOT/phaseD_cost_cubes_campaign250_R2_full_localdiag_20260705"
mkdir -p "$COST_DIR"
export SIAC_DUMP_COST_CUBE="$COST_DIR/${MID}.npz"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/phaseD_array_runner.py"

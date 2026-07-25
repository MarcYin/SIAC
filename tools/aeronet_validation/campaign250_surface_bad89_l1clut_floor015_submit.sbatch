#!/bin/bash
#SBATCH --job-name=b89l1cret
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-89%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/b89l1cret_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/b89l1cret_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8 MAIAC_AGG=median

# L1C MAIAC clean-day prior in SIAC LUT RT space, using the best full-62 H2
# recipe: direct L1C-LUT composite plus ExtraTrees visible predictor, no debias,
# uncertainty floor 0.015, solved with B02/B04 and staged MAIAC backstop.
export PHASE_D_CONFIG_PRESET=example
export PHASE_D_PREBUILT_DIR=l1c_seasonal_lut
export PHASE_D_SEASONAL_PREDICTOR=1
export PHASE_D_SEASONAL_DIR=l1c_seasonal_lut
export PHASE_D_SEASONAL_UNC_FLOOR=0.015
export PHASE_D_SEASONAL_DEBIAS=0
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_PREDICTOR_MODEL=et20
export PHASE_D_SOLVE_BANDS=B02,B04
export PHASE_D_SEASONAL_BANDS=B02,B04
export PHASE_D_SEASONAL_B01_UNC=0.015
export PHASE_D_MAIAC_DIR=maiac_qa
export PHASE_D_ALLOW_CLOUD_RETRIEVAL="${PHASE_D_ALLOW_CLOUD_RETRIEVAL:-0}"
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_MAX_WORKERS=1
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MIDS_FILE=campaign250_surface_bad89_mids.txt
export PHASE_D_RESULTS_DIR=phaseD_results_campaign250_surface_bad89_l1clut_floor015_20260705

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_surface_bad89_mids.txt"
PREBUILT="$ROOT/l1c_seasonal_lut"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

if [[ ! -f "$PREBUILT/${MID}.npz" ]]; then
  echo "Missing required L1C-LUT composite: $PREBUILT/${MID}.npz" >&2
  exit 3
fi

cd "$ROOT"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/phaseD_array_runner.py"

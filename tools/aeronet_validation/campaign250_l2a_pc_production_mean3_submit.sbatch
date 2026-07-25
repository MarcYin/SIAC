#!/bin/bash
#SBATCH --job-name=c250pcmean3
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-250%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250pcmean3_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250pcmean3_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8 MAIAC_AGG=median

# Generic fixed L2A seasonal-PC ensemble:
#   1. direct L2A seasonal PC composite solve,
#   2. ExtraTree-predicted visible-prior solve,
#   3. the same predictor with BOA uncertainty floor 0.03,
# fused by fixed arithmetic mean, then AOI/scene-mean pooled.
#
# This is the original simple recipe that reached 87% on full-62. The campaign
# runner below first verifies all 250 PC composites exist to avoid partial scores.
export PHASE_D_PREBUILT_DIR=l2a_seasonal_pc
export PHASE_D_SEASONAL_PREDICTOR=1
export PHASE_D_AOD_SELECTOR=production_l2a_mean3
export PHASE_D_AOD_SELECTOR_ALT_PREBUILT_DIR=none
export PHASE_D_AOD_SELECTOR_FALLBACK_SOLVE_BANDS=B02,B04
export PHASE_D_ENSEMBLE_PRED_PREBUILT_DIR=l2a_seasonal_pc
export PHASE_D_MAIAC_DIR=maiac_qa
export PHASE_D_ALLOW_CLOUD_RETRIEVAL="${PHASE_D_ALLOW_CLOUD_RETRIEVAL:-0}"
export PHASE_D_IGNORE_CLOUD_WATER=0
export PHASE_D_MAX_WORKERS=1
export PHASE_D_SERIAL_POOLS=1
export PHASE_D_RETRIEVAL_SPATIAL_POOL=scene_mean
export PHASE_D_MIDS_FILE=campaign250_mids.txt
export PHASE_D_RESULTS_DIR=phaseD_results_campaign250_l2a_pc_production_mean3_scene_mean

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
PREBUILT="$ROOT/l2a_seasonal_pc"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

if [[ ! -f "$PREBUILT/${MID}.npz" ]]; then
  echo "Missing required prebuilt L2A PC composite: $PREBUILT/${MID}.npz" >&2
  exit 3
fi

cd "$ROOT"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/phaseD_array_runner.py"

#!/bin/bash
#SBATCH --job-name=exttrain1000
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-1000%24
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/exttrain1000_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/exttrain1000_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/external_train1000_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

S2_L2A_AOT_OUTDIR=s2_l2a_aot_train1000 \
  "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/collect_s2_l2a_aot.py "$MID" || true

GEE_CAMS_AOD_OUTDIR=gee_cams_aod_train1000 \
GEE_MERRA_AOD_OUTDIR=gee_merra_aod_train1000 \
GEE_MODEL_MEAN_AOD_OUTDIR=gee_model_mean_aod_train1000 \
  "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/collect_gee_model_aod.py "$MID" || true

GEE_AEROSOL_CONTEXT_OUTDIR=gee_aerosol_context_train1000 \
  "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/collect_gee_aerosol_context.py "$MID" || true

VIIRS_AOD_OUTDIR_SUFFIX=train1000 \
  "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/collect_viirs_aod.py "$MID" || true

#!/bin/bash
#SBATCH --job-name=c250monpre
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-250%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250monpre_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250monpre_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export ENDPOINT=pc
export OUTDIR=l2a_seasonal_monthly

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

if [[ -f "$ROOT/l2a_seasonal_monthly/${MID}.npz" ]]; then
  echo "L2Amon ${MID} SKIP existing monthly composite"
  exit 0
fi

cd "$ROOT"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/build_l2a_monthly.py" "$MID"

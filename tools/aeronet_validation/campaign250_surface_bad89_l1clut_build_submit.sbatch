#!/bin/bash
#SBATCH --job-name=b89l1clut
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-89%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/b89l1clut_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/b89l1clut_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export AC_BACKEND=lut
export OUTDIR=l1c_seasonal_lut
export KEEP_FRAC=0.6

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_surface_bad89_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" "$ROOT/build_l1c_seasonal.py" "$MID"

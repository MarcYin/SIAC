#!/bin/bash
#SBATCH --job-name=c250l2apcfb
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-250%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250l2apcfb_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250l2apcfb_%A_%a.out
set -euo pipefail

export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export ENDPOINT=pc
export OUTDIR=l2a_seasonal_pc

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

if [[ -f "$ROOT/l2a_seasonal_pc/${MID}.npz" ]]; then
  echo "L2ApcFallback ${MID} SKIP existing seasonal PC composite"
  exit 0
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/build_l2a_monthly_no_maiac_fallback.py "$MID"

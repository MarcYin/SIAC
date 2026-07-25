#!/bin/bash
#SBATCH --job-name=c250monnomaiac
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-6%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250monnomaiac_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250monnomaiac_%A_%a.out
set -euo pipefail

export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export ENDPOINT=pc
export OUTDIR=l2a_seasonal_monthly

MIDS=/home/users/marcyin/SIAC/tools/aeronet_validation/campaign250_l2a_monthly_no_maiac_mids.txt
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")
if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/build_l2a_monthly_no_maiac_fallback.py "$MID"

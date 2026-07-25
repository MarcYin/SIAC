#!/bin/bash
#SBATCH --job-name=c250s2wvp
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-250%16
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250s2wvp_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250s2wvp_%A_%a.out
set -euo pipefail

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export S2_L2A_WVP_OUTDIR=s2_l2a_wvp_campaign250

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
for attempt in 1 2; do
  if "$PY" -u /home/users/marcyin/SIAC/tools/aeronet_validation/collect_s2_l2a_wvp.py "$MID"; then
    exit 0
  fi
  echo "Retrying $MID after attempt $attempt" >&2
  sleep $((attempt * 15))
done
exit 1

#!/bin/bash
#SBATCH --job-name=c250pq
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-250%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250pq_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250pq_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

OUT="$ROOT/prior_quality/${MID}.json"
if [[ -f "$OUT" && -z "${FORCE:-}" ]]; then
  echo "prior_quality exists: $OUT"
  exit 0
fi

cd /home/users/marcyin/SIAC
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
mkdir -p "$ROOT/prior_quality" "$ROOT/slurm"

T0=$SECONDS
pixi run python -u "$ROOT/diag_prior_quality.py" "$MID" "$ROOT/prior_quality"
if [[ ! -f "$OUT" ]]; then
  echo "Missing expected prior_quality output: $OUT" >&2
  exit 4
fi
echo "[$MID prior_quality elapsed $((SECONDS-T0))s]"

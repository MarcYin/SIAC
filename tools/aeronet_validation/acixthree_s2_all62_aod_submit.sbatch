#!/bin/bash
#SBATCH --job-name=acixs2aod62
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-61%4
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acixs2aod62_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acixs2aod62_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="${MIDS:-$ROOT/all62_mids.txt}"
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MIDS")
if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

cd /home/users/marcyin/SIAC
export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
mkdir -p "$ROOT/slurm" "$ROOT/acixthree_s2_all62_aod"

/home/users/marcyin/SIAC/.pixi/envs/default/bin/python -u \
  tools/aeronet_validation/acixthree_s2_replicate.py "$MID" \
  --hyper-backend none \
  --input-resolution-m 30 \
  --output-root "$ROOT/acixthree_s2_all62_aod" \
  --surface-dir "$ROOT/l1c_seasonal_species_lut" \
  --overwrite

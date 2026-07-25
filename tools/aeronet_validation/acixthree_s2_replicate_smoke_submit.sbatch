#!/bin/bash
#SBATCH --job-name=acixs2smoke
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acixs2smoke_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acixs2smoke_%j.out
set -euo pipefail

MID="${1:-Bambey-ISRA__T28PCB_20240610T113321}"
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor

cd /home/users/marcyin/SIAC
export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
mkdir -p "$ROOT/slurm" "$ROOT/acixthree_s2"

/home/users/marcyin/SIAC/.pixi/envs/default/bin/python -u \
  tools/aeronet_validation/acixthree_s2_replicate.py "$MID" \
  --hyper-backend 6s \
  --input-resolution-m 30 \
  --output-root "$ROOT/acixthree_s2" \
  --surface-dir "$ROOT/l1c_seasonal_species_lut" \
  --overwrite

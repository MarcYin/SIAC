#!/bin/bash
#SBATCH --job-name=l1cpqexp
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l1cpqexp_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l1cpqexp_%j.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
cd /home/users/marcyin/SIAC

/home/users/marcyin/SIAC/.pixi/envs/default/bin/python -u \
  tools/aeronet_validation/evaluate_l1c_prior_quality.py \
  --prior "baseline=$ROOT/l1c_seasonal_species_lut" \
  --prior "exact_wvp=$ROOT/l1c_seasonal_species_exactwvp_rel60_e1" \
  --prior "clean015_min3=$ROOT/l1c_seasonal_species_exactwvp_clean015min3_e2" \
  --output "$ROOT/l1c_prior_quality_e2_20260709.json"


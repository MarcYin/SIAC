#!/bin/bash
#SBATCH --job-name=surfrobust
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/surfrobust_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/surfrobust_%j.out
set -euo pipefail

cd /home/users/marcyin/SIAC
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

/home/users/marcyin/SIAC/.pixi/envs/default/bin/python -u \
  tools/aeronet_validation/evaluate_robust_surface_likelihood.py \
  --output /gws/ssde/j25a/nceo_isp/public/siac_refactor/robust_surface_likelihood_e3_20260709.json


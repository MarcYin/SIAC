#!/bin/bash
#SBATCH --job-name=l2atrain
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2atrain_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2atrain_%j.out
set -euo pipefail
umask 0022

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/train_l2a_l1c_harmonizer.py
"$PY" "$SCRIPT" --training-cutoff 2023-12-31

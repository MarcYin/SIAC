#!/bin/bash
#SBATCH --job-name=l2aterrainfit
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainfit_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainfit_%j.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
PAIRS=$ROOT/analysis/l2a_l1c_exact_pairs_terrain_mediumdev_20260714
ALL_OUTPUT=$ROOT/analysis/l2a_l1c_harmonizer_terrain_mediumdev_20260714
CLEAN_OUTPUT=$ROOT/analysis/l2a_l1c_harmonizer_terrain_clean015_mediumdev_20260714
FOLD_MAP=$ROOT/analysis/l2a_l1c_harmonizer_clean015_mediumdev_20260714/harmonizer.json

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/train_l2a_l1c_harmonizer.py

"$PY" "$SCRIPT" --pairs "$PAIRS" --output "$ALL_OUTPUT"
"$PY" "$SCRIPT" \
  --pairs "$PAIRS" \
  --output "$CLEAN_OUTPUT" \
  --training-maiac-aot-max 0.15 \
  --fold-map "$FOLD_MAP"

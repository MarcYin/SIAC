#!/bin/bash
#SBATCH --job-name=xrt_harm_fit
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
# The short QoS permits one CPU per node; the capped ridge fits are small.
#SBATCH --cpus-per-task=1
# Prior uncapped harmonizer fits peaked at 3.0 GB RSS.
#SBATCH --mem=16G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_harm_fit_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_harm_fit_%j.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
PAIRS=${CROSS_RT_PAIRS:-$ROOT/analysis/l2a_l1c_physical_pairs_lowcloud152_20260716}
CASES=${CROSS_RT_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
OUTPUT=${CROSS_RT_OUTPUT:-$ROOT/analysis/cross_rt_harmonizer_lowcloud152_20260716}

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/train_l2a_l1c_harmonizer.py

"$PY" "$SCRIPT" \
  --pairs "$PAIRS" \
  --cases "$CASES" \
  --output "$OUTPUT" \
  --training-cutoff 2023-12-31 \
  --evaluation-split-manifest \
  "$ROOT/analysis/l2a_harmonization_retrieval_score_20260714/summary.json" \
  --model-family cross-rt \
  --max-samples-per-scene 300

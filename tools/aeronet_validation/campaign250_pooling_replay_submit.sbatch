#!/bin/bash
#SBATCH --job-name=c250poolreplay
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-250%2
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250poolreplay_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250poolreplay_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
export POLARS_MAX_THREADS=1
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
OUT="$ROOT/campaign250_pooling_replay_20260711"
attempt=1
max_attempts=3
while true; do
  if "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/analyze_sparse_support_pooling.py \
    --comparison /home/users/marcyin/SIAC/notes/reports/masked-r2-regression-20260711/comparison.json \
    --rows-csv /home/users/marcyin/SIAC/notes/reports/masked-r2-regression-20260711/comparison-rows.csv \
    --cohort all \
    --output-dir "$OUT" \
    --index "$SLURM_ARRAY_TASK_ID"; then
    exit 0
  fi
  status=$?
  if (( attempt >= max_attempts )); then
    echo "Task failed after $attempt attempts: index=$SLURM_ARRAY_TASK_ID status=$status" >&2
    exit "$status"
  fi
  sleep_for=$((60 + (SLURM_ARRAY_TASK_ID % 60) + (attempt - 1) * 120))
  echo "Attempt $attempt failed with status=$status; retrying in ${sleep_for}s" >&2
  sleep "$sleep_for"
  attempt=$((attempt + 1))
done

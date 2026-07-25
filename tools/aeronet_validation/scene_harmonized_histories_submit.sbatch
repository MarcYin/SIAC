#!/bin/bash
#SBATCH --job-name=scenehist
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-36%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/scenehist_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/scenehist_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=$ROOT/reports/aod-medium-physics-20260713/downloads/cases.csv
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES" | cut -d, -f1)
if [[ -z "${MID:-}" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT=python-gee@gee-marc.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/home/users/marcyin/gee-service-account.json
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_scene_harmonized_l2a_histories.py
for attempt in 1 2 3; do
  if "$PY" "$SCRIPT" "$MID"; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((30 * attempt + SLURM_ARRAY_TASK_ID % 20))
  fi
done
echo "Scene harmonized history failed after three attempts: $MID" >&2
exit 1

#!/bin/bash
#SBATCH --job-name=xrt_nl_hist
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
# All numerical libraries are pinned to one thread below.
#SBATCH --cpus-per-task=1
# Ridge history builds peaked below 0.5 GB; the joblib fold models add ~1 GB.
#SBATCH --mem=6G
#SBATCH --array=1-152%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_nl_hist_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_nl_hist_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${CROSS_RT_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
PAIRS=${CROSS_RT_PAIRS:-$ROOT/analysis/l2a_l1c_physical_pairs_lowcloud152_20260716}
MODEL_ROOT=${CROSS_RT_NL_OUTPUT:-$ROOT/analysis/cross_rt_nonlinear_targetband_lowcloud152_20260716}
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES" | cut -d, -f1)
if [[ -z "${MID:-}" || "$MID" == "matchup_id" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID in $CASES" >&2
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
export GDAL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_harmonized_l2a_histories.py
for attempt in 1 2 3; do
  if "$PY" "$SCRIPT" "$MID" \
    --pairs "$PAIRS" \
    --model "$MODEL_ROOT/harmonizer.json" \
    --output "$MODEL_ROOT/daily_histories" \
    --cases "$CASES" \
    --history-cutoff 2023-12-31 \
    --model-scope crossfit \
    --candidate hgb_target_band:blue:0.030 \
    --candidate hgb_target_band:solver:0.015 \
    --candidate hgb_target_band:solver:0.030; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((45 * attempt + SLURM_ARRAY_TASK_ID % 30))
  fi
done
echo "Nonlinear cross-RT history construction failed after three attempts: $MID" >&2
exit 1

#!/bin/bash
#SBATCH --job-name=l2apairmap
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-36%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairmap_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairmap_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${PAIR_AFFINE_CASES:-$ROOT/reports/aod-medium-physics-20260713/downloads/cases.csv}
CASE_INDEX=${PAIR_AFFINE_CASE_INDEX:-$SLURM_ARRAY_TASK_ID}
MID=$(sed -n "$((CASE_INDEX + 1))p" "$CASES" | cut -d, -f1)
OUTPUT_REL=${PAIR_AFFINE_OUTPUT_REL:-analysis/l2a_l1c_per_acquisition_affine_20260715/daily_histories}
PAIRS_REL=${PAIR_AFFINE_PAIRS_REL:-analysis/l2a_l1c_exact_pairs_mediumdev_20260713}
PAIR_MAIAC_AOT_MAX=${PAIR_AFFINE_PAIR_MAIAC_AOT_MAX:-}
if [[ -z "${MID:-}" || "$MID" == "matchup_id" ]]; then
  echo "No matchup for task $CASE_INDEX in $CASES" >&2
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
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_per_acquisition_l2a_l1c_histories.py
PAIR_AOT_ARGS=()
if [[ -n "$PAIR_MAIAC_AOT_MAX" ]]; then
  PAIR_AOT_ARGS=(--pair-maiac-aot-max "$PAIR_MAIAC_AOT_MAX")
fi
for attempt in 1 2 3; do
  if "$PY" "$SCRIPT" "$MID" \
    --pairs "$ROOT/$PAIRS_REL" \
    "${PAIR_AOT_ARGS[@]}" \
    --output "$ROOT/$OUTPUT_REL" \
    --candidate solver:0.015 \
    --candidate solver:0.030 \
    --candidate visible:0.015 \
    --candidate visible:0.030 \
    --candidate all:0.030; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((30 * attempt + CASE_INDEX % 20))
  fi
done
echo "Same-acquisition affine history failed after three attempts: $MID" >&2
exit 1

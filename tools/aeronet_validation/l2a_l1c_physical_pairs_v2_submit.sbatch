#!/bin/bash
#SBATCH --job-name=l2apairv2
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=04:00:00
# The builder and LUT backend are explicitly single-threaded below.
#SBATCH --cpus-per-task=1
# Peak RSS in the 152-case physical archive build was 1.2 GB; retain ample headroom.
#SBATCH --mem=6G
#SBATCH --array=1-152%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairv2_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2apairv2_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${PHYSICAL_PAIR_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
CASE_INDEX=${PHYSICAL_PAIR_CASE_INDEX:-$SLURM_ARRAY_TASK_ID}
OUTPUT_REL=${PHYSICAL_PAIR_OUTPUT_REL:-analysis/l2a_l1c_physical_pairs_v2_lowcloud152_20260717}
SHARDS=${PHYSICAL_PAIR_SHARDS:-0}
MIDS=()
if (( SHARDS > 0 )); then
  # Optional sharding keeps the submitted array small on a congested account.
  # Every case index is assigned exactly once by index modulo shard count.
  CASE_COUNT=$(( $(wc -l < "$CASES") - 1 ))
  for (( index=SLURM_ARRAY_TASK_ID; index<=CASE_COUNT; index+=SHARDS )); do
    MIDS+=("$(sed -n "$((index + 1))p" "$CASES" | cut -d, -f1)")
  done
else
  MIDS+=("$(sed -n "$((CASE_INDEX + 1))p" "$CASES" | cut -d, -f1)")
fi
if (( ${#MIDS[@]} == 0 )) || [[ -z "${MIDS[0]:-}" || "${MIDS[0]}" == "matchup_id" ]]; then
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
export SIAC_LUT_SCENE_GEOMETRY_MODE=nearest

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_l2a_l1c_harmonization_pairs.py
for attempt in 1 2 3; do
  if "$PY" "$SCRIPT" "${MIDS[@]}" \
    --target-rt-mode physical_v2 \
    --include-terrain \
    --scene-day-max "${PHYSICAL_PAIR_SCENE_DAY_MAX:-2023-12-31}" \
    --output-dir "$ROOT/$OUTPUT_REL"; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((60 * attempt + CASE_INDEX % 30))
  fi
done
echo "Physical pair construction failed after three attempts: ${MIDS[*]}" >&2
exit 1

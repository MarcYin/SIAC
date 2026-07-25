#!/bin/bash
#SBATCH --job-name=lc20modisdb
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-152%2
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20modisdb_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/lc20modisdb_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
DATA_ROOT=${DATA_ROOT:-"$ROOT/experiments/lowcloud20_modis_monthly_database_20260711"}
APPROACH=${APPROACH:-monthly_database}
MIDS="$ROOT/campaign250_lowcloud20_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export SIAC_EDOWN_EXECUTABLE=/home/users/marcyin/SIAC/tools/edown_runtime
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
export POLARS_MAX_THREADS=1
export GDAL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
export TMPDIR="$ROOT/spx_tmp_lowcloud20_modis_monthly_database/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

echo "TASK $SLURM_ARRAY_TASK_ID -> mid=$MID approach=$APPROACH"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
RESULT_JSON="$DATA_ROOT/runs/$APPROACH/$MID/result.json"
for attempt in 1 2 3; do
  "$PY" -m tools.aeronet_validation.cli \
    --data-root "$DATA_ROOT" \
    run \
    --matchup-id "$MID" \
    --approaches "$APPROACH" \
    --max-workers 1 \
    --no-rasters

  status="missing"
  if [[ -f "$RESULT_JSON" ]]; then
    status=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$RESULT_JSON")
  fi
  if [[ "$status" == "ok" || "$status" == "no_valid_observation" ]]; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep_for=$((60 + SLURM_ARRAY_TASK_ID % 60 + (attempt - 1) * 120))
    echo "Attempt $attempt ended with status=$status; retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "Retrieval failed after three attempts: mid=$MID status=$status" >&2
exit 1

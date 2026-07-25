#!/bin/bash
#SBATCH --job-name=extmaiac
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --array=1-1000%12
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/extmaiac_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/extmaiac_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/${EXTERNAL_MAIAC_MIDS_FILE:-external_train1000_mids.txt}"
OUT_DIR="${EXTERNAL_MAIAC_OUT_DIR:-maiac_qa_external_balanced}"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

RESULT="$ROOT/$OUT_DIR/$MID.json"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python

if [[ -f "$RESULT" ]] && "$PY" -c \
  'import json,sys; rec=json.load(open(sys.argv[1])); assert "aot" in rec' \
  "$RESULT" 2>/dev/null; then
  echo "Existing terminal MAIAC record: $MID"
  exit 0
fi

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
for attempt in 1 2 3; do
  status=0
  "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/stage_maiac_prior.py \
    --root "$ROOT" \
    --output-dir "$ROOT/$OUT_DIR" \
    --cache-dir "$ROOT/cache/mcd19_native_stage" \
    --aggregation "${EXTERNAL_MAIAC_AGGREGATION:-median}" \
    "$MID" || status=$?
  if [[ -f "$RESULT" ]] && "$PY" -c \
    'import json,sys; rec=json.load(open(sys.argv[1])); assert "aot" in rec' \
    "$RESULT" 2>/dev/null; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep_for=$((30 + SLURM_ARRAY_TASK_ID % 30 + (attempt - 1) * 60))
    echo "Attempt $attempt failed (exit=$status); retrying in ${sleep_for}s" >&2
    sleep "$sleep_for"
  fi
done

echo "MAIAC staging failed after three attempts: $MID" >&2
exit 1

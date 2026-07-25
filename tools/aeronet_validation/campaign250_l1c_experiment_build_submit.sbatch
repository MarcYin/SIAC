#!/bin/bash
#SBATCH --job-name=c250priorx
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-250%4
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250priorx_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250priorx_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
TAG="${C250_PRIOR_TAG:-species_exactwvp_rel60}"
MIDS="${C250_PRIOR_MIDS:-$ROOT/campaign250_mids.txt}"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=8
export AC_BACKEND="${AC_BACKEND:-species_exact_wvp}"
export OUTDIR="${OUTDIR:-l1c_seasonal_$TAG}"
export WORKDIR="${WORKDIR:-l1c_test/seasonal}"
export KEEP_FRAC="${KEEP_FRAC:-0.6}"
export CLEAN_AOD_MAX="${CLEAN_AOD_MAX:-0}"
export CLEAN_MIN_SCENES="${CLEAN_MIN_SCENES:-3}"
export CLEAN_WINDOW_MAX_MEDIAN_AOD="${CLEAN_WINDOW_MAX_MEDIAN_AOD:-0}"
export CLEAN_WINDOW_MIN="${CLEAN_WINDOW_MIN:-0}"
export TMPDIR="$ROOT/spx_tmp_prior_${TAG}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$ROOT/$OUTDIR" "$ROOT/slurm" "$TMPDIR"

OUT="$ROOT/$OUTDIR/${MID}.npz"
if [[ -f "$OUT" && -z "${FORCE:-}" ]]; then
  echo "prior exists: $OUT"
  exit 0
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
T0=$SECONDS
"$PY" -u "$ROOT/build_l1c_seasonal.py" "$MID"
if [[ ! -f "$OUT" ]]; then
  echo "Missing expected prior: $OUT" >&2
  exit 4
fi
echo "[$MID $TAG elapsed $((SECONDS-T0))s]"

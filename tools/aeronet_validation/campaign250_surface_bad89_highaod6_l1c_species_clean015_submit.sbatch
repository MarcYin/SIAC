#!/bin/bash
#SBATCH --job-name=c250highaod6p
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --array=1-6%3
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6p_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6p_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_surface_bad89_highaod6_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")
OUTDIR="l1c_seasonal_species_lut_clean015_highaod6"
OUT="$ROOT/$OUTDIR/${MID}.npz"

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export TMPDIR="$ROOT/spx_tmp_clean015_highaod6"
export KEEP_FRAC=0.6
export CLEAN_AOD_MAX=0.15
export CLEAN_MIN_SCENES=3
export CLEAN_WINDOW_MIN=0
export AC_BACKEND=species_lut
export OUTDIR
export OMP_NUM_THREADS=8
export WORKDIR="l1c_test/seasonal_clean015_highaod6"

mkdir -p "$ROOT/$OUTDIR" "$ROOT/slurm" "$TMPDIR"

if [[ -f "$OUT" && -z "${FORCE:-}" ]]; then
  echo "clean015 prior exists: $OUT"
  exit 0
fi

cd "$ROOT"
PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" -u "$ROOT/build_l1c_seasonal.py" "$MID"
if [[ ! -f "$OUT" ]]; then
  echo "Missing expected clean015 prior: $OUT" >&2
  exit 4
fi
echo "[$MID clean015 prior done]"

#!/bin/bash
#SBATCH --job-name=sdspecprior
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=0-61%4
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/sdspecprior_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/sdspecprior_%A_%a.out
set -euo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
TAG="${REPRO_TAG:-20260708}"
MIDS="${MIDS:-$ROOT/all62_mids.txt}"
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

cd "$ROOT"
export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"
export TMPDIR="$ROOT/spx_tmp_sdspec_repro_$TAG"
export KEEP_FRAC=0.6
export AC_BACKEND=species_lut
export OUTDIR="l1c_seasonal_species_lut_repro_$TAG"
export WORKDIR="l1c_test/seasonal_sdspec_repro_$TAG"
export OMP_NUM_THREADS=8
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
mkdir -p "$ROOT/$OUTDIR" "$ROOT/slurm" "$TMPDIR"

OUT="$ROOT/$OUTDIR/${MID}.npz"
if [[ -f "$OUT" && -z "${FORCE:-}" ]]; then
  echo "sdspec repro prior exists: $OUT"
  exit 0
fi

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
T0=$SECONDS
"$PY" -u "$ROOT/build_l1c_seasonal.py" "$MID"
if [[ ! -f "$OUT" ]]; then
  echo "Missing expected sdspec repro prior: $OUT" >&2
  exit 4
fi
echo "[$MID sdspec prior repro elapsed $((SECONDS-T0))s]"

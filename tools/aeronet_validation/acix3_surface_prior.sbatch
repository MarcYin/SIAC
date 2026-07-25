#!/bin/bash
#SBATCH --job-name=acix3_prior
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-20%10
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acix3_prior_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acix3_prior_%A_%a.out
set -euo pipefail
umask 0022

# acixThree-style surface prior reimplementation (GEE qualityMosaic library +
# kNN predictor + self-prediction RMSE sigma), evaluated truth-referenced
# against the current et20 prior on the 20 calibration-dump example cases.
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS=$ROOT/calib_dump_examples_mids.txt
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")
if [[ -z "${MID:-}" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

export GEE_SERVICE_ACCOUNT=python-gee@gee-marc.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/home/users/marcyin/gee-service-account.json
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export TMPDIR="$ROOT/acix3_prior_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/acix3_surface_prior.py

mosaic_file="$ROOT/acix3_style_mosaics_20260719/${MID}.npz"
if [[ ! -f "$mosaic_file" ]]; then
  for attempt in 1 2 3; do
    if "$PY" "$SCRIPT" mosaics "$MID"; then
      break
    fi
    if (( attempt == 3 )); then
      echo "Mosaic build failed after three attempts: $MID" >&2
      exit 1
    fi
    sleep $((60 * attempt))
  done
fi
"$PY" "$SCRIPT" predict "$MID"

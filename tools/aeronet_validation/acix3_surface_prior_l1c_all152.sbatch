#!/bin/bash
#SBATCH --job-name=acix3l1c152
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1074,host1142,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --array=1-152%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acix3l1c152_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/acix3l1c152_%A_%a.out
set -euo pipefail
umask 0022

# Faithful acixThree L1C mosaics for the FULL lowcloud152 campaign. Task i maps
# to all-cases.csv row i+1 — the same mapping as cross_rt_calib_dump_all152 —
# so this array is submitted with --dependency=aftercorr:<dump job> and each
# task starts as soon as its site's calibration dump exists. %6 respects the
# GEE concurrency limit that killed Churchill in the Sen2Cor run.
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv
MID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES" | cut -d, -f1)
if [[ -z "${MID:-}" ]]; then
  echo "No matchup for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

DUMP=$ROOT/calib_dumps_crossrt_dt_20260719/${MID}.npz
if [[ ! -f "$DUMP" ]]; then
  echo "SKIP $MID: calibration dump missing (dump task failed?)" >&2
  exit 0
fi

export GEE_SERVICE_ACCOUNT=python-gee@gee-marc.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/home/users/marcyin/gee-service-account.json
export PYTHONPATH=/gws/ssde/j25a/nceo_isp/public/siac_env/bp_deps:/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export PHASE_D_LUT_PATH=/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip
export TMPDIR="$ROOT/acix3_prior_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/acix3_surface_prior.py

mosaic_file="$ROOT/acix3_style_mosaics_20260719/${MID}_l1c.npz"
if [[ ! -f "$mosaic_file" ]]; then
  for attempt in 1 2 3; do
    if "$PY" "$SCRIPT" mosaics_l1c "$MID"; then
      break
    fi
    if (( attempt == 3 )); then
      echo "L1C mosaic build failed after three attempts: $MID" >&2
      exit 1
    fi
    sleep $((90 * attempt))
  done
fi
"$PY" "$SCRIPT" predict_l1c "$MID" || echo "predict_l1c failed (non-fatal): $MID" >&2
exit 0

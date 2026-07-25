#!/bin/bash
#SBATCH --job-name=c250surfx
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-250%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250surfx_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250surfx_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
TAG="${SURFACE_EXP_TAG:-e1_l2awvp}"
MIDS="${SURFACE_EXP_MIDS:-$ROOT/campaign250_mids.txt}"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")
PRIOR_DIR="${SURFACE_EXP_PRIOR_DIR:-l1c_seasonal_species_lut}"

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
if [[ ! -f "$ROOT/$PRIOR_DIR/${MID}.npz" ]]; then
  echo "Missing experiment prior: $ROOT/$PRIOR_DIR/${MID}.npz" >&2
  exit 3
fi
if [[ ! -f "$ROOT/prior_quality/${MID}.json" ]]; then
  echo "Missing prior-quality sidecar: $ROOT/prior_quality/${MID}.json" >&2
  exit 3
fi

cd /home/users/marcyin/SIAC
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export TMPDIR="$ROOT/spx_tmp_surface_${TAG}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_surface_${TAG}"
export RETR_SEAS_DIR="$PRIOR_DIR"
export DEBIAS_CSV=debias_per_site_l1c.csv
export RETR_TAG="_${TAG}"
export AERO_PRIOR=maiac
export MAIAC_UNC=loose
export ANCHOR_ITERS=1
export ANCHOR_RT=1
export SURFDRV_ITERS="${SURFDRV_ITERS:-1}"
export SURFDRV_ANCHOR_INIT="${SURFDRV_ANCHOR_INIT:-cams}"
export ANCHOR_BAND_SET="${ANCHOR_BAND_SET:-b8a}"
export LOCALISER=meanvis
export SOLVER=surfdrv
export SPECIES_N="${SPECIES_N:-3}"
export AERO_SPECIES_MODEL=cci
export SPECIES_BACKEND=sixs
export REF="${REF:-composite}"
export PREDICTOR="${PREDICTOR:-tree}"
export SWIR_FEATS="${SWIR_FEATS:-B11,B12}"
export BACKSTOP="${BACKSTOP:-add}"
export GATE="${GATE:-0}"
export BACKSTOP_CAL=1
export COST_MODE="${COST_MODE:-auto2}"
export AOD_CLEAN=0.15
export SOLVE_BANDS="${SOLVE_BANDS:-B01,B02,B04}"
export AOD_HIGH=0.6
export AOD_AXIS_MAX="${AOD_AXIS_MAX:-4.0}"
export TARGET_TCWV_MODE="${TARGET_TCWV_MODE:-l2a}"
export TARGET_TCWV_FIXED="${TARGET_TCWV_FIXED:-2.0}"
export TARGET_TCWV_DIR="${TARGET_TCWV_DIR:-s2_l2a_wvp_campaign250}"
export PRIOR_MODEL_ERROR="${PRIOR_MODEL_ERROR:-0.0}"
export PRIOR_CONTAM_PENALTY="${PRIOR_CONTAM_PENALTY:-0.0}"
export ROBUST_CLIP="${ROBUST_CLIP:-0}"
export SURFDRV_AEROSOL_PROFILE="${SURFDRV_AEROSOL_PROFILE:-}"

mkdir -p "$ROOT/seasonal_val" "$ROOT/slurm" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT" "$SIAC_RT_RUN_CACHE_ROOT/rt6s" \
  "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache" 2>/dev/null || true

OUTJSON="$ROOT/seasonal_val/${MID}__seas_maiac_${TAG}.json"
if [[ -f "$OUTJSON" && -z "${FORCE:-}" ]]; then
  echo "surface experiment result exists: $OUTJSON"
  exit 0
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> $MID ($TAG prior=$PRIOR_DIR tcwv=$TARGET_TCWV_MODE)"
for attempt in 1 2; do
  pixi run -e rt6s python -u "$ROOT/diag_seasonal_retrieval.py" "$MID" 2>&1 \
    | grep -aE "^RET |surfdrv\]|seas\]|maiac\]|tcwv\]|prior-unc\]|EXC|Traceback|Error" \
    | grep -av excepthook \
    | tail -30
  if [[ -f "$OUTJSON" ]]; then
    echo "[ok $MID $TAG attempt=$attempt]"
    exit 0
  fi
  echo "[retry $MID $TAG attempt=$attempt]"
  sleep $((attempt * 20))
done

echo "[FAILED $MID $TAG]"
exit 1

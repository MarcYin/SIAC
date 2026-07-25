#!/bin/bash
#SBATCH --job-name=c250highaod6r
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-6%3
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6r_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250highaod6r_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_surface_bad89_highaod6_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export OMP_NUM_THREADS=4
export TMPDIR="$ROOT/spx_tmp_sdspec_highaod6_clean015"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_sdspec_highaod6_clean015"
export RETR_SEAS_DIR="l1c_seasonal_species_lut_clean015_highaod6"
export DEBIAS_CSV=debias_per_site_l1c.csv
export RETR_TAG="_sdspec_highaod6_clean015"
export AERO_PRIOR=maiac
export MAIAC_UNC=loose
export ANCHOR_ITERS=1
export ANCHOR_RT=1
export LOCALISER=meanvis
export SOLVER=surfdrv
export SPECIES_N=3
export AERO_SPECIES_MODEL=cci
export SPECIES_BACKEND=sixs
export REF=composite
export BACKSTOP=add
export GATE=0
export BACKSTOP_CAL=1
export COST_MODE=auto2
export AOD_CLEAN=0.15
export SOLVE_BANDS=B01,B02,B04
export AOD_HIGH=0.6
export AOD_AXIS_MAX=4.0
export TARGET_TCWV_MODE=l2a
export TARGET_TCWV_FIXED=2.0
export TARGET_TCWV_DIR=s2_l2a_wvp_campaign250
export PRIOR_MODEL_ERROR=0.0
export PRIOR_CONTAM_PENALTY=0.0

mkdir -p "$ROOT/seasonal_val" "$ROOT/slurm" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT" "$SIAC_RT_RUN_CACHE_ROOT/rt6s" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache" 2>/dev/null || true

OUTJSON="$ROOT/seasonal_val/${MID}__seas_maiac_sdspec_highaod6_clean015.json"
if [[ -f "$OUTJSON" && -z "${FORCE:-}" ]]; then
  echo "high-aod6 clean015 result exists: $OUTJSON"
  exit 0
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> $MID (acixThree-style clean015+L2A)"
for attempt in 1 2; do
  pixi run -e rt6s python -u "$ROOT/diag_seasonal_retrieval.py" "$MID" 2>&1 \
    | grep -aE "^RET |surfdrv\]|seas\]|maiac\]|tcwv\]|prior-unc\]|EXC|Traceback|Error" \
    | grep -av excepthook \
    | tail -20
  if [[ -f "$OUTJSON" ]]; then
    echo "[ok $MID highaod6 clean015 attempt=$attempt]"
    exit 0
  fi
  echo "[retry $MID highaod6 clean015 attempt=$attempt]"
  sleep $((attempt * 20))
done

echo "[FAILED $MID highaod6 clean015]"
exit 1

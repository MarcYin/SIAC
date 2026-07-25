#!/bin/bash
#SBATCH --job-name=c250sdspec
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-250%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250sdspec_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250sdspec_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
TAG="${C250_SDSPEC_TAG:-20260708}"
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

PRIOR_DIR="l1c_seasonal_species_lut"
PRIOR="$ROOT/$PRIOR_DIR/${MID}.npz"
PQ="$ROOT/prior_quality/${MID}.json"
if [[ ! -f "$PRIOR" ]]; then
  echo "Missing acixThree species prior: $PRIOR" >&2
  exit 3
fi
if [[ ! -f "$PQ" ]]; then
  echo "Missing prior_quality sidecar: $PQ" >&2
  exit 3
fi

cd /home/users/marcyin/SIAC
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export TMPDIR="$ROOT/spx_tmp_sdspec_c250_$TAG"
export SIAC_RT_RUN_CACHE_ROOT="$ROOT/rt_run_cache_sdspec_c250_$TAG"
export RETR_SEAS_DIR="$PRIOR_DIR"
export DEBIAS_CSV=debias_per_site_l1c.csv
export RETR_TAG="_sdspec_c250_$TAG"
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

mkdir -p "$ROOT/seasonal_val" "$ROOT/slurm" "$TMPDIR" "$SIAC_RT_RUN_CACHE_ROOT"
chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT" 2>/dev/null || true
if [[ -d "$SIAC_RT_RUN_CACHE_ROOT/rt6s" ]]; then
  chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT/rt6s" 2>/dev/null || true
fi
mkdir -p "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache"
chmod u+rwx,g+rx,o+rx "$SIAC_RT_RUN_CACHE_ROOT/rt6s" "$SIAC_RT_RUN_CACHE_ROOT/rt6s/run_cache" 2>/dev/null || true

OUTJSON="$ROOT/seasonal_val/${MID}__seas_maiac_sdspec_c250_${TAG}.json"
if [[ -f "$OUTJSON" && -z "${FORCE:-}" ]]; then
  echo "campaign250 sdspec result exists: $OUTJSON"
  exit 0
fi

echo "TASK ${SLURM_ARRAY_TASK_ID} -> $MID (campaign250 sdspec $TAG)"
for attempt in 1 2; do
  pixi run -e rt6s python -u "$ROOT/diag_seasonal_retrieval.py" "$MID" 2>&1 \
    | grep -aE "^RET |surfdrv\\]|seas\\]|maiac\\]|EXC|Traceback|Error" \
    | grep -av excepthook \
    | tail -20
  if [[ -f "$OUTJSON" ]]; then
    echo "[ok $MID campaign250 sdspec attempt=$attempt]"
    exit 0
  fi
  echo "[retry $MID campaign250 sdspec attempt=$attempt]"
  sleep $((attempt * 20))
done

echo "[FAILED $MID campaign250 sdspec]"
exit 1

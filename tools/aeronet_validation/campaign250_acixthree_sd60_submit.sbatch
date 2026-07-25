#!/bin/bash
#SBATCH --job-name=c250sd60
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --requeue
#SBATCH --array=1-250%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250sd60_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/c250sd60_%A_%a.out
set -uo pipefail

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

if [[ -z "$MID" ]]; then
  echo "No matchup_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

SEAS="$ROOT/l1c_seasonal_species_lut/${MID}.npz"
PQ="$ROOT/prior_quality/${MID}.json"
TAG="_sd60_c250_20260707"
OUTJSON="$ROOT/seasonal_val/${MID}__seas_maiac${TAG}.json"

if [[ ! -f "$SEAS" ]]; then
  echo "Missing required acixThree species prior: $SEAS" >&2
  exit 3
fi
if [[ ! -f "$PQ" ]]; then
  echo "Missing required acixThree prior_quality sidecar: $PQ" >&2
  exit 3
fi
if [[ -f "$OUTJSON" && -z "${FORCE:-}" ]]; then
  echo "acixThree sd60 result exists: $OUTJSON"
  exit 0
fi

cd /home/users/marcyin/SIAC
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:/gws/ssde/j25a/nceo_isp/public/siac_refactor
export TMPDIR="$ROOT/spx_tmp"
export GEE_SERVICE_ACCOUNT="python-gee@gee-marc.iam.gserviceaccount.com"
export GEE_SERVICE_ACCOUNT_KEY="/home/users/marcyin/gee-service-account.json"

# Faithful acixThree sd60 setup: species-corrected L1C seasonal prior, direct
# composite reference, MAIAC loose aerosol prior/backstop, B01/B02/B04 solve,
# and standalone surface-driven solver via diag_seasonal_retrieval.py.
export RETR_SEAS_DIR=l1c_seasonal_species_lut
export DEBIAS_CSV=debias_per_site_l1c.csv
export RETR_TAG="$TAG"
export SOLVE_RES=60
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

mkdir -p "$ROOT/seasonal_val" "$ROOT/slurm" "$TMPDIR"
echo "TASK ${SLURM_ARRAY_TASK_ID} -> $MID (acixThree sd60 campaign-250)"
for attempt in 1 2; do
  pixi run -e rt6s python -u "$ROOT/diag_seasonal_retrieval.py" "$MID" 2>&1 \
    | grep -aE "^RET |surfdrv\\]|seas\\]|maiac\\]|EXC|Traceback|Error" \
    | grep -av excepthook \
    | tail -20
  if [[ -f "$OUTJSON" ]]; then
    echo "[ok $MID acixThree sd60 attempt=$attempt]"
    exit 0
  fi
  echo "[retry $MID acixThree sd60 attempt=$attempt]"
  sleep $((attempt * 20))
done

echo "[FAILED $MID acixThree sd60]"
exit 1

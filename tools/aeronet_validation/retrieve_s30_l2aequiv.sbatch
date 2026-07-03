#!/bin/bash
#SBATCH --job-name=s30eq
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --requeue
#SBATCH --array=0-61%12
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/logs/%A_%a_s30eq.out
set -uo pipefail

D=/gws/ssde/j25a/nceo_isp/public/siac_refactor
REPO=/home/users/marcyin/SIAC

cd "$REPO"
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export PYTHONPATH="$REPO:$REPO/python${PYTHONPATH:+:$PYTHONPATH}"
export TMPDIR="${TMPDIR:-$D/spx_tmp}"
mkdir -p "$TMPDIR" "$D/slurm/logs"

# Same retrieval settings as the L2A 85.5% rc15 run; only the surface dictionary changes.
export RETR_SEAS_DIR="${RETR_SEAS_DIR:-s30maiac_l2aequiv20260701}"
export DEBIAS_CSV="${DEBIAS_CSV:-debias_per_site_l1c.csv}"
export RETR_TAG="${RETR_TAG:-_s30maiac_l2aequiv20260701}"
export AERO_PRIOR=maiac MAIAC_UNC=loose ANCHOR_ITERS=1 ANCHOR_RT=1 LOCALISER=meanvis
export SOLVER=surfdrv SPECIES_N=3 AERO_SPECIES_MODEL=cci SPECIES_BACKEND=sixs
export REF=composite BACKSTOP=add GATE=0 BACKSTOP_CAL=1 COST_MODE=auto2 AOD_CLEAN=0.15
export SOLVE_BANDS=B01,B02,B04 AOD_HIGH=0.6 ROBUST_CLIP=1.5
export MIE_CACHE=1 ELEV=0

MID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SITES:-$D/all62_mids.txt}")
[ -z "$MID" ] && exit 0

OUTJSON="$D/seasonal_val/${MID}__seas_maiac${RETR_TAG}.json"
echo "TASK $SLURM_ARRAY_TASK_ID -> $MID (seas=$RETR_SEAS_DIR tag=$RETR_TAG)"
for attempt in 1 2; do
  pixi run -e rt6s python -u "$D/diag_seasonal_retrieval.py" "$MID" 2>&1 \
    | grep -aE "^RET |surfdrv\]|seas\]|EXC|Traceback|Error|timing" | grep -av excepthook | tail -12
  [ -f "$OUTJSON" ] && { echo "[ok $MID s30eq $attempt]"; exit 0; }
  echo "[retry $MID s30eq $attempt]"
  sleep $((attempt * 20))
done
echo "[FAILED $MID s30eq]"
exit 1

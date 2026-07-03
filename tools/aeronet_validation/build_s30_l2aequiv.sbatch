#!/bin/bash
#SBATCH --job-name=bs30eq
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --requeue
#SBATCH --array=0-61%8
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/logs/%A_%a_bs30eq.out
set -uo pipefail

D=/gws/ssde/j25a/nceo_isp/public/siac_refactor
REPO=/home/users/marcyin/SIAC

cd "$D"
export TMPDIR="${TMPDIR:-$D/spx_tmp}"
export S30_OUT="${S30_OUT:-s30maiac_l2aequiv20260701}"
export KEEP_FRAC="${KEEP_FRAC:-0.6}"
export PYTHONPATH="$D:$REPO/python${PYTHONPATH:+:$PYTHONPATH}"
export GEE_SERVICE_ACCOUNT="${GEE_SERVICE_ACCOUNT:-python-gee@gee-marc.iam.gserviceaccount.com}"
export GEE_SERVICE_ACCOUNT_KEY="${GEE_SERVICE_ACCOUNT_KEY:-/home/users/marcyin/gee-service-account.json}"
mkdir -p "$TMPDIR" "$D/slurm/logs"

MID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SITES:-$D/all62_mids.txt}")
[ -z "$MID" ] && exit 0
export MID

echo "TASK $SLURM_ARRAY_TASK_ID -> $MID (S30_OUT=$S30_OUT KEEP_FRAC=$KEEP_FRAC)"
"$REPO/.pixi/envs/default/bin/python" -u - <<'PY'
import os
import sys

sys.path.insert(0, "/gws/ssde/j25a/nceo_isp/public/siac_refactor")

import build_s30_maiac_seasonal as m

m.EDOWN = "/home/users/marcyin/SIAC/tools/edown_runtime"
m.b.EDOWN = m.EDOWN
m.main([os.environ["MID"]])
PY

#!/bin/bash
# Historical surface-prior dictionary per scene, in the cci_climatology_exact RT space.
#SBATCH --job-name=corpus-dict
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_dict_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_dict_%A_%a.out

set -uo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
OUT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901
PY=$PROJECT/.pixi/envs/default/bin/python
SHARDS=${SHARDS:-200}
SHARD=${SLURM_ARRAY_TASK_ID:-0}

export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_dict_${SLURM_JOB_ID}_${SHARD}
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT/l1c_teacher20_cci_exact" 2>/dev/null || true
test -d "$TMPDIR"

# This shard's scenes: captured, and not already built. The builder skips
# existing output itself, but filtering here keeps the argument list short.
MIDS=$($PY - "$SHARD" "$SHARDS" <<'PY'
import sys
from pathlib import Path
shard, shards = int(sys.argv[1]), int(sys.argv[2])
out = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901")
captured = sorted(p.stem for p in (out / "current_l1c_multiscale").glob("*.npz"))
done = {p.stem for p in (out / "l1c_teacher20_cci_exact").glob("*.npz")}
mine = [m for i, m in enumerate(captured) if i % shards == shard and m not in done]
print(" ".join(mine))
PY
)
if [ -z "$MIDS" ]; then echo "shard $SHARD: nothing to do"; exit 0; fi
echo "shard $SHARD: $(echo $MIDS | wc -w) scenes"

# --aerosol-species is NOT the default and the teacher rejects a mismatched
# dictionary only after the whole build, so it is set explicitly here.
$PY -u -m experiments.current_date_toa_prior.build_expanded_l1c_teacher $MIDS \
  --current-root "$OUT/current_l1c_multiscale" \
  --index-root "$OUT/cloudscore_index20" \
  --output-root "$OUT/l1c_teacher20_cci_exact" \
  --aerosol-species cci_climatology_exact \
  --grid detail20 \
  --sixs-threads 1 \
  --fetch-workers 8

echo "shard $SHARD done"
rm -rf "$TMPDIR"

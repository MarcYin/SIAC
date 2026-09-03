#!/bin/bash
#SBATCH --job-name=corpus-reduce
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_reduce_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_reduce_%A_%a.out

set -uo pipefail
umask 0022
PROJECT=/home/users/marcyin/SIAC
OUT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901
PY=$PROJECT/.pixi/envs/default/bin/python
SHARDS=${SHARDS:-25}
SHARD=${SLURM_ARRAY_TASK_ID:-0}
export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_reduce_${SLURM_JOB_ID}_${SHARD}
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT/l1c_teacher20_cci_exact_reduced" 2>/dev/null || true
test -d "$TMPDIR"

MIDS=$($PY - "$SHARD" "$SHARDS" <<'PY'
import sys
from pathlib import Path
shard, shards = int(sys.argv[1]), int(sys.argv[2])
out = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901")
source = sorted(p.stem for p in (out / "l1c_teacher20_cci_exact").glob("*.npz"))
done = {p.stem for p in (out / "l1c_teacher20_cci_exact_reduced").glob("*.npz")}
print(" ".join(m for i, m in enumerate(source) if i % shards == shard and m not in done))
PY
)
if [ -z "$MIDS" ]; then echo "shard $SHARD: nothing to do"; exit 0; fi
echo "shard $SHARD: $(echo $MIDS | wc -w) scenes"
$PY -u -m experiments.current_date_toa_prior.prepare_reduced_teacher $MIDS \
  --teacher-root "$OUT/l1c_teacher20_cci_exact" \
  --output-root "$OUT/l1c_teacher20_cci_exact_reduced"
echo "shard $SHARD done"
rm -rf "$TMPDIR"

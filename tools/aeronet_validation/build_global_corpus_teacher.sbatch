#!/bin/bash
# Optimal M5 teacher: the training labels. 20 m surface, recurrent-library snow.
#SBATCH --job-name=corpus-teacher
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=20:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_teacher_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_teacher_%A_%a.out

set -uo pipefail
umask 0022
PROJECT=/home/users/marcyin/SIAC
OUT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901
PY=$PROJECT/.pixi/envs/default/bin/python
SHARDS=${SHARDS:-400}
SHARD=${SLURM_ARRAY_TASK_ID:-0}
export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_teach_${SLURM_JOB_ID}_${SHARD}
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT/optimal_m5_teacher20" 2>/dev/null || true
test -d "$TMPDIR"

# Only scenes with every input present: the teacher requires all four and
# raises FileNotFoundError otherwise, which would abort the shard's whole list.
MIDS=$($PY - "$SHARD" "$SHARDS" <<'PY'
import sys
from pathlib import Path
shard, shards = int(sys.argv[1]), int(sys.argv[2])
out = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901")
def stems(name):
    return {p.stem for p in (out / name).glob("*.npz")}
ready = sorted(
    stems("current_l1c_multiscale")
    & stems("current_l1c_fine20_visible")
    & stems("l1c_teacher20_cci_exact")
    & stems("l1c_teacher20_cci_exact_reduced")
)
done = stems("optimal_m5_teacher20")
print(" ".join(m for i, m in enumerate(ready) if i % shards == shard and m not in done))
PY
)
if [ -z "$MIDS" ]; then echo "shard $SHARD: nothing to do"; exit 0; fi
echo "shard $SHARD: $(echo $MIDS | wc -w) scenes"

# Two non-default settings. --surface-grid 20m gives the label genuine sub-60 m
# structure at bit-identical AOD. --snow-support-policy recurrent-library covers
# snow present in the scene but not the composites, or the reverse, which is a
# gross outlier the solver cannot otherwise detect. RT space comes from the
# dictionary and must not be passed here.
$PY -u -m experiments.current_date_toa_prior.build_optimal_m5_teacher $MIDS \
  --current-root "$OUT/current_l1c_multiscale" \
  --fine20-root "$OUT/current_l1c_fine20_visible" \
  --dictionary-root "$OUT/l1c_teacher20_cci_exact" \
  --reduced-root "$OUT/l1c_teacher20_cci_exact_reduced" \
  --output-root "$OUT/optimal_m5_teacher20" \
  --surface-grid 20m \
  --snow-support-policy recurrent-library

echo "shard $SHARD done"
rm -rf "$TMPDIR"

#!/bin/bash
# Attach per-image sun/view geometry to the locally composed winner indices.
#SBATCH --job-name=index-geom
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/index_geom_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/index_geom_%A_%a.out

set -euo pipefail
umask 0022
PROJECT=/home/users/marcyin/SIAC
OUT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/global_corpus_20260901
SHARDS=${SHARDS:-1}
SHARD=${SLURM_ARRAY_TASK_ID:-0}

export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_geom_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" 2>/dev/null || true
test -d "$TMPDIR"

"$PROJECT/.pixi/envs/default/bin/python" -u -m tools.aeronet_validation.add_index_geometry \
  --index-root "$OUT/cloudscore_index20" \
  --shard "$SHARD" --shard-count "$SHARDS" \
  --summary "$OUT/index_geom_shard${SHARD}.summary.json"

rm -rf "$TMPDIR"

#!/bin/bash
# Per-AOI: scout SCL -> targeted Cloud Score+ fetch -> local winner index.
#SBATCH --job-name=corpus-index
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_index_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/corpus_index_%A_%a.out

set -euo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
OUT=$ROOT/global_corpus_20260901
SHARDS=${SHARDS:-1}
SHARD=${SLURM_ARRAY_TASK_ID:-0}
LIMIT=${LIMIT:-0}
NAME=${NAME:-corpus_index}
MAIAC=${MAIAC:---maiac}

# Cloud Score+ rasters are a bulky reproducible intermediate: scratch.
CACHE=/work/scratch-pw5/marcyin/global_corpus/cloudscore
INDEX=$OUT/cloudscore_index20

export PYTHONPATH=$PROJECT:$PROJECT/python
export SIAC_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_cache
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_idx_${SLURM_JOB_ID}
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$INDEX" "$CACHE" 2>/dev/null || true
test -d "$TMPDIR" && test -d "$INDEX" && test -d "$CACHE"

test -s "$OUT/global_scene_catalog.csv"

"$PROJECT/.pixi/envs/default/bin/python" -u -m tools.aeronet_validation.build_global_corpus_index \
  --scene-catalog "$OUT/global_scene_catalog.csv" \
  --output-root "$INDEX" \
  --cache-root "$CACHE" \
  --shard "$SHARD" --shard-count "$SHARDS" --limit "$LIMIT" $MAIAC \
  --summary "$OUT/${NAME}_shard${SHARD}.summary.json"

rm -rf "$TMPDIR"

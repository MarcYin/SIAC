#!/bin/bash
# Current-date acquisition per catalogue AOI, from Earth Engine L1C listings.
#SBATCH --job-name=global-scenes
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/global_scenes_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/global_scenes_%A_%a.out

set -euo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
OUT=$ROOT/global_corpus_20260901
SHARDS=${SHARDS:-1}
SHARD=${SLURM_ARRAY_TASK_ID:-0}
LIMIT=${LIMIT:-0}
NAME=${NAME:-global_scene_catalog}

export PYTHONPATH=$PROJECT:$PROJECT/python:/gws/ssde/j25a/nceo_isp/public/siac_env/edown_runtime
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_scenes_${SLURM_JOB_ID}
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT"

test -s "$OUT/global_aoi_catalog.csv"

"$PROJECT/.pixi/envs/default/bin/python" -u -m tools.aeronet_validation.build_global_scene_catalog \
  --catalog "$OUT/global_aoi_catalog.csv" \
  --output "$OUT/${NAME}_shard${SHARD}.csv" \
  --shard "$SHARD" --shard-count "$SHARDS" --limit "$LIMIT"

rm -rf "$TMPDIR"

#!/bin/bash
# Global AOI catalogue: MCD12Q1 land-cover tiles -> class-balanced pool -> stratified draw.
#SBATCH --job-name=global-catalog
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=16:00:00
# SKIP_FETCH=1 reuses the cached land-cover tiles (all 540 present).
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/global_catalog_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/global_catalog_%j.out

set -euo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
OUT=$ROOT/global_corpus_20260901
# Land-cover tiles are a bulky, reproducible intermediate: scratch, not GWS.
CACHE=/work/scratch-pw5/marcyin/global_corpus/landcover_tiles

export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_$SLURM_JOB_ID
export MPLCONFIGDIR=$TMPDIR/matplotlib
# Array tasks race on the shared parent; mkdir -p can still report EEXIST
# on this filesystem, which set -e would turn into a dead shard.
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT" 2>/dev/null || true
test -d "$TMPDIR" && test -d "$OUT" "$CACHE"

"$PROJECT/.pixi/envs/default/bin/python" -u -m tools.aeronet_validation.build_global_catalog \
  --output "$OUT/global_aoi_catalog.csv" \
  --cache "$CACHE" \
  --total "${CATALOG_TOTAL:-5000}" \
  --tile-degrees 10.0 \
  --per-tile 400 \
  --seed 20260901 ${SKIP_FETCH:+--skip-fetch}

rm -rf "$TMPDIR"

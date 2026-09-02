#!/bin/bash
# Re-encode the Cloud Score+ cache to uint8 and move it from scratch to the GWS.
#SBATCH --job-name=cs-archive
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/cs_archive_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/cs_archive_%A_%a.out

set -euo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
SHARDS=${SHARDS:-1}
SHARD=${SLURM_ARRAY_TASK_ID:-0}

export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/global_corpus/tmp_arch_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$ROOT/global_corpus_20260901/cloudscore_cache_uint8" 2>/dev/null || true
test -d "$TMPDIR"

"$PROJECT/.pixi/envs/default/bin/python" -u -m tools.aeronet_validation.reencode_cloudscore_cache \
  --source /work/scratch-pw5/marcyin/global_corpus/cloudscore \
  --target "$ROOT/global_corpus_20260901/cloudscore_cache_uint8" \
  --shard "$SHARD" --shard-count "$SHARDS" \
  --summary "$ROOT/global_corpus_20260901/cs_archive_shard${SHARD}.summary.json"

rm -rf "$TMPDIR"

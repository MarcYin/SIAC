#!/bin/bash
#SBATCH --job-name=medaodpool
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-152%6
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodpool_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodpool_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ROOT/campaign250_lowcloud20_mids.txt")

if [[ -z "${MID:-}" ]]; then
  echo "No matchup_id for case ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
export MALLOC_CONF=background_thread:false,narenas:1

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
for attempt in 1 2 3; do
  if "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/replay_pooled_surface_predictor.py \
    --matchup-id "$MID"; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((30 * attempt + SLURM_ARRAY_TASK_ID % 20))
  fi
done

echo "Pooled surface replay failed after three attempts: $MID" >&2
exit 1

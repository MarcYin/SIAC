#!/bin/bash
#SBATCH --job-name=medaodoof
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --array=1-152%12
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodoof_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/medaodoof_%A_%a.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MIDS="$ROOT/campaign250_lowcloud20_mids.txt"
MID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MIDS")

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
  if "$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/replay_surface_oof_uncertainty.py \
    --matchup-id "$MID"; then
    exit 0
  fi
  if (( attempt < 3 )); then
    sleep $((30 * attempt + SLURM_ARRAY_TASK_ID % 20))
  fi
done

echo "Surface OOF replay failed after three attempts: $MID" >&2
exit 1

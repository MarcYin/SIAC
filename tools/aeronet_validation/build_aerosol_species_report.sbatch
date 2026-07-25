#!/bin/bash
#SBATCH --job-name=aod_species_report
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/aod_species_report_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/aod_species_report_%j.out
set -euo pipefail
umask 0022

export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MPLCONFIGDIR=/gws/ssde/j25a/nceo_isp/public/siac_refactor/mplconfig

mkdir -p "$MPLCONFIGDIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
"$PY" /home/users/marcyin/SIAC/tools/aeronet_validation/build_aerosol_species_report.py --force

REPORT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/reports/aod-aerosol-species-20260712
test -s "$REPORT/index.html"
test -s "$REPORT/data/report.json"
test "$(find "$REPORT/assets/spatial" -maxdepth 1 -name '*.png' | wc -l)" -eq 152
test "$(find "$REPORT/assets/diagnostic" -maxdepth 1 -name '*.png' | wc -l)" -eq 152

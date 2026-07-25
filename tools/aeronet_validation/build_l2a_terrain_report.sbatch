#!/bin/bash
#SBATCH --job-name=l2aterrainweb
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainweb_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainweb_%j.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT
export PIXI_CACHE_DIR=/gws/ssde/j25a/nceo_isp/public/siac_env/pixi_cache
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export MPLCONFIGDIR=$ROOT/mplconfig
mkdir -p "$MPLCONFIGDIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_l2a_harmonization_report.py
"$PY" "$SCRIPT" --reuse-figures

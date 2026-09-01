#!/bin/bash
# Quantify how the mosaic ordering changes when the Cloud Score+ band switches
# from cs_cdf (committed) to cs. Candidate days and daily scalars are held
# fixed from the reference index, so only the per-pixel argmax input varies.
#SBATCH --job-name=cs-band-order
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/cs_band_order_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/cs_band_order_%j.out
set -euo pipefail
PROJECT=/home/users/marcyin/SIAC
export PYTHONPATH=$PROJECT:$PROJECT/python
export PROJ_LIB=$PROJECT/.pixi/envs/default/share/proj PROJ_DATA=$PROJ_LIB
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT/.pixi/envs/default/bin/python" -u "$PROJECT/tools/aeronet_validation/compare_cloudscore_band_ordering.py" "$@"

#!/bin/bash
#SBATCH --job-name=highaodreport
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/highaod_case_report_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/highaod_case_report_%j.out
set -euo pipefail

export MPLCONFIGDIR=/gws/ssde/j25a/nceo_isp/public/siac_refactor/spx_tmp/highaod_report_mpl
mkdir -p "$MPLCONFIGDIR"

PY=/home/users/marcyin/SIAC/.pixi/envs/default/bin/python
SCRIPT=/home/users/marcyin/SIAC/tools/aeronet_validation/build_highaod_case_report.py
OUT=/gws/ssde/j25a/nceo_isp/public/siac_refactor/highaod6_surface_diagnostics_20260710

"$PY" "$SCRIPT" --output-dir "$OUT"

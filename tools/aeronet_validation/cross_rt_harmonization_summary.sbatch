#!/bin/bash
#SBATCH --job-name=xrt_harm_sum
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=nceo_isp
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_harm_sum_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/xrt_harm_sum_%j.out
set -euo pipefail
umask 0022

ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
MODEL_ROOT=${CROSS_RT_OUTPUT:-$ROOT/analysis/cross_rt_harmonizer_lowcloud152_20260716}
CASES=${CROSS_RT_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
COHORT_MANIFEST=${CROSS_RT_COHORT_MANIFEST:-$ROOT/campaign250_lowcloud20_mids.json}
export PYTHONPATH=/home/users/marcyin/SIAC:/home/users/marcyin/SIAC/python:$ROOT

/home/users/marcyin/SIAC/.pixi/envs/default/bin/python \
  /home/users/marcyin/SIAC/tools/aeronet_validation/summarize_cross_rt_harmonization.py \
  --model-root "$MODEL_ROOT" \
  --cases "$CASES" \
  --cohort-manifest "$COHORT_MANIFEST"

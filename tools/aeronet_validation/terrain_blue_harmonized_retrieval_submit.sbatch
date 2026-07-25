#!/bin/bash
#SBATCH --job-name=l2aterrainblu
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=nceo_isp
#SBATCH --exclude=host1042,host1057,host1153,host1201,host1205,host1206,host1243,host1261
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=28G
#SBATCH --array=1-72%4
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainblu_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/l2aterrainblu_%A_%a.out
set -euo pipefail
umask 0022

N_CASES=36
VARIANT_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / N_CASES ))
CASE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % N_CASES + 1 ))
case "$VARIANT_INDEX" in
  0)
    HISTORY_TAG=terrain_a100_blue_cap0p015
    RUN_TAG=terrain_a100_blue_cap0p015_direct_static
    ;;
  1)
    HISTORY_TAG=terrain_a100_blue_cap0p030
    RUN_TAG=terrain_a100_blue_cap0p030_direct_static
    ;;
  *)
    echo "Unknown terrain-blue variant index: $VARIANT_INDEX" >&2
    exit 2
    ;;
esac

export HARM_CASE_INDEX=$CASE_INDEX
export HARM_HISTORY_TAG=$HISTORY_TAG
export HARM_HISTORY_ROOT=analysis/l2a_l1c_harmonizer_terrain_mediumdev_20260714/daily_histories_blue
export HARM_DEBIAS=off
export HARM_DIRECT_HISTORY=1
export HARM_TAU_MODE=off
export HARM_BACKSTOP_CAL=1
export HARM_VISIBLE_FLOOR=0.006
export HARM_DUMP_COST=1
export HARM_DEM=none
export HARM_DEM_PRE_M3=0
export HARM_RUN_TAG=$RUN_TAG

exec bash /home/users/marcyin/SIAC/tools/aeronet_validation/harmonized_l2a_retrieval_submit.sbatch

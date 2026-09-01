#!/bin/bash
# Does scouting first change the winner index? Composed twice from one cache.
#SBATCH --job-name=scout-parity
#SBATCH --account=nceo_isp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=00:40:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/scout_parity_%j.out
#SBATCH --error=/gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/scout_parity_%j.out

set -euo pipefail
umask 0022

PROJECT=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CAMPAIGN=$ROOT/current_date_toa_surface_allera_stac_20260809
PY=$PROJECT/.pixi/envs/default/bin/python
MATCHUP=AAU_ET__T37PDL_20201224T074229
OUT=$ROOT/analysis/scout_shortlist_parity_20260901

export PYTHONPATH=$PROJECT:$PROJECT/python
export TMPDIR=/work/scratch-pw5/marcyin/scout_parity/$SLURM_JOB_ID
export MPLCONFIGDIR=$TMPDIR/matplotlib
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$OUT"

"$PY" -u -m tools.aeronet_validation.compare_scouted_shortlist_index \
  --matchup-id "$MATCHUP" \
  --teacher "$CAMPAIGN/l1c_teacher20_cci_exact/$MATCHUP.npz" \
  --cache "$ROOT/cloudscore_cache_full" \
  --reference "$CAMPAIGN/historical_l1c_cloudscore_index20_ablation/$MATCHUP.npz" \
  --library-years 2018 2019 2020 2021 2022 2023 \
  --top-k 5 --workers 8 \
  --output "$OUT/$MATCHUP.json"

rm -rf "$TMPDIR"

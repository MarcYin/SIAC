#!/bin/bash
set -euo pipefail

# Submit the complete leakage-controlled low-cloud cross-RT experiment. Pair
# cases with no exact atmospheric support terminate successfully and retain the
# same uncorrected L2A prior, so a normal afterok chain is sufficient.
REPO=/home/users/marcyin/SIAC
ROOT=/gws/ssde/j25a/nceo_isp/public/siac_refactor
CASES=${CROSS_RT_CASES:-$ROOT/reports/aod-final-performance-dashboard-20260713/data/all-cases.csv}
PAIR_OUTPUT_REL=${PHYSICAL_PAIR_OUTPUT_REL:-analysis/l2a_l1c_physical_pairs_lowcloud152_20260716}
PAIR_THROTTLE=${CROSS_RT_PAIR_THROTTLE:-10}
RETRIEVAL_THROTTLE=${CROSS_RT_RETRIEVAL_THROTTLE:-24}

case_count=$(( $(wc -l < "$CASES") - 1 ))
if (( case_count < 1 )); then
  echo "No experiment cases in $CASES" >&2
  exit 2
fi
if [[ ! "$PAIR_THROTTLE" =~ ^[1-9][0-9]*$ ]]; then
  echo "CROSS_RT_PAIR_THROTTLE must be a positive integer" >&2
  exit 2
fi
if [[ ! "$RETRIEVAL_THROTTLE" =~ ^[1-9][0-9]*$ ]]; then
  echo "CROSS_RT_RETRIEVAL_THROTTLE must be a positive integer" >&2
  exit 2
fi

pair_job=$(sbatch --parsable \
  --array="1-${case_count}%${PAIR_THROTTLE}" \
  --export="ALL,PHYSICAL_PAIR_CASES=$CASES,PHYSICAL_PAIR_OUTPUT_REL=$PAIR_OUTPUT_REL" \
  "$REPO/tools/aeronet_validation/l2a_l1c_physical_pairs_submit.sbatch")
pair_job=${pair_job%%;*}
fit_job=$(sbatch --parsable --dependency="afterok:$pair_job" \
  "$REPO/tools/aeronet_validation/cross_rt_harmonization_fit.sbatch")
fit_job=${fit_job%%;*}
history_job=$(sbatch --parsable --dependency="afterok:$fit_job" \
  "$REPO/tools/aeronet_validation/cross_rt_harmonization_histories.sbatch")
history_job=${history_job%%;*}
# Fix array concurrency at submission. On the deployed Slurm 25 controller,
# changing ArrayTaskThrottle while tasks were active cancelled a grouped set of
# otherwise healthy array elements with reason JobArrayTaskLimit.
retrieval_job=$(sbatch --parsable --dependency="afterok:$history_job" \
  --array="1-${case_count}%${RETRIEVAL_THROTTLE}" \
  "$REPO/tools/aeronet_validation/cross_rt_harmonization_retrievals.sbatch")
retrieval_job=${retrieval_job%%;*}
summary_job=$(sbatch --parsable --dependency="afterany:$retrieval_job" \
  "$REPO/tools/aeronet_validation/cross_rt_harmonization_summary.sbatch")
summary_job=${summary_job%%;*}

printf 'pair=%s\nfit=%s\nhistories=%s\nretrievals=%s\nsummary=%s\n' \
  "$pair_job" "$fit_job" "$history_job" "$retrieval_job" "$summary_job"

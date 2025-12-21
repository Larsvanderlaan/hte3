#!/usr/bin/env bash
set -euo pipefail

# Use the same R 4.4 Rscript that RStudio uses
RSCRIPT="/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript"
export R_LIBS_USER="/Users/larsvanderlaan/Library/R/arm64/4.4/library"

echo "Using: $RSCRIPT"
"$RSCRIPT" --version

RFILE="$HOME/repos/hte3/paper_EPlearner_experiments/SimulationRcode/EPlearnerCATE_wager.R"
LOGDIR="$HOME/repos/hte3/paper_EPlearner_experiments/results/logs"
mkdir -p "$LOGDIR"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for n in 500 1000 2000; do
for base_learner in rf; do
echo "Launching 4 parallel runs for n=$n, base_learner=$base_learner ..."
pids=(); tags=()
for sim_type in A C; do
tag="n${n}_learner${base_learner}_sim${sim_type}"
log="$LOGDIR/run_${tag}.out"
echo "  -> $tag"
"$RSCRIPT" "$RFILE" --args n="$n" d=5 sigma=1 base_learner="$base_learner" sim_type="$sim_type" \
> "$log" 2>&1 & pids+=($!); tags+=("$tag")
done
ok=1
for i in "${!pids[@]}"; do
if ! wait "${pids[$i]}"; then
echo "ERROR: ${tags[$i]} failed"; tail -n 50 "$LOGDIR/run_${tags[$i]}.out"; ok=0
fi
done
[[ $ok -eq 1 ]] || exit 1
echo "Completed n=$n, base_learner=$base_learner."
done
done
echo "All jobs done."

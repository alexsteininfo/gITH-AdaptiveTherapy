#!/usr/bin/env bash
# Recompute the statistics tables (snapshots.csv, samples.csv) of finished blocks from their raw
# tumour files, without simulating again — e.g. after a fix in a statistic. The summary tables
# (sfs_mean, tradeoff, trial, powerlaw, …) are not redone.
# Usage: bash code/01_datageneration/run_statistics.sh [block …]   (default: all blocks; WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
blocks=("$@")
[ ${#blocks[@]} -eq 0 ] && blocks=(01_dynamics 02_tradeoff 03_sensitivity 04_trial 05_metronomic 06_cost 07_validation)
args=()
for b in "${blocks[@]}"; do args+=("$(pblock "$b")"); done
"$JULIA" -p "$WORKERS" --project=. code/02_processing/02_statistics.jl "${args[@]}"
echo "# done $(date)"

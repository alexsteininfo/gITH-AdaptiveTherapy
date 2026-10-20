#!/usr/bin/env bash
# Part 3 — virtual trial (block 04_trial): simulations, tables, correlations and policies.
# Usage: bash code/01_datageneration/run_part3.sh   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
simulate_and_process 04_trial "$@"
summarise 06_trial.jl 04_trial
echo "# done $(date)"

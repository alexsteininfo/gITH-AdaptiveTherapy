#!/usr/bin/env bash
# Part 2 — trade-off and sensitivity (blocks 02_tradeoff, 03_sensitivity): the largest part.
# Usage: bash code/01_datageneration/run_part2.sh [cohort …]   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
simulate_and_process 02_tradeoff "$@"
simulate_and_process 03_sensitivity "$@"
summarise 05_tradeoff.jl 02_tradeoff 03_sensitivity
echo "# done $(date)"

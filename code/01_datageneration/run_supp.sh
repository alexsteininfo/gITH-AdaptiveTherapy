#!/usr/bin/env bash
# Supplement — metronomic (05), cost of resistance (06), validation (07).
# Usage: bash code/01_datageneration/run_supp.sh [cohort …]   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
simulate_and_process 05_metronomic "$@"
simulate_and_process 06_cost "$@"
simulate_and_process 07_validation "$@"
summarise 04_examples.jl 01_dynamics 05_metronomic
summarise 05_tradeoff.jl 05_metronomic 06_cost
summarise 08_powerlaw.jl 01_dynamics
summarise 07_validation.jl 07_validation
echo "# done $(date)"

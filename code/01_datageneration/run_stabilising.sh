#!/usr/bin/env bash
# Stabilising AT at fine time resolution (block 08_stabilising; Figure 4): simulations, tables,
# mean SFS and the neutral theory of the stabilising schedule.
# Usage: bash code/01_datageneration/run_stabilising.sh [cohort …]   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
simulate_and_process 08_stabilising "$@"
summarise 03_sfs.jl 08_stabilising
summarise 09_theory.jl 08_stabilising
echo "# done $(date)"

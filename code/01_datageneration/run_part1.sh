#!/usr/bin/env bash
# Part 1 — dynamics (block 01_dynamics): simulations, tables, mean SFS, Kaplan–Meier per arm, theory.
# Usage: bash code/01_datageneration/run_part1.sh [cohort …]   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
simulate_and_process 01_dynamics "$@"
summarise 03_sfs.jl 01_dynamics
summarise 05_tradeoff.jl 01_dynamics
summarise 09_theory.jl 01_dynamics
echo "# done $(date)"

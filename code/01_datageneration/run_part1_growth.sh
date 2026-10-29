#!/usr/bin/env bash
# Part 1 addendum — untreated records at 0.5 N_det for the finished block 01_dynamics (regrowth only,
# arms not re-run), then its statistics and mean SFS again (the mean SFS of the holidays without
# tumours that progress within the stored cycles).
# Usage: bash code/01_datageneration/run_part1_growth.sh [cohort …]   (WORKERS, JULIA, TEST=1)
source "$(dirname "$0")/common.sh"
"$JULIA" -p "$WORKERS" --project=. code/01_datageneration/01_dynamics_growth.jl "$@"
"$JULIA" -p "$WORKERS" --project=. code/02_processing/02_statistics.jl "$(pblock 01_dynamics)"
summarise 03_sfs.jl 01_dynamics
summarise 09_theory.jl 01_dynamics
echo "# done $(date)"

# Part 1: dynamics — neutral, full, partial and weak resistance under single round, MTD, holidays and modulation.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/01_dynamics/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/01_dynamics.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("01_dynamics"; parse_args(ARGS)...)

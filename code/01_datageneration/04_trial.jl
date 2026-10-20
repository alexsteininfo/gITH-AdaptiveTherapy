# Part 3: virtual trial — 500 Latin-hypercube patients (nu, alpha, N_det/K) under MTD, holidays and modulation.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/04_trial/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/04_trial.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("04_trial"; parse_args(ARGS)...)

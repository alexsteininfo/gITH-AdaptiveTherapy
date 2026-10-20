# Part 2: sensitivity — N_det/K = 0.5, 0.7, 0.9 and nu = 1e-5, 1e-3 under the containment scans.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/03_sensitivity/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/03_sensitivity.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("03_sensitivity"; parse_args(ARGS)...)

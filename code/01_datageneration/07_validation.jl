# Supplement: growth only, for checks against theory (SFS, Luria–Delbrück).
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/07_validation/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/07_validation.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("07_validation"; parse_args(ARGS)...)

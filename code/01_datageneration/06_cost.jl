# Supplement: random resistance levels (Beta(1,4)) with and without a cost of resistance.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/06_cost/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/06_cost.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("06_cost"; parse_args(ARGS)...)

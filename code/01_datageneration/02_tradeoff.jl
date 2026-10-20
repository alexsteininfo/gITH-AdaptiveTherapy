# Part 2: trade-off — full and partial resistance under MTD and the containment scans (holidays 0.5–0.9, modulation 0.5–1.0).
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/02_tradeoff/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/02_tradeoff.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("02_tradeoff"; parse_args(ARGS)...)

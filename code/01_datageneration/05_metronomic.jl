# Supplement: metronomic schedule on the neutral and full-resistance tumours of block 01.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/05_metronomic/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/05_metronomic.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("05_metronomic"; parse_args(ARGS)...)

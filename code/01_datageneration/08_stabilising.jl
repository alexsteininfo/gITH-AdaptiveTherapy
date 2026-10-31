# Stabilising AT at fine time resolution (block 08_stabilising; Figure 4): the first 200 tumours of
# the Part 1 cohorts under `mod` (neutral also `mod_0.5`, `mod_1.0`), distributions at the times
# HEAVY_FINE after detection and the untreated tumour at the stable size 0.8 N_det.
# Specification: code/01_datageneration/README.md (table "All planned simulations").
# Output: data/raw/08_stabilising/<cohort>/tNNNN.jld2 (one file per tumour), meta.toml, examples/.
# Usage : julia -p 12 --project=. code/01_datageneration/08_stabilising.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

run_block("08_stabilising"; parse_args(ARGS)...)

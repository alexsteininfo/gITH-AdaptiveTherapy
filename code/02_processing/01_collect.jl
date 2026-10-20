# Collect per-tumour outcomes and population counts of one or more blocks into tables.
# Input : data/raw/<block>/<cohort>/tNNNN.jld2      Output: data/processed/<block>/{outcomes.csv, counts.csv.gz}
# Usage : julia -p 8 --project=. code/02_processing/01_collect.jl <block> [<block> …]
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

isempty(ARGS) && error("usage: 01_collect.jl <block> [<block> …]")
foreach(collect_block, ARGS)

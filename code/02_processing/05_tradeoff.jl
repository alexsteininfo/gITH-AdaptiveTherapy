# Trade-off summaries (progression, Kaplan–Meier, diversity at the end, paired gain vs MTD).
# Input : data/processed/<block>/ tables and data/raw/<block>/ files     Output: data/processed/<block>/tradeoff.csv, km.csv
# Usage : julia -p 8 --project=. code/02_processing/05_tradeoff.jl <block> [<block> …]
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

isempty(ARGS) && error("usage: 05_tradeoff.jl <block> [<block> …]")
foreach(make_tradeoff, ARGS)

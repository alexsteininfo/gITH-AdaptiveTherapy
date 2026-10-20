# Validation tables (mean SFS vs theory; resistant clones at detection).
# Input : data/processed/<block>/ tables and data/raw/<block>/ files     Output: data/processed/<block>/validation_sfs.csv, validation_ld.csv
# Usage : julia -p 8 --project=. code/02_processing/07_validation.jl <block>
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

isempty(ARGS) && error("usage: 07_validation.jl <block>")
foreach(make_validation, ARGS)

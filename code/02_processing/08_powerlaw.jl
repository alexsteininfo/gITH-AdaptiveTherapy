# Power-law fits of the SFS tail over k_min, and the k_min chosen by the KS distance (Supp 6).
# Output: data/processed/<block>/powerlaw.csv
# Usage : julia -p 8 --project=. code/02_processing/08_powerlaw.jl <block>
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

length(ARGS) == 1 || error("usage: 08_powerlaw.jl <dynamics block>")
make_powerlaw(ARGS[1])

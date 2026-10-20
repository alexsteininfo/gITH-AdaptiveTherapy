# Statistics of all stored distributions (entire tree and subsamples) of one or more blocks.
# Input : data/raw/<block>/<cohort>/tNNNN.jld2
# Output: data/processed/<block>/{snapshots<suffix>.csv, samples<suffix>.csv}
# Usage : julia -p 8 --project=. code/02_processing/02_statistics.jl <block> [<block> …] [--kmin-frac 0.01] [--suffix _x]
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

function main(args)
    kmin_frac = 0.01; suffix = ""; blocks = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--kmin-frac"; kmin_frac = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--suffix"; suffix = args[i+1]; i += 2
        else push!(blocks, args[i]); i += 1
        end
    end
    isempty(blocks) && error("usage: 02_statistics.jl <block> [<block> …] [--kmin-frac x] [--suffix s]")
    foreach(b -> statistics_block(b; kmin_frac, suffix), blocks)
end
main(ARGS)

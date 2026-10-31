# Expected neutral SFS and Hill numbers under the holidays (deterministic logistic cycle; master
# formula and harmonic-mean law, .docs/theory/math/ANALYSIS.md §3.1, §4.6), the single round (master
# formula) and, for a stabilising block (08_stabilising), the stabilising schedule (master formula on
# growth, reduction and hold; .docs/theory/STABILISING.md). No raw data needed.
# Output: data/processed/<block>/theory_hill.csv, theory_sfs.csv, theory_single_hill.csv, theory_single_sfs.csv
#         (dynamics block) or theory_stab_hill.csv, theory_stab_sfs.csv (stabilising block)
# Usage : julia --project=. code/02_processing/09_theory.jl <dynamics or stabilising block>
include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

length(ARGS) == 1 || error("usage: 09_theory.jl <block>")
if occursin("08_stabilising", ARGS[1])
    make_theory_stabilising(ARGS[1])
else
    make_theory(ARGS[1])
    make_theory_single(ARGS[1])
end

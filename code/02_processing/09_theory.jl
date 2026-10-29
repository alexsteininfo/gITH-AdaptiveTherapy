# Expected neutral SFS and Hill numbers under the holidays (deterministic logistic cycle; master
# formula and harmonic-mean law, .docs/theory/math/ANALYSIS.md §3.1, §4.6) and the single round
# (master formula). No raw data needed.
# Output: data/processed/<block>/theory_hill.csv, theory_sfs.csv, theory_single_hill.csv, theory_single_sfs.csv
# Usage : julia --project=. code/02_processing/09_theory.jl <dynamics block>
include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

length(ARGS) == 1 || error("usage: 09_theory.jl <dynamics block>")
make_theory(ARGS[1])
make_theory_single(ARGS[1])

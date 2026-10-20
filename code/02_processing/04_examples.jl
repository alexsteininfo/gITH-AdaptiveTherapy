# Example tumour of cohort `full` (closest to the median resistant fraction at detection, among
# the tumours with stored clone histories): trajectories and resistant-clone sizes for Figure 1.
# Output: data/processed/<block01>/{example_traj.csv, example_clones.csv}
# Usage : julia --project=. code/02_processing/04_examples.jl <block01> <block05>
include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

length(ARGS) == 2 || error("usage: 04_examples.jl <dynamics block> <metronomic block>")
make_examples(ARGS[1], ARGS[2])

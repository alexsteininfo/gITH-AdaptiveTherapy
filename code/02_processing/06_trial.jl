# Virtual trial: per-patient table, correlations of diversity with the treatment course and AT
# benefit (within one parameter set, across patients, partial), treatment-selection policies.
# Output: data/processed/<block>/{trial.csv, trial_cor.csv, trial_policy.csv, trial_km.csv}
# Usage : julia --project=. code/02_processing/06_trial.jl <trial block> [<dynamics block>]
include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

isempty(ARGS) && error("usage: 06_trial.jl <trial block> [<dynamics block>]")
make_trial(ARGS[1]; within_block = length(ARGS) > 1 ? ARGS[2] : replace(ARGS[1], "04_trial" => "01_dynamics"))

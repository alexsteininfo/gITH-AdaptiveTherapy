# Part 1 addendum: untreated records at every size of `growth_heavy` (0.1 and 0.5 N_det) for tumours
# of block 01_dynamics simulated with fewer sizes. Each tumour is regrown from its seed (the same
# tumour) and its growth records are replaced; the arms are not re-run. The regrown detection record
# and the growth records stored before must equal the stored ones, otherwise the file is left unchanged
# and the run fails at the end, listing those tumours. Resumable.
# Output: data/raw/01_dynamics/<cohort>/tNNNN.jld2 (growth records replaced in place, atomically).
# Usage : julia -p 12 --project=. code/01_datageneration/01_dynamics_growth.jl [cohort …]    (TEST=1: smoke test)
using Distributed
@everywhere include(joinpath(@__DIR__, "..", "X_helpers", "helpers.jl"))

regrow_block("01_dynamics"; parse_args(ARGS)...)

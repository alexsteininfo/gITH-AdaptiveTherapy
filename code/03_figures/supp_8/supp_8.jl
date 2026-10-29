# Supplementary Figure 8 — as Figure 3d–k at the lower threshold 0.5 N_det (end of the treatment phase
# of every cycle), relative to the same tumour untreated at 0.5 N_det (growth record). Only tumours
# that complete the first six cycles without progressing; theory as in Figure 3 (dotted).
#   a–d   passenger Hill numbers q = 0, 1, 2, ∞ (SFS; q = 0: number of mutations S)
#   e–h   resistant-clone Hill numbers q = 0, 1, 2, ∞ (clone size distribution), weak, partial, full
#   legend
# Data: 01_dynamics/snapshots.csv, theory_hill.csv
# Run : julia --project=. code/03_figures/supp_8/supp_8.jl
include(joinpath(@__DIR__, "..", "style.jl"))
include(joinpath(@__DIR__, "..", "fig_3", "cycles.jl"))

const NAME = "supp_8"

hol = cycles()
hill_panels(hol, "lower", NAME, ["a", "b", "c", "d", "e", "f", "g", "h"])
hill_legend(NAME)

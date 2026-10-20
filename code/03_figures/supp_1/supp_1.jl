# Supplementary Figure 1 — validation of the simulator against theory (no competition, no drug).
#   a  mean SFS at N = 10⁵ (pure birth, d = 0) and the expectation 2μN Σₗ (d/b)ˡ/((k+l)(k+l+1))
#   b  same for the birth–death process (b = 1, d = 0.5)
#   c  number of resistant clones at detection (ν = 10⁻³); dotted: expectation
#      ∫ 2νb N(s) P(clone born at s alive at T) ds, N(s) = e^{rs}, P = r/(b − d e^{−r(T−s)}), r = b − d
#   d  number of resistant cells at detection (log scale); dotted: expected mean 2ν(b/r) N ln N
# Data: 07_validation/validation_sfs.csv, validation_ld.csv, snapshots.csv
# Run : julia --project=. code/03_figures/supp_1/supp_1.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_1"
const B, D, NU = 1.0, 0.5, 1e-3      # cohort `ld` (code/X_helpers/params.jl)

sfs = readproc("07_validation", "validation_sfs.csv")
ld = readproc("07_validation", "validation_ld.csv")

for (p, cohort) in (("a", "sfs_pb"), ("b", "sfs_bd"))
    d = sort(sfs[sfs.cohort .== cohort, :], :k)
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
    lines!(ax, d.k, d.S_mean; color = COL_MEAN, linewidth = 2, label = "Simulations")
    lines!(ax, d.k, d.S_pred; color = COL_PRED, linestyle = LS_PRED, label = "Theory")
    axislegend(ax; position = :rt)
    savepanel(fig, NAME, p)
end

snap = readproc("07_validation", "snapshots.csv")
const NDET_LD = first(snap.N[snap.cohort .== "ld"])
r = B - D; T = log(NDET_LD) / r
s = range(0, T; length = 4001)
alive = [2NU * B * exp(r * x) * r / (B - D * exp(-r * (T - x))) for x in s]
expected_clones = sum(alive) * step(s)
expected_cells = 2NU * (B / r) * NDET_LD * log(NDET_LD)

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:resclones], ylabel = "Tumours")
hist!(ax, ld.n_res_clones; bins = 20, color = (COL_MEAN, 0.6))
vlines!(ax, [expected_clones]; color = COL_PRED, linestyle = LS_PRED)
savepanel(fig, NAME, "c")

fig = panel(); ax = Axis(fig[1, 1]; xlabel = "Resistant cells", ylabel = "Tumours", LOGX...)
x = filter(>(0), ld.n_res_cells)
edges = 10 .^ range(floor(log10(minimum(x))), ceil(log10(maximum(x))); length = 16)
hist!(ax, x; bins = edges, color = (COL_MEAN, 0.6))
vlines!(ax, [expected_cells]; color = COL_PRED, linestyle = LS_PRED)
vlines!(ax, [sum(ld.n_res_cells) / nrow(ld)]; color = COL_MEAN)
savepanel(fig, NAME, "d")

# Supplementary Figure 9 — SFS under cyclic adaptive therapy for every resistance cohort, at the lower
# (0.5 N_det) and upper (N_det) threshold of every cycle, paler with every cycle (as Figure 3a, b).
# Only tumours that complete the first six cycles without progressing (03_sfs.jl).
#   a–d   mean SFS, cohorts none, weak, partial, full
#   e–h   compensated SFS k(k+1)·S_k (without the clonal bin)
#   legend
# Data: 01_dynamics/sfs_mean.csv
# Run : julia --project=. code/03_figures/supp_9/supp_9.jl
include(joinpath(@__DIR__, "..", "style.jl"))
include(joinpath(@__DIR__, "..", "fig_3", "cycles.jl"))

const NAME = "supp_9"

for (row, compensated) in ((["a", "b", "c", "d"], false), (["e", "f", "g", "h"], true))
    curves = Dict((c, th, cyc) => sfs_curve(c, th, cyc; compensated)
                  for c in COHORTS for th in THRESHOLDS for cyc in 0:NCYCLES if !(th == "lower" && cyc == 0))
    ys = vcat([last(v) for v in values(curves)]...); xs = vcat([first(v) for v in values(curves)]...)
    for (p, c) in zip(row, COHORTS)
        fig = panel()
        ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = compensated ? "𝑘(𝑘+1)·𝑆ₖ" : LABELS[:Sk], LOGX..., LOGY...)
        for th in THRESHOLDS, cyc in NCYCLES:-1:0                  # pale (late) first, dark on top
            haskey(curves, (c, th, cyc)) || continue
            lines!(ax, curves[(c, th, cyc)]...; color = cyclecolour(th, cyc), linewidth = 0.8)
        end
        limits!(ax, minimum(xs) / 1.5, maximum(xs) * 1.5, minimum(ys) / 1.3, maximum(ys) * 1.3)
        savepanel(fig, NAME, p)
    end
end

fig = panel(; w = 180, h = 8)
Legend(fig[1, 1], [LineElement(color = cyclecolour(th, c), linewidth = 1.5) for th in THRESHOLDS for c in (th == "lower" ? 1 : 0, NCYCLES)],
       ["$(THRESHOLD_NAME[th]), cycle $c" for th in THRESHOLDS for c in (th == "lower" ? 1 : 0, NCYCLES)];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

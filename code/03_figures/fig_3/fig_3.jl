# Figure 3 — genetic diversity increases with every cycle of cyclic adaptive therapy (arm `hol`: full
# dose until 0.5 N_det, off until regrowth to N_det), here at the upper threshold N_det (detection =
# cycle 0); the lower threshold is Supp 8, the SFS of the resistance cohorts Supp 9. Only tumours that
# complete the first six cycles without progressing. Theory: deterministic logistic cycle,
# growth-phase mutations by the master formula, mutations born during cycling by the harmonic-mean
# law (02_processing/09_theory.jl).
#   a     mean SFS, no-resistance cohort, at detection and the upper threshold of every cycle (paler
#         with every cycle)
#   b     compensated SFS k(k+1)·S_k of a (without the clonal bin), theory dashed
#   c     peak position of b per cycle, simulation and theory; dotted: b_eff·τ, the cutoff of the
#         new mutations' 1/k law
#   d–g   passenger Hill numbers q = 0, 1, 2, ∞ (SFS; q = 0: number of mutations S), treated / the
#         same tumour at detection, median and interquartile range per cohort, theory dotted
#   h–k   resistant-clone Hill numbers q = 0, 1, 2, ∞ (clone size distribution), weak, partial, full
#   legend (d–k), legend_sfs (a–c)
# Data: 01_dynamics/snapshots.csv, sfs_mean.csv, theory_hill.csv, theory_sfs.csv
# Run : julia --project=. code/03_figures/fig_3/fig_3.jl
include(joinpath(@__DIR__, "..", "style.jl"))
include(joinpath(@__DIR__, "cycles.jl"))

const NAME = "fig_3"
const SFS_COHORT = "neutral"

# ── a, b: SFS per cycle ───────────────────────────────────────────────────────
for (p, compensated) in (("a", false), ("b", true))
    fig = panel()
    ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = compensated ? "𝑘(𝑘+1)·𝑆ₖ" : LABELS[:Sk], LOGX..., LOGY...)
    for cyc in NCYCLES:-1:0                                     # pale (late) first, dark on top
        lines!(ax, sfs_curve(SFS_COHORT, "upper", cyc; compensated)...; color = cyclecolour("upper", cyc), linewidth = 1.0)
    end
    compensated && for cyc in 0:NCYCLES                         # theory above all simulated curves
        lines!(ax, theory_curve("upper", cyc; compensated)...; color = :grey30, linestyle = :dash, linewidth = 0.6)
    end
    savepanel(fig, NAME, p)
end

# ── c: peak of the compensated SFS ────────────────────────────────────────────
"Position of the maximum of `y` over `k` below 10⁴ (vertex of a parabola in log–log through the
largest bin and its neighbours); NaN without data."
function peakpos(k, y)
    keep = k .< 1e4; x = log.(k[keep]); z = log.(y[keep])
    isempty(z) && return NaN                                    # cycle not reached (test data)
    i = argmax(z)
    (i == 1 || i == length(z)) && return exp(x[i])
    (x1, x2, x3), (y1, y2, y3) = x[i-1:i+1], z[i-1:i+1]
    A = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1)
    B = (y2 - y1) / (x2 - x1) - A * (x1 + x2)
    return exp(-B / 2A)
end

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:cycle], ylabel = "Peak of 𝑘(𝑘+1)·𝑆ₖ", xticks = 0:NCYCLES, LOGY...)
th = theory_hill[(theory_hill.model .== THEORY_MODEL) .& (theory_hill.threshold .== "upper") .& (theory_hill.kind .== "treated"), :]
lines!(ax, th.cycle, th.kc; color = COL_PRED, linestyle = LS_PRED, linewidth = 1, label = rich("𝑏", subscript("eff"), "𝜏"))
lines!(ax, 0:NCYCLES, [peakpos(theory_curve("upper", c; compensated = true)...) for c in 0:NCYCLES];
       color = :grey30, linestyle = :dash, linewidth = 1, label = "Theory")
scatter!(ax, 0:NCYCLES, [peakpos(sfs_curve(SFS_COHORT, "upper", c; compensated = true)...) for c in 0:NCYCLES];
         color = [cyclecolour("upper", c) for c in 0:NCYCLES], strokecolor = COL_THRESHOLD["upper"], strokewidth = 0.6,
         markersize = 6, label = "Simulation")
xlims!(ax, -0.5, NCYCLES + 0.5)
axislegend(ax; position = :rb)
savepanel(fig, NAME, "c")

# ── d–k: Hill numbers per cycle ───────────────────────────────────────────────
hol = cycles()
hill_panels(hol, "upper", NAME, ["d", "e", "f", "g", "h", "i", "j", "k"])
hill_legend(NAME)

# ── legend of a–c ─────────────────────────────────────────────────────────────
fig = panel()                                                   # fills the free panel slot of the first row
Legend(fig[1, 1], [[LineElement(color = cyclecolour("upper", c), linewidth = 1.5) for c in 0:NCYCLES];
                   [LineElement(color = :grey30, linestyle = :dash, linewidth = 0.8)]],
       [["Cycle $c" for c in 0:NCYCLES]; "Theory"]; nbanks = 2, labelsize = 7, patchsize = (12, 5),
       rowgap = 4, colgap = 8)
savepanel(fig, NAME, "legend_sfs")

# Supplementary Figure 7 — the SFS after regrowth depends on the depth of treatment. As Figure 2D,
# but for a 50% decrease (a, d), and compensated versions of Figure 2A and 2D (b, c). 50% decrease:
# full dose until 0.5 N_det, then off until regrowth to N_det (first cycle of arm `hol`, end of block 2); compared with the untreated tumour at detection. Tumours that never
# shrink to 0.5 N_det have no second block and are thereby excluded (1 of 500 with full resistance),
# as are tumours that progress within the first six cycles (mean SFS of the holidays, 03_sfs.jl).
#   a  mean SFS after regrowth (N): untreated (no-resistance cohort) and treated per cohort
#   b, c  compensated SFS k(k+1)·S_k of Figure 2A (after treatment, 0.1N) and 2D (after regrowth, N);
#         flat for the neutral expectation 2μN/(k(k+1)) of a growing population; without the clonal bin
#   d  same for the 50% decrease (panel a)
#   legend
# Data: 01_dynamics/sfs_mean.csv
# Run : julia --project=. code/03_figures/supp_7/supp_7.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_7"
const COHORTS = ["neutral", "weak", "partial", "full"]
const UNTREATED = "neutral"     # cohort of the untreated SFS (identical in law across cohorts)

sfs = readproc("01_dynamics", "sfs_mean.csv")

# ── a: mean SFS after regrowth from 0.5 N ─────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
u = sfs[(sfs.cohort .== UNTREATED) .& (sfs.tag .== "detection") .& (sfs.S_mean .> 0), :]
lines!(ax, u.k, u.S_mean; LINE_UNTREATED...)
for c in COHORTS
    d = sfs[(sfs.cohort .== c) .& (sfs.arm .== "hol") .& (sfs.tag .== "blockend") .& (sfs.idx .== 2) .&
            (sfs.S_mean .> 0) .& (sfs.k .<= maximum(u.k)), :]
    lines!(ax, d.k, d.S_mean; color = COL_COHORT[c], LINE_COHORT[c]...)
end
savepanel(fig, NAME, "a")

# ── b–d: compensated SFS ────────────────────────────────────────────────────
"Compensated mean SFS k(k+1)·S_k, without the clonal bin: untreated record `untreated = (tag, idx)`
(style, k, y) first, then the treated records `(arm, idx)` per cohort."
function compensated(arm, idx, untreated)
    u = sfs[(sfs.cohort .== UNTREATED) .& (sfs.tag .== untreated[1]) .& (sfs.idx .== untreated[2]) .& (sfs.S_mean .> 0), :]
    kmax = 0.9maximum(u.k)
    curves = Any[(LINE_UNTREATED, u.k[u.k .< kmax], (u.k .* (u.k .+ 1) .* u.S_mean)[u.k .< kmax])]
    for c in COHORTS
        d = sfs[(sfs.cohort .== c) .& (sfs.arm .== arm) .& (sfs.tag .== "blockend") .& (sfs.idx .== idx) .&
                (sfs.S_mean .> 0) .& (sfs.k .< kmax), :]
        push!(curves, ((; color = COL_COHORT[c], LINE_COHORT[c]...), d.k, d.k .* (d.k .+ 1) .* d.S_mean))
    end
    return curves
end

"Ticks of a log axis spanning less than a decade or two: 1, 2, 3, 4, 5, 7 per decade."
fineticks(lo, hi) = [m * 10.0^e for e in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 3, 4, 5, 7)
                     if lo <= m * 10.0^e <= hi]

function compensated_panel(curves; limits = extrema(vcat([c[3] for c in curves]...)))
    fig = panel()
    ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = "𝑘(𝑘+1)·𝑆ₖ", LOGX..., yscale = log10, yticks = fineticks,
              ytickformat = logformat)
    for (style, k, y) in curves
        lines!(ax, k, y; style...)
    end
    ylims!(ax, limits[1] / 1.05, limits[2] * 1.05)
    return fig
end

savepanel(compensated_panel(compensated("single", 1, ("growth", 10))), NAME, "b")
regrown = (compensated("single", 2, ("detection", 0)), compensated("hol", 2, ("detection", 0)))
shared = extrema(vcat([c[3] for cs in regrown for c in cs]...))       # c and d on the same axis
savepanel(compensated_panel(regrown[1]; limits = shared), NAME, "c")
savepanel(compensated_panel(regrown[2]; limits = shared), NAME, "d")

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 160, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)]; [LineElement(; color = COL_COHORT[c], LINE_COHORT[c]...) for c in COHORTS]],
       ["Untreated"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

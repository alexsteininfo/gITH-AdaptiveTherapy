# Supplementary Figure 7 — the SFS after regrowth depends on the depth of treatment. As Figure 2E, F,
# but for a 50% decrease: full dose until 0.5 N_det, then off until regrowth to N_det (first cycle of
# arm `hol`, end of block 2); compared with the untreated tumour at detection. Tumours that never
# shrink to 0.5 N_det have no second block and are thereby excluded (1 of 500 with full resistance),
# as are tumours that progress within the first six cycles (mean SFS of the holidays, 03_sfs.jl).
# Theory: deterministic logistic cycle, master formula (theory_sfs.csv, cycle 1, upper threshold).
#   a  mean SFS after regrowth (N): untreated (no-resistance cohort) and treated per cohort, theory dashed
#   b  compensated SFS k(k+1)·S_k (without the clonal bin), theory dashed; same y-axis as Figure 2F
#   legend
# Data: 01_dynamics/sfs_mean.csv, theory_sfs.csv, theory_single_sfs.csv (y-axis of Figure 2F)
# Run : julia --project=. code/03_figures/supp_7/supp_7.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_7"
const COHORTS = ["neutral", "weak", "partial", "full"]
const UNTREATED = "neutral"     # cohort of the untreated SFS (identical in law across cohorts)
const THEORY = (color = :grey30, linestyle = :dash, linewidth = 0.6)

sfs = readproc("01_dynamics", "sfs_mean.csv")
theory_hol = readproc("01_dynamics", "theory_sfs.csv")
theory_single = readproc("01_dynamics", "theory_single_sfs.csv")

"Curves `(style, k, y)` after regrowth to N_det from `arm` (`hol`: 0.5 N_det, `single`: 0.1 N_det):
untreated, treated per cohort, then theory (untreated, treated); compensated without the clonal bin."
function sfs_curves(arm; compensated)
    u = sfs[(sfs.cohort .== UNTREATED) .& (sfs.tag .== "detection") .& (sfs.S_mean .> 0), :]
    kmax = maximum(u.k)
    keep(k) = compensated ? k .< 0.9kmax : k .<= kmax
    y(d) = compensated ? d.k .* (d.k .+ 1) .* d.S_mean : d.S_mean
    curves = Any[(LINE_UNTREATED, u.k[keep(u.k)], y(u)[keep(u.k)])]
    for c in COHORTS
        d = sfs[(sfs.cohort .== c) .& (sfs.arm .== arm) .& (sfs.tag .== "blockend") .& (sfs.idx .== 2) .& (sfs.S_mean .> 0), :]
        push!(curves, ((; color = COL_COHORT[c], LINE_COHORT[c]...), d.k[keep(d.k)], y(d)[keep(d.k)]))
    end
    th = arm == "hol" ?
         [theory_hol[(theory_hol.model .== "master") .& (theory_hol.threshold .== "upper") .& (theory_hol.kind .== kind) .&
                     (theory_hol.cycle .== cyc), :] for (kind, cyc) in (("untreated", 0), ("treated", 1))] :
         [theory_single[(theory_single.threshold .== "upper") .& (theory_single.kind .== kind), :] for kind in ("untreated", "treated")]
    for d in th
        d = d[d.S_mean .> 0, :]
        push!(curves, (THEORY, d.k[keep(d.k)], y(d)[keep(d.k)]))
    end
    return curves
end

# ── a, b: SFS and compensated SFS after regrowth from 0.5 N ───────────────────
for (p, compensated) in (("a", false), ("b", true))
    fig = panel()
    ax = compensated ?
         Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = "𝑘(𝑘+1)·𝑆ₖ", LOGX..., yscale = log10, yticks = fineticks,
              ytickformat = logformat) :
         Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
    for (style, k, y) in sfs_curves("hol"; compensated)
        lines!(ax, k, y; style...)
    end
    if compensated                                      # the y-axis of Figure 2F
        v = [y for (_, _, ys) in sfs_curves("single"; compensated) for y in ys]
        ylims!(ax, minimum(v) / 1.05, maximum(v) * 1.05)
    end
    savepanel(fig, NAME, p)
end

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 180, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)]; [LineElement(; color = COL_COHORT[c], LINE_COHORT[c]...) for c in COHORTS];
                   [LineElement(; THEORY..., linewidth = 0.8)]],
       ["Untreated"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]; "Theory"];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

# Figure 2 — one treatment round leaves high diversity behind. Arm `single`: full dose until
# 0.1 N_det ("after treatment"), then off until regrowth to N_det ("relapse"); compared with the
# same tumour untreated at these sizes (growth record at 0.1 N_det, detection). Only tumours that
# respond (reach 0.1 N_det); the others had high resistance at detection and never left the
# treatment phase. Cohorts: no resistance (ν = 0), weak, partial, full.
# Top row after treatment (0.1N), bottom row after regrowth (N):
#   a, d  mean SFS: untreated (no-resistance cohort) and treated per cohort
#   b, e  passenger Hill numbers ᵠD (q = 0: number of mutations S), treated / untreated, same tumour
#   c, f  resistant-clone Hill numbers (clone size distribution), treated / untreated
#         (weak, partial, full; none in the no-resistance cohort)
#   legend
# Data: 01_dynamics/sfs_mean.csv, snapshots.csv
# Run : julia --project=. code/03_figures/fig_2/fig_2.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_2"
const COHORTS = ["neutral", "weak", "partial", "full"]
const RES_COHORTS = ["weak", "partial", "full"]
const UNTREATED = "neutral"     # cohort of the untreated SFS (identical in law across cohorts)
# SFS lines, drawn in this order: thick to thin, so that overlapping curves stay visible
const LINE_UNTREATED = (color = (COL_ARM["untreated"], 0.6), linewidth = 3.0)
const LINE = Dict("neutral" => (linestyle = :solid, linewidth = 1.8), "weak" => (linestyle = :dash, linewidth = 1.4),
                  "partial" => (linestyle = :dot, linewidth = 1.2), "full" => (linestyle = :dashdot, linewidth = 1.0))
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]
# matched sizes: single-round block (1: end of treatment, 2: regrown) and untreated record
const SIZES = ((block = 1, untreated = "growth"), (block = 2, untreated = "detection"))

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")         # single round: responders only

"Single-round record of `block` (1: end of treatment, 2: regrown) of every responding tumour."
function treated(block)
    t = snap[(snap.arm .== "single") .& (snap.tag .== "blockend"), :]
    det = snap[snap.tag .== "detection", [:cohort, :tumour, :N]]
    b1 = innerjoin(t[t.idx .== 1, [:cohort, :tumour, :N]], det; on = [:cohort, :tumour], renamecols = "" => "_det")
    resp = b1[b1.N .<= b1.N_det ./ 10, [:cohort, :tumour]]
    return semijoin(t[t.idx .== block, :], resp; on = [:cohort, :tumour])
end

# ── a, b: mean SFS ────────────────────────────────────────────────────────────
for (p, s) in zip(("a", "d"), SIZES)
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
    u = sfs[(sfs.cohort .== UNTREATED) .& (sfs.tag .== s.untreated) .& (sfs.S_mean .> 0), :]
    lines!(ax, u.k, u.S_mean; LINE_UNTREATED...)
    for c in COHORTS                                     # up to the tumour size, as the untreated SFS
        d = sfs[(sfs.cohort .== c) .& (sfs.arm .== "single") .& (sfs.tag .== "blockend") .& (sfs.idx .== s.block) .&
                (sfs.S_mean .> 0) .& (sfs.k .<= maximum(u.k)), :]
        lines!(ax, d.k, d.S_mean; color = COL_COHORT[c], LINE[c]...)
    end
    savepanel(fig, NAME, p)
end

# ── c–f: treated / untreated at matched size ──────────────────────────────────
"Per-tumour ratios treated / untreated of each column for `cohort` at size `s`."
function ratios(cohort, s, cols)
    t = treated(s.block); t = t[t.cohort .== cohort, [:tumour; cols]]
    u = snap[(snap.cohort .== cohort) .& (snap.tag .== s.untreated), [:tumour; cols]]
    j = innerjoin(t, u; on = :tumour, renamecols = "_t" => "_u")
    return [filter(x -> isfinite(x) && x > 0, Float64.(j[!, "$(c)_t"]) ./ Float64.(j[!, "$(c)_u"])) for c in cols]
end

"Two panels (0.1N, N) of median and interquartile range of the ratios per q and cohort, same y-axis;
`source` (SFS, CSD) labels the distribution the Hill numbers come from."
function ratio_panels(cohorts, cols, names, source)
    r = Dict((c, i) => ratios(c, s, cols) for c in cohorts for (i, s) in enumerate(SIZES))
    v = [f(vals) for rs in values(r) for vals in rs if !isempty(vals) for f in (q25, q75)]
    lims = isempty(v) ? (0.5, 2.0) : (min(minimum(v), 1.0) / 1.15, max(maximum(v), 1.0) * 1.15)
    for (i, p) in enumerate(names)
        fig = panel()
        ylabel = rich(superscript("𝑞"), "𝐷", subscript(source), " (Treated/Untreated)")
        ax = Axis(fig[1, 1]; xlabel = LABELS[:q], ylabel, xticks = (1:length(cols), Q_LABELS), LOGY...)
        hlines!(ax, [1.0]; LINE_UNTREATED...)
        for (k, c) in enumerate(cohorts)
            x = Float64[]; med = Float64[]; lo = Float64[]; hi = Float64[]
            for (q, vals) in enumerate(r[(c, i)])
                isempty(vals) && continue
                push!(x, q + (k - (length(cohorts) + 1) / 2) * 0.15)
                push!(med, median(vals)); push!(lo, q25(vals)); push!(hi, q75(vals))
            end
            isempty(x) && continue
            errorbars!(ax, x, med, med .- lo, hi .- med; color = COL_COHORT[c])
            scatter!(ax, x, med; color = COL_COHORT[c])
        end
        xlims!(ax, 0.5, length(cols) + 0.5); ylims!(ax, lims...)
        savepanel(fig, NAME, p)
    end
end

ratio_panels(COHORTS, PASSENGER_Q, ("b", "e"), "SFS")
ratio_panels(RES_COHORTS, RESISTANT_Q, ("c", "f"), "CSD")

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 160, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)];
                   [[LineElement(; color = COL_COHORT[c], LINE[c]...), MarkerElement(color = COL_COHORT[c], marker = :circle, markersize = 4)]
                    for c in COHORTS]],
       ["Untreated"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

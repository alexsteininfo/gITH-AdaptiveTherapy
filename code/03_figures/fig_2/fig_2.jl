# Figure 2 — one treatment round leaves high diversity behind. Arm `single`: full dose until
# 0.1 N_det ("after treatment"), then off until regrowth to N_det ("relapse"); compared with the
# same tumour untreated at these sizes (growth record at 0.1 N_det, detection). Only tumours that
# respond (reach 0.1 N_det); the others had high resistance at detection and never left the
# treatment phase. Cohorts: no resistance (ν = 0), weak, partial, full. Theory: deterministic
# logistic trajectory, master formula (02_processing/09_theory.jl, theory_single_*).
# Top row after treatment (0.1N), bottom row after regrowth (N):
#   a, e  mean SFS: untreated (no-resistance cohort) and treated per cohort, theory dashed
#   b, f  compensated SFS k(k+1)·S_k (without the clonal bin), theory dashed
#   c, g  passenger Hill numbers ᵠD (q = 0: number of mutations S), treated / untreated, same tumour;
#         theory as bars
#   d, h  resistant-clone Hill numbers (clone size distribution), treated / untreated
#         (weak, partial, full; none in the no-resistance cohort)
#   legend
# Data: 01_dynamics/sfs_mean.csv, snapshots.csv, theory_single_sfs.csv, theory_single_hill.csv
# Run : julia --project=. code/03_figures/fig_2/fig_2.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_2"
const COHORTS = ["neutral", "weak", "partial", "full"]
const RES_COHORTS = ["weak", "partial", "full"]
const UNTREATED = "neutral"     # cohort of the untreated SFS (identical in law across cohorts)
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]
# matched sizes: single-round block (1: end of treatment, 2: regrown), untreated record (tag, idx)
# and the threshold of the theory tables
const SIZES = ((block = 1, untreated = "growth", uidx = 10, threshold = "lower"),
               (block = 2, untreated = "detection", uidx = 0, threshold = "upper"))
const THEORY = (color = :grey30, linestyle = :dash, linewidth = 0.6)

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")         # single round: responders only
theory_sfs = readproc("01_dynamics", "theory_single_sfs.csv")
theory_hill = readproc("01_dynamics", "theory_single_hill.csv")

"Single-round record of `block` (1: end of treatment, 2: regrown) of every responding tumour."
function treated(block)
    t = snap[(snap.arm .== "single") .& (snap.tag .== "blockend"), :]
    det = snap[snap.tag .== "detection", [:cohort, :tumour, :N]]
    b1 = innerjoin(t[t.idx .== 1, [:cohort, :tumour, :N]], det; on = [:cohort, :tumour], renamecols = "" => "_det")
    resp = b1[b1.N .<= b1.N_det ./ 10, [:cohort, :tumour]]
    return semijoin(t[t.idx .== block, :], resp; on = [:cohort, :tumour])
end

# ── a, b, e, f: mean SFS and compensated SFS ──────────────────────────────────
"Curves `(style, k, y)` at size `s`: untreated, treated per cohort, then theory (untreated, treated);
up to the tumour size (compensated: without the clonal bin)."
function sfs_curves(s; compensated)
    u = sfs[(sfs.cohort .== UNTREATED) .& (sfs.tag .== s.untreated) .& (sfs.idx .== s.uidx) .& (sfs.S_mean .> 0), :]
    kmax = maximum(u.k)
    keep(k) = compensated ? k .< 0.9kmax : k .<= kmax
    y(d) = compensated ? d.k .* (d.k .+ 1) .* d.S_mean : d.S_mean
    curves = Any[(LINE_UNTREATED, u.k[keep(u.k)], y(u)[keep(u.k)])]
    for c in COHORTS
        d = sfs[(sfs.cohort .== c) .& (sfs.arm .== "single") .& (sfs.tag .== "blockend") .& (sfs.idx .== s.block) .&
                (sfs.S_mean .> 0), :]
        push!(curves, ((; color = COL_COHORT[c], LINE_COHORT[c]...), d.k[keep(d.k)], y(d)[keep(d.k)]))
    end
    for kind in ("untreated", "treated")
        d = theory_sfs[(theory_sfs.kind .== kind) .& (theory_sfs.threshold .== s.threshold) .& (theory_sfs.S_mean .> 0), :]
        push!(curves, (THEORY, d.k[keep(d.k)], y(d)[keep(d.k)]))
    end
    return curves
end

for (s, (p, pc)) in zip(SIZES, (("a", "b"), ("e", "f")))
    for (name, compensated) in ((p, false), (pc, true))
        fig = panel()
        ax = compensated ?
             Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = "𝑘(𝑘+1)·𝑆ₖ", LOGX..., yscale = log10, yticks = fineticks,
                  ytickformat = logformat) :
             Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
        curves = sfs_curves(s; compensated)
        for (style, k, y) in curves
            lines!(ax, k, y; style...)
        end
        if compensated                                  # explicit, so that Supp 7B can share the y-axis of f
            v = [y for (_, _, ys) in curves for y in ys]
            ylims!(ax, minimum(v) / 1.05, maximum(v) * 1.05)
        end
        savepanel(fig, NAME, name)
    end
end

# ── c, d, g, h: treated / untreated at matched size ───────────────────────────
"Per-tumour ratios treated / untreated of each column for `cohort` at size `s`."
function ratios(cohort, s, cols)
    t = treated(s.block); t = t[t.cohort .== cohort, [:tumour; cols]]
    u = snap[(snap.cohort .== cohort) .& (snap.tag .== s.untreated) .& (snap.idx .== s.uidx), [:tumour; cols]]
    j = innerjoin(t, u; on = :tumour, renamecols = "_t" => "_u")
    return [filter(x -> isfinite(x) && x > 0, Float64.(j[!, "$(c)_t"]) ./ Float64.(j[!, "$(c)_u"])) for c in cols]
end

"Theory of the passenger fold changes (q = 0, 1, 2, ∞) at size `s`."
theory_ratios(s) = (r = only(eachrow(theory_hill[(theory_hill.kind .== "treated") .& (theory_hill.threshold .== s.threshold), :]));
                    [r.S_rel, r.D1n_rel, r.D2n_rel, r.Dinfn_rel])

"Two panels (0.1N, N) of median and interquartile range of the ratios per q and cohort, same y-axis;
`source` (SFS, CSD) labels the distribution the Hill numbers come from; `theory`: bars at the
predicted fold changes."
function ratio_panels(cohorts, cols, names, source; theory = false)
    r = Dict((c, i) => ratios(c, s, cols) for c in cohorts for (i, s) in enumerate(SIZES))
    v = [f(vals) for rs in values(r) for vals in rs if !isempty(vals) for f in (q25, q75)]
    theory && append!(v, vcat(theory_ratios.(SIZES)...))
    lims = isempty(v) ? (0.5, 2.0) : (min(minimum(v), 1.0) / 1.15, max(maximum(v), 1.0) * 1.15)
    for (i, p) in enumerate(names)
        fig = panel()
        ylabel = rich("Fold change in ", superscript("𝑞"), "𝐷", subscript(source))
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
        theory && for (q, y) in enumerate(theory_ratios(SIZES[i]))      # bar across the cohorts of each q, on top
            lines!(ax, [q - 0.34, q + 0.34], [y, y]; THEORY..., linewidth = 1.0)
        end
        xlims!(ax, 0.5, length(cols) + 0.5); ylims!(ax, lims...)
        savepanel(fig, NAME, p)
    end
end

ratio_panels(COHORTS, PASSENGER_Q, ("c", "g"), "SFS"; theory = true)
ratio_panels(RES_COHORTS, RESISTANT_Q, ("d", "h"), "CSD")

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 180, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)];
                   [[LineElement(; color = COL_COHORT[c], LINE_COHORT[c]...), MarkerElement(color = COL_COHORT[c], marker = :circle, markersize = 4)]
                    for c in COHORTS];
                   [LineElement(; THEORY..., linewidth = 0.8)]],
       ["Untreated"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]; "Theory"];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

# Figure 3 — diversity grows with every cycle of cyclic adaptive therapy (holidays: full dose until
# 0.5 N_det, off until regrowth to N_det). Only tumours that complete the first six cycles (all
# stored distributions) without progressing; detection is the upper threshold of cycle 0.
# Exploratory grids, one per threshold (lower 0.5 N_det, upper N_det), median and interquartile
# range over tumours per resistance cohort:
#   grid_<threshold>      Hill numbers treated / untreated: relative to the same tumour untreated at
#                         the same size (growth record at 0.5 N_det, detection); top row passengers
#                         (SFS, q = 0: number of mutations S), bottom row resistant clones (CSD)
#   grid_<threshold>_raw  the same, absolute values
#   grid_both             treated / untreated, both thresholds in one panel (upper solid/filled, lower dashed/open)
#   sfs                   mean SFS at the lower and upper threshold of every cycle, one panel per
#                         cohort, paler with every cycle
#   c  modulation: total mutational burden B, split into mutations on lineages born before
#      (pre-existing, B_pre) and after detection (new, B_new); dotted: 2μN ln(1 + t), new mutations at constant size N
#   d  modulation: mean SFS at detection, on reaching the stable size and later
# Data: 01_dynamics/snapshots.csv, sfs_mean.csv
# Run : julia --project=. code/03_figures/fig_3/fig_3.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_3"
const MU = 1.0      # neutral mutations per daughter cell (code/X_helpers/params.jl)
const NCYCLES = 6   # holiday cycles with stored distributions (HEAVY_BLOCKS = 12 in params.jl)
const COHORTS = ["neutral", "weak", "partial", "full"]
const RES_COHORTS = ["weak", "partial", "full"]
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]
const THRESHOLDS = ("lower", "upper")
const COL_THRESHOLD = Dict("lower" => PAIRED[8], "upper" => PAIRED[2])
const THRESHOLD_NAME = Dict("lower" => "Lower (0.5𝑁)", "upper" => "Upper (𝑁)")

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")
any((snap.tag .== "growth") .& (snap.idx .== 50)) ||
    error("no untreated records at 0.5 N_det: run code/01_datageneration/run_part1_growth.sh")

# ── Hill numbers per cycle ────────────────────────────────────────────────────
"""
Holiday threshold records of the tumours that complete all `NCYCLES` cycles without progressing
(a block ending in progression ends above N_det), plus their detection record (upper threshold of
cycle 0), with `cycle`, `threshold` and every Hill number relative to the same tumour untreated at
the threshold size (`<col>_rel`; lower: growth record at 0.5 N_det, upper: detection).
"""
function cycles()
    det = snap[snap.tag .== "detection", :]
    hol = snap[(snap.arm .== "hol") .& (snap.tag .== "blockend") .& (snap.idx .<= 2NCYCLES), :]
    ndet = Dict((r.cohort, r.tumour) => r.N for r in eachrow(det))
    ok = combine(groupby(hol, [:cohort, :tumour]),
                 [:cohort, :tumour, :N] => ((c, t, n) -> length(n) == 2NCYCLES && all(n .<= ndet[(c[1], t[1])])) => :ok)
    ok = ok[ok.ok, [:cohort, :tumour]]
    d = semijoin(vcat(det, hol), ok; on = [:cohort, :tumour])
    d.cycle = [r.tag == "detection" ? 0 : (r.idx + 1) ÷ 2 for r in eachrow(d)]
    d.threshold = [r.tag == "detection" || iseven(r.idx) ? "upper" : "lower" for r in eachrow(d)]
    cols = [PASSENGER_Q; RESISTANT_Q]
    untreated = snap[(snap.tag .== "detection") .| ((snap.tag .== "growth") .& (snap.idx .== 50)), :]
    untreated.threshold = [r.tag == "detection" ? "upper" : "lower" for r in eachrow(untreated)]
    d = innerjoin(d, untreated[:, [:cohort; :tumour; :threshold; cols]];
                  on = [:cohort, :tumour, :threshold], renamecols = "" => "_untreated")
    for c in cols
        d[!, "$(c)_rel"] = Float64.(d[!, c]) ./ Float64.(d[!, "$(c)_untreated"])
    end
    println("tumours completing $NCYCLES cycles: ", Dict(r.cohort => r.n for r in eachrow(combine(groupby(ok, :cohort), nrow => :n))))
    return d
end

"Median and quartiles of `col` (finite, positive values) per cycle of `cohort` at `threshold`."
function percycle(d, cohort, threshold, col)
    s = d[(d.cohort .== cohort) .& (d.threshold .== threshold), :]
    g = combine(groupby(s, :cycle), col => (v -> (f = filter(x -> isfinite(x) && x > 0, Float64.(v));
                                                   isempty(f) ? (med = NaN, lo = NaN, hi = NaN, n = 0) :
                                                   (med = median(f), lo = q25(f), hi = q75(f), n = length(f)))) => AsTable)
    return sort!(g[g.n .> 0, :], :cycle)
end

"Grid of Hill numbers at `thresholds` (one, or both: upper solid with filled markers, lower dashed
with open markers): passengers (top) and resistant clones (bottom) per order q, treated / untreated
at the same size (`rel`) or absolute; y axes shared within a row when relative."
function hill_grid(d, thresholds; rel)
    both = length(thresholds) > 1
    style(th) = both && th == "lower" ? (linestyle = :dash, open = true) : (linestyle = :solid, open = false)
    fig = panel(; w = 4PANEL_MM[1], h = 2PANEL_MM[2] + 8)
    rows = ((PASSENGER_Q, COHORTS, "SFS", "Passengers"), (RESISTANT_Q, RES_COHORTS, "CSD", "Resistant clones"))
    for (i, (cols, cohorts, source, what)) in enumerate(rows)
        ylabel = rel ? rich(superscript("𝑞"), "𝐷", subscript(source), " (Treated/Untreated)") :
                       rich(superscript("𝑞"), "𝐷", subscript(source))
        Label(fig[i, 0], what; rotation = π / 2, fontsize = 8, font = :bold, tellheight = false)
        axs = Axis[]
        for (j, col) in enumerate(cols)
            ax = Axis(fig[i, j]; xlabel = i == 2 ? LABELS[:cycle] : "", ylabel = j == 1 ? ylabel : "",
                      title = i == 1 ? "𝑞 = $(Q_LABELS[j])" : "", titlesize = 8, titlefont = :regular,
                      xticks = 0:NCYCLES, LOGY...)
            rel && hlines!(ax, [1.0]; LINE_UNTREATED...)
            y = rel ? "$(col)_rel" : col
            for th in thresholds, (k, c) in enumerate(cohorts)
                g = percycle(d, c, th, y)
                isempty(g) && continue
                st = style(th)
                x = g.cycle .+ (k - (length(cohorts) + 1) / 2) * 0.15 .+ (both ? (th == "lower" ? -0.05 : 0.05) : 0.0)
                lines!(ax, x, g.med; color = COL_COHORT[c], linewidth = 0.8, linestyle = st.linestyle)
                errorbars!(ax, x, g.med, g.med .- g.lo, g.hi .- g.med; color = COL_COHORT[c])
                scatter!(ax, x, g.med; color = st.open ? :white : COL_COHORT[c], strokecolor = COL_COHORT[c],
                         strokewidth = st.open ? 0.8 : 0)
            end
            xlims!(ax, -0.5, NCYCLES + 0.5)
            push!(axs, ax)
        end
        rel && linkyaxes!(axs...)
    end
    elements = [[LineElement(color = COL_COHORT[c], linewidth = 0.8), MarkerElement(color = COL_COHORT[c], marker = :circle, markersize = 4)]
                for c in COHORTS]
    labels = [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]
    if both
        for th in ("upper", "lower")
            st = style(th)
            push!(elements, [LineElement(color = :black, linewidth = 0.8, linestyle = st.linestyle),
                             MarkerElement(color = st.open ? :white : :black, strokecolor = :black,
                                           strokewidth = st.open ? 0.8 : 0, marker = :circle, markersize = 4)])
            push!(labels, "$(THRESHOLD_NAME[th]) threshold")
        end
    end
    Legend(fig[3, 1:4], elements, labels; orientation = :horizontal, tellwidth = false, tellheight = true)
    Label(fig[0, 1:4], both ? "Both thresholds" : "$(THRESHOLD_NAME[only(thresholds)]) threshold"; fontsize = 8, font = :bold)
    return fig
end

hol = cycles()
for th in THRESHOLDS
    savepanel(hill_grid(hol, (th,); rel = true), NAME, "grid_$th")
    savepanel(hill_grid(hol, (th,); rel = false), NAME, "grid_$(th)_raw")
end
savepanel(hill_grid(hol, THRESHOLDS; rel = true), NAME, "grid_both")

# ── Mean SFS per cycle ────────────────────────────────────────────────────────
"Colour `c` blended towards white by fraction `f`."
pale(c, f) = (c = Makie.to_color(c); Makie.RGBAf(c.r + f * (1 - c.r), c.g + f * (1 - c.g), c.b + f * (1 - c.b), 1))

const NDET = maximum(sfs.k[sfs.arm .== "detection"])

"Mean SFS of every cycle at both thresholds, one panel per cohort; `compensated`: k(k+1)·S_k without
the clonal bin (flat for the neutral expectation 2μN/(k(k+1)) of a growing population)."
function sfs_cycles(; compensated)
    fig = panel(; w = 4PANEL_MM[1], h = PANEL_MM[2] + 6)
    ylabel = compensated ? "𝑘(𝑘+1)·𝑆ₖ" : LABELS[:Sk]
    axs = Axis[]
    for (j, c) in enumerate(COHORTS)
        ax = Axis(fig[1, j]; xlabel = LABELS[:k], ylabel = j == 1 ? ylabel : "", LOGX..., LOGY...,
                  title = c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance", titlesize = 8, titlefont = :regular)
        s = sfs[(sfs.cohort .== c) .& (sfs.S_mean .> 0), :]
        for th in THRESHOLDS, cyc in NCYCLES:-1:0                       # pale (late) first, dark on top
            (th == "lower" && cyc == 0) && continue
            r = cyc == 0 ? s[s.arm .== "detection", :] :
                s[(s.arm .== "hol") .& (s.tag .== "blockend") .& (s.idx .== 2cyc - (th == "lower")), :]
            kmax = (th == "lower" ? 0.5 : 1.0) * NDET                   # progressed tumours are larger
            r = r[compensated ? (r.k .< 0.9kmax) : (r.k .<= kmax), :]    # compensated: without the clonal bin
            y = compensated ? r.k .* (r.k .+ 1) .* r.S_mean : r.S_mean
            lines!(ax, r.k, y; color = pale(COL_THRESHOLD[th], 0.85cyc / NCYCLES), linewidth = 0.8)
        end
        push!(axs, ax)
    end
    linkaxes!(axs...)
    Legend(fig[0, 1:4], [LineElement(color = pale(COL_THRESHOLD[th], f)) for th in THRESHOLDS for f in (0.0, 0.85)],
           ["$(THRESHOLD_NAME[th]), cycle $(cyc)" for th in THRESHOLDS for cyc in (th == "lower" ? 1 : 0, NCYCLES)];
           orientation = :horizontal, tellwidth = false, tellheight = true)
    return fig
end

savepanel(sfs_cycles(; compensated = false), NAME, "sfs")
savepanel(sfs_cycles(; compensated = true), NAME, "sfs_compensated")

# ── c, d: modulation ──────────────────────────────────────────────────────────
snap = snap[snap.cohort .== "neutral", :]; sfs = sfs[sfs.cohort .== "neutral", :]

"Modulation records over time: detection, stable size reached (block end 1), grid times, end."
function modulation()
    d = snap[(snap.tag .== "detection") .| ((snap.arm .== "mod") .& ((snap.tag .== "grid") .| (snap.tag .== "end") .|
             ((snap.tag .== "blockend") .& (snap.idx .== 1)))), :]
    d.point = [r.tag == "detection" ? "detection" : "$(r.tag)_$(r.idx)" for r in eachrow(d)]
    return d
end

mods = modulation()
times = sort(combine(groupby(mods, :point), :t_rel => median => :t, :N => median => :N), :t)

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:B])
for (y, col, label) in ((:B, COL_ALL, "Total"), (:B_pre, COL_PRE, "Pre-existing"), (:B_new, COL_NEW, "New"))
    g = sort(combine(groupby(mods, :point), :t_rel => median => :t, y => median => :med, y => q25 => :lo,
                     y => q75 => :hi), :t)
    medianband!(ax, g.t, g.med, g.lo, g.hi; color = col, label)
end
Ns = times.N[findfirst(!=("detection"), times.point)]
tt = range(0, maximum(times.t); length = 200)
lines!(ax, tt, 2MU * Ns .* log.(1 .+ tt); color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
toplegend(fig, ax)
savepanel(fig, NAME, "c")

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
sel = times.point[unique(round.(Int, range(1, nrow(times); length = min(4, nrow(times)))))]
cols = shades(COL_ARM["mod"], length(sel))
for (p, col) in zip(sel, cols)
    d = p == "detection" ? sfs[sfs.arm .== "detection", :] :
        sfs[(sfs.arm .== "mod") .& (string.(sfs.tag, "_", sfs.idx) .== p), :]
    d = d[d.S_mean .> 0, :]
    t = times.t[times.point .== p][1]
    label = p == "detection" ? "Detection" : p == "blockend_1" ? "Stable size" : "t = $(round(Int, t))"
    lines!(ax, d.k, d.S_mean; color = col, label)
end
axislegend(ax; position = :rt)
savepanel(fig, NAME, "d")

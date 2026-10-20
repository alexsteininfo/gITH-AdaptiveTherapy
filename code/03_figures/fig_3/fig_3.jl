# Figure 3 — holidays accumulate diversity, modulation removes it transiently (neutral cohort).
#   a  holidays: number of mutations S at the lower (0.5 N_det) and upper (N_det) threshold of each
#      cycle; median and interquartile range over tumours
#   b  holidays: tail exponent γ̂ of the SFS at the same times (k_min = 0.01 N)
#   c  modulation: total mutational burden B, split into mutations on lineages born before
#      (pre-existing, B_pre) and after detection (new, B_new); dotted: 2μN ln(1 + t), new mutations at constant size N
#   d  modulation: mean SFS at detection, on reaching the stable size and later
# Data: 01_dynamics/snapshots.csv, sfs_mean.csv
# Run : julia --project=. code/03_figures/fig_3/fig_3.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_3"
const MU = 1.0      # neutral mutations per daughter cell (code/X_helpers/params.jl)

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")
snap = snap[snap.cohort .== "neutral", :]; sfs = sfs[sfs.cohort .== "neutral", :]

# ── a, b: holiday cycles ──────────────────────────────────────────────────────
"Holiday records with cycle number and threshold (detection counts as the upper threshold of cycle 0)."
function cycles()
    det = snap[snap.tag .== "detection", :]
    hol = snap[(snap.arm .== "hol") .& (snap.tag .== "blockend"), :]
    d = vcat(det, hol)
    d.cycle = [r.tag == "detection" ? 0 : (r.idx + 1) ÷ 2 for r in eachrow(d)]
    d.threshold = [r.tag == "detection" || iseven(r.idx) ? "up" : "low" for r in eachrow(d)]
    return d
end

function cycle_panel(y; ylabel)
    d = cycles()
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:cycle], ylabel)
    for (th, col, label) in (("low", PAIRED[1], "Lower"), ("up", PAIRED[2], "Upper"))
        g = combine(groupby(d[d.threshold .== th, :], :cycle), y => (v -> median(filter(isfinite, v))) => :med,
                    y => (v -> q25(filter(isfinite, v))) => :lo, y => (v -> q75(filter(isfinite, v))) => :hi)
        sort!(g, :cycle)
        lines!(ax, g.cycle, g.med; color = col)
        errorbars!(ax, g.cycle, g.med, g.med .- g.lo, g.hi .- g.med; color = col)
        scatter!(ax, g.cycle, g.med; color = col, label)
    end
    return fig, ax
end

fig, ax = cycle_panel(:S; ylabel = LABELS[:S]); toplegend(fig, ax); savepanel(fig, NAME, "a")
fig, ax = cycle_panel(:gamma; ylabel = LABELS[:gamma]); savepanel(fig, NAME, "b")

# ── c, d: modulation ──────────────────────────────────────────────────────────
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

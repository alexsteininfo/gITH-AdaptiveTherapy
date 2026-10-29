# Stabilising adaptive therapy (arm `mod`: full dose until 0.8 N_det, then the dose that keeps the
# size constant), no-resistance cohort; moved out of Figure 3, placement in the paper open.
#   a  total mutational burden B, split into mutations on lineages born before (pre-existing, B_pre)
#      and after detection (new, B_new); dotted: 2μN ln(1 + t), new mutations at constant size N
#   b  mean SFS at detection, on reaching the stable size and later
# Data: 01_dynamics/snapshots.csv, sfs_mean.csv
# Run : julia --project=. code/03_figures/stabilising/stabilising.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "stabilising"
const MU = 1.0      # neutral mutations per daughter cell (code/X_helpers/params.jl)

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")
# ── a, b ──────────────────────────────────────────────────────────────────────
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
savepanel(fig, NAME, "a")

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
savepanel(fig, NAME, "b")

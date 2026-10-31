# Supplementary Figure 10 — the balance of the number of mutations under stabilising adaptive therapy
# (Figure 4B, C): mutations from before and after detection, and the singletons that set the loss,
# dE[B]/dτ = b(2μN* − E S₁) at constant size. Neutral cohort of block 08_stabilising (200 tumours),
# arms mod_0.5, mod and mod_1.0 (stable sizes 0.5, 0.8, 1.0 N_det); time τ since stabilisation (end
# of the full-dose phase, median over the arm's tumours), first 10 time units. Theory: master formula
# on growth → reduction → hold (theory_stab_hill.csv: S_old born during growth, S_red + S_new after).
#   a, b  number of mutations relative to the same tumour at detection, B / B_det, at 0.5 (a) and
#         0.8 N_det (b): total, pre-existing (born before detection, B_pre) and new (born after it,
#         B_new, including the full-dose phase), labelled on the curves; median and interquartile
#         range, theory dotted
#   c     loss / gain of mutations, S₁ / 2μN*: mean singletons over tumours (sfs_mean.csv) at the
#         three stable sizes, theory dotted; B decreases while the ratio is above 1 (dashed)
#   d     the same ratio at stabilisation vs stable size: theory (dotted), simulated sizes (dots);
#         vertical line where it crosses 1 (≈ 0.72 N_det)
#   legend
# Data: 08_stabilising/snapshots.csv, sfs_mean.csv, theory_stab_hill.csv, theory_stab_levels.csv
# Run : julia --project=. code/03_figures/supp_10/supp_10.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_10"
const BLOCK = "08_stabilising"
const LEVEL_ARMS = ["mod_0.5" => 0.5, "mod" => 0.8, "mod_1.0" => 1.0]
const COL_LEVEL = Dict(0.5 => PAIRED[8], 0.8 => PAIRED[4], 1.0 => PAIRED[10])     # as in Figure 4
const LEVEL_NAME = Dict(0.5 => rich("𝑁* = 0.5 𝑁", subscript("det")), 0.8 => rich("𝑁* = 0.8 𝑁", subscript("det")),
                        1.0 => rich("𝑁* = 𝑁", subscript("det")))
const PARTS = [:B => "Total", :B_pre => "Pre-existing", :B_new => "New"]
# curve labels of a, b: time and offset from the theory curve (above > 0, below < 0)
const LABEL_AT = Dict(0.5 => Dict(:B => (8.0, 0.08), :B_pre => (7.0, 0.09), :B_new => (8.0, 0.08)),
                      0.8 => Dict(:B => (8.0, -0.09), :B_pre => (7.0, 0.09), :B_new => (8.0, -0.09)))
const MU = 1.0          # neutral mutations per daughter cell (code/X_helpers/params.jl)
const T_EARLY = 10.0
const LABEL_T = "Time since stabilisation, 𝜏"
const LABEL_RATIO = rich("Loss / gain, 𝑆", subscript("1"), " / 2𝜇𝑁*")
const LINE_BALANCE = (color = :black, linestyle = :dash, linewidth = 0.8)

snap = readproc(BLOCK, "snapshots.csv")
sfs = readproc(BLOCK, "sfs_mean.csv")
th = readproc(BLOCK, "theory_stab_hill.csv")
thl = sort(readproc(BLOCK, "theory_stab_levels.csv"), :level)
const NDET = maximum(snap.N[snap.tag .== "detection"])

"""
Records of the neutral cohort under `arm` from stabilisation on (end of the full-dose phase, fine
records) with `tau` (time since stabilisation; `tred` = median end of the full-dose phase) and
`B`, `B_pre`, `B_new` relative to the same tumour's `B` at detection (`_rel`).
"""
function records(arm)
    d = snap[(snap.cohort .== "neutral") .& (snap.arm .== arm) .&
             (((snap.tag .== "blockend") .& (snap.idx .== 1)) .| (snap.tag .== "fine")), :]
    tred = median(d.t_rel[d.tag .== "blockend"])
    d.tau = [r.tag == "fine" ? r.idx / 100 - tred : 0.0 for r in eachrow(d)]
    det = snap[(snap.cohort .== "neutral") .& (snap.tag .== "detection"), [:tumour, :B]]
    d = innerjoin(d, rename(det, :B => :B_det); on = :tumour)
    for c in (:B, :B_pre, :B_new)
        d[!, "$(c)_rel"] = d[!, c] ./ d.B_det
    end
    return (d = d[d.tau .<= T_EARLY, :], tred = tred)
end

"Median and quartiles of `col` per time since stabilisation."
quartiles(d, col) = sort!(combine(groupby(d, :tau), col => median => :med, col => q25 => :lo, col => q75 => :hi), :tau)

th_tr(level) = sort(th[(th.level .== level) .& (th.kind .== "treated"), :], :tau)
th_det(level) = only(th[(th.level .== level) .& (th.kind .== "detection"), :S])

"Singletons S₁ / 2μN* of the neutral cohort's mean SFS under `arm`: at stabilisation (`which = :reduced`) or at the fine record `which` (time since detection)."
function singletons(arm, L, which)
    s = sfs[(sfs.cohort .== "neutral") .& (sfs.arm .== arm) .& (sfs.k .== 1), :]
    r = which === :reduced ? s[(s.tag .== "blockend") .& (s.idx .== 1), :] : s[(s.tag .== "fine") .& (s.idx .== round(Int, 100which)), :]
    return nrow(r) == 1 ? r.S_mean[1] / (2MU * L * NDET) : NaN
end

"Stable size at which `r` (vs `L`) crosses 1, by linear interpolation; NaN if it does not."
function crossing(L, r)
    i = findfirst(j -> (r[j] - 1) * (r[j+1] - 1) <= 0, 1:length(r)-1)
    i === nothing && return NaN
    return L[i] + (1 - r[i]) * (L[i+1] - L[i]) / (r[i+1] - r[i])
end

recs = Dict(a => records(a) for (a, _) in LEVEL_ARMS)

# ── a, b: pre-existing and new mutations ─────────────────────────────────────
const SPLIT = (("a", "mod_0.5", 0.5), ("b", "mod", 0.8))
ymax = 1.05 * maximum(maximum(quartiles(recs[a].d, "$(c)_rel").hi) for (_, a, _) in SPLIT for (c, _) in PARTS)
for (p, arm, L) in SPLIT
    local fig = panel()
    local ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = rich("Number of mutations, 𝐵 / 𝐵", subscript("det")))
    hlines!(ax, [1.0]; LINE_UNTREATED...)
    t = th_tr(L); keep = t.tau .<= T_EARLY
    theory = Dict(:B => t.S, :B_pre => t.S_old, :B_new => t.S_red .+ t.S_new)
    for (col, label) in PARTS
        g = quartiles(recs[arm].d, "$(col)_rel")
        medianband!(ax, g.tau, g.med, g.lo, g.hi; color = COL_LEVEL[L])
        y = theory[col][keep] ./ th_det(L)
        lines!(ax, t.tau[keep], y; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
        τl, dy = LABEL_AT[L][col]
        text!(ax, τl, y[argmin(abs.(t.tau[keep] .- τl))] + dy; text = label, fontsize = 7, align = (:center, :center))
    end
    xlims!(ax, -0.3, T_EARLY); ylims!(ax, 0, ymax)
    savepanel(fig, NAME, p)
end

# ── c: loss / gain over time ─────────────────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = LABEL_RATIO)
hlines!(ax, [1.0]; LINE_BALANCE...)
for (arm, L) in LEVEL_ARMS
    tred = recs[arm].tred
    times = sort(unique(sfs.idx[(sfs.arm .== arm) .& (sfs.tag .== "fine")])) ./ 100
    whichs = [[:reduced]; [t for t in times if t - tred <= T_EARLY]]
    ts = [w === :reduced ? 0.0 : w - tred for w in whichs]
    scatterlines!(ax, ts, [singletons(arm, L, w) for w in whichs]; color = COL_LEVEL[L], markersize = 3, linewidth = 0.6)
    t = th_tr(L); keep = t.tau .<= T_EARLY
    lines!(ax, t.tau[keep], t.S1[keep] ./ (2MU * L * NDET); color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
end
xlims!(ax, -0.3, T_EARLY)
savepanel(fig, NAME, "c")

# ── d: loss / gain at stabilisation vs stable size ───────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = rich("Stable size, 𝑁* / 𝑁", subscript("det")), ylabel = LABEL_RATIO)
hlines!(ax, [1.0]; LINE_BALANCE...)
lc = crossing(thl.level, thl.r0)
isfinite(lc) && vlines!(ax, [lc]; color = :grey60, linewidth = 0.6)
isfinite(lc) && text!(ax, lc + 0.01, 0.95maximum(thl.r0); text = string(round(lc; digits = 2)), fontsize = 7, align = (:left, :center))
lines!(ax, thl.level, thl.r0; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
for (arm, L) in LEVEL_ARMS
    scatter!(ax, [L], [singletons(arm, L, :reduced)]; color = COL_LEVEL[L], markersize = 6)
end
xlims!(ax, 0.25, 1.05)
savepanel(fig, NAME, "d")

# ── legend ───────────────────────────────────────────────────────────────────
fig = panel(; w = 180, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)]; [LineElement(; LINE_BALANCE...)];
                   [[LineElement(color = COL_LEVEL[L], linewidth = 1.2), MarkerElement(color = COL_LEVEL[L], marker = :circle, markersize = 4)]
                    for (_, L) in LEVEL_ARMS];
                   [LineElement(color = COL_PRED, linestyle = LS_PRED, linewidth = 1)]],
       ["Detection"; "Loss = gain"; [LEVEL_NAME[L] for (_, L) in LEVEL_ARMS]; "Theory"];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

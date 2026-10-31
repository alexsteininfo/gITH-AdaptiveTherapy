# Figure 4 — stabilising adaptive therapy (arm `mod`: full dose until 0.8 N_det, then the dose that
# holds the size): the transient of the mutation load, the moving front of the SFS and the Hill
# numbers over time. Block 08_stabilising: 200 tumours per cohort, distributions every 0.25 time
# units at first; the neutral cohort also at the stable sizes 0.5 and 1.0 N_det. Theory: master
# formula on the deterministic path growth → reduction → hold (code/X_helpers/theory.jl,
# .docs/theory/STABILISING.md).
# Time τ = time since stabilisation (the end of the full-dose reduction, when the tumour first
# reaches the stable size; median over the arm's tumours: 0.058 after detection at 0.8 N_det, 0.183
# at 0.5 N_det, 0 at 1.0 N_det — the spread between tumours is ~10⁻⁴). Panels over time start at
# stabilisation (τ = 0; the full-dose phase is not shown) but stay relative to detection; the fine
# records, stored at fixed times t since detection, lie at τ = t − t_red.
#   a        ⁰D_SFS treated / detection over the first 10 time units of the hold (from stabilisation), no-resistance
#            cohort at the stable sizes 0.5, 0.8, 1.0 N_det (median, interquartile range), theory dotted
#   b        singletons S₁ / 2μN* over the same time (the mechanism of a: the load falls while S₁ is
#            above its stationary value 1, dashed), simulation points and theory
#   c_tmin   time of the minimum of ⁰D (since stabilisation, τ_min) vs stable size: theory, simulated levels
#   c_smin   ⁰D at its minimum (dotted) and at stabilisation (dashed) relative to detection vs stable
#            size: theory, simulated levels
#   d_sfs    mean SFS, no-resistance cohort, at detection, at stabilisation and τ ≈ 1, 2, 5, 10, 20,
#            40 (records at t = 1, 2, 5, … after detection, i.e. τ = 0.94, 1.94, …; paler with time)
#   d        compensated SFS k(k+1)·S_k of d_sfs (without the clonal bin), theory dashed: the front
#   e        position of the front (peak of d) over time per cohort, theory dashed, bτ dotted
#   f–i      passenger ⁰D, ¹D, ²D, ∞D treated / detection over time per cohort, theory dotted
#   f_size–i_size  the same relative to the untreated tumour at the stable size (Fig 3's convention)
#   j–m      resistant-clone ⁰D, ¹D, ²D, ∞D treated / detection over time: weak, partial, full
#   dose     dose over time per cohort (log): the rise that ends the hold under full/partial resistance
#   legend (cohorts, f–m), legend_levels (a–c; vertical, transparent, for the free corner of c_tmin), legend_sfs (d)
# Hill-number axes show the fold change (treated / detection, or / untreated at N*). No and weak
# resistance overlap, and so do partial and full: weak and full are dashed and drawn on top.
# Data: 08_stabilising/snapshots.csv, sfs_mean.csv, counts.csv.gz, theory_stab_hill.csv,
#       theory_stab_sfs.csv, theory_stab_levels.csv
# Run : julia --project=. code/03_figures/fig_4/fig_4.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_4"
const BLOCK = "08_stabilising"
const COHORTS = ["neutral", "weak", "partial", "full"]
const RES_COHORTS = ["weak", "partial", "full"]
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]
const LEVEL_ARMS = ["mod_0.5" => 0.5, "mod" => 0.8, "mod_1.0" => 1.0]
const COL_LEVEL = Dict(0.5 => PAIRED[8], 0.8 => PAIRED[4], 1.0 => PAIRED[10])
const LEVEL_NAME = Dict(0.5 => rich("𝑁* = 0.5 𝑁", subscript("det")), 0.8 => rich("𝑁* = 0.8 𝑁", subscript("det")),
                        1.0 => rich("𝑁* = 𝑁", subscript("det")))
const LS_COHORT = Dict("neutral" => :solid, "weak" => :dash, "partial" => :solid, "full" => :dash)   # overlapping pairs: dashed on top
const LABEL_LEVEL = rich("Stable size, 𝑁* / 𝑁", subscript("det"))
const MU = 1.0          # neutral mutations per daughter cell (code/X_helpers/params.jl)
const B0 = 1.0          # birth rate
const LABEL_T = "Time since stabilisation, 𝜏"
const T_EARLY = 10.0    # window of a, b (time since stabilisation)

"Axis label of a Hill-number fold change, e.g. Fold change in ⁰𝐷_SFS."
foldlabel(q, source) = rich("Fold change in ", superscript(q), "𝐷", subscript(source))

snap = readproc(BLOCK, "snapshots.csv")
sfs = readproc(BLOCK, "sfs_mean.csv")
counts = readproc(BLOCK, "counts.csv.gz")
th = readproc(BLOCK, "theory_stab_hill.csv")
ths = readproc(BLOCK, "theory_stab_sfs.csv")
thl = readproc(BLOCK, "theory_stab_levels.csv")
const NDET = maximum(snap.N[snap.tag .== "detection"])
const TEND = maximum(snap.idx[snap.tag .== "fine"]) / 100            # last fine record (40; 8 on test data)
const SFS_TIMES = [[:detection, :reduced]; [t for t in (1, 2, 5, 10, 20, 40) if t <= TEND]]

# ── Records over time with per-tumour baselines ───────────────────────────────
"""
Snapshot records of arm `armid` from stabilisation on (end of the reduction, fine records) with the
time since detection `t` and since stabilisation `tau = t − tred` (`tred`: median end of the
reduction over the arm's tumours), and every Hill number relative to the same tumour at detection
(`_rel_det`) and, for passengers, to the same tumour untreated at the stable size 0.8 N_det
(`_rel_size`, growth record idx 80).
"""
function records(armid)
    det = snap[snap.tag .== "detection", :]
    d = snap[(snap.arm .== armid) .& (((snap.tag .== "blockend") .& (snap.idx .== 1)) .| (snap.tag .== "fine")), :]
    tred = median(d.t_rel[d.tag .== "blockend"])                        # end of the reduction: one common time
    d.t = [r.tag == "fine" ? r.idx / 100 : tred for r in eachrow(d)]
    d.tred .= tred; d.tau = d.t .- tred
    cols = [PASSENGER_Q; RESISTANT_Q]
    d = innerjoin(d, det[:, [:cohort; :tumour; cols]]; on = [:cohort, :tumour], renamecols = "" => "_det")
    unt = snap[(snap.tag .== "growth") .& (snap.idx .== 80), :]
    d = leftjoin(d, unt[:, [:cohort; :tumour; PASSENGER_Q]]; on = [:cohort, :tumour], renamecols = "" => "_size")
    for c in cols
        d[!, "$(c)_rel_det"] = Float64.(d[!, c]) ./ Float64.(d[!, "$(c)_det"])
    end
    for c in PASSENGER_Q
        d[!, "$(c)_rel_size"] = Float64.(d[!, c]) ./ Float64.(coalesce.(d[!, "$(c)_size"], NaN))
    end
    return sort!(d, [:cohort, :tumour, :t])
end

"Median and quartiles of `col` (finite, positive) per time `tau` of `cohort`, while at least half of its tumours have a record."
function overtime(d, cohort, col)
    s = d[d.cohort .== cohort, :]
    ntotal = length(unique(s.tumour))
    g = combine(groupby(s, :tau), col => (v -> (f = filter(x -> isfinite(x) && x > 0, Float64.(v));
                                                 isempty(f) ? (med = NaN, lo = NaN, hi = NaN, n = 0) :
                                                 (med = median(f), lo = q25(f), hi = q75(f), n = length(f)))) => AsTable)
    return sort!(g[g.n .>= ntotal / 2, :], :tau)
end

# ── Theory ────────────────────────────────────────────────────────────────────
th_tr(level) = sort(th[(th.level .== level) .& (th.kind .== "treated"), :], :tau)
"Duration of the reduction phase on the theory path (detection at τ = −t_red)."
th_tred(level) = (t = th_tr(level); t.t_rel[1] - t.tau[1])

"Theory of `col` relative to detection (`_rel_det`) or to the untreated tumour at N* (`_rel_size`) vs time since stabilisation."
function theory_time(level, col; rel = "_rel_det")
    t = th_tr(level)
    return t.tau, t[!, "$(col)$(rel)"]
end

# ── Mean SFS ──────────────────────────────────────────────────────────────────
"Mean SFS of `cohort` under `arm`: `which` = :detection, :reduced (end of the reduction) or a time (fine record); `compensated`: k(k+1)·S_k without the clonal bin."
function sfs_rec(cohort, which; arm = "mod", level = 0.8, compensated = false)
    s = sfs[(sfs.cohort .== cohort) .& (sfs.S_mean .> 0), :]
    r = which === :detection ? s[s.arm .== "detection", :] :
        which === :reduced ? s[(s.arm .== arm) .& (s.tag .== "blockend") .& (s.idx .== 1), :] :
        s[(s.arm .== arm) .& (s.tag .== "fine") .& (s.idx .== round(Int, 100which)), :]
    kmax = (which === :detection ? 1.0 : level) * NDET
    r = r[compensated ? (r.k .< 0.9kmax) : (r.k .<= kmax), :]
    return r.k, compensated ? r.k .* (r.k .+ 1) .* r.S_mean : r.S_mean
end

"Theory of the mean SFS at `level` (same keys and cuts as `sfs_rec`; a time `which` is a fine record's time since detection, evaluated at τ = which − t_red; `tau` gives a hold time directly)."
function theory_sfs(level, which; compensated = false, tau = nothing)
    s = ths[(ths.level .== level) .& (ths.S_mean .> 0), :]
    target = tau !== nothing ? Float64(tau) : which === :reduced ? 0.0 : which === :detection ? NaN : Float64(which) - th_tred(level)
    s = which === :detection && tau === nothing ? s[s.kind .== "detection", :] :
        s[(s.kind .== "treated") .& isapprox.(s.tau, target; atol = 1e-6), :]
    kmax = (which === :detection ? 1.0 : level) * NDET
    s = s[compensated ? (s.k .< 0.9kmax) : (s.k .<= kmax), :]
    return s.k, compensated ? s.k .* (s.k .+ 1) .* s.S_mean : s.S_mean
end

"Position of the maximum of `y` over `k` below 10⁴ (vertex of a parabola in log–log through the largest bin and its neighbours); NaN without data."
function peakpos(k, y)
    keep = (k .< 1e4) .& (y .> 0); x = log.(k[keep]); z = log.(y[keep])
    isempty(z) && return NaN
    i = argmax(z)
    (i == 1 || i == length(z)) && return exp(x[i])
    (x1, x2, x3), (y1, y2, y3) = x[i-1:i+1], z[i-1:i+1]
    A = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1)
    B = (y2 - y1) / (x2 - x1) - A * (x1 + x2)
    return exp(-B / 2A)
end

"Singletons S₁ / 2μN* of the mean SFS record (NaN without data)."
function singletons(cohort, which; arm, level)
    k, y = sfs_rec(cohort, which; arm, level)
    i = findfirst(==(1), k)
    return i === nothing ? NaN : y[i] / (2MU * level * NDET)
end

recs = Dict(a => records(a) for (a, _) in LEVEL_ARMS)
mod = recs["mod"]
timelabel(which) = which === :detection ? "Detection" : which === :reduced ? "Stabilisation" : "𝜏 ≈ $(which)"

# ── a: the richness transient by stable size ─────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = foldlabel("0", "SFS"))
hlines!(ax, [1.0]; LINE_UNTREATED...)
for (a, L) in LEVEL_ARMS
    g = overtime(recs[a], "neutral", "S_rel_det"); g = g[g.tau .<= T_EARLY, :]
    medianband!(ax, g.tau, g.med, g.lo, g.hi; color = COL_LEVEL[L])
    x, y = theory_time(L, :S); keep = x .<= T_EARLY
    lines!(ax, x[keep], y[keep]; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
end
xlims!(ax, -0.3, T_EARLY)
savepanel(fig, NAME, "a")

# ── b: the singleton balance ──────────────────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = "Singletons, 𝑆₁ / 2𝜇𝑁*")
hlines!(ax, [1.0]; color = :black, linestyle = :dash, linewidth = 0.8)
for (a, L) in LEVEL_ARMS
    tred = recs[a].tred[1]
    whichs = [[:reduced]; [t for t in 0.25:0.25:(T_EARLY + 1) if t <= TEND && t - tred <= T_EARLY &&
                                       any((sfs.arm .== a) .& (sfs.tag .== "fine") .& (sfs.idx .== round(Int, 100t)))]]
    ts = [w === :reduced ? 0.0 : Float64(w) - tred for w in whichs]
    ys = [singletons("neutral", w; arm = a, level = L) for w in whichs]
    scatterlines!(ax, ts, ys; color = COL_LEVEL[L], markersize = 3, linewidth = 0.6)
    t = th_tr(L); keep = t.tau .<= T_EARLY
    lines!(ax, t.tau[keep], t.S1[keep] ./ (2MU * L * NDET); color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
end
xlims!(ax, -0.3, T_EARLY)
savepanel(fig, NAME, "b")

# ── c: level dependence of the minimum ───────────────────────────────────────
"Time since stabilisation and depth of the minimum of the median ⁰D trajectory of the neutral cohort under `arm` (first 10 time units of the hold)."
function sim_minimum(arm)
    g = overtime(recs[arm], "neutral", :S); g = g[g.tau .<= T_EARLY, :]
    isempty(g) && return (tau = NaN, smin = NaN, s0 = NaN)
    i = argmin(g.med); sdet = median(recs[arm].S_det[recs[arm].cohort .== "neutral"])
    red = recs[arm][(recs[arm].tag .== "blockend") .& (recs[arm].cohort .== "neutral"), :S]
    return (tau = g.tau[i], smin = g.med[i] / sdet, s0 = median(red) / sdet)
end
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_LEVEL, ylabel = rich("Time of the minimum, 𝜏", subscript("min")))
lines!(ax, thl.level, thl.tau_min; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
for (a, L) in LEVEL_ARMS
    scatter!(ax, [L], [sim_minimum(a).tau]; color = COL_LEVEL[L], markersize = 6)
end
xlims!(ax, 0.25, 1.05)
savepanel(fig, NAME, "c_tmin")

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_LEVEL, ylabel = foldlabel("0", "SFS"))
hlines!(ax, [1.0]; LINE_UNTREATED...)
lines!(ax, thl.level, thl.Smin_rel; color = COL_PRED, linestyle = LS_PRED, linewidth = 1, label = "Minimum")
lines!(ax, thl.level, thl.S0_rel; color = COL_PRED, linestyle = :dash, linewidth = 0.8, label = "At stabilisation")
for (a, L) in LEVEL_ARMS
    m = sim_minimum(a)
    scatter!(ax, [L], [m.smin]; color = COL_LEVEL[L], markersize = 6)
    scatter!(ax, [L], [m.s0]; color = COL_LEVEL[L], marker = :utriangle, markersize = 6)
end
xlims!(ax, 0.25, 1.05); axislegend(ax; position = :rb)
savepanel(fig, NAME, "c_smin")

# ── d: SFS and compensated SFS over time (no resistance) ─────────────────────
sfscolour(i) = pale(COL_LEVEL[0.8], 0.85 * (i - 1) / max(1, length(SFS_TIMES) - 1))
"Colour `c` blended towards white by fraction `f`."
pale(c, f) = (c = Makie.to_color(c); Makie.RGBAf(c.r + f * (1 - c.r), c.g + f * (1 - c.g), c.b + f * (1 - c.b), 1))
for (p, compensated) in (("d_sfs", false), ("d", true))
    local fig = panel()
    local ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = compensated ? "𝑘(𝑘+1)·𝑆ₖ" : LABELS[:Sk], LOGX..., LOGY...)
    for (i, w) in Iterators.reverse(collect(enumerate(SFS_TIMES)))      # pale (late) first, dark on top
        k, y = sfs_rec("neutral", w; compensated); isempty(k) && continue
        lines!(ax, k, y; color = sfscolour(i), linewidth = 1.0)
    end
    compensated && for w in SFS_TIMES
        k, y = theory_sfs(0.8, w; compensated); isempty(k) && continue
        lines!(ax, k, y; color = :grey30, linestyle = :dash, linewidth = 0.6)
    end
    savepanel(fig, NAME, p)
end

# ── e: position of the front over time ───────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = "Peak of 𝑘(𝑘+1)·𝑆ₖ")
tt = range(0, TEND; length = 50)
lines!(ax, tt, B0 .* tt; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
let taus = sort(unique(ths.tau[(ths.level .== 0.8) .& (ths.kind .== "treated")]))
    lines!(ax, taus, [peakpos(theory_sfs(0.8, nothing; compensated = true, tau = τ)...) for τ in taus]; color = :grey30, linestyle = :dash, linewidth = 1)
end
fine_times = sort(unique(sfs.idx[(sfs.tag .== "fine") .& (sfs.arm .== "mod")])) ./ 100
for c in COHORTS
    ts = [t for t in fine_times if t >= 1]
    ys = [peakpos(sfs_rec(c, t; compensated = true)...) for t in ts]
    keep = isfinite.(ys) .& (ys .< 100)                                 # the peak jumps to the resistant bump after the takeover
    scatterlines!(ax, ts[keep] .- mod.tred[1], ys[keep]; color = COL_COHORT[c], markersize = 3, linewidth = 0.6)
end
xlims!(ax, -1, TEND + 1); ylims!(ax, 0, nothing)
savepanel(fig, NAME, "e")

# ── f–m: Hill numbers over time ──────────────────────────────────────────────
"y limits (log axis) of the panel of `col` over `cohorts`, including 1 and the theory."
function hill_limits(d, col, cohorts; rel, theory)
    v = [1.0]
    for c in cohorts
        g = overtime(d, c, "$(col)$(rel)"); append!(v, g.lo); append!(v, g.hi)
    end
    theory && append!(v, theory_time(0.8, col; rel)[2])
    v = filter(x -> isfinite(x) && x > 0, v)
    return (minimum(v) / 1.1, maximum(v) * 1.1)
end

"Ticks of a log axis with at least three inside [lo, hi]: 1–2–5 per decade, finer for narrow ranges."
function hillticks(lo, hi)
    for ms in ((1, 2, 5), (1, 2, 3, 5, 7), (1, 1.5, 2, 3, 4, 5, 7), (1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5, 6, 7, 8, 9))
        t = [m * 10.0^e for e in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in ms if lo <= m * 10.0^e <= hi]
        length(t) >= 3 && return t
    end
    return [lo, hi]
end

"Tick labels with two significant digits (1.2, 0.9, 20)."
hilllabels(vs) = [(r = round(v; sigdigits = 2); isinteger(r) ? string(Int(r)) : string(r)) for v in vs]

"""
Panel of the fold change of `col` (order `qlabel`, distribution `source`) relative to detection or to
the untreated tumour at N*: median and interquartile band per cohort (bands first, then the lines in
cohort order, so that the dashed weak and full lines lie on top of no and partial resistance), theory dotted.
"""
function hill_panel(d, col; cohorts, qlabel, source, rel, lims, theory)
    fig = panel()
    ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = foldlabel(qlabel, source), LOGY..., yticks = hillticks(lims...), ytickformat = hilllabels)
    hlines!(ax, [1.0]; LINE_UNTREATED...)
    theory && lines!(ax, theory_time(0.8, col; rel)...; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
    gs = [(c, overtime(d, c, "$(col)$(rel)")) for c in cohorts]
    for (c, g) in gs
        isempty(g) || band!(ax, g.tau, g.lo, g.hi; color = (COL_COHORT[c], 0.2))
    end
    for (c, g) in gs
        isempty(g) || lines!(ax, g.tau, g.med; color = COL_COHORT[c], linestyle = LS_COHORT[c])
    end
    xlims!(ax, -1, TEND + 1); ylims!(ax, lims...)
    return fig
end

for (names, rel) in ((["f", "g", "h", "i"], "_rel_det"), (["f_size", "g_size", "h_size", "i_size"], "_rel_size"))
    for (q, col, p) in zip(Q_LABELS, PASSENGER_Q, names)
        lims = hill_limits(mod, col, COHORTS; rel, theory = true)
        savepanel(hill_panel(mod, col; cohorts = COHORTS, qlabel = q, source = "SFS", rel, lims, theory = true), NAME, p)
    end
end
for (q, col, p) in zip(Q_LABELS, RESISTANT_Q, ["j", "k", "l", "m"])
    lims = hill_limits(mod, col, RES_COHORTS; rel = "_rel_det", theory = false)
    savepanel(hill_panel(mod, col; cohorts = RES_COHORTS, qlabel = q, source = "CSD", rel = "_rel_det", lims, theory = false), NAME, p)
end

# ── dose ─────────────────────────────────────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABEL_T, ylabel = LABELS[:dose], LOGY...)
for c in COHORTS
    d = counts[(counts.cohort .== c) .& (counts.arm .== "mod"), :]; nrow(d) == 0 && continue
    q = running_quantiles(d, :m; ntotal = length(unique(d.tumour)))
    medianband!(ax, q.t .- mod.tred[1], max.(q.med, 1e-3), max.(q.lo, 1e-3), max.(q.hi, 1e-3); color = COL_COHORT[c], linestyle = LS_COHORT[c])
end
xlims!(ax, -1, TEND + 1)
savepanel(fig, NAME, "dose")

# ── legends ──────────────────────────────────────────────────────────────────
fig = panel(; w = 180, h = 8)
Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)]; [LineElement(color = COL_COHORT[c], linewidth = 1.2, linestyle = LS_COHORT[c]) for c in COHORTS];
                   [LineElement(color = COL_PRED, linestyle = LS_PRED, linewidth = 1)]],
       ["Detection"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]; "Theory"];
       orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
savepanel(fig, NAME, "legend")

# stable sizes of a–c, top to bottom as the curves of a; transparent, to lie in the empty upper right of c_tmin
fig = Figure(size = (24MM_PT, 13MM_PT), backgroundcolor = :transparent)
Legend(fig[1, 1], [[LineElement(color = COL_LEVEL[L], linewidth = 1.2), MarkerElement(color = COL_LEVEL[L], marker = :circle, markersize = 5)]
                   for (_, L) in reverse(LEVEL_ARMS)],
       [LEVEL_NAME[L] for (_, L) in reverse(LEVEL_ARMS)]; labelsize = 7, patchsize = (14, 5), rowgap = 2, backgroundcolor = :transparent)
savepanel(fig, NAME, "legend_levels")

fig = panel()
Legend(fig[1, 1], [[LineElement(color = sfscolour(i), linewidth = 1.5) for i in eachindex(SFS_TIMES)]; [LineElement(color = :grey30, linestyle = :dash, linewidth = 0.8)]],
       [[timelabel(w) for w in SFS_TIMES]; "Theory"]; nbanks = 2, labelsize = 7, patchsize = (12, 5), rowgap = 4, colgap = 8)
savepanel(fig, NAME, "legend_sfs")

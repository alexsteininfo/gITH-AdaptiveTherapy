# Shared by Figure 3, Supp 8 and Supp 9 (cyclic adaptive therapy = arm `hol`): tumours that complete
# the first six cycles without progressing, Hill numbers per cycle relative to the same tumour
# untreated at the same size, the theory (theory_hill.csv, theory_sfs.csv), Hill-number panels and
# SFS curves per cycle. Include after style.jl.

const NCYCLES = 6   # holiday cycles with stored distributions (HEAVY_BLOCKS = 12 in params.jl)
const COHORTS = ["neutral", "weak", "partial", "full"]
const RES_COHORTS = ["weak", "partial", "full"]
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]
const THRESHOLDS = ("lower", "upper")
const COL_THRESHOLD = Dict("lower" => PAIRED[8], "upper" => PAIRED[2])
const THRESHOLD_NAME = Dict("lower" => "Lower (0.5𝑁)", "upper" => "Upper (𝑁)")
const THEORY_MODEL = "harmonic"     # mutations born during cycling by the harmonic-mean law

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")
theory_hill = readproc("01_dynamics", "theory_hill.csv")
theory_sfs = readproc("01_dynamics", "theory_sfs.csv")
any((snap.tag .== "growth") .& (snap.idx .== 50)) ||
    error("no untreated records at 0.5 N_det: run code/01_datageneration/run_part1_growth.sh")
const NDET = maximum(sfs.k[sfs.arm .== "detection"])

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

"Theory of Hill number `col` relative to untreated per cycle at `threshold` (upper: from detection)."
function theory_percycle(col, threshold)
    t = theory_hill[(theory_hill.model .== THEORY_MODEL) .& (theory_hill.threshold .== threshold) .&
                    ((theory_hill.kind .== "treated") .| (theory_hill.threshold .== "upper")), :]
    return sort(t, :cycle).cycle, t[sortperm(t.cycle), "$(col)_rel"]
end

"Common y limits (log axis) of the Hill-number panels of `cols` at `threshold`, including 1 and the theory."
function hill_limits(d, cols, cohorts, threshold; theory)
    v = [1.0]
    for col in cols, c in cohorts
        g = percycle(d, c, threshold, "$(col)_rel"); append!(v, g.lo); append!(v, g.hi)
    end
    theory && for col in cols
        append!(v, theory_percycle(col, threshold)[2])
    end
    v = filter(isfinite, v)
    return (minimum(v) / 1.1, maximum(v) * 1.1)
end

"""
Panel of Hill number `col` (order `qlabel`, from distribution `source`) treated / untreated per cycle
at `threshold`: median and interquartile range per cohort, theory (neutral model) dotted.
"""
function hill_panel(d, col, threshold; cohorts, qlabel, source, lims, theory)
    fig = panel()
    ylabel = rich(superscript(qlabel), "𝐷", subscript(source), " (Treated/Untreated)")
    first = threshold == "lower" ? 1 : 0                 # the lower threshold has no cycle 0
    ax = Axis(fig[1, 1]; xlabel = LABELS[:cycle], ylabel, xticks = first:NCYCLES, LOGY...)
    hlines!(ax, [1.0]; LINE_UNTREATED...)
    if theory
        x, y = theory_percycle(col, threshold)
        lines!(ax, x, y; color = COL_PRED, linestyle = LS_PRED, linewidth = 1)
    end
    for (k, c) in enumerate(cohorts)
        g = percycle(d, c, threshold, "$(col)_rel")
        isempty(g) && continue
        x = g.cycle .+ (k - (length(cohorts) + 1) / 2) * 0.15
        lines!(ax, x, g.med; color = COL_COHORT[c], linewidth = 0.6)
        errorbars!(ax, x, g.med, g.med .- g.lo, g.hi .- g.med; color = COL_COHORT[c])
        scatter!(ax, x, g.med; color = COL_COHORT[c])
    end
    xlims!(ax, first - 0.5, NCYCLES + 0.5); ylims!(ax, lims...)
    return fig
end

"""
The eight Hill-number panels of `threshold` (passengers q = 0, 1, 2, ∞ with theory, then resistant
clones), saved as `names[1:8]` of figure `name`; y axes shared within each group of four.
"""
function hill_panels(d, threshold, name, names)
    groups = ((PASSENGER_Q, COHORTS, "SFS", true), (RESISTANT_Q, RES_COHORTS, "CSD", false))
    i = 0
    for (cols, cohorts, source, theory) in groups
        lims = hill_limits(d, cols, cohorts, threshold; theory)
        for (q, col) in zip(Q_LABELS, cols)
            i += 1
            savepanel(hill_panel(d, col, threshold; cohorts, qlabel = q, source, lims, theory), name, names[i])
        end
    end
end

"Legend of the Hill-number panels: cohorts, untreated reference and theory."
function hill_legend(name)
    fig = panel(; w = 180, h = 8)
    Legend(fig[1, 1], [[LineElement(; LINE_UNTREATED...)];
                       [[LineElement(color = COL_COHORT[c], linewidth = 0.6), MarkerElement(color = COL_COHORT[c], marker = :circle, markersize = 4)]
                        for c in COHORTS];
                       [LineElement(color = COL_PRED, linestyle = LS_PRED, linewidth = 1)]],
           ["Untreated, same size"; [c == "neutral" ? COHORT_NAME[c] : "$(COHORT_NAME[c]) resistance" for c in COHORTS]; "Theory"];
           orientation = :horizontal, labelsize = 7, patchsize = (18, 5))
    savepanel(fig, name, "legend")
end

# ── SFS per cycle ─────────────────────────────────────────────────────────────
"Colour `c` blended towards white by fraction `f`."
pale(c, f) = (c = Makie.to_color(c); Makie.RGBAf(c.r + f * (1 - c.r), c.g + f * (1 - c.g), c.b + f * (1 - c.b), 1))

"Shade of `threshold`'s colour for `cycle`: dark at the first record, paler with every cycle."
cyclecolour(threshold, cycle) = pale(COL_THRESHOLD[threshold], 0.85cycle / NCYCLES)

"""
Mean SFS of `cohort` at `threshold` of `cycle` (cycle 0, upper: detection), up to the threshold
size; `compensated`: k(k+1)·S_k without the clonal bin. Returns `(k, y)`.
"""
function sfs_curve(cohort, threshold, cycle; compensated)
    s = sfs[(sfs.cohort .== cohort) .& (sfs.S_mean .> 0), :]
    r = cycle == 0 ? s[s.arm .== "detection", :] :
        s[(s.arm .== "hol") .& (s.tag .== "blockend") .& (s.idx .== 2cycle - (threshold == "lower")), :]
    kmax = (threshold == "lower" ? 0.5 : 1.0) * NDET
    r = r[compensated ? (r.k .< 0.9kmax) : (r.k .<= kmax), :]
    return r.k, compensated ? r.k .* (r.k .+ 1) .* r.S_mean : r.S_mean
end

"Theory of the mean SFS at `threshold` of `cycle` (same bins and cuts as `sfs_curve`)."
function theory_curve(threshold, cycle; compensated)
    s = theory_sfs[(theory_sfs.model .== THEORY_MODEL) .& (theory_sfs.threshold .== threshold) .&
                   (theory_sfs.cycle .== cycle) .& (theory_sfs.kind .== (cycle == 0 ? "untreated" : "treated")) .&
                   (theory_sfs.S_mean .> 0), :]
    kmax = (threshold == "lower" ? 0.5 : 1.0) * NDET
    s = s[compensated ? (s.k .< 0.9kmax) : (s.k .<= kmax), :]
    return s.k, compensated ? s.k .* (s.k .+ 1) .* s.S_mean : s.S_mean
end

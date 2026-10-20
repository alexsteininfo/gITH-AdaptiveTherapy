# Supplementary Figure 5 — how well does a sample measure the tumour? Virtual-trial patients at
# detection: statistics of uniform samples of n cells (mutations in fewer than f·n cells removed)
# against the same statistics of the entire tumour.
#   a  Spearman ρ between sample and entire-tumour values across patients vs n, f = 1% (default)
#   b  same, f = 0 (all mutations)
#   c  same, f = 5%
#   d  median relative error of the exponent γ̂ vs n, per f (γ is scale-free; S, ²D, … are not)
# Data: 04_trial/samples.csv, snapshots.csv
# Run : julia --project=. code/03_figures/supp_5/supp_5.jl
include(joinpath(@__DIR__, "..", "style.jl"))
using StatsBase: corspearman

const NAME = "supp_5"
const STATS = ["S", "D2n", "H1n", "gamma", "mean_burden", "var_burden"]
const STAT_NAME = Dict("S" => "S", "D2n" => "²D", "H1n" => "H₁", "gamma" => "γ", "mean_burden" => "Mean burden",
                       "var_burden" => "Burden variance")
const COL_STAT = Dict(zip(STATS, PAIRED[2:2:12]))

sm = readproc("04_trial", "samples.csv")
snap = readproc("04_trial", "snapshots.csv")
det = snap[snap.tag .== "detection", :]
sm = sm[sm.tag .== "detection", :]

"Sample and entire-tumour value of `y` per patient, for sample size n and threshold f."
function paired(y, n, f)
    d = sm[(sm.n .== n) .& (sm.fmin .== f), [:tumour, Symbol(y)]]
    j = innerjoin(rename(d, Symbol(y) => :sample), select(det, :tumour, Symbol(y) => :tumour_value); on = :tumour)
    return j[isfinite.(j.sample) .& isfinite.(j.tumour_value), :]
end

ns = sort(unique(sm.n))
for (p, f) in (("a", 0.01), ("b", 0.0), ("c", 0.05))
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:n], ylabel = LABELS[:rho], LOGX...)
    for y in STATS
        r = [(j = paired(y, n, f); nrow(j) >= 3 ? corspearman(j.sample, j.tumour_value) : NaN) for n in ns]
        lines!(ax, ns, r; color = COL_STAT[y], label = STAT_NAME[y])
        scatter!(ax, ns, r; color = COL_STAT[y])
    end
    p == "a" && toplegend(fig, ax; nbanks = 2)
    savepanel(fig, NAME, p)
end

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:n], ylabel = "Relative error of γ", LOGX...)
for (f, col) in zip(sort(unique(sm.fmin)), shades(:black, length(unique(sm.fmin))))
    e = [(j = paired("gamma", n, f); isempty(j.sample) ? NaN : median(abs.(j.sample .- j.tumour_value) ./ j.tumour_value)) for n in ns]
    lines!(ax, ns, e; color = col, label = "f = $(round(Int, 100f))%")
    scatter!(ax, ns, e; color = col)
end
toplegend(fig, ax)
savepanel(fig, NAME, "d")

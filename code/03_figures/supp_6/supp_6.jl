# Supplementary Figure 6 — power-law fit of the SFS tail (discrete MLE, Clauset et al. 2009) over the
# lower cut-off k_min, at detection (neutral and full resistance) and after the last stored holiday
# cycle (neutral). Median and interquartile range over tumours.
#   a  exponent γ̂ vs k_min / N; dotted: the median k_min selected by the Kolmogorov–Smirnov distance
#   b  standard error of γ̂, (γ̂ − 1)/√n_tail
#   c  KS-selected k_min / N over tumours
# Data: 01_dynamics/powerlaw.csv
# Run : julia --project=. code/03_figures/supp_6/supp_6.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_6"

pl = readproc("01_dynamics", "powerlaw.csv")
const GROUPS = [("neutral", "detection", PAIRED[2], "Neutral"),
                ("neutral", "blockend", PAIRED[1], "After holidays"),
                ("full", "detection", PAIRED[6], "Full resistance")]
group(c, tag) = pl[(pl.cohort .== c) .& (pl.tag .== tag), :]
"KS-selected k_min / N, one value per tumour."
selected(d) = (u = unique(d[:, [:tumour, :ks_kmin, :N]]); u.ks_kmin ./ u.N)

for (p, y, ylabel, log) in (("a", :gamma, LABELS[:gamma], false), ("b", :se, "Standard error of γ", true))
    fig = panel()
    ax = log ? Axis(fig[1, 1]; xlabel = "Lower cut-off, kₘᵢₙ / N", ylabel, LOGX..., LOGY...) :
               Axis(fig[1, 1]; xlabel = "Lower cut-off, kₘᵢₙ / N", ylabel, LOGX...)
    for (c, tag, col, label) in GROUPS
        d = group(c, tag); nrow(d) == 0 && continue
        g = sort(combine(groupby(d[isfinite.(d[!, y]), :], :kmin_frac), y => median => :med, y => q25 => :lo,
                         y => q75 => :hi), :kmin_frac)
        medianband!(ax, g.kmin_frac, g.med, g.lo, g.hi; color = col, label)
        p == "a" && vlines!(ax, [median(selected(d))]; color = col, linestyle = LS_PRED, linewidth = 0.8)
    end
    p == "a" && toplegend(fig, ax; nbanks = 2)
    savepanel(fig, NAME, p)
end

fig = panel(); ax = Axis(fig[1, 1]; xlabel = "Selected kₘᵢₙ / N", ylabel = "Tumours", LOGX...)
for (c, tag, col, _) in GROUPS
    x = selected(group(c, tag)); isempty(x) && continue
    edges = 10 .^ range(-4.0625, -0.9375; length = 26)
    hist!(ax, x; bins = edges, color = (col, 0.5))
end
savepanel(fig, NAME, "c")

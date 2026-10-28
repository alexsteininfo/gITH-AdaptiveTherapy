# Shared style of all figure panels: theme, colours, labels, panel size, reading tables, saving
# panels and summaries over tumours. Figure scripts read only data/processed/. With FIG_DATA=test
# they read the TEST=1 tables (data/processed/00_test_<block>) and write to figures/00_test/.

using CairoMakie, CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "..", "X_helpers", "paths.jl"))

const TEST_DATA = get(ENV, "FIG_DATA", "") == "test"

# ── Size and theme ────────────────────────────────────────────────────────────
const MM_PT = 72 / 25.4             # CairoMakie: 1 unit = 1 pt
const PANEL_MM = (45, 45)
const PX_PER_UNIT = 300 / 72        # 300 dpi

set_theme!(Theme(
    fontsize = 7, figure_padding = 4, rowgap = 3, colgap = 3,
    Axis = (xlabelsize = 8, ylabelsize = 8, xticklabelsize = 7, yticklabelsize = 7,
            xgridvisible = false, ygridvisible = false, xtickalign = 1, ytickalign = 1,
            xticksize = 3, yticksize = 3, spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8,
            xminorticksvisible = false, yminorticksvisible = false, xlabelpadding = 2, ylabelpadding = 2),
    Legend = (framevisible = false, labelsize = 6, patchsize = (8, 5), patchlabelgap = 2, rowgap = 0, colgap = 4,
              padding = (2, 2, 2, 2), margin = (0, 0, 0, 0)),
    Lines = (linewidth = 1.2,), Scatter = (markersize = 4,), Errorbars = (whiskerwidth = 3, linewidth = 0.8),
    Colorbar = (labelsize = 7, ticklabelsize = 6, size = 5, tickalign = 1, ticksize = 2)))

"Empty figure of `w × h` mm."
panel(; w = PANEL_MM[1], h = PANEL_MM[2]) = Figure(size = (w * MM_PT, h * MM_PT))

"Save panel `p` of figure `name` as figures/<name>/<p>.png (300 dpi); returns the path."
function savepanel(fig, name::AbstractString, p::AbstractString)
    path = joinpath(figdir(TEST_DATA ? joinpath("00_test", name) : name), "$p.png")
    save(path, fig; px_per_unit = PX_PER_UNIT)
    return path
end

"Processed table `file` of `block` (its 00_test_ version with FIG_DATA=test)."
readproc(block::AbstractString, file::AbstractString) =
    CSV.read(joinpath(ROOT, "data", "processed", TEST_DATA ? "00_test_$block" : block, file), DataFrame)

# ── Colours ───────────────────────────────────────────────────────────────────
const PAIRED = Makie.to_colormap(:Paired_12)    # ColorBrewer Paired, as in Stein & Werner (2025)

# trajectories (BirthDeathMutation_AT Figure 1)
const COL_TOTAL = :black
const COL_SENSITIVE = (:royalblue, 0.9)
const COL_RESISTANT = (:firebrick, 0.9)
const COL_SUBCLONE = (:firebrick, 0.2)
const COL_SHADE = :grey50                       # treatment shading, opacity 0.3·m

# schedules
const COL_ARM = Dict("single" => PAIRED[8], "mtd" => PAIRED[6], "hol" => PAIRED[2], "mod" => PAIRED[4],
                     "met" => PAIRED[10], "untreated" => :grey60)
const ARM_NAME = Dict("single" => "Single round", "mtd" => "MTD", "hol" => "Holidays", "mod" => "Modulation",
                      "met" => "Metronomic", "untreated" => "Untreated")

# simulations vs theory
const COL_SIM = PAIRED[9]       # single simulation
const COL_MEAN = PAIRED[2]      # average of simulations
const COL_PRED = :black         # prediction, dotted
const LS_PRED = :dot

# resistance cohorts (none → full; checked for colour-vision deficiency); untreated in COL_ARM
const COL_COHORT = Dict(c => Makie.to_color(h) for (c, h) in ("neutral" => "#e08214", "weak" => "#d7191c",
                                                               "partial" => "#7b3294", "full" => "#2c7bb6"))
const COHORT_NAME = Dict("neutral" => "No resistance", "weak" => "Weak", "partial" => "Partial", "full" => "Full")
# SFS lines per cohort, drawn in this order: thick to thin, so that overlapping curves stay visible
const LINE_UNTREATED = (color = (COL_ARM["untreated"], 0.6), linewidth = 3.0)
const LINE_COHORT = Dict("neutral" => (linestyle = :solid, linewidth = 1.8), "weak" => (linestyle = :dash, linewidth = 1.4),
                         "partial" => (linestyle = :dot, linewidth = 1.2), "full" => (linestyle = :dashdot, linewidth = 1.0))

# burden split (Stein & Werner 2025: pre-existing orange, new green, total blue)
const COL_PRE = PAIRED[8]
const COL_NEW = PAIRED[4]
const COL_ALL = PAIRED[2]

"`n` shades of colour `c`, light to dark (opacity 0.35–1)."
shades(c, n) = n == 1 ? [(c, 1.0)] : [(c, a) for a in range(0.35, 1.0; length = n)]

# ── Labels ────────────────────────────────────────────────────────────────────
const LABELS = Dict(
    :time => "Time", :cells => "Number of cells", :k => "Number of cells, 𝑘",
    :Sk => "Number of mutations, 𝑆ₖ", :S => "Number of mutations, 𝑆", :B => "Total mutational burden, 𝘉",
    :clones => "Number of clones", :resclones => "Resistant clones",
    :D => rich("Diversity, ", superscript("𝑞"), "𝐷"), :D2 => rich("Diversity, ", superscript("2"), "𝐷"),
    :H1 => "Entropy, 𝐻₁",
    :ttp => "Time to progression", :gain => "Gain in time to progression", :dose => "Dose",
    :freq => "Frequency", :resfrac => "Resistant fraction", :largest => "Largest resistant clone",
    :gamma => "Exponent, γ", :cycle => "Cycle", :rho => "Correlation, ρ", :absrho => "|Correlation|, |ρ|",
    :n => "Sample size", :density => "Density", :pfs => "Fraction progression-free", :q => "Order, 𝑞",
    :ratio => "Treated / untreated", :level => "Containment level", :burden => "Mean burden",
    :relative => "Relative to detection")

const Q_LABELS = ["0", "1", "2", "∞"]

# ── Log axes ──────────────────────────────────────────────────────────────────
"Ticks of a log axis: powers of ten over wide ranges, 1–3 or 1–2–5 steps over narrow ones."
function logticks(lo, hi)
    a, b = floor(Int, log10(lo)), ceil(Int, log10(hi))
    for mult in ([1.0], [1.0, 3.0], [1.0, 2.0, 5.0], 1.0:9.0)
        t = [m * 10.0^e for e in a:b for m in mult if lo <= m * 10.0^e <= hi]
        length(t) >= 3 && return t
    end
    return [lo, hi]
end

"One tick label: 10ᵏ, m×10ᵏ from 10⁴ on, a plain number otherwise."
function ticklabel(v)
    e = floor(Int, log10(v) + 1e-9); m = round(v / 10.0^e; sigdigits = 2)
    (v < 1e4 && v >= 0.01) && return string(v >= 1 ? round(Int, v) : round(v; sigdigits = 2))
    m == 1 && return rich("10", superscript(string(e)))
    return rich(string(isinteger(m) ? Int(m) : m), "×10", superscript(string(e)))
end

"Labels of a log axis: all as 10ᵏ when every tick is a power of ten over ≥ 3 decades."
function logformat(vs)
    pow = length(vs) >= 3 && all(v -> isapprox(v, 10.0^round(log10(v))), vs) && log10(vs[end] / vs[1]) >= 3
    return pow ? [rich("10", superscript(string(round(Int, log10(v))))) for v in vs] : [ticklabel(v) for v in vs]
end

const LOGX = (xscale = log10, xticks = logticks, xtickformat = logformat)
const LOGY = (yscale = log10, yticks = logticks, ytickformat = logformat)

"Legend above an axis (horizontal), for panels whose data fill the axis."
toplegend(fig, ax; kw...) = Legend(fig[0, 1], ax; orientation = :horizontal, tellwidth = false, tellheight = true, kw...)

# ── Summaries over tumours ────────────────────────────────────────────────────
q25(v) = quantile(v, 0.25)
q75(v) = quantile(v, 0.75)

"Restricted time to progression of an outcome row: progression time, else the arm's follow-up."
restricted_ttp(r) = r.outcome == "progressed" ? Float64(r.ttp) : Float64(r.tmax)

"Kind of an arm id (`hol_0.6` → `hol`)."
armkind(a::AbstractString) = String(first(split(a, "_")))

"""
    running_quantiles(df, y; ntotal, t = :t_rel)

Median and quartiles of column `y` per time point over tumours, while at least half of the
`ntotal` tumours are still running (have a row at that time; NaN values count as running but are
left out of the quantiles). Returns `t, med, lo, hi, n`.
"""
function running_quantiles(df::AbstractDataFrame, y; ntotal::Integer, t = :t_rel)
    d = DataFrame(t = round.(Float64.(df[!, t]); digits = 2), y = Float64.(coalesce.(df[!, y], NaN)))
    g = combine(groupby(d, :t), :y => (v -> (f = filter(isfinite, v); isempty(f) ? NaN : median(f))) => :med,
                :y => (v -> (f = filter(isfinite, v); isempty(f) ? NaN : q25(f))) => :lo,
                :y => (v -> (f = filter(isfinite, v); isempty(f) ? NaN : q75(f))) => :hi, nrow => :n)
    sort!(g, :t)
    stop = findfirst(<(ntotal / 2), g.n)
    return stop === nothing ? g : g[1:stop - 1, :]
end

"Median line with an interquartile band."
function medianband!(ax, x, med, lo, hi; color, linestyle = :solid, label = nothing, alpha = 0.2)
    ok = isfinite.(med) .& isfinite.(lo) .& isfinite.(hi)
    any(ok) || return nothing
    band!(ax, x[ok], lo[ok], hi[ok]; color = (Makie.to_color(color), alpha))
    lines!(ax, x[ok], med[ok]; color, linestyle, label)
end

"Kaplan–Meier step curve from a table with columns `t, S`."
kmcurve!(ax, d; kw...) = stairs!(ax, d.t, d.S; step = :post, kw...)

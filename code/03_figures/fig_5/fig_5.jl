# Figure 5 — the competition–diversity trade-off. Each schedule is scanned over its containment
# level (holidays: lower threshold 0.5–0.9 N_det; modulation: stable size 0.5–1.0 N_det); points
# get darker with the level; MTD is the single reference point. Time to progression = restricted
# mean up to the follow-up (Kaplan–Meier).
#   a  time to progression vs resistant clones at the end of follow-up (median, IQR), full resistance
#   b  same, partial resistance
#   c  time to progression vs total mutational burden B at the end (∝ risk of resistance to a
#      second drug), full resistance
#   d_hol, d_mod  paired gain in time to progression vs MTD (median, 95% CI) by level, for
#      N_det/K = 0.5, 0.7, 0.8, 0.9
#   e_hol, e_mod  same for ν = 10⁻⁵, 10⁻⁴, 10⁻³
#   legend  schedules
# Data: 02_tradeoff/tradeoff.csv, 03_sensitivity/tradeoff.csv
# Run : julia --project=. code/03_figures/fig_5/fig_5.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_5"
const KINDS = ["hol", "mod"]

tr = vcat(readproc("02_tradeoff", "tradeoff.csv"), readproc("03_sensitivity", "tradeoff.csv"); cols = :union)

# ── a–c: trade-off curves ─────────────────────────────────────────────────────
function tradeoff_panel(cohort, y, lo, hi; ylabel)
    d = tr[tr.cohort .== cohort, :]
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:ttp], ylabel)
    for kind in KINDS
        g = sort(d[d.kind .== kind, :], :level); nrow(g) == 0 && continue
        cols = shades(COL_ARM[kind], nrow(g))
        lines!(ax, g.ttp_rmst, g[!, y]; color = COL_ARM[kind])
        errorbars!(ax, g.ttp_rmst, g[!, y], g[!, y] .- g[!, lo], g[!, hi] .- g[!, y]; color = cols)
        scatter!(ax, g.ttp_rmst, g[!, y]; color = cols, strokecolor = COL_ARM[kind], strokewidth = 0.5)
    end
    m = d[d.kind .== "mtd", :]
    if nrow(m) > 0
        errorbars!(ax, m.ttp_rmst, m[!, y], m[!, y] .- m[!, lo], m[!, hi] .- m[!, y]; color = COL_ARM["mtd"])
        scatter!(ax, m.ttp_rmst, m[!, y]; color = COL_ARM["mtd"], marker = :diamond, markersize = 6)
    end
    return fig
end

savepanel(tradeoff_panel("full", :resD0_med, :resD0_q25, :resD0_q75; ylabel = LABELS[:resclones]), NAME, "a")
savepanel(tradeoff_panel("partial", :resD0_med, :resD0_q25, :resD0_q75; ylabel = LABELS[:resclones]), NAME, "b")
savepanel(tradeoff_panel("full", :B_med, :B_q25, :B_q75; ylabel = LABELS[:B]), NAME, "c")

# ── d, e: gain vs level across competition and resistance rates ──────────────
function gain_panel(kind, cohorts, labels)
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:level], ylabel = LABELS[:gain])
    hlines!(ax, [0.0]; color = :grey70, linestyle = :dot, linewidth = 0.8)
    cols = shades(COL_ARM[kind], length(cohorts))
    for (k, (c, label)) in enumerate(zip(cohorts, labels))
        g = sort(tr[(tr.cohort .== c) .& (tr.kind .== kind), :], :level); nrow(g) == 0 && continue
        x = g.level .+ (k - (length(cohorts) + 1) / 2) * 0.008
        lines!(ax, x, g.gain_med; color = cols[k], label)
        errorbars!(ax, x, g.gain_med, g.gain_med .- g.gain_lo, g.gain_hi .- g.gain_med; color = cols[k])
        scatter!(ax, x, g.gain_med; color = cols[k])
    end
    toplegend(fig, ax; nbanks = 2)
    return fig
end

const FK_COHORTS = (["fK0.5", "fK0.7", "full", "fK0.9"], [rich("N/K = ", s) for s in ("0.5", "0.7", "0.8", "0.9")])
const NU_COHORTS = (["nu1e-5", "full", "nu1e-3"], [rich("ν = 10", superscript(e)) for e in ("−5", "−4", "−3")])

for kind in KINDS
    savepanel(gain_panel(kind, FK_COHORTS...), NAME, "d_$kind")
    savepanel(gain_panel(kind, NU_COHORTS...), NAME, "e_$kind")
end

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 90, h = 8)
Legend(fig[1, 1], [MarkerElement(marker = :diamond, color = COL_ARM["mtd"], markersize = 6);
                   [[LineElement(color = COL_ARM[k]), MarkerElement(marker = :circle, color = COL_ARM[k], markersize = 4)]
                    for k in KINDS]],
       [ARM_NAME["mtd"]; [ARM_NAME[k] for k in KINDS]]; orientation = :horizontal, labelsize = 7)
savepanel(fig, NAME, "legend")

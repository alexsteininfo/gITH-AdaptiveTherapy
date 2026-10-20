# Supplementary Figure 2 — metronomic schedule (full dose for 1 of every 3 time units), on the same
# tumours as Part 1 (block 05 = cohorts of block 01).
#   a  mean mutational burden per cell over time, neutral cohort: metronomic vs the Part 1 arms
#      (distributions are stored only at detection and the end for metronomic, counts every 2)
#   b  resistant clones over time, full resistance
#   c  Kaplan–Meier curves, full resistance
#   d  restricted time to progression of each tumour: metronomic vs MTD (diagonal: equal)
# Data: 05_metronomic/counts.csv.gz, km.csv, outcomes.csv; 01_dynamics/counts.csv.gz, km.csv, outcomes.csv
# Run : julia --project=. code/03_figures/supp_2/supp_2.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_2"

counts = vcat(readproc("01_dynamics", "counts.csv.gz"), readproc("05_metronomic", "counts.csv.gz"); cols = :union)
km = vcat(readproc("01_dynamics", "km.csv"), readproc("05_metronomic", "km.csv"); cols = :union)
out = vcat(readproc("01_dynamics", "outcomes.csv"), readproc("05_metronomic", "outcomes.csv"); cols = :union)

function time_panel(cohort, y, arms; ylabel, log = false)
    fig = panel()
    ax = log ? Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel, LOGY...) : Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel)
    for arm in arms
        d = counts[(counts.cohort .== cohort) .& (counts.arm .== arm), :]; nrow(d) == 0 && continue
        q = running_quantiles(d, y; ntotal = length(unique(d.tumour)))
        med = log ? [x > 0 ? x : NaN for x in q.med] : q.med
        lines!(ax, q.t, med; color = COL_ARM[arm], label = ARM_NAME[arm])
    end
    return fig, ax
end

fig, ax = time_panel("neutral", :mean_burden, ["single", "hol", "mod", "met"]; ylabel = LABELS[:burden])
toplegend(fig, ax; nbanks = 2)
savepanel(fig, NAME, "a")

fig, ax = time_panel("full", :n_res, ["mtd", "hol", "mod", "met"]; ylabel = LABELS[:resclones], log = true)
toplegend(fig, ax; nbanks = 2)
savepanel(fig, NAME, "b")

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:pfs])
for arm in ["mtd", "hol", "mod", "met"]
    d = km[(km.cohort .== "full") .& (km.arm .== arm), :]; nrow(d) == 0 && continue
    kmcurve!(ax, d; color = COL_ARM[arm])
end
ylims!(ax, -0.02, 1.02)
savepanel(fig, NAME, "c")

arm_ttp(arm) = select(out[(out.cohort .== "full") .& (out.arm .== arm), :], :tumour, AsTable(:) => ByRow(restricted_ttp) => arm)
j = innerjoin(arm_ttp("mtd"), arm_ttp("met"); on = :tumour)
fig = panel(); ax = Axis(fig[1, 1]; xlabel = "TTP, MTD", ylabel = "TTP, metronomic")
hi = maximum([j.mtd; j.met; 1.0])
lines!(ax, [0, hi], [0, hi]; color = :grey70, linestyle = :dot, linewidth = 0.8)
scatter!(ax, j.mtd, j.met; color = (COL_ARM["met"], 0.6), markersize = 3)
savepanel(fig, NAME, "d")

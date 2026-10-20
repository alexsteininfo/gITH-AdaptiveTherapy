# Figure 4 — with resistance: composition under the schedules (cohorts `full` and `partial`; one
# panel per cohort, suffix _full / _partial). Median and interquartile range over tumours, shown
# while at least half of the tumours are still running.
#   a  number of resistant clones
#   b  resistant lineages (resistance arisen in a sensitive cell) arisen before (solid) and after
#      detection (dashed, de novo); medians
#   c  resistant-clone diversity ²D (entire tree, every 10 time units)
#   d  resistant fraction of the tumour
#   e  Kaplan–Meier curves of progression
#   legend  schedules and line styles
# Data: 01_dynamics/counts.csv.gz, snapshots.csv, km.csv
# Run : julia --project=. code/03_figures/fig_4/fig_4.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_4"
const ARMS = ["mtd", "hol", "mod"]
const COHORTS = ["full", "partial"]

counts = readproc("01_dynamics", "counts.csv.gz")
snap = readproc("01_dynamics", "snapshots.csv")
km = readproc("01_dynamics", "km.csv")
counts.resfrac = counts.nR ./ max.(counts.N, 1)

rows(df, cohort, arm) = df[(df.cohort .== cohort) .& (df.arm .== arm), :]

positive(v) = [x > 0 ? Float64(x) : NaN for x in v]

function time_panel(cohort, y; ylabel, log = false, data = counts)
    fig = panel()
    ax = log ? Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel, LOGY...) : Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel)
    for arm in ARMS
        d = rows(data, cohort, arm); nrow(d) == 0 && continue
        q = running_quantiles(d, y; ntotal = length(unique(d.tumour)))
        f = log ? positive : identity
        medianband!(ax, q.t, f(q.med), f(q.lo), f(q.hi); color = COL_ARM[arm])
    end
    return fig
end

for c in COHORTS
    any(counts.cohort .== c) || continue
    savepanel(time_panel(c, :n_res; ylabel = LABELS[:resclones], log = true), NAME, "a_$c")

    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:resclones], LOGY...)
    for arm in ARMS, (y, ls) in ((:res_pre, :solid), (:res_new, :dash))
        d = rows(counts, c, arm); nrow(d) == 0 && continue
        q = running_quantiles(d, y; ntotal = length(unique(d.tumour)))
        lines!(ax, q.t, positive(q.med); color = COL_ARM[arm], linestyle = ls)
    end
    savepanel(fig, NAME, "b_$c")

    # ²D of resistant clones: detection and every 10 time units
    det = snap[(snap.cohort .== c) .& (snap.tag .== "detection"), :]
    grid = vcat([vcat(transform(det, :arm => ByRow(_ -> a) => :arm),
                      snap[(snap.cohort .== c) .& (snap.arm .== a) .& (snap.tag .== "grid"), :]) for a in ARMS]...)
    grid.t_rel = [r.tag == "grid" ? Float64(r.idx) : r.t_rel for r in eachrow(grid)]    # nominal grid times
    savepanel(time_panel(c, :res_D2; ylabel = LABELS[:D2], data = grid), NAME, "c_$c")

    savepanel(time_panel(c, :resfrac; ylabel = LABELS[:resfrac]), NAME, "d_$c")

    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:pfs])
    for arm in ARMS
        d = km[(km.cohort .== c) .& (km.arm .== arm), :]; nrow(d) == 0 && continue
        kmcurve!(ax, d; color = COL_ARM[arm])
    end
    ylims!(ax, -0.02, 1.02)
    savepanel(fig, NAME, "e_$c")
end

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 120, h = 8)
Legend(fig[1, 1], [[LineElement(color = COL_ARM[a]) for a in ARMS]; LineElement(color = :black);
                   LineElement(color = :black, linestyle = :dash)],
       [[ARM_NAME[a] for a in ARMS]; "Pre-existing"; "New"]; orientation = :horizontal, labelsize = 7)
savepanel(fig, NAME, "legend")

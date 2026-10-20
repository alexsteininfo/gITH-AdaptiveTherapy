# Figure 1 — model and schedules.
#   schematic  model schematic (panels a–c; Figure 1 of BirthDeathMutation_AT.jl, copied as .ai/.png)
#   d–g        example tumour of cohort `full` under MTD, metronomic, holidays and modulation:
#              total, sensitive and resistant cells, resistant clones (thin), dose (grey ∝ m)
#   legend     legend of d–g
# Data: 01_dynamics/example_traj.csv, example_clones.csv
# Run : julia --project=. code/03_figures/fig_1/fig_1.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_1"
const SCHEMATIC = joinpath(ROOT, "..", "BirthDeatMutation_AT.jl", "plots", "adobeillustrator")

# ── a–c: schematic ────────────────────────────────────────────────────────────
for ext in ("ai", "png")
    src = joinpath(SCHEMATIC, "Figure_1.$ext")
    if isfile(src)
        cp(src, joinpath(figdir(TEST_DATA ? joinpath("00_test", NAME) : NAME), "schematic.$ext"); force = true)
    else
        @warn "schematic not found: $src"
    end
end

# ── d–g: example trajectories ─────────────────────────────────────────────────
traj = readproc("01_dynamics", "example_traj.csv")
clones = readproc("01_dynamics", "example_clones.csv")

"Intervals (t₁, t₂, m) of constant dose m > 0 along a trajectory."
function dose_intervals(t, m)
    out = Tuple{Float64,Float64,Float64}[]
    i = 1
    while i < length(t)
        j = i
        while j < length(t) - 1 && m[j + 1] == m[i]; j += 1; end
        m[i] > 0 && push!(out, (t[i], t[j + 1], m[i]))
        i = j + 1
    end
    return out
end

function trajectory_panel(arm)
    tr = sort(traj[traj.arm .== arm, :], :t_rel)
    cl = clones[clones.arm .== arm, :]
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:cells], LOGY...)
    iv = dose_intervals(tr.t_rel, tr.m)
    isempty(iv) || vspan!(ax, first.(iv), getindex.(iv, 2); color = [(COL_SHADE, 0.3m) for (_, _, m) in iv])
    for g in groupby(sort(cl, :t_rel), :clone)
        t = collect(g.t_rel); s = Float64.(g.size)
        t[1] > 0 && (pushfirst!(t, t[1] - 1); pushfirst!(s, 1.0))   # born between two records
        lines!(ax, t, s; color = COL_SUBCLONE, linewidth = 0.5)
    end
    gone(v) = [x == 0 ? NaN : Float64(x) for x in v]                 # lines end at extinction
    lines!(ax, tr.t_rel, gone(tr.nR); color = COL_RESISTANT)
    lines!(ax, tr.t_rel, gone(tr.nS); color = COL_SENSITIVE)
    lines!(ax, tr.t_rel, tr.N; color = COL_TOTAL, linewidth = 1.6)
    xlims!(ax, 0, maximum(tr.t_rel)); ylims!(ax, 1, 3 * maximum(tr.N))
    return fig
end

for (p, arm) in zip(("d", "e", "f", "g"), ("mtd", "met", "hol", "mod"))
    any(traj.arm .== arm) ? savepanel(trajectory_panel(arm), NAME, p) : @warn "no example trajectory for $arm"
end

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 120, h = 8)
Legend(fig[1, 1],
       [LineElement(color = COL_TOTAL, linewidth = 1.6), LineElement(color = COL_SENSITIVE),
        LineElement(color = COL_RESISTANT), LineElement(color = (:firebrick, 0.4), linewidth = 0.8),
        PolyElement(color = (COL_SHADE, 0.3))],
       ["Total", "Sensitive", "Resistant", "Resistant clone", "Treatment"];
       orientation = :horizontal, labelsize = 7)
savepanel(fig, NAME, "legend")

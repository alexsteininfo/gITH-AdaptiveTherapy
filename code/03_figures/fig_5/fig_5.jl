# Figure 5 — virtual trial: adaptive therapy trades time to progression for genetic diversity.
# 500 virtual patients (Latin hypercube over ν, α and N_det/K); the tumour of every patient is
# treated with MTD, holidays and modulation (paired arms, each patient their own control).
#   a     trial design: schematic, drawn in Illustrator (not made here)
#   b     Kaplan–Meier curves per arm; median time to progression in the legend
#   c–f   per patient: gain in time to progression (AT − MTD) vs fold change of the passenger
#         diversity ᵠD_SFS at progression (AT / MTD; q = 0, 1, 2, ∞), holidays and modulation;
#         large markers: medians
#   g–j   same for the diversity of resistant clones ᵠD_CSD
#   k     partial rank correlation (PRCC) of the parameters with time to progression and gain
#   l_<arm>     PRCC of the parameters with ᵠD at progression under MTD, holidays, modulation
#   l_fc_<arm>  PRCC with the fold change of ᵠD (holidays / modulation vs MTD)
#   legend  arms
# Data: 04_trial/trial.csv, trial_km.csv, trial_prcc.csv
# Run : julia --project=. code/03_figures/fig_5/fig_5.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_5"
const ARMS = ["mtd", "hol", "mod"]
const AT = ["hol", "mod"]
const SFS_STATS = ["S", "D1n", "D2n", "Dinfn"]
const CSD_STATS = ["res_D0", "res_D1", "res_D2", "res_Dinf"]
const PARAMS = ["log10_nu", "alpha", "fK"]
const PARAM_NAME = Dict("log10_nu" => "ν", "alpha" => "α", "fK" => rich("𝑁", subscript("det"), "/𝐾"))

pat = readproc("04_trial", "trial.csv")
km = readproc("04_trial", "trial_km.csv")
prcc = readproc("04_trial", "trial_prcc.csv")

hill_label(q, source) = rich(superscript(q), "𝐷", subscript(source))

# ── b: Kaplan–Meier ───────────────────────────────────────────────────────────
"Median of a Kaplan–Meier table (first time with S ≤ ½); NaN if never reached."
km_med(d) = (i = findfirst(<=(0.5), d.S); i === nothing ? NaN : d.t[i])

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:pfs])
for a in ARMS
    d = sort(km[km.arm .== a, :], :t); nrow(d) == 0 && continue
    kmcurve!(ax, d; color = COL_ARM[a], label = "$(ARM_NAME[a]) ($(round(km_med(d); digits = 1)))")
end
xlims!(ax, 0, maximum(pat.ttp_mod) * 1.05); ylims!(ax, 0, 1.02)
axislegend(ax; position = :rt)
savepanel(fig, NAME, "b")

# ── c–j: gain vs fold change of diversity at progression ─────────────────────
function tradeoff_panel(stat, q, source)
    fig = panel()
    ax = Axis(fig[1, 1]; xlabel = LABELS[:gain], ylabel = rich("Fold change in ", hill_label(q, source)))
    hlines!(ax, [1.0]; color = :grey70, linestyle = :dot, linewidth = 0.8)
    vlines!(ax, [0.0]; color = :grey70, linestyle = :dot, linewidth = 0.8)
    for a in AT
        x = Float64.(pat[!, "benefit_$a"]); y = Float64.(pat[!, "fc_$(stat)_$a"])
        ok = isfinite.(x) .& isfinite.(y) .& (y .> 0); any(ok) || continue
        scatter!(ax, x[ok], y[ok]; color = (COL_ARM[a], 0.35), markersize = 2.5)
    end
    for a in AT                                           # medians on top of both clouds
        x = Float64.(pat[!, "benefit_$a"]); y = Float64.(pat[!, "fc_$(stat)_$a"])
        ok = isfinite.(x) .& isfinite.(y) .& (y .> 0); any(ok) || continue
        scatter!(ax, [median(x[ok])], [median(y[ok])]; color = COL_ARM[a], markersize = 7,
                 strokecolor = :black, strokewidth = 0.6)
    end
    return fig
end

for (stats, source, letters) in ((SFS_STATS, "SFS", ["c", "d", "e", "f"]), (CSD_STATS, "CSD", ["g", "h", "i", "j"]))
    for (stat, q, p) in zip(stats, Q_LABELS, letters)
        savepanel(tradeoff_panel(stat, q, source), NAME, p)
    end
end

# ── k, l: parameter importance (PRCC) ─────────────────────────────────────────
"Heatmap of the PRCC of the parameters (columns) with `targets` (rows), values written in the cells."
function prcc_panel(targets, labels; w = PANEL_MM[1], h = 45, split = nothing)
    M = [(d = prcc[(prcc.parameter .== p) .& (prcc.target .== t), :]; nrow(d) == 0 ? NaN : Float64(first(d.rho)))
         for p in PARAMS, t in targets]
    fig = panel(; w, h)
    ax = Axis(fig[1, 1]; xticks = (1:length(PARAMS), [PARAM_NAME[p] for p in PARAMS]),
              yticks = (1:length(targets), labels), xticksize = 0, yticksize = 0)
    hm = heatmap!(ax, 1:length(PARAMS), 1:length(targets), M; colormap = :RdBu, colorrange = (-1, 1), nan_color = :white)
    for i in eachindex(PARAMS), j in eachindex(targets)
        isfinite(M[i, j]) || continue
        text!(ax, i, j; text = string(round(M[i, j]; digits = 2)), align = (:center, :center), fontsize = 6,
              color = abs(M[i, j]) > 0.6 ? :white : :black)
    end
    split === nothing || hlines!(ax, [split + 0.5]; color = :black, linewidth = 0.8)
    ylims!(ax, length(targets) + 0.5, 0.5)                   # reversed: first target on top
    Colorbar(fig[1, 2], hm; label = "PRCC")
    return fig
end

savepanel(prcc_panel(["ttp_mtd", "ttp_hol", "ttp_mod", "benefit_hol", "benefit_mod"],
                     ["TTP, MTD", "TTP, holidays", "TTP, modulation", "Gain, holidays", "Gain, modulation"];
                     w = 60, h = 35, split = 3), NAME, "k")

const DIV_LABELS = [[hill_label(q, "SFS") for q in Q_LABELS]; [hill_label(q, "CSD") for q in Q_LABELS]]
for a in ARMS
    savepanel(prcc_panel(["end_$(s)_$a" for s in [SFS_STATS; CSD_STATS]], DIV_LABELS; split = 4), NAME, "l_$a")
end
for a in AT
    savepanel(prcc_panel(["fc_$(s)_$a" for s in [SFS_STATS; CSD_STATS]], DIV_LABELS; split = 4), NAME, "l_fc_$a")
end

# ── legend ────────────────────────────────────────────────────────────────────
fig = panel(; w = 90, h = 8)
Legend(fig[1, 1], [MarkerElement(marker = :circle, color = COL_ARM[a], markersize = 5) for a in ARMS],
       [ARM_NAME[a] for a in ARMS]; orientation = :horizontal, labelsize = 7)
savepanel(fig, NAME, "legend")

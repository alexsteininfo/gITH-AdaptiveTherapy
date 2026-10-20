# Figure 2 — one treatment round leaves high diversity behind. Arm `single`: full dose until
# 0.1 N_det ("treated, 0.1N"), then regrowth to N_det ("treated, N"); compared with the same tumour
# untreated at these sizes (growth record at 0.1 N_det, detection).
#   a  mean SFS at the four states (neutral)
#   b  passenger Hill numbers ᵠD (q = 0: number of mutations S), treated / untreated, same tumour,
#      at 0.1N and at N (neutral)
#   c  same for weak resistance (α = 0.8)
#   d  resistant-clone Hill numbers at the four states (weak)
#   e  frequency of the largest resistant clone at the four states (weak)
# Data: 01_dynamics/sfs_mean.csv, snapshots.csv
# Run : julia --project=. code/03_figures/fig_2/fig_2.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_2"
const STATES = ["growth", "bottom", "detection", "regrown"]
const PASSENGER_Q = [:S, :D1n, :D2n, :Dinfn]
const RESISTANT_Q = [:res_D0, :res_D1, :res_D2, :res_Dinf]

snap = readproc("01_dynamics", "snapshots.csv")
sfs = readproc("01_dynamics", "sfs_mean.csv")

"State of a record: untreated at 0.1N / N, treated (single round) at 0.1N / N; nothing otherwise."
function state(arm, tag, idx)
    tag == "growth" && return "growth"
    tag == "detection" && return "detection"
    arm == "single" && tag == "blockend" && return idx == 1 ? "bottom" : idx == 2 ? "regrown" : nothing
    return nothing
end

snap.state = [state(r.arm, r.tag, r.idx) for r in eachrow(snap)]
sfs.state = [state(r.arm, r.tag, r.idx) for r in eachrow(sfs)]
snap = snap[.!isnothing.(snap.state), :]

# ── a: mean SFS ───────────────────────────────────────────────────────────────
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:k], ylabel = LABELS[:Sk], LOGX..., LOGY...)
for s in STATES
    d = sfs[(sfs.cohort .== "neutral") .& (sfs.state .=== s) .& (sfs.S_mean .> 0), :]
    lines!(ax, d.k, d.S_mean; color = COL_STATE[s], label = STATE_NAME[s])
end
axislegend(ax; position = :rt)
savepanel(fig, NAME, "a")

# ── b, c: treated / untreated at matched size ────────────────────────────────
"Per-tumour ratio treated / untreated of each column, at 0.1N (bottom/growth) and N (regrown/detection)."
function ratios(cohort, cols)
    d = snap[snap.cohort .== cohort, :]
    out = Dict{Tuple{String,Int},Vector{Float64}}()
    for (size, treated, untreated) in (("0.1N", "bottom", "growth"), ("N", "regrown", "detection"))
        j = innerjoin(d[d.state .== treated, [:tumour; cols]], d[d.state .== untreated, [:tumour; cols]];
                      on = :tumour, renamecols = "_t" => "_u")
        for (i, c) in enumerate(cols)
            r = Float64.(j[!, "$(c)_t"]) ./ Float64.(j[!, "$(c)_u"])
            out[(size, i)] = filter(isfinite, r)
        end
    end
    return out
end

function ratio_panel(cohort)
    r = ratios(cohort, PASSENGER_Q)
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:q], ylabel = LABELS[:ratio], xticks = (1:4, Q_LABELS))
    hlines!(ax, [1.0]; color = :grey70, linestyle = :dot, linewidth = 0.8)
    for (size, dx, col) in (("0.1N", -0.12, COL_STATE["bottom"]), ("N", 0.12, COL_STATE["regrown"]))
        x = Float64[]; med = Float64[]; lo = Float64[]; hi = Float64[]
        for i in 1:4
            v = r[(size, i)]; isempty(v) && continue
            push!(x, i + dx); push!(med, median(v)); push!(lo, q25(v)); push!(hi, q75(v))
        end
        isempty(x) && continue
        errorbars!(ax, x, med, med .- lo, hi .- med; color = col)
        scatter!(ax, x, med; color = col, label = size == "0.1N" ? "At 0.1N" : "At N")
    end
    xlims!(ax, 0.5, 4.5)
    toplegend(fig, ax)
    return fig
end

savepanel(ratio_panel("neutral"), NAME, "b")
savepanel(ratio_panel("weak"), NAME, "c")

# ── d, e: resistant clones of the weak cohort at the four states ─────────────
function states_panel(cols; ylabel, log = false, xticks = (1:length(cols), Q_LABELS[1:length(cols)]),
                      xlabel = LABELS[:q])
    d = snap[snap.cohort .== "weak", :]
    fig = panel(); ax = log ? Axis(fig[1, 1]; xlabel, ylabel, xticks, LOGY...) : Axis(fig[1, 1]; xlabel, ylabel, xticks)
    for (k, s) in enumerate(STATES)
        x = Float64[]; med = Float64[]; lo = Float64[]; hi = Float64[]
        for (i, c) in enumerate(cols)
            v = filter(x -> isfinite(x) && (!log || x > 0), Float64.(d[d.state .== s, c]))
            isempty(v) && continue
            push!(x, i + (k - 2.5) * 0.15); push!(med, median(v)); push!(lo, q25(v)); push!(hi, q75(v))
        end
        isempty(x) && continue
        errorbars!(ax, x, med, med .- lo, hi .- med; color = COL_STATE[s])
        scatter!(ax, x, med; color = COL_STATE[s], label = STATE_NAME[s])
    end
    xlims!(ax, 0.5, length(cols) + 0.5)
    return fig, ax
end

fig, ax = states_panel(RESISTANT_Q; ylabel = LABELS[:D], log = true)
Legend(fig[0, 1], ax; orientation = :horizontal, nbanks = 2, tellwidth = false)
savepanel(fig, NAME, "d")

fig, ax = states_panel([:largest]; ylabel = LABELS[:largest], xticks = ([1], [""]), xlabel = "")
savepanel(fig, NAME, "e")

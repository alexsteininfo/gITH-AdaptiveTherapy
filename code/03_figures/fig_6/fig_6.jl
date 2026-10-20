# Figure 6 — does heterogeneity at diagnosis predict the course and the benefit of adaptive therapy?
# Virtual trial: per-patient ν, α and N_det/K; every patient treated with MTD, holidays and
# modulation. Measured statistics come from a sample of 1000 cells (mutations in ≥ 1% of the
# sample); oracle statistics (*) from the entire tumour and not measurable.
#   a  Spearman ρ of each statistic with the time to progression per arm and with the gain of
#      holidays / modulation over MTD — within one parameter set (Part 1, cohort `full`), across
#      the trial, and partial (across, controlling for ν, α, N_det/K); "early" = end of the first
#      full-dose phase
#   b  gain of holidays vs the statistic that predicts it best, coloured by N_det/K
#   c  |ρ| of that statistic with the gain vs sample size, per detection threshold; dashed: oracle
#   d  treatment policies on the held-out half of the patients: mean time to progression vs median
#      cumulative dose; always MTD / holidays / modulation (coloured) and "adaptive if statistic
#      below (above) a threshold, else MTD" (grey)
# Data: 04_trial/trial.csv, trial_cor.csv, trial_policy.csv
# Run : julia --project=. code/03_figures/fig_6/fig_6.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "fig_6"
const STATS = ["S", "D2n", "H1n", "gamma", "mean_burden", "var_burden"]
const ORACLE = ["det_res_frac", "det_n_res", "det_largest"]
const STAT_NAME = Dict("S" => "Mutations, S", "D2n" => "Diversity, ²D", "H1n" => "Entropy, H₁",
                       "gamma" => "Exponent, γ", "mean_burden" => "Mean burden", "var_burden" => "Burden variance",
                       "det_res_frac" => "Resistant fraction*", "det_n_res" => "Resistant clones*",
                       "det_largest" => "Largest resistant clone*")
const TARGETS = ["ttp_mtd", "ttp_hol", "ttp_mod", "benefit_hol", "benefit_mod"]
const TARGET_NAME = Dict("ttp_mtd" => "TTP MTD", "ttp_hol" => "TTP holidays", "ttp_mod" => "TTP modulation",
                         "benefit_hol" => "Gain holidays", "benefit_mod" => "Gain modulation")
const SETTINGS = ["within", "across", "partial"]

pat = readproc("04_trial", "trial.csv")
cor = readproc("04_trial", "trial_cor.csv")
pol = readproc("04_trial", "trial_policy.csv")

# ── a: correlation heatmap ────────────────────────────────────────────────────
"ρ of one predictor row (source, time) with one target in one setting; NaN if not computed."
function rho(pred, source, time, target, setting)
    d = cor[(cor.predictor .== pred) .& (cor.source .== source) .& (cor.time .== time) .& (cor.target .== target) .&
            (cor.setting .== setting) .& ((source == "oracle") .| ((cor.n .== 1000) .& (cor.fmin .== 0.01))), :]
    return nrow(d) == 0 ? NaN : Float64(first(d.rho))
end

rowspec = [[(p, "sample", "detection", STAT_NAME[p]) for p in STATS];
           [(p, "sample", "early", STAT_NAME[p] * ", early") for p in STATS];
           [(p, "oracle", "detection", STAT_NAME[p]) for p in ORACLE]]
colspec = [(t, s) for s in SETTINGS for t in TARGETS]
M = [rho(p, src, tm, t, s) for (t, s) in colspec, (p, src, tm, _) in rowspec]   # x = columns, y = rows

fig = panel(; w = 110, h = 95)
ax = Axis(fig[1, 1]; xticks = (1:length(colspec), [TARGET_NAME[t] for (t, _) in colspec]), xticklabelrotation = π / 2,
          yticks = (1:length(rowspec), [r[4] for r in rowspec]), xticksize = 0, yticksize = 0)
hm = heatmap!(ax, 1:length(colspec), 1:length(rowspec), M; colormap = :RdBu, colorrange = (-1, 1), nan_color = :white)
vlines!(ax, [5.5, 10.5]; color = :black, linewidth = 0.8)
hlines!(ax, [6.5, 12.5]; color = :black, linewidth = 0.8)
for (i, s) in enumerate(SETTINGS)
    text!(ax, 5i - 2, -0.05; text = uppercasefirst(s), align = (:center, :center), fontsize = 7)
end
ylims!(ax, length(rowspec) + 0.5, -0.6)                     # reversed: first row on top
Colorbar(fig[1, 2], hm; label = LABELS[:rho])
savepanel(fig, NAME, "a")

# ── b, c: the best predictor of the holiday gain ─────────────────────────────
cand = cor[(cor.source .== "sample") .& (cor.time .== "detection") .& (cor.n .== 1000) .& (cor.fmin .== 0.01) .&
           (cor.setting .== "across") .& (cor.target .== "benefit_hol") .& isfinite.(cor.rho), :]
best = nrow(cand) == 0 ? "S" : String(cand.predictor[argmax(abs.(cand.rho))])

fig = panel()
ax = Axis(fig[1, 1]; xlabel = STAT_NAME[best], ylabel = "Gain holidays")
hlines!(ax, [0.0]; color = :grey70, linestyle = :dot, linewidth = 0.8)
sc = scatter!(ax, pat[!, "sample_$best"], pat.benefit_hol; color = pat.fK, colormap = :viridis, markersize = 3)
Colorbar(fig[1, 2], sc; label = "N/K")
savepanel(fig, NAME, "b")

fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:n], ylabel = LABELS[:absrho], LOGX...)
d = cor[(cor.predictor .== best) .& (cor.source .== "sample") .& (cor.time .== "detection") .&
        (cor.setting .== "across") .& (cor.target .== "benefit_hol"), :]
fmins = sort(unique(d.fmin))
for (fm, col) in zip(fmins, shades(:black, length(fmins)))
    g = sort(d[d.fmin .== fm, :], :n)
    lines!(ax, g.n, abs.(g.rho); color = col, label = "f ≥ $(round(Int, 100fm))%")
    scatter!(ax, g.n, abs.(g.rho); color = col)
end
o = rho(best, "oracle", "detection", "benefit_hol", "across")
isfinite(o) && hlines!(ax, [abs(o)]; color = :black, linestyle = :dash, linewidth = 0.8)
ylims!(ax, 0, 1)
toplegend(fig, ax)
savepanel(fig, NAME, "c")

# ── d: policies on the held-out half ──────────────────────────────────────────
test = pol[pol.half .== "test", :]
fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:dose], ylabel = LABELS[:ttp])
rules = test[startswith.(test.policy, "rule_"), :]
scatter!(ax, rules.dose_med, rules.ttp_mean; color = :grey60, markersize = 3, label = "Rules")
for a in ("mtd", "hol", "mod")
    g = test[test.policy .== "always_$a", :]; nrow(g) == 0 && continue
    scatter!(ax, g.dose_med, g.ttp_mean; color = COL_ARM[a], markersize = 6, label = ARM_NAME[a])
end
toplegend(fig, ax; nbanks = 2)
savepanel(fig, NAME, "d")

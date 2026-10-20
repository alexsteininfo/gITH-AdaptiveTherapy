# Supplementary Figure 3 — random resistance levels with a cost: α ~ Beta(1, 4) per resistance
# mutation, cost Δs = −c(1 − α) on the birth rate, c = 0 (cost0) and c = 0.02 (cost0.02); one panel
# per cohort (suffix).
#   a  resistant clones over time per arm (median; while ≥ 50% of tumours run)
#   b  Kaplan–Meier curves per arm
#   legend  arms (modulation at N_det dashed)
# Data: 06_cost/counts.csv.gz, km.csv
# Run : julia --project=. code/03_figures/supp_3/supp_3.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_3"
const ARMS = ["single", "mtd", "hol", "mod", "mod_1.0"]
arm_style(a) = (color = COL_ARM[armkind(a)], linestyle = a == "mod_1.0" ? :dash : :solid)
arm_label(a) = a == "mod_1.0" ? "Modulation at N" : ARM_NAME[a]

counts = readproc("06_cost", "counts.csv.gz")
km = readproc("06_cost", "km.csv")

for c in ("cost0", "cost0.02")
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:resclones], LOGY...)
    for a in ARMS
        d = counts[(counts.cohort .== c) .& (counts.arm .== a), :]; nrow(d) == 0 && continue
        q = running_quantiles(d, :n_res; ntotal = length(unique(d.tumour)))
        lines!(ax, q.t, [x > 0 ? x : NaN for x in q.med]; arm_style(a)...)
    end
    savepanel(fig, NAME, "a_$c")

    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:pfs])
    for a in ARMS
        d = km[(km.cohort .== c) .& (km.arm .== a), :]; nrow(d) == 0 && continue
        kmcurve!(ax, d; arm_style(a)...)
    end
    ylims!(ax, -0.02, 1.02)
    savepanel(fig, NAME, "b_$c")
end

fig = panel(; w = 120, h = 8)
Legend(fig[1, 1], [LineElement(; arm_style(a)...) for a in ARMS], arm_label.(ARMS); orientation = :horizontal, labelsize = 7)
savepanel(fig, NAME, "legend")

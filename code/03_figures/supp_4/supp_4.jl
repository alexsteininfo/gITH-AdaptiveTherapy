# Supplementary Figure 4 — two definitions of passenger diversity: Hill numbers of the normalised
# SFS (main text; solid) and of cell fractions (v1 eq. 10; dashed), relative to the same tumour at
# detection, under holidays and modulation (neutral cohort; median over tumours).
#   a  Shannon entropy H₁
#   b  ²D
# Data: 01_dynamics/snapshots.csv
# Run : julia --project=. code/03_figures/supp_4/supp_4.jl
include(joinpath(@__DIR__, "..", "style.jl"))

const NAME = "supp_4"
const ARMS = ["hol", "mod"]

snap = readproc("01_dynamics", "snapshots.csv")
snap = snap[snap.cohort .== "neutral", :]
det = snap[snap.tag .== "detection", :]

"Per-tumour value of `y` relative to detection, median over tumours per stored time point."
function relative(arm, y)
    d = snap[(snap.arm .== arm) .& in.(snap.tag, Ref(["blockend", "grid", "end"])), :]
    d = innerjoin(d, select(det, :tumour, y => :ref); on = :tumour)
    d.rel = d[!, y] ./ d.ref
    g = combine(groupby(d, [:tag, :idx]), :t_rel => median => :t, :rel => (v -> median(filter(isfinite, v))) => :rel)
    g = sort!(vcat(DataFrame(tag = "detection", idx = 0, t = 0.0, rel = 1.0), g), :t)
    return unique(g, :t)
end

for (p, normalised, cellfrac) in (("a", :H1n, :H1cf), ("b", :D2n, :D2cf))
    fig = panel(); ax = Axis(fig[1, 1]; xlabel = LABELS[:time], ylabel = LABELS[:relative])
    for arm in ARMS, (y, ls) in ((normalised, :solid), (cellfrac, :dash))
        g = relative(arm, y)
        lines!(ax, g.t, g.rel; color = COL_ARM[arm], linestyle = ls)
    end
    p == "a" && Legend(fig[0, 1], [[LineElement(color = COL_ARM[a]) for a in ARMS]; LineElement(color = :black);
                                   LineElement(color = :black, linestyle = :dash)],
                       [[ARM_NAME[a] for a in ARMS]; "Normalised"; "Cell fraction"];
                       orientation = :horizontal, nbanks = 2, tellwidth = false)
    savepanel(fig, NAME, p)
end

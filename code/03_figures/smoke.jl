# Smoke test of the figure scripts on the TEST=1 data: every script runs and writes its panels
# (non-empty PNGs in figures/00_test/<figure>/). Missing test tables are made first with the
# TEST=1 bash scripts (a few minutes each).
# Run: julia --project=. code/03_figures/smoke.jl [regex]     e.g. fig_2 or "fig_[45]"
using Test

ENV["FIG_DATA"] = "test"
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FILTER = isempty(ARGS) ? "" : ARGS[1]

const P1 = "run_part1.sh" => ["00_test_01_dynamics/sfs_mean.csv", "00_test_01_dynamics/km.csv",
                             "00_test_01_dynamics/theory_hill.csv", "00_test_01_dynamics/theory_single_sfs.csv"]
const P2 = "run_part2.sh" => ["00_test_02_tradeoff/tradeoff.csv", "00_test_03_sensitivity/tradeoff.csv"]
const P3 = "run_part3.sh" => ["00_test_04_trial/trial_cor.csv", "00_test_04_trial/trial_policy.csv",
                             "00_test_04_trial/trial_prcc.csv"]
const PS = "run_supp.sh" => ["00_test_01_dynamics/example_traj.csv", "00_test_01_dynamics/powerlaw.csv",
                             "00_test_05_metronomic/km.csv", "00_test_06_cost/km.csv",
                             "00_test_07_validation/validation_ld.csv"]
const P8 = "run_stabilising.sh" => ["00_test_08_stabilising/sfs_mean.csv", "00_test_08_stabilising/counts.csv.gz",
                                    "00_test_08_stabilising/theory_stab_hill.csv", "00_test_08_stabilising/theory_stab_levels.csv"]

# figure => (panels it must write, test tables it reads: bash script => files)
const FIGURES = [
    "fig_1" => (["schematic", "d", "e", "f", "g", "legend"], [P1, PS]),
    "fig_2" => (["a", "b", "c", "d", "e", "f", "g", "h", "legend"], [P1]),
    "fig_3" => (["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "legend", "legend_sfs"], [P1]),
    "fig_4" => (["a", "b", "c_tmin", "c_smin", "d_sfs", "d", "e", "f", "g", "h", "i", "f_size", "g_size", "h_size", "i_size",
                 "j", "k", "l", "m", "dose", "legend", "legend_levels", "legend_sfs"], [P8]),
    "fig_5" => (["b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l_mtd", "l_hol", "l_mod", "l_fc_hol", "l_fc_mod",
                 "legend"], [P3]),
    "supp_1" => (["a", "b", "c", "d"], [PS]),
    "supp_2" => (["a", "b", "c", "d"], [P1, PS]),
    "supp_3" => (["a_cost0", "b_cost0", "a_cost0.02", "b_cost0.02", "legend"], [PS]),
    "supp_4" => (["a", "b"], [P1]),
    "supp_5" => (["a", "b", "c", "d"], [P3]),
    "supp_6" => (["a", "b", "c"], [PS]),
    "supp_7" => (["a", "b", "legend"], [P1]),
    "supp_8" => (["a", "b", "c", "d", "e", "f", "g", "h", "legend"], [P1]),
    "supp_9" => (["a", "b", "c", "d", "e", "f", "g", "h", "legend"], [P1]),
    "supp_10" => (["a", "b", "c", "d", "legend"], [P8]),
    "supp_11" => (["a", "b", "c", "d_hol", "d_mod", "e_hol", "e_mod", "legend"], [P2]),
    "supp_12" => (["a", "b", "c", "d"], [P1, P3]),
]

selected = [f for f in FIGURES if occursin(Regex(FILTER), first(f))]
for (script, files) in unique(vcat([last(last(f)) for f in selected]...))
    all(f -> isfile(joinpath(ROOT, "data", "processed", f)), files) && continue
    println("making test tables: TEST=1 bash code/01_datageneration/$script")
    run(addenv(`bash $(joinpath(ROOT, "code", "01_datageneration", script))`, "TEST" => "1", "WORKERS" => "4"))
end

@testset "figures smoke" begin
    for (name, (panels, _)) in selected
        @testset "$name" begin
            dir = joinpath(ROOT, "figures", "00_test", name)
            isdir(dir) && rm(dir; recursive = true)
            script = joinpath(ROOT, "code", "03_figures", name, "$name.jl")
            @test isfile(script)
            isfile(script) || continue
            m = Module(Symbol(name))                       # each script in its own namespace
            Core.eval(m, :(include(p) = Base.include($m, p)))
            t = @elapsed Base.include(m, script)
            println("$name: $(round(t; digits = 1)) s")
            for p in panels
                f = joinpath(dir, "$p.png")
                @test isfile(f) && filesize(f) > 1000
            end
        end
    end
end

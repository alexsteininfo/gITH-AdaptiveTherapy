# Tests of the shared helpers. Run: julia --project=. code/X_helpers/test/runtests.jl [filter]
include(joinpath(@__DIR__, "..", "helpers.jl"))
using Test, CSV, DataFrames

"Remove raw and processed data of a test block (leftovers of smoke tests would change counts)."
clean_test(b) = for d in (joinpath(ROOT, "data", "raw", b), joinpath(ROOT, "data", "processed", b))
    isdir(d) && rm(d; recursive = true)
end

const FILTER = isempty(ARGS) ? "" : ARGS[1]
runs(name) = isempty(FILTER) || occursin(FILTER, name)

runs("paths") && @testset "paths" begin
    @test isfile(joinpath(ROOT, "Project.toml"))
    @test rawdir("00_test", "c") == joinpath(ROOT, "data", "raw", "00_test", "c") && isdir(rawdir("00_test", "c"))
    @test endswith(procdir("00_test"), joinpath("data", "processed", "00_test")) && endswith(figdir("fig_1"), joinpath("figures", "fig_1"))
    rm(joinpath(ROOT, "data", "raw", "00_test"); recursive = true); rm(procdir("00_test"); recursive = true)
end

runs("registry") && @testset "registry" begin
    @test Set(keys(BLOCKS)) == Set(["01_dynamics", "02_tradeoff", "03_sensitivity", "04_trial", "05_metronomic", "06_cost", "07_validation"])
    f1 = only(filter(c -> c.name == "full", BLOCKS["01_dynamics"])); f2 = only(filter(c -> c.name == "full", BLOCKS["02_tradeoff"]))
    @test f1.key == f2.key && tumour_seed(f1, 7) == tumour_seed(f2, 7) && f2.n == 300 && f1.fK == 0.8
    @test [a.id for a in f1.arms] == ["single", "mtd", "hol", "mod"]
    @test [a.id for a in f2.arms] == ["mtd", "hol", "hol_0.6", "hol_0.7", "hol_0.8", "hol_0.9", "mod_0.5", "mod_0.6", "mod_0.7", "mod", "mod_0.9", "mod_1.0"]
    @test only(filter(c -> c.name == "partial", BLOCKS["01_dynamics"])).alpha == 0.006
    @test nprog(10^6, 0.8) == 1_125_000
    d = lhs_design(500, UInt64(1)); ν = getfield.(d, :nu)
    @test length(d) == 500 && all(1e-5 .<= ν .<= 1e-3) && all(0.5 .<= getfield.(d, :fK) .<= 0.9) && all(0 .<= getfield.(d, :alpha) .<= 0.006)
    @test sort(floor.(Int, 500 .* (log10.(ν) .+ 5) ./ 2)) == 0:499
    tr = only(BLOCKS["04_trial"]); @test params(tr, 3) == params(tr, 3) && params(tr, 3).cost == 0
    @test [c.name for c in BLOCKS["06_cost"]] == ["cost0", "cost0.02"] && [c.cost for c in BLOCKS["06_cost"]] == [0.0, 0.02]
    # a fully resistant clone paying the cost must be able to reach N_prog: K(1 − 2c) > N_prog
    @test all((B0 * (1 - c.cost) - D0) / ((B0 - D0) * c.fK / c.ndet) > nprog(c.ndet, c.fK) for c in BLOCKS["06_cost"])
end

runs("schedules") && @testset "schedules" begin
    p = (nu = 1e-3, alpha = 0.0, fK = 0.8, cost = 0.0, d0 = D0); nd = 10^4
    pop = simulate!(initialize_population(1), growth_block(p, nd), Xoshiro(1)); @test popsize(pop) == nd
    for (k, lab) in ((:mtd, "MTD"), (:hol, "AT on"), (:mod, "stabilising: reduce"), (:single, "single on"), (:met, "metronomic on"))
        q = simulate!(copy(pop), schedule(arm(k; tmax = 20.0), p, nd), Xoshiro(2)); @test q.log[2].label == lab
    end
    q = simulate!(copy(pop), schedule(arm(:mod; level = 1.0, tmax = 5.0), p, nd), Xoshiro(3)); @test q.log[2].ticks == 0
end

runs("stats") && @testset "stats" begin
    h = resistant_hill([1.0, 0.0, 0.0], [100, 30, 10]); @test h.D0 == 2 && h.Dinf ≈ 40 / 30
    @test apply_fmin([5, 4, 3, 2, 1, 0, 0, 0, 0, 1], 0.3) == [0, 0, 3, 2, 1, 0, 0, 0, 0, 1]
    sfs = [round(Int, 1e7 / (k * (k + 1))) for k in 1:10^4]
    f = powerlaw_fit(sfs; kmin = 50); @test isapprox(f.gamma, 2.0; atol = 0.05) && f.se > 0
    @test 2 <= kmin_ks(sfs; kmin_range = 2:200) <= 200
    st = sfs_stats([4, 0, 3]); @test st.S == 7 && st.D2n ≈ 169 / 31 && st.D2cf ≈ 9 / 31
    @test burden_stats(Dict(1 => 2, 3 => 2)) == (mean = 2.0, var = 4 / 3)
    lb = logbin(Float64[10, 5, 3, 2, 1, 1, 1, 1, 0, 0]); @test issorted(lb.k) && all(lb.S .>= 0)
    ps = paired_summary([1.0, 2.0, 3.0, missing], [2.0, 4.0, 5.0, 1.0]); @test ps.median_diff == 2.0 && ps.lo <= 2.0 <= ps.hi
end

runs("Kaplan–Meier with censoring") && @testset "Kaplan–Meier with censoring" begin   # Review Focus 1
    @test km_median([1.0, 2.0, 3.0, 4.0], [true, true, false, false]) == 2.0
    @test ismissing(km_median([5.0, 5.0], [false, false])) && rmst([5.0, 5.0], [false, false], 5.0) == 5.0
    @test rmst([1.0, 3.0], [true, true], 4.0) ≈ 2.0    # S = 1 on [0,1), 0.5 on [1,3), 0 after
end

runs("resistant Hill numbers of an empty clone set") && @testset "resistant Hill numbers of an empty clone set" begin   # Review Focus 2
    h = resistant_hill([1.0], [10]); @test h.D0 == 0 && isnan(h.D2)
end

runs("power-law fit with a short tail") && @testset "power-law fit with a short tail" begin   # Review Focus 3
    f = powerlaw_fit([100, 10, 1, 0, 0, 0, 0, 0, 1, 0]; kmin = 5); @test isnan(f.gamma) && f.n == 1   # k = N is clonal, not tail
end

runs("restricted time of extinct arms") && @testset "restricted time of extinct arms" begin   # Review Focus 1
    o = DataFrame(outcome = ["progressed", "extinct", "censored"], ttp = [5.0, missing, missing],
                  t_end = [5.0, 3.0, 50.000000001], tmax = [50.0, 50.0, 50.0])
    tt, ev = restricted(o)
    @test tt == [5.0, 50.0, 50.0] && ev == [true, false, false]   # cured = progression-free to t_max
end

runs("measure") && @testset "measure" begin
    p = (nu = 1e-3, alpha = 0.0, fK = 0.8, cost = 0.0, d0 = D0)
    pop = simulate!(initialize_population(1), growth_block(p, 3000), Xoshiro(4))
    ctx = (t_det = pop.t, ids = (block = "t", cohort = "c", tumour = 1, arm = "detection"))
    h = measure_heavy(pop, ctx, UInt64(2); tag = "detection", idx = 0)
    @test h.B == total_mutational_burden(pop) && h.B_new == 0 && h.N == 3000
    @test sum(h.tree.sfs.c) == h.B && sum(h.tree.scmb.c) == 3000 && sum(h.tree.clones.size) == 3000
    @test sort(collect(keys(h.samples))) == [100, 1000]                      # 10^4 > N
    s = h.samples[1000]; @test sum(Int.(s.scmb.j) .* s.scmb.c) == sum(Int.(s.sfs.k) .* s.sfs.c)
    c = count_row(pop, ctx); @test c.nS + c.nR == 3000 && c.res_pre <= c.n_res && c.res_new == 0
    @test h.clonal == clonal_mutations(pop)                          # clonal count read from the SFS
end

runs("tail fit ignores clonal mutations") && @testset "tail fit ignores clonal mutations" begin
    v = [round(Int, 1e7 / (k * (k + 1))) for k in 1:10_000]; v[end] = 0
    w = copy(v); w[end] = 2000                                          # 2000 mutations in all N cells
    @test powerlaw_fit(w; kmin = 100).gamma ≈ powerlaw_fit(v; kmin = 100).gamma
    @test abs(powerlaw_fit(w; kmin = 100).gamma - 2) < 0.1          # 1/k² spectrum, truncated where < ½
    @test kmin_ks(w; kmin_range = [10, 50, 100, 500]) == kmin_ks(v; kmin_range = [10, 50, 100, 500])
end

runs("tail fit of thresholded samples") && @testset "tail fit of thresholded samples" begin
    v = [round(Int, 1e6 / (k * (k + 1))) for k in 1:1000]; v[end] = 0
    s = (N = 1000, sfs = sparse_counts(v), scmb = (j = Int32[1, 2], c = Int64[500, 500]),
         clones = (alpha = Float64[], size = Int[]))
    r = (block = "b", cohort = "c", tumour = 1, arm = "a", tag = "detection", idx = 0, t_rel = 0.0, N = 10^6,
         samples = Dict(1000 => s))
    rows = sample_rows(r, 0.0; kmin_frac = 0.01)
    g(f) = only(x.gamma for x in rows if x.fmin == f)
    @test abs(g(0.05) - g(0.0)) < 0.2       # fit starts at the threshold (from 0.01·n it gives ≈ 1.4)
end

runs("grid records at exact times") && @testset "grid records at exact times" begin
    tslip = 25.23                     # adding 10.0 eleven times overshoots t_det + 110·1.0 (count time)
    @test foldl((x, _) -> x + 10.0, 1:11; init = tslip) > tslip + 110 * 1.0
    p = (nu = 0.0, alpha = 1.0, fK = 0.0, cost = 0.0, d0 = D0)       # no competition: never progresses
    pop = simulate!(initialize_population(1), growth_block(p, 1000), Xoshiro(1)); pop.t = tslip
    rec = run_arm(pop, Arm("mod", :mod, 0.8, 115.0, :full), p, 1000, UInt64(7), tslip,
                  (block = "x", cohort = "c", tumour = 1); arm_index = 1)
    g = [h for h in rec["heavy"] if h.tag == "grid"]
    @test length(g) == 11 && all(abs(h.t_rel - h.idx) < 1e-6 for h in g)
end

runs("de novo resistance counts new lineages") && @testset "de novo resistance counts new lineages" begin
    p = (nu = 2e-2, alpha = 0.0, fK = 0.8, cost = 0.0, d0 = D0)
    pop = simulate!(initialize_population(1), growth_block(p, 2000), Xoshiro(3)); t_det = pop.t
    simulate!(pop, remake(base_block(p, 2000); dose = FixedDose(0.0), tmax = t_det + 3), Xoshiro(4))
    ctx = (t_det = t_det, ids = (block = "t", cohort = "c", tumour = 1, arm = "x"))
    ct = pop.clones; res = [k for k in ct.alive if ct.α[k] < 1]
    nested = count(k -> ct.origin_time[k] > t_det && ct.parent[k] != 1, res)
    @test nested > 0                                                  # second hits inside resistant lineages exist
    c = count_row(pop, ctx)
    @test c.res_new == count(k -> ct.origin_time[k] > t_det && ct.parent[k] == 1, res)
    @test c.res_pre == count(k -> ct.origin_time[k] <= t_det && ct.parent[k] == 1, res)
    cs = clone_stats(clone_record(pop), popsize(pop), t_det)
    @test cs.res_new == c.res_new && cs.res_pre == c.res_pre && cs.res_D0 == length(res)
end

runs("sample larger than the population") && @testset "sample larger than the population" begin   # Review Focus 5
    pop = initialize_population(50); ctx = (t_det = 0.0, ids = (block = "t", cohort = "c", tumour = 1, arm = "x"))
    @test isempty(measure_heavy(pop, ctx, UInt64(1); tag = "x", idx = 0).samples)
end

runs("runner") && @testset "runner" begin
    c = test_cohort(only(filter(c -> c.name == "full", BLOCKS["01_dynamics"])))
    rec = run_tumour(c, 1; block = "t", test = true)
    @test Set(keys(rec["arms"])) == Set(["single", "mtd", "hol", "mod"])
    o = rec["arms"]["mtd"]["outcome"]; @test (o.outcome == "progressed") == !ismissing(o.ttp)
    @test rec["detection"].tag == "detection" && !isempty(rec["detection"].samples)
    @test run_tumour(c, 1; block = "t", test = true)["arms"]["hol"]["outcome"] == rec["arms"]["hol"]["outcome"]
    tr = test_cohort(only(BLOCKS["04_trial"])); r = run_tumour(tr, 1; block = "t", test = true)
    @test any(h -> h.tag == "early", r["arms"]["hol"]["heavy"])
end

runs("outcome of extinct and censored arms") && @testset "outcome of extinct and censored arms" begin   # Review Focus 1
    w = test_cohort(only(filter(c -> c.name == "weak", BLOCKS["01_dynamics"]))); r = run_tumour(w, 1; block = "t", test = true)
    @test r["arms"]["mtd"]["outcome"].outcome == "extinct" && ismissing(r["arms"]["mtd"]["outcome"].ttp)
    @test r["arms"]["mod"]["outcome"].outcome == "censored"
end

runs("single round ends regrown") && @testset "single round ends regrown" begin
    c = test_cohort(only(filter(c -> c.name == "full", BLOCKS["01_dynamics"])))
    o = run_tumour(c, 1; block = "t", test = true)["arms"]["single"]["outcome"]
    @test o.outcome in ("regrown", "extinct", "censored") && ismissing(o.ttp)
    d = DataFrame(cohort = "c", arm = ["single", "mtd"], kind = ["single", "mtd"], outcome = "progressed",
                  ttp = [3.0, 4.0], t_end = [3.0, 4.0], tmax = 10.0)
    @test Set(km_table(d).arm) == Set(["mtd"])                        # no survival curve for the single round
end

runs("policies scored on all patients") && @testset "policies scored on all patients" begin
    rng = StableRNG(5); n = 40
    pat = DataFrame(tumour = 1:n)
    for a in ("mtd", "hol", "mod")
        pat[!, "ttp_$a"] = 10 .* rand(rng, n); pat[!, "event_$a"] = trues(n); pat[!, "dose_$a"] = rand(rng, n)
    end
    for p in PREDICTORS; pat[!, "sample_$p"] = rand(rng, n); end
    pat.sample_gamma[1:2:n] .= NaN                                    # undefined for half of the patients
    pol = trial_policies(pat)
    for h in ("fit", "test")
        d = pol[pol.half .== h, :]; @test all(d.n .== d.n[1])
    end
end

runs("early predictors after a completed first phase") && @testset "early predictors after a completed first phase" begin
    sm = DataFrame(tag = "early", arm = "hol", n = 1000, fmin = 0.01, tumour = [1, 2], t_rel = [0.2, 5.0], S = [7.0, 9.0])
    o = DataFrame(arm = "hol", tumour = [1, 2], t_end = [30.0, 5.0])  # tumour 2 progressed in its first phase
    @test isequal(early_values(sm, o, "hol", :S, [1, 2]), [7.0, NaN])
end

runs("resume skips finished, redoes partial") && @testset "resume skips finished, redoes partial" begin   # Review Focus 4
    clean_test("00_test_01_dynamics")
    d = rawdir("00_test_01_dynamics", "full"); write(joinpath(d, "t0002.jld2.tmp"), "partial")
    run_block("01_dynamics"; cohorts = ["full"], test = true)
    @test isfile(joinpath(d, "t0002.jld2")) && !any(endswith(".tmp"), readdir(d))
    m = mtime(joinpath(d, "t0001.jld2")); run_block("01_dynamics"; cohorts = ["full"], test = true)
    @test mtime(joinpath(d, "t0001.jld2")) == m
    rm(joinpath(ROOT, "data", "raw", "00_test_01_dynamics"); recursive = true)
end

runs("theory") && @testset "theory" begin
    p = cyclic_path(; ndet = 10^4, fK = 0.8, level = 0.5, ncycles = 2)
    @test path_N(p, path_time(p, 0, "upper")) ≈ 10^4 && path_N(p, path_time(p, 2, "lower")) ≈ 5000 &&
          path_N(p, path_time(p, 2, "upper")) ≈ 10^4 && path_N(p, path_growth_time(p, 5000)) ≈ 5000
    t = path_time(p, 2, "upper"); lo, hi = fine_bins(10^4)
    lw, hw = fine_bins(10^6)                                             # no truncation at N
    C = sfs_master(p, t, lw, hw; s0 = 0.0, s1 = t)
    @test sum(C .* (lw .+ hw) ./ 2) ≈ 2MU * 10^4 * t rtol = 0.01        # mean burden 2μ∫b for any trajectory
    h = sfs_harmonic(p, t, lo, hi)
    @test sum(h.C) ≈ 2MU * h.N_harm * log(1 + h.b_eff * h.tau) rtol = 0.01
    @test sum(h.C) ≈ sum(sfs_master(p, t, lo, hi; s0 = path_time(p, 0, "upper"), s1 = t)) rtol = 0.15
    r = logbin_ranges(1000); @test r.k ≈ logbin(ones(1000)).k && r.lo[1] == 1 && r.hi[end] == 1000
    q = cyclic_path(; ndet = 10^4, fK = 0.8, level = 0.1, ncycles = 1)          # single round
    @test path_N(q, path_time(q, 1, "lower")) ≈ 1000 && path_N(q, path_time(q, 1, "upper")) ≈ 10^4
end

runs("regrowth adds growth records") && @testset "regrowth adds growth records" begin
    c = test_cohort(only(filter(c -> c.name == "full", BLOCKS["01_dynamics"])))
    old = Cohort(c.name, c.key, c.nu, c.alpha, c.fK, c.cost, c.d0, c.n, Arm[], c.ndet, [0.1])   # stored with 0.1 N_det only
    dir = mktempdir(); p1 = tumour_file(dir, 1); p2 = tumour_file(dir, 2)
    save_tumour(p1, run_tumour(old, 1; block = "00_test_01_dynamics")); cp(p1, p2)
    @test regrow_tumour!(p1, c, 1; block = "00_test_01_dynamics") == :done
    g = load(p1, "growth")
    @test [x.idx for x in g] == [10, 50] && [x.N for x in g] == [c.ndet ÷ 10, c.ndet ÷ 2]
    @test regrow_tumour!(p1, c, 1; block = "00_test_01_dynamics") == :skipped
    m = mtime(p2)                                       # tumour 1's file regrown as tumour 2: another tumour
    @test regrow_tumour!(p2, c, 2; block = "00_test_01_dynamics") == :mismatch && mtime(p2) == m
end

runs("processing on test data") && @testset "processing on test data" begin
    b = "00_test_01_dynamics"; clean_test(b)
    run_block("01_dynamics"; cohorts = ["neutral", "full", "weak"], test = true)
    collect_block(b); statistics_block(b)
    o = CSV.read(joinpath(procdir(b), "outcomes.csv"), DataFrame)
    @test nrow(o) == 3 * (3 + 4 + 4) && all((o.outcome .== "progressed") .== .!ismissing.(o.ttp))
    @test issubset(["nu", "alpha", "fK", "kind", "level", "tmax", "det_N", "det_res_frac", "det_B"], names(o))
    @test all(o.tmax[o.cohort .== "neutral"] .== 40.0 / 5) && all(o.tmax[o.cohort .== "full"] .== TMAX / 5)
    cnt = CSV.read(joinpath(procdir(b), "counts.csv.gz"), DataFrame); @test issubset(["t_rel", "N", "n_res", "res_new"], names(cnt))
    s = CSV.read(joinpath(procdir(b), "snapshots.csv"), DataFrame)
    @test Set(s.tag) ⊇ Set(["growth", "detection", "blockend", "end"]) && issubset(["res_D0", "S", "D2n", "gamma", "var_burden"], names(s))
    @test Set(s.idx[s.tag .== "growth"]) == Set([10, 50])
    sm = CSV.read(joinpath(procdir(b), "samples.csv"), DataFrame)
    @test all(sm.n .<= sm.N) && Set(sm.fmin) == Set(FMIN)
    statistics_block(b; kmin_frac = 0.05, suffix = "_kmin0.05")
    s2 = CSV.read(joinpath(procdir(b), "snapshots_kmin0.05.csv"), DataFrame)
    @test nrow(s2) == nrow(s) && !isequal(s2.gamma, s.gamma)
    rm(joinpath(ROOT, "data", "raw", b); recursive = true); rm(procdir(b); recursive = true)
end

runs("block without arms") && @testset "block without arms" begin
    b = "00_test_07_validation"; clean_test(b)
    run_block("07_validation"; cohorts = ["sfs_bd"], test = true)
    collect_block(b); statistics_block(b)
    s = CSV.read(joinpath(procdir(b), "snapshots.csv"), DataFrame)
    @test nrow(s) == 3 && all(s.tag .== "detection")
    rm(joinpath(ROOT, "data", "raw", b); recursive = true); rm(procdir(b); recursive = true)
end

runs("summaries on test data") && @testset "summaries on test data" begin
    for b in ("01_dynamics", "02_tradeoff", "04_trial", "05_metronomic", "07_validation"); clean_test("00_test_$b"); end
    run_block("01_dynamics"; cohorts = ["full"], test = true); run_block("02_tradeoff"; cohorts = ["full"], test = true)
    run_block("04_trial"; test = true); run_block("05_metronomic"; cohorts = ["full"], test = true)
    run_block("07_validation"; test = true)
    for b in ("00_test_01_dynamics", "00_test_02_tradeoff", "00_test_04_trial", "00_test_05_metronomic", "00_test_07_validation")
        collect_block(b); statistics_block(b)
    end
    make_sfs("00_test_01_dynamics")
    s = CSV.read(joinpath(procdir("00_test_01_dynamics"), "sfs_mean.csv"), DataFrame)
    @test all(issorted(g.k) for g in groupby(s, [:cohort, :arm, :tag, :idx])) && all(s.S_mean .>= 0)
    make_examples("00_test_01_dynamics", "00_test_05_metronomic")
    ex = CSV.read(joinpath(procdir("00_test_01_dynamics"), "example_traj.csv"), DataFrame)
    @test Set(ex.arm) == Set(["mtd", "hol", "mod", "met"]) && length(unique(ex.tumour)) == 1
    make_tradeoff("00_test_02_tradeoff")                          # also checks determinism against block 01
    tr = CSV.read(joinpath(procdir("00_test_02_tradeoff"), "tradeoff.csv"), DataFrame)
    @test nrow(tr) == 1 + length(LEVELS_HOL) + length(LEVELS_MOD) && all(tr.n .== 3)
    km = CSV.read(joinpath(procdir("00_test_02_tradeoff"), "km.csv"), DataFrame)
    @test all(issorted(g.S; rev = true) for g in groupby(km, [:cohort, :arm]))
    @test all(g.t[end] == TMAX / 5 for g in groupby(km, [:cohort, :arm]))     # curves run to the follow-up
    make_trial("00_test_04_trial"; within_block = "00_test_01_dynamics")
    t = CSV.read(joinpath(procdir("00_test_04_trial"), "trial.csv"), DataFrame); @test nrow(t) == 3
    cr = CSV.read(joinpath(procdir("00_test_04_trial"), "trial_cor.csv"), DataFrame)
    @test all(x -> isnan(x) || -1 <= x <= 1, cr.rho) && Set(cr.setting) ⊆ Set(["within", "across", "partial"])
    @test isfile(joinpath(procdir("00_test_04_trial"), "trial_policy.csv"))
    tk = CSV.read(joinpath(procdir("00_test_04_trial"), "trial_km.csv"), DataFrame)
    @test all(g.t[end] == TMAX / 5 for g in groupby(tk, :arm))
    make_powerlaw("00_test_01_dynamics")
    pl = CSV.read(joinpath(procdir("00_test_01_dynamics"), "powerlaw.csv"), DataFrame)
    @test issubset(["cohort", "tumour", "arm", "tag", "kmin", "kmin_frac", "gamma", "se", "n_tail", "ks_kmin"], names(pl))
    @test Set(pl.tag) == Set(["detection", "blockend"]) && all(in.(pl.ks_kmin, Ref(Set(pl.kmin))))
    @test count(isfinite, pl.gamma) > nrow(pl) / 2 && all(pl.se[isfinite.(pl.se)] .> 0)
    make_validation("00_test_07_validation")
    v = CSV.read(joinpath(procdir("00_test_07_validation"), "validation_sfs.csv"), DataFrame)
    @test Set(v.cohort) == Set(["sfs_pb", "sfs_bd"]) && all(v.S_pred .> 0)
    for b in ("01_dynamics", "02_tradeoff", "04_trial", "05_metronomic", "07_validation"); clean_test("00_test_$b"); end
end

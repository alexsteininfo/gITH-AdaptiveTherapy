# Summaries across tumours, written next to the block's tables in data/processed/<block>/:
#   make_sfs         sfs_mean.csv           log-binned mean SFS per cohort × arm × time point
#   make_examples    example_traj.csv, example_clones.csv   (Figure 1)
#   make_tradeoff    tradeoff.csv, km.csv   (containment scans; paired gains against MTD)
#   make_trial       trial.csv, trial_cor.csv, trial_policy.csv, trial_km.csv
#   make_powerlaw    powerlaw.csv           γ̂ over k_min and the KS choice (supplement)
#   make_validation  validation_sfs.csv, validation_ld.csv

using StatsBase: corspearman, tiedrank, quantile
using Distributions: TDist, ccdf

readtable(block, name) = CSV.read(joinpath(ROOT, "data", "processed", block, name), DataFrame)

# ── Mean SFS ───────────────────────────────────────────────────────────────────
function _sfs_tumour(path)
    out = Dict{Tuple{String,String,String,Int},NamedTuple}()
    jldopen(path, "r") do f
        recs = Any[f["growth"]...; f["detection"]]
        ndet = f["detection"].tree.N
        for id in arm_ids(f)
            heavy = f["arms/$id/heavy"]
            # single round: only tumours that respond (first block ends at 0.1 N_det)
            id == "single" && !any(r -> r.tag == "blockend" && r.idx == 1 && r.tree.N <= ndet / 10, heavy) && continue
            append!(recs, heavy)
        end
        for r in recs
            lb = logbin(expand_counts(r.tree.sfs, r.tree.N))
            out[(r.cohort, r.arm, r.tag, r.idx)] = lb
        end
    end
    return out
end

"Log-binned SFS averaged over tumours for every cohort × arm × stored time point."
function make_sfs(block::AbstractString)
    parts = Distributed.pmap(x -> _sfs_tumour(x[2]), tumour_files(block))
    groups = Dict{Tuple{String,String,String,Int},Vector{NamedTuple}}()
    for p in parts, (k, v) in p
        push!(get!(groups, k, NamedTuple[]), v)
    end
    rows = NamedTuple[]
    for ((cohort, arm, tag, idx), lbs) in groups
        L = maximum(length(lb.k) for lb in lbs)
        ref = lbs[argmax([length(lb.k) for lb in lbs])]
        for i in 1:L
            v = Float64[lb.S[i] for lb in lbs if length(lb.S) >= i]
            push!(rows, (cohort = cohort, arm = arm, tag = tag, idx = idx, k = ref.k[i], S_mean = mean(v),
                         S_q25 = quantile(v, 0.25), S_q75 = quantile(v, 0.75), n_tumours = length(v)))
        end
    end
    df = sort!(DataFrame(rows), [:cohort, :arm, :tag, :idx, :k])
    CSV.write(joinpath(procdir(block), "sfs_mean.csv"), df)
    return df
end

# ── Example trajectories (Figure 1) ───────────────────────────────────────────
"""
    make_examples(block01, block05)

Example tumour of cohort `full`: among the tumours with stored clone histories, the one whose
resistant fraction at detection is closest to the cohort median. Writes its trajectories and
resistant-clone sizes for `mtd`, `hol`, `mod` (block 01) and `met` (block 05).
"""
function make_examples(block01::AbstractString, block05::AbstractString)
    o = readtable(block01, "outcomes.csv")
    full = unique(o[o.cohort .== "full", [:tumour, :det_res_frac]])
    med = median(full.det_res_frac)
    cand = full[full.tumour .<= N_EXAMPLES, :]
    t = cand.tumour[argmin(abs.(cand.det_res_frac .- med))]
    traj = NamedTuple[]; clones = NamedTuple[]
    for (block, arms) in ((block01, ("mtd", "hol", "mod")), (block05, ("met",)))
        path = tumour_file(joinpath(ROOT, "data", "raw", block, "full"), t)
        isfile(path) || continue
        jldopen(path, "r") do f
            for a in arms
                haskey(f, "arms/$a") || continue
                tr = f["arms/$a/traj"]
                for i in eachindex(tr.t_rel)
                    push!(traj, (cohort = "full", tumour = t, arm = a, t_rel = tr.t_rel[i], N = tr.N[i],
                                 nS = tr.nS[i], nR = tr.nR[i], m = tr.m[i]))
                end
                haskey(f, "arms/$a/clones") || continue
                for snap in f["arms/$a/clones"], j in eachindex(snap.clones.id)
                    (snap.clones.alpha[j] < 1 && snap.clones.size[j] > 0) || continue
                    push!(clones, (cohort = "full", tumour = t, arm = a, t_rel = snap.t_rel, clone = snap.clones.id[j],
                                   alpha = snap.clones.alpha[j], size = snap.clones.size[j]))
                end
            end
        end
    end
    CSV.write(joinpath(procdir(block01), "example_traj.csv"), DataFrame(traj))
    CSV.write(joinpath(procdir(block01), "example_clones.csv"), DataFrame(clones))
    return t
end

# ── Trade-off ─────────────────────────────────────────────────────────────────
"""
Time to progression restricted to the follow-up, and the event flag. Tumours that did not progress
count as progression-free until the arm's follow-up `tmax` — censored ones and extinct (cured) ones.
"""
restricted(o) = ([r.outcome == "progressed" ? Float64(r.ttp) : Float64(r.tmax) for r in eachrow(o)],
                 o.outcome .== "progressed")

"Kaplan–Meier curve as points `(t, S)` from (0, 1) to the end of follow-up `tmax`."
function km_points(tt, ev, tmax)
    k = km(tt, ev)
    t = [0.0; k.t]; S = [1.0; k.S]
    t[end] < tmax && (push!(t, tmax); push!(S, S[end]))
    return (t = t, S = S)
end

"Kaplan–Meier curves per cohort × arm, from (0, 1) to the arm's follow-up (not for the single round)."
function km_table(o)
    rows = NamedTuple[]
    for g in groupby(o[o.kind .!= "single", :], [:cohort, :arm])
        tt, ev = restricted(g)
        k = km_points(tt, ev, maximum(g.tmax))
        append!(rows, [(cohort = g.cohort[1], arm = g.arm[1], t = t, S = s) for (t, s) in zip(k.t, k.S)])
    end
    return DataFrame(rows)
end

"Check that shared tumours give identical outcomes in two blocks (determinism of the pairing)."
function check_determinism(block, other)
    isfile(joinpath(ROOT, "data", "processed", other, "outcomes.csv")) || return nothing
    a = readtable(block, "outcomes.csv"); b = readtable(other, "outcomes.csv")
    j = innerjoin(a[:, [:cohort, :tumour, :arm, :outcome, :t_end]], b[:, [:cohort, :tumour, :arm, :outcome, :t_end]],
                  on = [:cohort, :tumour, :arm], makeunique = true)
    bad = j[(j.outcome .!= j.outcome_1) .| (j.t_end .!= j.t_end_1), :]
    nrow(bad) == 0 || error("determinism check failed: $(nrow(bad)) arms differ between $block and $other")
    return nrow(j)
end

"""
    make_tradeoff(block)

Per cohort × arm: progression, Kaplan–Meier median and restricted mean time to progression,
resistant diversity and burden at the end of the arm and resistant clones at t = 50, and the
paired gain in restricted time to progression against MTD.
"""
function make_tradeoff(block::AbstractString)
    o = readtable(block, "outcomes.csv"); s = readtable(block, "snapshots.csv"); c = readtable(block, "counts.csv.gz")
    o = o[o.kind .!= "single", :]                      # the single round ends on regrowth, not progression
    ends = s[s.tag .== "end", :]
    rows = NamedTuple[]
    for g in groupby(o, [:cohort, :arm])
        cohort, a = g.cohort[1], g.arm[1]
        tt, ev = restricted(g)
        tmax = maximum(g.tmax)
        e = innerjoin(ends[(ends.cohort .== cohort) .& (ends.arm .== a), :], g[:, [:tumour]], on = :tumour)
        e50 = combine(groupby(c[(c.cohort .== cohort) .& (c.arm .== a) .& (abs.(c.t_rel .- 50) .<= 1), :], :tumour),
                      :n_res => (x -> isempty(x) ? missing : first(x)) => :n_res)
        mt = sort(o[(o.cohort .== cohort) .& (o.arm .== "mtd"), :], :tumour)
        gs = sort(g, :tumour)
        gain = nrow(mt) == nrow(gs) ? paired_summary(restricted(mt)[1], restricted(gs)[1]) :
               (median_diff = NaN, lo = NaN, hi = NaN, p = NaN)
        q(v, p) = isempty(v) ? NaN : quantile(collect(skipmissing(v)), p)
        push!(rows, (cohort = cohort, arm = a, kind = g.kind[1], level = g.level[1], n = nrow(g),
                     n_progressed = count(ev), n_censored = count(.!ev),
                     ttp_km_median = coalesce(km_median(tt, ev), NaN), ttp_rmst = rmst(tt, ev, tmax),
                     resD0_med = q(e.res_D0, 0.5), resD0_q25 = q(e.res_D0, 0.25), resD0_q75 = q(e.res_D0, 0.75),
                     resD2_med = q(filter(!isnan, e.res_D2), 0.5), B_med = q(e.B, 0.5), B_q25 = q(e.B, 0.25),
                     B_q75 = q(e.B, 0.75), resD0_t50_med = q(e50.n_res, 0.5),
                     gain_med = gain.median_diff, gain_lo = gain.lo, gain_hi = gain.hi, gain_p = gain.p))
    end
    df = sort!(DataFrame(rows), [:cohort, :kind, :level])
    CSV.write(joinpath(procdir(block), "tradeoff.csv"), df)
    CSV.write(joinpath(procdir(block), "km.csv"), km_table(o))
    startswith(registry_block(block), "02") && check_determinism(block, replace(block, "02_tradeoff" => "01_dynamics"))
    return df
end

# ── Virtual trial ─────────────────────────────────────────────────────────────
"""
    early_values(sm, o, a, p, tumours)

Sample statistic `p` (n = 1000, f_min = 1%) at the end of the first full-dose phase of arm `a`, per
tumour; NaN if missing or if the arm ended in that phase (the record would be taken at progression
or extinction and carry the outcome).
"""
function early_values(sm, o, a, p, tumours)
    e = sm[(sm.tag .== "early") .& (sm.arm .== a) .& (sm.n .== 1000) .& (sm.fmin .== 0.01), :]
    ends = Dict(zip(o.tumour[o.arm .== a], o.t_end[o.arm .== a]))
    return [(v = e[e.tumour .== t, :]; (isempty(v) || first(v.t_rel) >= get(ends, t, Inf) - 1e-9) ? NaN :
             Float64(first(v[!, p]))) for t in tumours]
end

const PREDICTORS = [:S, :D2n, :H1n, :gamma, :mean_burden, :var_burden]
const ORACLE_ONLY = [:det_res_frac, :det_n_res, :det_largest]

"Spearman ρ with t-approximation p and bootstrap 95% CI (NaN if fewer than 6 complete pairs)."
function spearman(x, y; rng = StableRNG(2))
    keep = [!ismissing(a) && !ismissing(b) && isfinite(a) && isfinite(b) for (a, b) in zip(x, y)]
    x = Float64.(x[keep]); y = Float64.(y[keep]); n = length(x)
    n < 6 && return (rho = NaN, p = NaN, lo = NaN, hi = NaN, n = n)
    ρ = corspearman(x, y)
    tstat = ρ * sqrt((n - 2) / max(1e-12, 1 - ρ^2))
    boots = [(i = rand(rng, 1:n, n); corspearman(x[i], y[i])) for _ in 1:500]
    boots = filter(isfinite, boots)
    return (rho = ρ, p = 2ccdf(TDist(n - 2), abs(tstat)), lo = quantile(boots, 0.025), hi = quantile(boots, 0.975), n = n)
end

"Partial Spearman ρ: correlation of rank residuals after regressing on the ranks of `Z`'s columns."
function partial_spearman(x, y, Z::AbstractMatrix)
    keep = [isfinite(a) && isfinite(b) for (a, b) in zip(x, y)]
    sum(keep) < size(Z, 2) + 6 && return (rho = NaN, p = NaN, lo = NaN, hi = NaN, n = sum(keep))
    X = hcat(ones(sum(keep)), [tiedrank(Z[keep, j]) for j in 1:size(Z, 2)]...)
    res(v) = (r = tiedrank(v[keep]); r .- X * (X \ r))
    return spearman(res(x), res(y))
end

"""
    make_trial(block; within_block = "01_dynamics")

Per-patient table, correlations of diversity at diagnosis (and early on treatment) with time to
progression and AT benefit — within one parameter set (`within_block`, cohort `full`), across the
trial, and partial (controlling for log ν, α, fK) — plus simple treatment-selection policies.
"""
function make_trial(block::AbstractString; within_block::AbstractString = "01_dynamics")
    o = readtable(block, "outcomes.csv"); s = readtable(block, "snapshots.csv"); sm = readtable(block, "samples.csv")
    arms = ["mtd", "hol", "mod"]
    pat = unique(o[:, [:tumour, :nu, :alpha, :fK, :det_res_frac, :det_n_res, :det_largest]])
    for a in arms
        g = sort(o[o.arm .== a, :], :tumour); tt, ev = restricted(g)
        pat[!, "ttp_$a"] = tt; pat[!, "event_$a"] = ev; pat[!, "dose_$a"] = g.cum_dose
    end
    pat.benefit_hol = pat.ttp_hol .- pat.ttp_mtd; pat.benefit_mod = pat.ttp_mod .- pat.ttp_mtd
    det = s[s.tag .== "detection", :]; detm = sm[(sm.tag .== "detection") .& (sm.n .== 1000) .& (sm.fmin .== 0.01), :]
    for p in PREDICTORS
        pat[!, "oracle_$p"] = [first(det[det.tumour .== t, p]) for t in pat.tumour]
        pat[!, "sample_$p"] = [(v = detm[detm.tumour .== t, p]; isempty(v) ? NaN : first(v)) for t in pat.tumour]
    end
    CSV.write(joinpath(procdir(block), "trial.csv"), pat)

    targets = ["ttp_mtd", "ttp_hol", "ttp_mod", "benefit_hol", "benefit_mod"]
    Z = hcat(log10.(pat.nu), pat.alpha, pat.fK)
    rows = NamedTuple[]
    addrow!(pred, src, time, n, fmin, target, setting, r) =
        push!(rows, (predictor = String(pred), source = src, time = time, n = n, fmin = fmin, target = target,
                     setting = setting, rho = r.rho, p = r.p, lo = r.lo, hi = r.hi, n_patients = r.n))
    for target in targets
        y = Float64.(pat[!, target])
        for p in PREDICTORS, src in ("oracle", "sample")
            x = Float64.(pat[!, "$(src)_$p"])
            addrow!(p, src, "detection", src == "sample" ? 1000 : 0, src == "sample" ? 0.01 : 0.0, target, "across", spearman(x, y))
            addrow!(p, src, "detection", src == "sample" ? 1000 : 0, src == "sample" ? 0.01 : 0.0, target, "partial", partial_spearman(x, y, Z))
        end
        for p in ORACLE_ONLY
            x = Float64.(pat[!, p])
            addrow!(p, "oracle", "detection", 0, 0.0, target, "across", spearman(x, y))
            addrow!(p, "oracle", "detection", 0, 0.0, target, "partial", partial_spearman(x, y, Z))
        end
        # every sample size and threshold (across)
        for n in SAMPLE_N, fm in FMIN, p in PREDICTORS
            d = sm[(sm.tag .== "detection") .& (sm.n .== n) .& (sm.fmin .== fm), :]
            x = [(v = d[d.tumour .== t, p]; isempty(v) ? NaN : first(v)) for t in pat.tumour]
            addrow!(p, "sample", "detection", n, fm, target, "across", spearman(x, y))
        end
    end
    # early on treatment (end of the first full-dose phase), for the arm's own course
    for a in ("hol", "mod"), p in PREDICTORS
        x = early_values(sm, o, a, p, pat.tumour)
        for target in ("ttp_$a", "benefit_$a")
            addrow!(p, "sample", "early", 1000, 0.01, target, "across", spearman(x, Float64.(pat[!, target])))
            addrow!(p, "sample", "early", 1000, 0.01, target, "partial", partial_spearman(x, Float64.(pat[!, target]), Z))
        end
    end
    # within one parameter set (Part 1, cohort full): stochastic variation only
    if isfile(joinpath(ROOT, "data", "processed", within_block, "outcomes.csv"))
        wo = readtable(within_block, "outcomes.csv"); ws = readtable(within_block, "samples.csv")
        wo = wo[wo.cohort .== "full", :]
        wdet = ws[(ws.cohort .== "full") .& (ws.tag .== "detection") .& (ws.n .== 1000) .& (ws.fmin .== 0.01), :]
        tum = sort(unique(wo.tumour))
        wt = Dict(a => (g = sort(wo[wo.arm .== a, :], :tumour); restricted(g)[1]) for a in arms)
        wy = Dict("ttp_mtd" => wt["mtd"], "ttp_hol" => wt["hol"], "ttp_mod" => wt["mod"],
                  "benefit_hol" => wt["hol"] .- wt["mtd"], "benefit_mod" => wt["mod"] .- wt["mtd"])
        for p in PREDICTORS
            x = [(v = wdet[wdet.tumour .== t, p]; isempty(v) ? NaN : first(v)) for t in tum]
            for target in targets
                addrow!(p, "sample", "detection", 1000, 0.01, target, "within", spearman(x, wy[target]))
            end
        end
    end
    CSV.write(joinpath(procdir(block), "trial_cor.csv"), DataFrame(rows))
    CSV.write(joinpath(procdir(block), "trial_policy.csv"), trial_policies(pat))
    km_rows = NamedTuple[]
    for a in arms
        k = km_points(pat[!, "ttp_$a"], pat[!, "event_$a"], maximum(o.tmax[o.arm .== a]))
        append!(km_rows, [(arm = a, t = t, S = v) for (t, v) in zip(k.t, k.S)])
    end
    CSV.write(joinpath(procdir(block), "trial_km.csv"), DataFrame(km_rows))
    return pat
end

"""
    trial_policies(pat)

Always-MTD / always-holidays / always-modulation, and rules "adaptive (holidays or modulation) if a
sample statistic at diagnosis is below (or above) a threshold, otherwise MTD" (also MTD when the
statistic is undefined, so every policy is scored on the same patients). The threshold and
direction are fitted on a random half of the patients (largest mean restricted TTP) and evaluated
on the other half.
"""
function trial_policies(pat)
    n = nrow(pat); rng = StableRNG(3)
    test = falses(n); test[randperm(rng, n)[1:fld(n, 2)]] .= true
    rows = NamedTuple[]
    summ(policy, thr, ttp, ev, dose, half) = push!(rows, (policy = policy, threshold = thr,
        ttp_mean = mean(ttp), ttp_km_median = coalesce(km_median(ttp, ev), NaN), dose_med = median(dose),
        n = length(ttp), half = half))
    for a in ("mtd", "hol", "mod"), (half, idx) in (("fit", .!test), ("test", test))
        any(idx) && summ("always_$a", NaN, pat[idx, "ttp_$a"], pat[idx, "event_$a"], pat[idx, "dose_$a"], half)
    end
    for at in ("hol", "mod"), p in PREDICTORS
        x = pat[!, "sample_$p"]; ok = isfinite.(x)
        sum(ok .& .!test) < 4 && continue
        rule(thr, dir) = ok .& (dir .* x .< dir * thr)          # adaptive if the rule holds, MTD otherwise
        best = (score = -Inf, thr = NaN, dir = 1)
        for thr in quantile(x[ok .& .!test], 0.1:0.1:0.9), dir in (1, -1)
            choose = rule(thr, dir); idx = .!test
            sc = mean(ifelse.(choose[idx], pat[idx, "ttp_$at"], pat[idx, "ttp_mtd"]))
            sc > best.score && (best = (score = sc, thr = thr, dir = dir))
        end
        for (half, idx) in (("fit", .!test), ("test", test))
            any(idx) || continue
            choose = rule(best.thr, best.dir)
            ttp = ifelse.(choose[idx], pat[idx, "ttp_$at"], pat[idx, "ttp_mtd"])
            ev = ifelse.(choose[idx], pat[idx, "event_$at"], pat[idx, "event_mtd"])
            dose = ifelse.(choose[idx], pat[idx, "dose_$at"], pat[idx, "dose_mtd"])
            summ("rule_$(at)_$(p)_$(best.dir == 1 ? "below" : "above")", best.thr, ttp, ev, dose, half)
        end
    end
    return DataFrame(rows)
end

# ── Power law over k_min (supplement) ─────────────────────────────────────────
const KMIN_FRACS = 10 .^ range(-4, -1; length = 25)    # k_min / N, 8 per decade

function _powerlaw_tumour(path, cohorts)
    rows = NamedTuple[]
    jldopen(path, "r") do f
        m = f["meta"]
        m.cohort in cohorts || return
        recs = Any[f["detection"]]
        if haskey(f, "arms/hol")
            hol = [r for r in f["arms/hol/heavy"] if r.tag == "blockend"]
            isempty(hol) || push!(recs, hol[argmax([r.idx for r in hol])])
        end
        for r in recs
            N = r.tree.N; sfs = expand_counts(r.tree.sfs, N)
            kmins = unique(max.(2, round.(Int, KMIN_FRACS .* N)))
            ks = kmin_ks(sfs; kmin_range = kmins)
            for k in kmins
                fit = powerlaw_fit(sfs; kmin = k)
                push!(rows, (cohort = m.cohort, tumour = m.tumour, arm = r.arm, tag = r.tag, idx = r.idx, N = N,
                             kmin = k, kmin_frac = k / N, gamma = fit.gamma, se = fit.se, n_tail = fit.n,
                             ks_kmin = ks))
            end
        end
    end
    return rows
end

"""
    make_powerlaw(block; cohorts = ["neutral", "full"])

Tail exponent γ̂ with standard error over a grid of k_min (10⁻⁴–10⁻¹ N) and the k_min chosen by
the Kolmogorov–Smirnov distance, for the entire tree at detection and after the last stored
holiday cycle. Writes `powerlaw.csv`.
"""
function make_powerlaw(block::AbstractString; cohorts = ["neutral", "full"])
    files = [p for (c, p) in tumour_files(block) if c in cohorts]
    rows = reduce(vcat, Distributed.pmap(p -> _powerlaw_tumour(p, cohorts), files); init = NamedTuple[])
    df = DataFrame(rows)
    CSV.write(joinpath(procdir(block), "powerlaw.csv"), df)
    return df
end

# ── Validation ────────────────────────────────────────────────────────────────
"Expected SFS of a birth–death process at size N (Stein & Werner 2025, eq. 9; μ per daughter)."
sfs_pred(k, N, b, d; μ = MU) = 2μ * N * sum((d / b)^l / ((k + l) * (k + l + 1)) for l in 0:2000)

"Mean SFS (k ≤ 200) against theory for `sfs_pb`/`sfs_bd`; resistant clones and cells at detection for `ld`."
function make_validation(block::AbstractString)
    rows = NamedTuple[]; ld = NamedTuple[]
    for cohort in ("sfs_pb", "sfs_bd", "ld")
        files = [p for (c, p) in tumour_files(block) if c == cohort]
        isempty(files) && continue
        c = only(filter(x -> x.name == cohort, BLOCKS[registry_block(block)]))
        if cohort == "ld"
            for p in files
                jldopen(p, "r") do f
                    cl = f["detection"].tree.clones; res = cl.alpha .< 1
                    push!(ld, (tumour = f["meta"].tumour, n_res_clones = count(res), n_res_cells = sum(cl.size[res]; init = 0)))
                end
            end
        else
            sfs = [jldopen(p -> expand_counts(p["detection"].tree.sfs, p["detection"].tree.N)[1:200], f, "r") for f in files]
            N = jldopen(f -> f["detection"].tree.N, files[1], "r")
            for k in 1:200
                push!(rows, (cohort = cohort, k = k, S_mean = mean(s[k] for s in sfs),
                             S_pred = sfs_pred(k, N, B0, c.d0)))
            end
        end
    end
    CSV.write(joinpath(procdir(block), "validation_sfs.csv"), DataFrame(rows))
    CSV.write(joinpath(procdir(block), "validation_ld.csv"), DataFrame(ld))
    return nothing
end

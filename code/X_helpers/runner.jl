# Running tumours: growth to detection, every arm on a copy of the detection population with the
# same post-detection seed, atomic per-tumour files, resumable blocks on worker processes.

using Distributed: pmap, myid

"Plain NamedTuple of a block summary (no package types in the stored files)."
log_row(b) = (label = b.label, t_start = b.t_start, t_end = b.t_end, end_reason = String(b.end_reason),
              dose = b.dose, m_end = b.m_end)

"Cumulative dose ∫ m dt from trajectory points `(t, m)` (trapezoid rule)."
function cum_dose(traj)
    length(traj) < 2 && return 0.0
    return sum((traj[i].t - traj[i-1].t) * (traj[i].m + traj[i-1].m) / 2 for i in 2:length(traj))
end

"Number of completed holiday cycles: full-dose blocks that ended at the lower threshold."
n_cycles(log) = count(b -> b.label == "AT on" && b.end_reason == :below, log)

"""
    run_arm(pop0, a, p, ndet, seed, t_det, ids; arm_index, example = false)

Treat a copy of the detection population with arm `a`. Returns outcome row, block log,
trajectory, count rows, distribution records (`heavy`) and, for example tumours, clone records.
"""
function run_arm(pop0, a::Arm, p, ndet, seed, t_det, ids; arm_index::Integer, example::Bool = false)
    pop = copy(pop0); nlog0 = length(pop.log)
    sch = schedule(a, p, ndet)
    rng = Xoshiro(derive_seed(seed, 1))
    ctx = (t_det = t_det, ids = merge(ids, (arm = a.id,)))
    heavy = Any[]; counter = Ref(0)
    record!(q, tag, idx) = (counter[] += 1;
        push!(heavy, measure_heavy(q, ctx, derive_seed(seed, 2, arm_index, counter[]); tag, idx)))
    # grid points t_det + j·HEAVY_EVERY, computed like the count times t_det + k·COUNT_EVERY
    # (repeated addition drifts by an ulp and moves a record to the next count)
    j = Ref(1)
    at_grid(q) = q.t >= t_det + j[] * HEAVY_EVERY - 1e-9
    function grid_heavy(q)
        a.heavy == :full || return nothing
        if at_grid(q)
            record!(q, "grid", round(Int, j[] * HEAVY_EVERY))
            while at_grid(q); j[] += 1; end
        end
        return nothing
    end
    stats = Any[:counts => q -> count_row(q, ctx), :heavy => grid_heavy]
    example && push!(stats, :clones => q -> clone_record(q))
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = TRAJ_DT,
        triggers = [Every(a.heavy == :end ? COUNT_EVERY_SCAN : COUNT_EVERY)], stats = stats))
    t0 = pop.t; i = 0
    while (b = next_block(sch, pop, i, t0)) !== nothing
        simulate!(pop, b, rng; accumulator = acc)
        i += 1
        a.heavy == :full && i <= HEAVY_BLOCKS && record!(pop, "blockend", i)
        a.heavy == :early && i == 1 && record!(pop, "early", 1)
    end
    record!(pop, "end", i)
    m = finalize_measurements(acc)
    out = schedule_outcome(pop)
    # the single round stops on regrowth to N_det, which is not progression (N_prog)
    outcome = out == :progressed ? (a.kind == :single ? "regrown" : "progressed") :
              out == :extinct ? "extinct" : "censored"
    alllog = pop.log[nlog0+1:end]
    traj = [(t = q.t, m = q.m) for q in m.trajectory]
    rec = Dict{String,Any}(
        "outcome" => (; ctx.ids..., outcome = outcome,
                      ttp = outcome == "progressed" ? pop.t - t_det : missing,
                      t_end = pop.t - t_det, cum_dose = cum_dose(traj), n_cycles = n_cycles(alllog),
                      n_blocks = length(alllog)),
        "log" => [log_row(b) for b in alllog],
        "traj" => (t_rel = [q.t - t_det for q in m.trajectory], N = [q.N for q in m.trajectory],
                   m = [q.m for q in m.trajectory], nS = [q.n_sensitive for q in m.trajectory],
                   nR = [q.n_resistant for q in m.trajectory]),
        "counts" => [s.stats.counts for s in m.snapshots],
        "heavy" => heavy)
    example && (rec["clones"] = [(t_rel = s.t - t_det, clones = s.stats.clones) for s in m.snapshots])
    return rec
end

"""
    grow_tumour(c, t; block)

Grow tumour `t` of cohort `c` to detection (its own seed, so the same tumour every time). Returns
the detection population, the growth records at the sizes `c.growth_heavy` (fractions of N_det;
`idx` = size in percent of N_det, subsample seed by position, so that adding sizes leaves the
earlier records unchanged) and the detection record.
"""
function grow_tumour(c::Cohort, t::Integer; block::AbstractString)
    seed = tumour_seed(c, t); p = params(c, t); ndet = c.ndet
    ids = (block = String(block), cohort = c.name, tumour = Int(t))
    pop = initialize_population(1)
    growth = Any[]
    if isempty(c.growth_heavy)
        simulate!(pop, growth_block(p, ndet), Xoshiro(seed))
    else
        gctx = (t_det = 0.0, ids = merge(ids, (arm = "growth",)))
        j = Ref(0)
        acc = MeasurementAccumulator(MeasurementSpec(
            triggers = [AtSize(round(Int, f * ndet)) for f in sort(c.growth_heavy)],
            stats = [:heavy => q -> (j[] += 1; measure_heavy(q, gctx, derive_seed(seed, 2, 0, j[]); tag = "growth",
                                                             idx = round(Int, 100 * popsize(q) / ndet)))]))
        simulate!(pop, growth_block(p, ndet), Xoshiro(seed); accumulator = acc)
        growth = [s.stats.heavy for s in finalize_measurements(acc).snapshots]
    end
    t_det = pop.t
    growth = [merge(g, (t_rel = g.t_rel - t_det, B_pre = g.B, B_new = 0)) for g in growth]
    detection = measure_heavy(pop, (t_det = t_det, ids = merge(ids, (arm = "detection",))),
                              derive_seed(seed, 2, 0, 0); tag = "detection", idx = 0)
    return pop, growth, detection
end

"""
    run_tumour(c, t; block, test = false, example_dir = nothing)

Grow tumour `t` of cohort `c` to detection, record distributions during growth and at detection,
and run every arm. Example tumours (t ≤ N_EXAMPLES) keep clone records over time and, if
`example_dir` is given, their detection population is saved there.
"""
function run_tumour(c::Cohort, t::Integer; block::AbstractString, test::Bool = false, example_dir = nothing)
    seed = tumour_seed(c, t); p = params(c, t); ndet = c.ndet
    ids = (block = String(block), cohort = c.name, tumour = Int(t))
    pop, growth, detection = grow_tumour(c, t; block)
    t_det = pop.t
    example = t <= N_EXAMPLES
    example && example_dir !== nothing && save_population(joinpath(example_dir, @sprintf("t%04d.cepop", t)), pop)
    arms = Dict{String,Any}(a.id => run_arm(pop, a, p, ndet, seed, t_det, ids; arm_index = j, example)
                            for (j, a) in enumerate(c.arms))
    meta = (cohort = c.name, key = c.key, tumour = Int(t), seed = seed, nu = p.nu, alpha = string(p.alpha),
            fK = p.fK, cost = p.cost, d0 = p.d0, ndet = ndet, t_det = t_det,
            package = string(pkgversion(ClonalEvolutionTx)), julia = string(VERSION), test = test)
    return Dict{String,Any}("meta" => meta, "growth" => growth, "detection" => detection, "arms" => arms)
end

# Nested dictionaries become JLD2 groups, so a reader can load e.g. `arms/mtd/outcome` alone.
function _write!(f, prefix, v::AbstractDict)
    for (k, x) in v
        _write!(f, isempty(prefix) ? string(k) : prefix * "/" * string(k), x)
    end
end
_write!(f, path, v) = (f[path] = v)

"Write a tumour record atomically (`.tmp`, then rename); nested dictionaries become groups."
function save_tumour(path::AbstractString, rec::AbstractDict)
    tmp = path * ".tmp"
    jldopen(tmp, "w"; compress = true) do f
        _write!(f, "", rec)
    end
    mv(tmp, path; force = true)
    return path
end

tumour_file(dir, t) = joinpath(dir, @sprintf("t%04d.jld2", t))

"""
    run_block(block; cohorts = nothing, test = false)

Simulate all (or the named) cohorts of a block on the available worker processes, skipping
tumours whose file exists. `test = true` runs 3 small tumours per cohort into `00_test_<block>`.
"""
function run_block(block::AbstractString; cohorts = nothing, test::Bool = false)
    bname = test ? "00_test_$block" : String(block)
    list = [c for c in BLOCKS[block] if cohorts === nothing || c.name in cohorts]
    for c0 in list
        c = test ? test_cohort(c0) : c0
        dir = rawdir(bname, c.name); exdir = mkpath(joinpath(dir, "examples"))
        open(joinpath(dir, "meta.toml"), "w") do io
            TOML.print(io, Dict("block" => block, "cohort" => c.name, "key" => c.key, "n" => c.n,
                "nu" => string(c.nu), "alpha" => string(c.alpha), "fK" => string(c.fK), "cost" => c.cost,
                "d0" => c.d0, "ndet" => c.ndet, "arms" => [a.id for a in c.arms],
                "arm_tmax" => [a.tmax for a in c.arms], "base_seed" => string(BASE_SEED),
                "package" => string(pkgversion(ClonalEvolutionTx)), "julia" => string(VERSION),
                "date" => string(now()), "test" => test))
        end
        todo = [t for t in 1:c.n if !isfile(tumour_file(dir, t))]
        println("[$bname/$(c.name)] $(length(todo)) of $(c.n) tumours to run on $(max(1, nworkers_safe())) workers")
        t_start = time()
        pmap(todo) do t
            t1 = time()
            rec = run_tumour(c, t; block = bname, test, example_dir = exdir)
            save_tumour(tumour_file(dir, t), rec)
            println(@sprintf("[%s/%s] t%04d done in %.1f s (worker %d)", bname, c.name, t, time() - t1, myid()))
            t
        end
        println(@sprintf("[%s/%s] finished in %.1f min", bname, c.name, (time() - t_start) / 60))
    end
    return nothing
end

nworkers_safe() = try; Distributed.nworkers(); catch; 1; end

"""
    regrow_tumour!(path, c, t; block)

Replace the growth records of the stored tumour `t` (file `path`) by those at the current sizes
`c.growth_heavy`, regrowing it from its seed; the treatment arms are kept. The regrown detection
record and every growth record stored before (apart from its `idx`) must equal the stored ones,
otherwise the file is left unchanged. Returns `:done`, `:skipped` (records already present) or
`:mismatch`.
"""
function regrow_tumour!(path::AbstractString, c::Cohort, t::Integer; block::AbstractString)
    want = sort([round(Int, 100f) for f in c.growth_heavy])
    old = jldopen(f -> (growth = f["growth"], detection = f["detection"]), path, "r")
    sort([g.idx for g in old.growth]) == want && return :skipped
    _, growth, detection = grow_tumour(c, t; block)
    isequal(detection, old.detection) || return :mismatch
    for g in old.growth
        i = findfirst(h -> h.N == g.N, growth)
        (i !== nothing && isequal(merge(g, (idx = growth[i].idx,)), growth[i])) || return :mismatch
    end
    rec = load(path)                                    # every dataset, keyed by its full path
    rec["growth"] = growth
    save_tumour(path, rec)
    return :done
end

"""
    regrow_block(block; cohorts = nothing, test = false)

`regrow_tumour!` for every stored tumour of the (named) cohorts of a block with growth records,
on the available worker processes; resumable. Throws if any tumour is not reproduced.
"""
function regrow_block(block::AbstractString; cohorts = nothing, test::Bool = false)
    bname = test ? "00_test_$block" : String(block)
    list = [c for c in BLOCKS[block] if (cohorts === nothing || c.name in cohorts) && !isempty(c.growth_heavy)]
    bad = String[]
    for c0 in list
        c = test ? test_cohort(c0) : c0
        dir = rawdir(bname, c.name)
        todo = [t for t in 1:c.n if isfile(tumour_file(dir, t))]
        println("[$bname/$(c.name)] regrowing $(length(todo)) of $(c.n) tumours (sizes $(c.growth_heavy) N_det) on $(max(1, nworkers_safe())) workers")
        t_start = time()
        res = pmap(todo) do t
            r = regrow_tumour!(tumour_file(dir, t), c, t; block = bname)
            println(@sprintf("[%s/%s] t%04d %s (worker %d)", bname, c.name, t, r, myid()))
            r
        end
        append!(bad, ["$(c.name)/t$(lpad(t, 4, '0'))" for (t, r) in zip(todo, res) if r == :mismatch])
        println(@sprintf("[%s/%s] %d done, %d skipped, %d not reproduced in %.1f min", bname, c.name,
                         count(==(:done), res), count(==(:skipped), res), count(==(:mismatch), res),
                         (time() - t_start) / 60))
    end
    isempty(bad) || error("not reproduced (files unchanged): ", join(bad, ", "))
    return nothing
end

"Cohort names from the command line; `TEST=1` in the environment selects the smoke test."
parse_args(args) = (cohorts = isempty(args) ? nothing : collect(String, args), test = get(ENV, "TEST", "0") == "1")

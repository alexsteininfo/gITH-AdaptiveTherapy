# Processing of raw tumour files into tables (data/processed/<block>/):
#   collect_block     outcomes.csv (one row per tumour × arm) and counts.csv.gz (time series)
#   statistics_block  snapshots.csv (entire tree) and samples.csv (subsamples × f_min), computed
#                     from the stored distributions; re-runnable with other k_min.

using CSV, DataFrames

"Registry block of a (possibly test) block name."
registry_block(b::AbstractString) = replace(b, r"^00_test_" => "")

"All tumour files of a block as `(cohort, path)`."
function tumour_files(block::AbstractString)
    root = joinpath(ROOT, "data", "raw", block)
    out = Tuple{String,String}[]
    for c in sort(readdir(root))
        dir = joinpath(root, c); isdir(dir) || continue
        for f in sort(readdir(dir))
            occursin(r"^t\d{4}\.jld2$", f) && push!(out, (c, joinpath(dir, f)))
        end
    end
    return out
end

"Arm ids stored in an open tumour file (none for growth-only blocks)."
arm_ids(f) = haskey(f, "arms") ? collect(keys(f["arms"])) : String[]

"vcat of data frames that also works for an empty list."
vcat_frames(dfs) = isempty(dfs) ? DataFrame() : vcat(dfs...; cols = :union)

"Kind, level and follow-up `tmax` of arm `id` of `cohort` in a (possibly test) block."
function arm_info(block, cohort, id)
    c = only(filter(x -> x.name == cohort, BLOCKS[registry_block(block)]))
    startswith(block, "00_test_") && (c = test_cohort(c))
    a = only(filter(x -> x.id == id, c.arms))
    return (kind = String(a.kind), level = a.level, tmax = a.tmax)
end

"""
Resistant cells, Hill numbers, largest clone, and resistant lineages (resistance arisen in a
sensitive cell) that arose before (`res_pre`) and after detection (`res_new`), from a clone record.
"""
function clone_stats(clones, N, t_det)
    res = clones.alpha .< 1
    lineage = res .& (clones.parent .== FOUNDER)
    h = resistant_hill(clones.alpha, clones.size)
    pre = count(i -> lineage[i] && clones.origin_time[i] <= t_det, eachindex(res))
    new = count(i -> lineage[i] && clones.origin_time[i] > t_det, eachindex(res))
    return (res_frac = N == 0 ? 0.0 : sum(clones.size[res]; init = 0) / N, res_D0 = h.D0, res_D1 = h.D1,
            res_D2 = h.D2, res_Dinf = h.Dinf, res_pre = pre, res_new = new,
            largest = largest_frac(clones.alpha, clones.size, N))
end

hist_dict(h) = Dict{Int,Int}(Int(j) => Int(c) for (j, c) in zip(h.j, h.c))

function _collect_tumour(block, cohort, path)
    outs = DataFrame[]; cnts = DataFrame[]
    jldopen(path, "r") do f
        meta = f["meta"]; det = f["detection"]
        cs = clone_stats(det.tree.clones, det.N, meta.t_det)
        detf = (det_N = det.N, det_res_frac = cs.res_frac, det_n_res = cs.res_D0, det_largest = cs.largest,
                det_B = det.B, det_mean_burden = burden_stats(hist_dict(det.tree.scmb)).mean)
        for id in arm_ids(f)
            a = arm_info(block, cohort, id)
            o = f["arms/$id/outcome"]
            push!(outs, DataFrame([merge(o, (nu = meta.nu, alpha = meta.alpha, fK = meta.fK, cost = meta.cost,
                                              kind = a.kind, level = a.level, tmax = a.tmax), detf)]))
            c = f["arms/$id/counts"]
            isempty(c) || push!(cnts, DataFrame(c))
        end
    end
    return (vcat_frames(outs), vcat_frames(cnts))
end

"Collect outcomes and counts of all tumours of `block` into `outcomes.csv` and `counts.csv.gz`."
function collect_block(block::AbstractString)
    files = tumour_files(block)
    res = Distributed.pmap(x -> _collect_tumour(block, x...), files)
    dir = procdir(block)
    CSV.write(joinpath(dir, "outcomes.csv"), vcat_frames(first.(res)))
    CSV.write(joinpath(dir, "counts.csv.gz"), vcat_frames(last.(res)); compress = true)
    println("[$block] collected $(length(files)) tumours")
    return nothing
end

"One snapshot row (entire tree) from a stored distribution record."
function snapshot_row(r, t_det; kmin_frac)
    tree = r.tree; N = tree.N
    return (block = r.block, cohort = r.cohort, tumour = r.tumour, arm = r.arm, tag = r.tag, idx = r.idx,
            t_rel = r.t_rel, N = N, m = r.m, B = r.B, B_pre = r.B_pre, B_new = r.B_new, clonal = r.clonal,
            clone_stats(tree.clones, N, t_det)...,
            mean_burden = burden_stats(hist_dict(tree.scmb)).mean, var_burden = burden_stats(hist_dict(tree.scmb)).var,
            sfs_stats(expand_counts(tree.sfs, N); kmin_frac)...)
end

"Rows for every subsample of a record × detection threshold f_min (tail fit from k_min ≥ f_min·n)."
function sample_rows(r, t_det; kmin_frac, fmin = FMIN)
    rows = NamedTuple[]
    for n in sort(collect(keys(r.samples)))
        s = r.samples[n]
        sfs = expand_counts(s.sfs, n)
        b = burden_stats(hist_dict(s.scmb)); h = resistant_hill(s.clones.alpha, s.clones.size)
        for fm in fmin
            push!(rows, (block = r.block, cohort = r.cohort, tumour = r.tumour, arm = r.arm, tag = r.tag,
                         idx = r.idx, t_rel = r.t_rel, N = r.N, n = n, fmin = fm,
                         sfs_stats(apply_fmin(sfs, fm); kmin_frac = max(kmin_frac, fm))..., mean_burden = b.mean,
                         var_burden = b.var, res_D0 = h.D0, res_D2 = h.D2, res_Dinf = h.Dinf))
        end
    end
    return rows
end

function _statistics_tumour(path; kmin_frac)
    snaps = NamedTuple[]; samps = NamedTuple[]
    jldopen(path, "r") do f
        t_det = f["meta"].t_det
        recs = Any[f["growth"]...; f["detection"]]
        for id in arm_ids(f)
            append!(recs, f["arms/$id/heavy"])
        end
        for r in recs
            push!(snaps, snapshot_row(r, t_det; kmin_frac))
            append!(samps, sample_rows(r, t_det; kmin_frac))
        end
    end
    return (DataFrame(snaps), vcat_frames([DataFrame([s]) for s in samps]))
end

"""
    statistics_block(block; kmin_frac = 0.01, suffix = "")

Statistics of every stored distribution of `block`: `snapshots<suffix>.csv` (entire tree) and
`samples<suffix>.csv` (subsamples × f_min). Use `suffix` for experiments, e.g. other `kmin_frac`.
"""
function statistics_block(block::AbstractString; kmin_frac::Real = 0.01, suffix::AbstractString = "")
    files = tumour_files(block)
    res = Distributed.pmap(x -> _statistics_tumour(x[2]; kmin_frac), files)
    dir = procdir(block)
    CSV.write(joinpath(dir, "snapshots$suffix.csv"), vcat_frames(first.(res)))
    CSV.write(joinpath(dir, "samples$suffix.csv"), vcat_frames(last.(res)))
    println("[$block] statistics of $(length(files)) tumours (kmin_frac = $kmin_frac)")
    return nothing
end

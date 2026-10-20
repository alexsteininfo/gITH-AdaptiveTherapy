# What a simulation stores about a population: population counts (often) and distributions of the
# entire tree and of nested subsamples (at key time points). No statistics are computed here; see
# stats.jl and code/02_processing/02_statistics.jl.

const KIND_CODE = Dict(:founder => 0x00, :driver => 0x01, :resistance => 0x02, :custom => 0x03)

# Id of the founder clone. Resistance is the only mutation type in every cohort, so the founder is
# the only sensitive clone: a resistant clone whose parent is the founder starts a new resistant
# lineage; any other resistant clone is a second hit inside a resistant lineage.
const FOUNDER = 1

"Sparse representation `(k, c)` of the nonzero entries of a count vector indexed by k = 1, 2, …"
function sparse_counts(v::AbstractVector{<:Integer})
    k = Int32[]; c = Int64[]
    for (i, x) in enumerate(v)
        x == 0 && continue
        push!(k, i); push!(c, x)
    end
    return (k = k, c = c)
end

"Dense vector of length `N` from a sparse `(k, c)` record."
function expand_counts(sp, N::Integer)
    v = zeros(Int, N)
    for (k, c) in zip(sp.k, sp.c)
        v[k] = c
    end
    return v
end

"Burden histogram `(j, c)`: number of cells `c` carrying `j` neutral mutations."
function burden_hist(burdens::AbstractVector{<:Integer})
    h = Dict{Int,Int}()
    for b in burdens
        h[b] = get(h, b, 0) + 1
    end
    j = sort!(collect(keys(h)))
    return (j = Int32.(j), c = Int64[h[x] for x in j])
end

"Living clones: id, parent, α, s, size, origin time and kind (0 founder, 1 driver, 2 resistance, 3 custom)."
function clone_record(pop)
    ct = pop.clones; ids = ct.alive
    return (id = Int32.(ids), parent = ct.parent[ids], alpha = ct.α[ids], s = ct.s[ids],
            size = ct.size[ids], origin_time = ct.origin_time[ids],
            kind = UInt8[KIND_CODE[ct.kind[i]] for i in ids])
end

"Distributions of one population: SFS, scMB histogram and clone sizes."
dist_record(pop::TreePopulation) = (N = popsize(pop), sfs = sparse_counts(site_frequency_spectrum(pop)),
                                    scmb = burden_hist(mutations_per_cell(pop)), clones = clone_record(pop))

"""
    burden_split(pop, t_det)

Total mutational burden `B`, split into mutations on lineages born up to detection (`B_pre`,
including clonal mutations of the founder) and after it (`B_new`).
"""
function burden_split(pop::TreePopulation, t_det::Real)
    tree = pop.tree; pre = pop.clonal_mutations; new = 0
    for i in eachindex(tree.parent)
        tree.parent[i] == ClonalEvolutionTx.REMOVED && continue
        tree.birthtime[i] <= t_det ? (pre += tree.neutral[i]) : (new += tree.neutral[i])
    end
    return (B = pre + new, B_pre = pre, B_new = new)
end

"Mutations carried by all N cells: the SFS entry at k = N (0 if absent)."
clonal_from_sfs(sfs, N) = (i = findfirst(==(N), sfs.k); i === nothing ? 0 : Int(sfs.c[i]))

"""
    count_row(pop, ctx)

Population counts at one time: `N`, dose `m`, sensitive/resistant cells, number of resistant clones
`n_res`, resistant lineages arisen before (`res_pre`) and after detection (`res_new`, de novo),
largest resistant clone (fraction of N), mean burden.
"""
function count_row(pop, ctx)
    ct = pop.clones; N = popsize(pop)
    n_res = 0; pre = 0; new = 0; largest = 0
    for k in ct.alive
        ct.α[k] < 1 || continue
        n_res += 1
        largest = max(largest, ct.size[k])
        ct.parent[k] == FOUNDER || continue                 # second hit inside a resistant lineage
        ct.origin_time[k] <= ctx.t_det ? (pre += 1) : (new += 1)
    end
    return (; ctx.ids..., t_rel = pop.t - ctx.t_det, N = N, m = pop.m, nS = n_sensitive(pop),
            nR = n_resistant(pop), n_res = n_res, res_pre = pre, res_new = new,
            largest = N == 0 ? 0.0 : largest / N, mean_burden = N == 0 ? NaN : mean_burden(pop))
end

"""
    measure_heavy(pop, ctx, seed; tag, idx)

Distributions of the entire tree and of nested subsamples of `SAMPLE_N` cells (one seed for all
sizes; sizes larger than the population are skipped), plus the burden split.
"""
function measure_heavy(pop::TreePopulation, ctx, seed::Integer; tag::AbstractString, idx::Integer)
    samples = Dict{Int,Any}()
    for n in SAMPLE_N
        n <= popsize(pop) || continue
        samples[n] = dist_record(sample_population(pop, n; seed = seed))
    end
    tree = dist_record(pop)
    return (; ctx.ids..., tag = String(tag), idx = Int(idx), t_rel = pop.t - ctx.t_det,
            N = popsize(pop), m = pop.m, burden_split(pop, ctx.t_det)...,
            clonal = clonal_from_sfs(tree.sfs, tree.N), tree = tree, samples = samples)
end

# Expected neutral SFS and its Hill numbers under cyclic adaptive therapy (.docs/theory/math/
# ANALYSIS.md): deterministic logistic trajectory, master formula (3.1) for mutations born during
# growth, harmonic-mean law (§4.6) or the master formula for mutations born during cycling.
#   cyclic_path      growth from one cell to N_det, then ON (full dose) to N_low, OFF to N_det, …
#   sfs_master       expected counts per k-bin of the mutations born in [s0, s1] (master formula)
#   sfs_harmonic     the same for the mutations born during cycling by the harmonic-mean law
#   bin_hill         Hill numbers (q = 0, 1, 2, ∞; normalised frequencies) of an expected SFS
#   make_theory      theory_hill.csv, theory_sfs.csv of a dynamics block
#   make_theory_single  theory_single_hill.csv, theory_single_sfs.csv (single round, Figure 2)

# ── Deterministic trajectory ──────────────────────────────────────────────────
# Every phase is dN/dt = N(a − lN), so x = 1/N obeys x' = −a x + l and is explicit.
_phase_x(x0, a, l, u) = l / a + (x0 - l / a) * exp(-a * u)
_phase_intx(x0, a, l, u) = l / a * u + (x0 - l / a) * (-expm1(-a * u)) / a        # ∫₀ᵘ dv/N
_phase_time(N0, N1, a, l) = log((1 / N1 - l / a) / (1 / N0 - l / a)) / (-a)       # N0 → N1

"""
    cyclic_path(; ndet, fK, level, ncycles, b = B0, d0 = D0, ds = DS)

Mean-field trajectory of a tumour under holidays: logistic growth (death d0 + lN, l = (b − d0)fK/N_det)
from one cell to N_det, then `ncycles` cycles of full dose (death + ds) down to `level`·N_det and
regrowth to N_det. Phases `(t0, N0, a, T)` with net rate a − lN; time 0 is the founder cell.
"""
function cyclic_path(; ndet::Real, fK::Real, level::Real, ncycles::Integer, b::Real = B0, d0::Real = D0, ds::Real = DS)
    l = (b - d0) * fK / ndet; nlow = level * ndet
    ph = NTuple{4,Float64}[]; t = 0.0
    add!(N0, N1, a) = (T = _phase_time(N0, N1, a, l); push!(ph, (t, N0, a, T)); t += T)
    add!(1.0, ndet, b - d0)
    for _ in 1:ncycles
        add!(ndet, nlow, b - d0 - ds); add!(nlow, ndet, b - d0)
    end
    return (b = Float64(b), l = l, ndet = Float64(ndet), nlow = nlow, phases = ph)
end

"Size of the path at time `s`."
path_N(p, s) = (ph = p.phases[something(findlast(q -> q[1] <= s, p.phases), 1)];
                1 / _phase_x(1 / ph[2], ph[3], p.l, s - ph[1]))

"∫ₛᵗ du/N(u) along the path."
function path_intinv(p, s, t)
    tot = 0.0
    for (t0, N0, a, T) in p.phases
        lo = max(s, t0); hi = min(t, t0 + T)
        hi > lo && (tot += _phase_intx(1 / N0, a, p.l, hi - t0) - _phase_intx(1 / N0, a, p.l, lo - t0))
    end
    return tot
end

"Detection time; end of the ON (`lower`) or OFF (`upper`) phase of `cycle` (cycle 0, upper: detection)."
path_time(p, cycle, threshold) =
    cycle == 0 ? p.phases[1][4] : (ph = p.phases[threshold == "lower" ? 2cycle : 2cycle + 1]; ph[1] + ph[4])

"Time at which the growth phase reaches size `N`."
path_growth_time(p, N) = _phase_time(1.0, N, p.phases[1][3], p.l)

# ── Expected SFS on k-bins ────────────────────────────────────────────────────
"Bins `(lo, hi)` covering 1..K: single k up to 2000, then 200 per decade (for Hill numbers)."
function fine_bins(K::Integer)
    e = unique([collect(1:min(K, 2000)); round.(Int, 10 .^ (log10(2001):1/200:log10(K)))])
    hi = [e[2:end] .- 1; K]
    keep = e .<= hi
    return e[keep], hi[keep]
end

"Σ_{k=lo}^{hi} (1 − 1/W)^{k−1}, stable for W ≫ 1."
_geosum(W, lo, hi) = W <= 1 + 1e-12 ? Float64(lo == 1) :
    (L = log1p(-1 / W); exp((lo - 1) * L) * (-expm1((hi - lo + 1) * L)) * W)

"""
    sfs_master(p, t, lo, hi; s0, s1, m = MU, n = 3000)

Expected number of mutations per bin `lo[j]:hi[j]` at time `t`, of those born in [s0, s1]
(ANALYSIS eq. 3.1 written in W: Σ_k S_k = 2m ∫ N(s(W)) W⁻² (1 − 1/W)^{k−1} dW, with
W(s, t) = 1 + bN(t)∫ₛᵗ du/N); quadrature on `n` points uniform in ln W. The bins end at the
population size K = N(t); the geometric lineage sizes beyond K (the oldest, clonal mutations) are
counted in the last bin, so that the burden Σ k S_k = 2m N(t) ∫ b ds holds.
"""
function sfs_master(p, t, lo, hi; s0::Real, s1::Real, m::Real = MU, n::Integer = 3000)
    Nt = path_N(p, t); K = hi[end]
    Wof(s) = 1 + p.b * Nt * path_intinv(p, s, t)
    lw = range(log(Wof(s1)), log(Wof(s0)); length = n); dl = step(lw)
    C = zeros(length(lo))
    for (i, L) in enumerate(lw)
        a, z = Float64(s0), Float64(s1)              # s with W(s) = e^L (W decreases in s)
        for _ in 1:60
            mid = (a + z) / 2
            Wof(mid) > exp(L) ? (a = mid) : (z = mid)
        end
        W = exp(L)
        c = 2m * path_N(p, (a + z) / 2) / W * dl * ((i == 1 || i == n) ? 0.5 : 1.0)
        @inbounds for j in eachindex(lo)
            C[j] += c * _geosum(W, lo[j], hi[j])
        end
        W > 1 + 1e-12 && (C[end] += c * W * exp(K * log1p(-1 / W)))        # Σ_{k>K}(1 − 1/W)^{k−1}
    end
    return C
end

"""
    sfs_harmonic(p, t, lo, hi; m = MU)

Expected counts per bin of the mutations born during cycling (detection to `t`) by the
harmonic-mean law (ANALYSIS §4.6): S_k = 2m N_harm k⁻¹ (b_eff τ/(1 + b_eff τ))^k with τ the cycling
time, N_harm = τ / ∫ du/N and b_eff = b N(t)/N_harm. Returns `(C, tau, N_harm, b_eff)`.
"""
function sfs_harmonic(p, t, lo, hi; m::Real = MU)
    tdet = path_time(p, 0, "upper"); τ = t - tdet
    C = zeros(length(lo))
    τ <= 0 && return (C = C, tau = 0.0, N_harm = NaN, b_eff = NaN)
    Nh = τ / path_intinv(p, tdet, t); beff = p.b * path_N(p, t) / Nh; x = beff * τ / (1 + beff * τ)
    j = 1
    for k in 1:hi[end]
        while hi[j] < k; j += 1; end
        v = 2m * Nh / k * x^k
        v < 1e-12 && break
        C[j] += v
    end
    return (C = C, tau = τ, N_harm = Nh, b_eff = beff)
end

"""
    bin_hill(C, lo, hi, K)

Hill numbers of an expected SFS given as counts per bin (bin centre for k), with the
normalised frequencies of `sfs_hill_number` (p = k/Σ jS_j; ∞D with the largest class K).
"""
function bin_hill(C, lo, hi, K)
    kb = (lo .+ hi) ./ 2; M = sum(C .* kb); p = kb ./ M; nz = C .> 0
    return (S = sum(C), D1n = exp(-sum(C[nz] .* p[nz] .* log.(p[nz]))), D2n = 1 / sum(C .* p .^ 2), Dinfn = M / K)
end

# ── Tables ────────────────────────────────────────────────────────────────────
"""
    make_theory(block)

Expected neutral SFS and Hill numbers under the holiday arm of a dynamics block (its `neutral`
cohort; test blocks use the test sizes): at the untreated sizes (detection; growth at the lower
threshold) and at both thresholds of every stored cycle, for mutations born during cycling by the
harmonic-mean law (`model = harmonic`) or the master formula (`master`). Writes
`theory_hill.csv` (`kind` untreated / treated, `threshold, cycle, t_rel, N, tau, N_harm, b_eff, kc`
= b_eff τ, `model`, `S, D1n, D2n, Dinfn` and `*_rel` = treated / untreated at the same size) and
`theory_sfs.csv` (same keys, log-binned as `sfs_mean.csv`: `k, S_old, S_new, S_mean`).
"""
function make_theory(block::AbstractString)
    c = only(filter(x -> x.name == "neutral", BLOCKS[registry_block(block)]))
    startswith(block, "00_test_") && (c = test_cohort(c))
    a = only(filter(x -> x.kind == :hol, c.arms))
    ncyc = HEAVY_BLOCKS ÷ 2
    p = cyclic_path(; ndet = c.ndet, fK = c.fK, level = a.level, ncycles = ncyc, d0 = c.d0)
    tdet = path_time(p, 0, "upper")
    points = vcat([(kind = "untreated", threshold = "upper", cycle = 0, t = tdet),
                   (kind = "untreated", threshold = "lower", cycle = 0, t = path_growth_time(p, p.nlow))],
                  [(kind = "treated", threshold = th, cycle = n, t = path_time(p, n, th)) for n in 1:ncyc for th in ("lower", "upper")])
    hills = NamedTuple[]; sfs = NamedTuple[]
    for pt in points, model in ("harmonic", "master")
        K = round(Int, path_N(p, pt.t))
        born = min(pt.t, tdet)                          # untreated: everything is "old"
        flo, fhi = fine_bins(K)
        r = logbin_ranges(K)
        old_f = sfs_master(p, pt.t, flo, fhi; s0 = 0.0, s1 = born); old_b = sfs_master(p, pt.t, r.lo, r.hi; s0 = 0.0, s1 = born)
        if model == "harmonic"
            hf = sfs_harmonic(p, pt.t, flo, fhi); new_f = hf.C; new_b = sfs_harmonic(p, pt.t, r.lo, r.hi).C
        else
            hf = sfs_harmonic(p, pt.t, flo, fhi)        # τ, N_harm, b_eff only
            new_f, new_b = pt.t > tdet ? (sfs_master(p, pt.t, flo, fhi; s0 = tdet, s1 = pt.t),
                                          sfs_master(p, pt.t, r.lo, r.hi; s0 = tdet, s1 = pt.t)) :
                                         (zeros(length(flo)), zeros(length(r.lo)))
        end
        h = bin_hill(old_f .+ new_f, flo, fhi, K)
        push!(hills, (; pt.kind, pt.threshold, pt.cycle, t_rel = pt.t - tdet, N = K, hf.tau, hf.N_harm, hf.b_eff,
                      kc = hf.b_eff * hf.tau, model, h...))
        n = r.hi .- r.lo .+ 1
        for j in eachindex(r.lo)
            push!(sfs, (; pt.kind, pt.threshold, pt.cycle, model, k = r.k[j], S_old = old_b[j] / n[j],
                        S_new = new_b[j] / n[j], S_mean = (old_b[j] + new_b[j]) / n[j]))
        end
    end
    h = DataFrame(hills)
    for col in (:S, :D1n, :D2n, :Dinfn)
        base = Dict((r.threshold, r.model) => r[col] for r in eachrow(h) if r.kind == "untreated")
        h[!, Symbol(col, "_rel")] = [r[col] / base[(r.threshold, r.model)] for r in eachrow(h)]
    end
    CSV.write(joinpath(procdir(block), "theory_hill.csv"), h)
    CSV.write(joinpath(procdir(block), "theory_sfs.csv"), DataFrame(sfs))
    return h
end

"""
    make_theory_single(block)

Expected neutral SFS and Hill numbers for the single round of a dynamics block (its `neutral` cohort;
test blocks use the test sizes): full dose from detection down to 0.1 N_det (the single-round
threshold in `build.jl`), then regrowth to N_det, all by the master formula. Untreated at the same
sizes (growth at 0.1 N_det, detection) and treated after treatment and after regrowth. Writes
`theory_single_hill.csv` (`kind` untreated / treated, `threshold` lower (0.1 N_det) / upper (N_det),
`t_rel, N, S, D1n, D2n, Dinfn` and `*_rel` = treated / untreated at the same size) and
`theory_single_sfs.csv` (same keys, log-binned as `sfs_mean.csv`: `k, S_old, S_new, S_mean`).
"""
function make_theory_single(block::AbstractString)
    c = only(filter(x -> x.name == "neutral", BLOCKS[registry_block(block)]))
    startswith(block, "00_test_") && (c = test_cohort(c))
    p = cyclic_path(; ndet = c.ndet, fK = c.fK, level = 0.1, ncycles = 1, d0 = c.d0)
    tdet = path_time(p, 0, "upper")
    points = [(kind = "untreated", threshold = "lower", t = path_growth_time(p, p.nlow)),
              (kind = "untreated", threshold = "upper", t = tdet),
              (kind = "treated", threshold = "lower", t = path_time(p, 1, "lower")),
              (kind = "treated", threshold = "upper", t = path_time(p, 1, "upper"))]
    hills = NamedTuple[]; sfs = NamedTuple[]
    for pt in points
        K = round(Int, path_N(p, pt.t)); born = min(pt.t, tdet)
        flo, fhi = fine_bins(K); r = logbin_ranges(K)
        master(lo, hi, s0, s1) = s1 > s0 ? sfs_master(p, pt.t, lo, hi; s0, s1) : zeros(length(lo))
        old_f = master(flo, fhi, 0.0, born); new_f = master(flo, fhi, tdet, pt.t)
        old_b = master(r.lo, r.hi, 0.0, born); new_b = master(r.lo, r.hi, tdet, pt.t)
        push!(hills, (; pt.kind, pt.threshold, t_rel = pt.t - tdet, N = K, bin_hill(old_f .+ new_f, flo, fhi, K)...))
        n = r.hi .- r.lo .+ 1
        for j in eachindex(r.lo)
            push!(sfs, (; pt.kind, pt.threshold, k = r.k[j], S_old = old_b[j] / n[j], S_new = new_b[j] / n[j],
                        S_mean = (old_b[j] + new_b[j]) / n[j]))
        end
    end
    h = DataFrame(hills)
    for col in (:S, :D1n, :D2n, :Dinfn)
        base = Dict(r.threshold => r[col] for r in eachrow(h) if r.kind == "untreated")
        h[!, Symbol(col, "_rel")] = [r[col] / base[r.threshold] for r in eachrow(h)]
    end
    CSV.write(joinpath(procdir(block), "theory_single_hill.csv"), h)
    CSV.write(joinpath(procdir(block), "theory_single_sfs.csv"), DataFrame(sfs))
    return h
end

# Expected neutral SFS and its Hill numbers under cyclic adaptive therapy (.docs/theory/math/
# ANALYSIS.md): deterministic logistic trajectory, master formula (3.1) for mutations born during
# growth, harmonic-mean law (§4.6) or the master formula for mutations born during cycling.
#   cyclic_path      growth from one cell to N_det, then ON (full dose) to N_low, OFF to N_det, …
#   sfs_master       expected counts per k-bin of the mutations born in [s0, s1] (master formula)
#   sfs_harmonic     the same for the mutations born during cycling by the harmonic-mean law
#   bin_hill         Hill numbers (q = 0, 1, 2, ∞; normalised frequencies) of an expected SFS
#   make_theory      theory_hill.csv, theory_sfs.csv of a dynamics block
#   make_theory_single  theory_single_hill.csv, theory_single_sfs.csv (single round, Figure 2)
#   stabilising_path growth from one cell to N_det, full dose to N* = level·N_det, hold at N*
#   stab_richness    richness and singletons at the hold for one stable level (cheap exact bins)
#   make_theory_stabilising  theory_stab_hill.csv, theory_stab_sfs.csv, theory_stab_levels.csv of a stabilising block (Figure 4)

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

# ── Stabilising adaptive therapy ──────────────────────────────────────────────
const STAB_TAUS = [0:0.05:1; 1.1:0.1:3; 3.25:0.25:10; 10.5:0.5:40]     # hold times of theory_stab_hill.csv
const STAB_SFS_TAUS = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0]   # … of theory_stab_sfs.csv
const STAB_SFS_DET_TIMES = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0]              # … and at these times since detection (hold time t − t_red; the fine records of Fig 4)
const STAB_LEVELS = 0.3:0.05:1.0                                           # stable sizes of theory_stab_levels.csv
const STAB_LEVEL_TAUS = [0:0.1:3; 3.5:0.5:10; 12:2:40]                     # hold times of the level scan

"""
    stabilising_path(; ndet, fK, level, tau, b = B0, d0 = D0, ds = DS)

Mean-field trajectory under stabilising AT: logistic growth from one cell to N_det, full dose
(death + ds) down to the stable size N* = `level`·N_det (no reduction phase for level 1), then
the hold at N* (net rate 0: the dose makes the death rate b) for `tau` time units. Same fields as
`cyclic_path` (`nlow` = N*).
"""
function stabilising_path(; ndet::Real, fK::Real, level::Real, tau::Real, b::Real = B0, d0::Real = D0, ds::Real = DS)
    l = (b - d0) * fK / ndet; nstab = level * ndet
    ph = NTuple{4,Float64}[]; t = 0.0
    add!(N0, N1, a) = (T = _phase_time(N0, N1, a, l); push!(ph, (t, N0, a, T)); t += T)
    add!(1.0, ndet, b - d0)
    level < 1 && add!(ndet, nstab, b - d0 - ds)
    push!(ph, (t, nstab, l * nstab, Float64(tau)))          # a − lN* = 0: constant size
    return (b = Float64(b), l = l, ndet = Float64(ndet), nlow = nstab, phases = ph)
end

"Detection time and start of the hold of a stabilising path."
stab_times(p) = (t_det = p.phases[1][4], t_hold = p.phases[end][1])

"""
    stab_richness(level, taus; ndet, fK, d0 = D0)

Richness `S` and singletons `S1` of a neutral tumour at hold times `taus` for the stable level
`level` (growth to N_det, full dose to N* = level·N_det, hold), split into lineages from before
the hold (`_pre`) and born during it (`_new`), with the doubletons `S2_pre`, the detection values
`S_det, S1_det` and the duration of the reduction `t_red`; exact per-bin sums of the master formula
on the log bins (cheap).
"""
function stab_richness(level, taus; ndet, fK, d0 = D0)
    p = stabilising_path(; ndet, fK, level, tau = maximum(taus) + 1, d0); st = stab_times(p)
    K = round(Int, level * ndet); r = logbin_ranges(K)
    rd = logbin_ranges(ndet); det = sfs_master(p, st.t_det, rd.lo, rd.hi; s0 = 0.0, s1 = st.t_det)
    rows = NamedTuple[]
    for τ in taus
        t = st.t_hold + τ
        pre = sfs_master(p, t, r.lo, r.hi; s0 = 0.0, s1 = st.t_hold)
        new = τ > 0 ? sfs_master(p, t, r.lo, r.hi; s0 = st.t_hold, s1 = t) : zeros(length(r.lo))
        push!(rows, (level = level, tau = τ, S = sum(pre) + sum(new), S_pre = sum(pre), S_new = sum(new),
                     S1 = pre[1] + new[1], S1_pre = pre[1], S1_new = new[1], S2_pre = pre[2], S_det = sum(det), S1_det = det[1],
                     t_red = st.t_hold - st.t_det))
    end
    return DataFrame(rows)
end

"""
    stab_level_scan(; ndet, fK, d0 = D0, levels = STAB_LEVELS, taus = STAB_LEVEL_TAUS)

Summary of the richness transient per stable level (`.docs/theory/STABILISING.md` §2.3): `r0` =
S₁/2μN* at the start of the hold (> 1: a minimum at the hold), `S0_rel` = S after the reduction /
S at detection, `tau_min` (hold time of the minimum, parabolic refinement), `t_min` = t_red + tau_min
(since detection), `Smin_rel` = S_min / S_det, `tau_approx` = (S₁(0) − 2μN*)/(b(S₁(0) − 2S₂(0))),
`tau_recover` (S back at S_det; NaN if not within the grid), `t_red`.
"""
function stab_level_scan(; ndet, fK, d0 = D0, levels = STAB_LEVELS, taus = STAB_LEVEL_TAUS, m = MU)
    rows = NamedTuple[]
    for L in levels
        d = stab_richness(L, Float64.(taus); ndet, fK, d0)
        ex = 2m * L * ndet; r0 = d.S1[1] / ex
        i = argmin(d.S); τm = d.tau[i]
        if 1 < i < nrow(d)
            x1, x2, x3 = d.tau[i-1:i+1]; y1, y2, y3 = d.S[i-1:i+1]
            A = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1); Bq = (y2 - y1) / (x2 - x1) - A * (x1 + x2)
            A > 0 && (τm = -Bq / 2A)
        end
        j = findfirst(d.S .>= d.S_det[1])
        τrec = j === nothing ? NaN : j == 1 ? 0.0 : d.tau[j-1] + (d.S_det[1] - d.S[j-1]) * (d.tau[j] - d.tau[j-1]) / (d.S[j] - d.S[j-1])
        push!(rows, (level = L, r0 = r0, S0_rel = d.S[1] / d.S_det[1], tau_min = τm, t_min = d.t_red[1] + τm,
                     Smin_rel = d.S[i] / d.S_det[1], tau_approx = r0 > 1 ? (d.S1_pre[1] - ex) / (B0 * (d.S1_pre[1] - 2d.S2_pre[1])) : 0.0,
                     tau_recover = τrec, t_red = d.t_red[1]))
    end
    return DataFrame(rows)
end

"""
    make_theory_stabilising(block; taus = STAB_TAUS, sfs_taus = STAB_SFS_TAUS)

Expected neutral SFS and Hill numbers under stabilising AT (its `neutral` cohort, one set of rows
per `mod` level; test blocks use the test sizes): the master formula on the deterministic path
(growth → full dose to N* → hold), mutations split by the phase they were born in. Writes
`theory_stab_hill.csv` (`level`; `kind` detection / untreated (growth at N*) / treated; `tau` time
at the stable size; `t_rel, N`; `S, S_old, S_red, S_new` (born during growth / reduction / hold);
singletons `S1, S1_pre, S1_new` (lineages from before the hold / born during it); `D1n, D2n,
Dinfn`; `*_rel_det` = treated / detection and `*_rel_size` = treated / untreated at N*) and
`theory_stab_sfs.csv` (`level, kind, tau`, log-binned `k, S_old, S_new, S_mean` as `sfs_mean.csv`;
`S_new` = born after detection; at the hold times `STAB_SFS_TAUS` and at the times since detection
`STAB_SFS_DET_TIMES`) and `theory_stab_levels.csv` (`stab_level_scan` over `STAB_LEVELS`).
"""
function make_theory_stabilising(block::AbstractString; taus = STAB_TAUS, sfs_taus = STAB_SFS_TAUS)
    c = only(filter(x -> x.name == "neutral", BLOCKS[registry_block(block)]))
    startswith(block, "00_test_") && (c = test_cohort(c))
    levels = unique([a.level for a in c.arms if a.kind == :mod])
    tmax = maximum(a.tmax for a in c.arms if a.kind == :mod)
    taus = filter(<=(tmax), Float64.(taus)); sfs_taus = filter(<=(tmax), Float64.(sfs_taus))
    hills = NamedTuple[]; sfs = NamedTuple[]
    for L in levels
        p = stabilising_path(; ndet = c.ndet, fK = c.fK, level = L, tau = tmax, d0 = c.d0)
        st = stab_times(p)
        extra = [t - (st.t_hold - st.t_det) for t in STAB_SFS_DET_TIMES if 0 <= t - (st.t_hold - st.t_det) <= tmax]
        ltaus = sort(unique([taus; extra])); lsfs = Set([sfs_taus; extra])
        points = [(kind = "detection", tau = NaN, t = st.t_det),
                  (kind = "untreated", tau = NaN, t = path_growth_time(p, p.nlow))]
        append!(points, [(kind = "treated", tau = τ, t = st.t_hold + τ) for τ in ltaus])
        for pt in points
            K = round(Int, path_N(p, pt.t))
            master(lo, hi, s0, s1) = s1 > s0 ? sfs_master(p, pt.t, lo, hi; s0, s1) : zeros(length(lo))
            treated = pt.kind == "treated"
            split(lo, hi) = treated ? (master(lo, hi, 0.0, st.t_det), master(lo, hi, st.t_det, st.t_hold), master(lo, hi, st.t_hold, pt.t)) :
                                      (master(lo, hi, 0.0, pt.t), zeros(length(lo)), zeros(length(lo)))
            flo, fhi = fine_bins(K)
            old, red, new = split(flo, fhi); tot = old .+ red .+ new
            h = bin_hill(tot, flo, fhi, K)
            push!(hills, (; level = L, pt.kind, pt.tau, t_rel = pt.t - st.t_det, N = K, S = h.S, S_old = sum(old),
                          S_red = sum(red), S_new = sum(new), S1 = tot[1], S1_pre = old[1] + red[1], S1_new = new[1],
                          h.D1n, h.D2n, h.Dinfn))
            if !treated || pt.tau in lsfs
                r = logbin_ranges(K); n = r.hi .- r.lo .+ 1
                old_b, red_b, new_b = split(r.lo, r.hi)
                for j in eachindex(r.lo)
                    push!(sfs, (; level = L, pt.kind, pt.tau, k = r.k[j], S_old = old_b[j] / n[j],
                                S_new = (red_b[j] + new_b[j]) / n[j], S_mean = (old_b[j] + red_b[j] + new_b[j]) / n[j]))
                end
            end
        end
    end
    h = DataFrame(hills)
    for col in (:S, :D1n, :D2n, :Dinfn)
        det = Dict(r.level => r[col] for r in eachrow(h) if r.kind == "detection")
        unt = Dict(r.level => r[col] for r in eachrow(h) if r.kind == "untreated")
        h[!, Symbol(col, "_rel_det")] = [r[col] / det[r.level] for r in eachrow(h)]
        h[!, Symbol(col, "_rel_size")] = [r[col] / unt[r.level] for r in eachrow(h)]
    end
    CSV.write(joinpath(procdir(block), "theory_stab_hill.csv"), h)
    CSV.write(joinpath(procdir(block), "theory_stab_sfs.csv"), DataFrame(sfs))
    CSV.write(joinpath(procdir(block), "theory_stab_levels.csv"), stab_level_scan(; ndet = c.ndet, fK = c.fK, d0 = c.d0))
    return h
end

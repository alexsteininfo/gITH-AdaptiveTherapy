# Statistics on plain vectors: diversity of resistant clones and passenger mutations, power-law
# exponents, detection thresholds, burdens, log-binning, paired tests and survival summaries.
# Used by the processing step; nothing here touches a population.

using StatsBase: quantile
using HypothesisTests: SignedRankTest, pvalue

const QS = (0, 1, 2, Inf)

"Hill numbers ^qD (q = 0, 1, 2, ∞) of the resistant clones (α < 1); `D0 = 0` and NaN if none."
function resistant_hill(alpha::AbstractVector, size::AbstractVector)
    sizes = [s for (a, s) in zip(alpha, size) if a < 1 && s > 0]
    isempty(sizes) && return (D0 = 0.0, D1 = NaN, D2 = NaN, Dinf = NaN)
    return (D0 = Float64(length(sizes)), D1 = hill_number(sizes, 1), D2 = hill_number(sizes, 2),
            Dinf = hill_number(sizes, Inf))
end

"Frequency of the largest resistant clone (0 if none)."
largest_frac(alpha, size, N) =
    N == 0 ? 0.0 : maximum((s for (a, s) in zip(alpha, size) if a < 1); init = 0) / N

"Detection threshold: drop mutations carried by fewer than `fmin` of the cells."
function apply_fmin(sfs::AbstractVector{<:Integer}, fmin::Real)
    out = collect(sfs)
    kcut = ceil(Int, fmin * length(sfs))
    out[1:min(kcut - 1, length(out))] .= 0
    return out
end

"""
    powerlaw_fit(sfs; kmin)

Approximate discrete MLE of the tail exponent (Clauset et al. 2009, eq. 3.7):
γ̂ = 1 + n [Σᵢ ln(kᵢ/(kmin − ½))]⁻¹ over all mutations with kmin ≤ kᵢ < N, se = (γ̂ − 1)/√n.
Mutations carried by all N cells (clonal) are not part of the subclonal tail and are left out.
NaN if fewer than 10 mutations lie in the tail.
"""
function powerlaw_fit(sfs::AbstractVector{<:Integer}; kmin::Integer)
    n = 0; s = 0.0
    for k in kmin:length(sfs) - 1
        c = sfs[k]
        c > 0 || continue
        n += c; s += c * log(k / (kmin - 0.5))
    end
    (n < 10 || s <= 0) && return (gamma = NaN, se = NaN, n = n)
    γ = 1 + n / s
    return (gamma = γ, se = (γ - 1) / sqrt(n), n = n)
end

"""
    kmin_ks(sfs; kmin_range, n_min = 5)

k_min minimising the Kolmogorov–Smirnov distance between the empirical tail and the fitted
discrete power law truncated at N − 1 (clonal mutations left out, as in `powerlaw_fit`;
Clauset et al. 2009 §3.3; port of `_estimate_kmin_ks` from BirthDeathMutation_AT.jl).
"""
function kmin_ks(sfs::AbstractVector{<:Integer}; kmin_range, n_min::Integer = 5)
    K = length(sfs) - 1
    suffix = zeros(Int, K + 1)
    for k in K:-1:1
        suffix[k] = suffix[k + 1] + sfs[k]
    end
    best, bestD = first(kmin_range), Inf
    for kmin in kmin_range
        kmin > K && break
        n = suffix[kmin]
        n < n_min && continue
        denom = sum(sfs[k] * log(k / (kmin - 0.5)) for k in kmin:K if sfs[k] > 0; init = 0.0)
        denom <= 0 && continue
        α = 1 + n / denom
        Z = sum(Float64(k)^(-α) for k in kmin:K)
        D = 0.0; S = 0.0; P = 0.0
        for k in kmin:K
            S += sfs[k] / n
            P += Float64(k)^(-α) / Z
            D = max(D, abs(S - P))
        end
        D < bestD && ((best, bestD) = (kmin, D))
    end
    return best
end

"""
    sfs_stats(sfs; kmin_frac = 0.01)

Number of mutations `S`; normalised Hill numbers `D1n, D2n, Dinfn` and Shannon entropy `H1n`;
cell-fraction (v1 eq. 10) entropy `H1cf` and `D2cf`; tail exponent with SE at
kmin = max(2, ⌈kmin_frac·N⌉).
"""
function sfs_stats(sfs::AbstractVector{<:Integer}; kmin_frac::Real = 0.01)
    S = sum(sfs)
    nan = (S = S, D1n = NaN, D2n = NaN, Dinfn = NaN, H1n = NaN, H1cf = NaN, D2cf = NaN,
           gamma = NaN, gamma_se = NaN, n_tail = 0)
    S == 0 && return nan
    f = powerlaw_fit(sfs; kmin = max(2, ceil(Int, kmin_frac * length(sfs))))
    D1 = sfs_hill_number(sfs, 1)
    return (S = S, D1n = D1, D2n = sfs_hill_number(sfs, 2), Dinfn = sfs_hill_number(sfs, Inf),
            H1n = log(D1), H1cf = sfs_renyi_entropy(sfs, 1; frequencies = :cell_fraction),
            D2cf = sfs_hill_number(sfs, 2; frequencies = :cell_fraction),
            gamma = f.gamma, gamma_se = f.se, n_tail = f.n)
end

"Mean and (sample) variance of a burden histogram `j => number of cells`."
function burden_stats(hist::AbstractDict)
    n = sum(values(hist); init = 0)
    n == 0 && return (mean = NaN, var = NaN)
    m = sum(j * c for (j, c) in hist) / n
    v = n > 1 ? sum(c * (j - m)^2 for (j, c) in hist) / (n - 1) : NaN
    return (mean = m, var = v)
end

"""
    logbin(y; per_decade = 8)

Geometric bins of y[k], k = 1…: fixed edges 10^(i/per_decade) (the same for every length, so
spectra of different populations can be averaged bin by bin); mean y per integer k in each bin.
"""
function logbin(y::AbstractVector{<:Real}; per_decade::Integer = 8)
    K = length(y)
    edges = unique(round.(Int, 10 .^ (0:1/per_decade:log10(K) + 1/per_decade)))
    ks = Float64[]; Ss = Float64[]
    for i in 1:length(edges) - 1
        lo, hi = edges[i], min(edges[i + 1] - 1, K)
        lo > hi && continue
        push!(ks, sqrt(lo * hi)); push!(Ss, sum(@view y[lo:hi]) / (hi - lo + 1))
    end
    return (k = ks, S = Ss)
end

"Median paired difference y − x, 95% bootstrap CI and Wilcoxon signed-rank p (missing pairs dropped)."
function paired_summary(x::AbstractVector, y::AbstractVector)
    keep = [!ismissing(a) && !ismissing(b) for (a, b) in zip(x, y)]
    d = Float64[y[i] - x[i] for i in eachindex(x) if keep[i]]
    isempty(d) && return (median_diff = NaN, lo = NaN, hi = NaN, p = NaN)
    rng = StableRNG(1)
    boots = [median(rand(rng, d, length(d))) for _ in 1:2000]
    p = all(iszero, d) ? 1.0 : pvalue(SignedRankTest(d))
    return (median_diff = median(d), lo = quantile(boots, 0.025), hi = quantile(boots, 0.975), p = p)
end

"Kaplan–Meier curve: event times `t` and survival `S` just after each."
function km(times::AbstractVector{<:Real}, event::AbstractVector{Bool})
    order = sortperm(times)
    ts = Float64[]; Ss = Float64[]
    S = 1.0; atrisk = length(times); i = 1
    while i <= length(order)
        t = times[order[i]]; d = 0; c = 0
        while i <= length(order) && times[order[i]] == t
            event[order[i]] ? (d += 1) : (c += 1); i += 1
        end
        if d > 0
            S *= 1 - d / atrisk
            push!(ts, t); push!(Ss, S)
        end
        atrisk -= d + c
    end
    return (t = ts, S = Ss)
end

"Kaplan–Meier median (first time S ≤ ½); missing if never reached."
function km_median(times, event)
    k = km(times, event)
    i = findfirst(<=(0.5), k.S)
    return i === nothing ? missing : k.t[i]
end

"Restricted mean survival time up to τ (area under the Kaplan–Meier curve)."
function rmst(times, event, τ::Real)
    k = km(times, event)
    area = 0.0; tprev = 0.0; S = 1.0
    for (t, s) in zip(k.t, k.S)
        t >= τ && break
        area += S * (t - tprev); tprev = t; S = s
    end
    return area + S * (τ - tprev)
end

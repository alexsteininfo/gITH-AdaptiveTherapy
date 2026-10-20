# Registry of all simulations: constants, treatment arms, cohorts, seeds and the Latin-hypercube
# design of the virtual trial. Single source of truth in code for code/01_datageneration/README.md.

# ── Constants ──────────────────────────────────────────────────────────────────
const B0 = 1.0                  # birth rate (time unit = one cell cycle)
const D0 = 0.5                  # intrinsic death rate
const DS = 4.0                  # drug kill at full dose
const MU = 1.0                  # neutral mutations per daughter (Poisson mean)
const NDET = 10^6               # detection size
const FK = 0.8                  # default detection size relative to carrying capacity
const TMAX = 250.0              # follow-up after detection
const BASE_SEED = 0x6A1B2C3D
const SAMPLE_N = [100, 1000, 10_000]
const FMIN = [0.0, 0.01, 0.05]
const HEAVY_BLOCKS = 12         # distributions at the end of the first 12 treatment blocks
const HEAVY_EVERY = 10.0        # … and every 10 time units
const COUNT_EVERY = 1.0         # population counts every time unit
const COUNT_EVERY_SCAN = 2.0    # … every 2 time units for scan arms
const TRAJ_DT = 0.25
const LEVEL_DEFAULT = Dict(:hol => 0.5, :mod => 0.8)
const LEVELS_HOL = [0.5, 0.6, 0.7, 0.8, 0.9]
const LEVELS_MOD = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
const N_EXAMPLES = 10           # tumours per cohort with full trees and clone histories

"Progression size: halfway between detection size and carrying capacity (no competition: never)."
nprog(ndet, fK) = fK == 0 ? typemax(Int) : round(Int, (ndet + ndet / fK) / 2)

# ── Arms ───────────────────────────────────────────────────────────────────────
"""
    Arm(id, kind, level, tmax, heavy)

One treatment arm. `kind ∈ (:single, :mtd, :hol, :mod, :met)`; `level` is the containment level as
a fraction of the detection size (holidays: lower threshold; modulation: stable size); `heavy`
selects when distributions are stored: `:full` (detection, first block ends, every
`HEAVY_EVERY`, end), `:early` (detection, end of the first full-dose block, end), `:end`.
"""
struct Arm
    id::String
    kind::Symbol
    level::Float64
    tmax::Float64
    heavy::Symbol
end

function arm(kind::Symbol; level = get(LEVEL_DEFAULT, kind, NaN), tmax = TMAX, heavy = :full)
    kind in (:single, :mtd, :hol, :mod, :met) || throw(ArgumentError("unknown arm kind $kind"))
    heavy in (:full, :early, :end) || throw(ArgumentError("unknown heavy policy $heavy"))
    default = get(LEVEL_DEFAULT, kind, NaN)
    id = (isnan(level) || level == default) ? string(kind) : "$(kind)_$(level)"
    return Arm(id, kind, Float64(level), Float64(tmax), heavy)
end

scan_arms(; tmax = TMAX, heavy = :end) = [arm(:mtd; tmax, heavy);
                                          [arm(:hol; level = L, tmax, heavy) for L in LEVELS_HOL];
                                          [arm(:mod; level = L, tmax, heavy) for L in LEVELS_MOD]]

# ── Cohorts ────────────────────────────────────────────────────────────────────
"""
    Cohort

A set of tumours with common parameters. `nu`, `alpha`, `fK` are numbers, a `Distribution`
(`alpha` in block 06), or `:lhs` (per-patient values of the virtual trial). Cohorts with the same
`key` share tumour seeds, so their detection populations are identical.
"""
struct Cohort
    name::String
    key::Int
    nu::Any
    alpha::Any
    fK::Any
    cost::Float64
    d0::Float64
    n::Int
    arms::Vector{Arm}
    ndet::Int
    growth_heavy::Vector{Float64}
end

Cohort(name, key; nu = 1e-4, alpha = 0.0, fK = FK, cost = 0.0, d0 = D0, n, arms, ndet = NDET,
       growth_heavy = Float64[]) = Cohort(name, key, nu, alpha, fK, cost, d0, n, arms, ndet, growth_heavy)

part1_arms(tmax) = [arm(:single; tmax), arm(:mtd; tmax), arm(:hol; tmax), arm(:mod; tmax)]

const BLOCKS = Dict{String,Vector{Cohort}}(
    "01_dynamics" => [
        Cohort("neutral", 1; nu = 0.0, alpha = 1.0, n = 500, growth_heavy = [0.1],
               arms = [arm(:single; tmax = 40.0), arm(:hol; tmax = 40.0), arm(:mod; tmax = 40.0)]),
        Cohort("full", 2; alpha = 0.0, n = 500, growth_heavy = [0.1], arms = part1_arms(TMAX)),
        Cohort("partial", 3; alpha = 0.006, n = 500, growth_heavy = [0.1], arms = part1_arms(TMAX)),
        Cohort("weak", 4; alpha = 0.8, n = 500, growth_heavy = [0.1], arms = part1_arms(60.0)),
    ],
    "02_tradeoff" => [
        Cohort("full", 2; alpha = 0.0, n = 300, arms = scan_arms()),
        Cohort("partial", 3; alpha = 0.006, n = 300, arms = scan_arms()),
    ],
    "03_sensitivity" => [
        Cohort("fK0.5", 5; fK = 0.5, n = 200, arms = scan_arms()),
        Cohort("fK0.7", 6; fK = 0.7, n = 200, arms = scan_arms()),
        Cohort("fK0.9", 7; fK = 0.9, n = 200, arms = scan_arms()),
        Cohort("nu1e-5", 8; nu = 1e-5, n = 200, arms = scan_arms()),
        Cohort("nu1e-3", 9; nu = 1e-3, n = 200, arms = scan_arms()),
    ],
    "04_trial" => [
        Cohort("trial", 10; nu = :lhs, alpha = :lhs, fK = :lhs, n = 500,
               arms = [arm(:mtd; heavy = :early), arm(:hol; heavy = :early), arm(:mod; heavy = :early)]),
    ],
    "05_metronomic" => [
        Cohort("neutral", 1; nu = 0.0, alpha = 1.0, n = 500, arms = [arm(:met; tmax = 40.0, heavy = :end)]),
        Cohort("full", 2; alpha = 0.0, n = 500, arms = [arm(:met; heavy = :end)]),
    ],
    "06_cost" => [
        Cohort(name, key; alpha = Beta(1, 4), cost = c, n = 300,
               arms = [arm(:single; heavy = :end), arm(:mtd; heavy = :end), arm(:hol; heavy = :end),
                       arm(:mod; heavy = :end), arm(:mod; level = 1.0, heavy = :end)])
        for (name, key, c) in (("cost0", 11, 0.0), ("cost0.02", 12, 0.02))
    ],
    "07_validation" => [
        Cohort("sfs_pb", 13; nu = 0.0, alpha = 1.0, fK = 0.0, d0 = 0.0, n = 200, ndet = 10^5, arms = Arm[]),
        Cohort("sfs_bd", 14; nu = 0.0, alpha = 1.0, fK = 0.0, n = 200, ndet = 10^5, arms = Arm[]),
        Cohort("ld", 15; nu = 1e-3, alpha = 0.0, fK = 0.0, n = 200, ndet = 10^5, arms = Arm[]),
    ],
)

"Smoke-test version of a cohort: 3 tumours of 10⁴ cells, follow-up divided by 5."
test_cohort(c::Cohort) = Cohort(c.name, c.key, c.nu, c.alpha, c.fK, c.cost, c.d0, 3,
                                [Arm(a.id, a.kind, a.level, a.tmax / 5, a.heavy) for a in c.arms],
                                min(c.ndet, 10^4), c.growth_heavy)

# ── Seeds and parameters ──────────────────────────────────────────────────────
tumour_seed(c::Cohort, t::Integer)::UInt64 = derive_seed(BASE_SEED, c.key, t)

"""
    lhs_design(n, seed)

Latin hypercube of `n` virtual patients over log10 ν ∈ [−5, −3], α ∈ [0, 0.006] and
fK ∈ [0.5, 0.9]: one stratum per patient in every dimension, uniform jitter within strata.
"""
function lhs_design(n::Integer, seed::Integer)
    rng = StableRNG(seed)
    u() = (randperm(rng, n) .- 1 .+ rand(rng, n)) ./ n
    a, b, c = u(), u(), u()
    return [(nu = 10.0^(-5 + 2a[i]), alpha = 0.006b[i], fK = 0.5 + 0.4c[i]) for i in 1:n]
end

"Parameters of tumour `t` of cohort `c` (resolves the trial's Latin hypercube)."
function params(c::Cohort, t::Integer)
    if c.nu === :lhs
        d = lhs_design(c.n, derive_seed(BASE_SEED, c.key))[t]
        return (nu = d.nu, alpha = d.alpha, fK = d.fK, cost = c.cost, d0 = c.d0)
    end
    return (nu = c.nu, alpha = c.alpha, fK = c.fK, cost = c.cost, d0 = c.d0)
end

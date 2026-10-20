# Blocks and treatment schedules built from a parameter set `p = (nu, alpha, fK, cost, d0)`.

"Birth–death block with competition, linear drug response and resistance mutations (no stop rule)."
function base_block(p, ndet::Integer)
    density = p.fK == 0 ? Returns(0.0) : (N -> (B0 - p.d0) * p.fK / ndet * N)
    Δs = p.cost == 0 ? 0.0 : (α -> -p.cost * (1 - α))
    mutations = p.nu == 0 ? () : (Resistance(ν = p.nu, α = p.alpha, Δs = Δs),)
    return BirthDeathBlock(birth = s -> B0 * (1 + s), death = s -> p.d0, death_density = density,
                           response = LinearResponse(DS), μ = MU, mutations = mutations)
end

"Untreated growth from one cell to the detection size (restart on extinction)."
growth_block(p, ndet::Integer) =
    remake(base_block(p, ndet); stop_above = ndet, restart_on_extinction = true, label = "growth")

"The schedule of arm `a`, starting at detection."
function schedule(a::Arm, p, ndet::Integer)
    base = base_block(p, ndet)
    N_prog = nprog(ndet, p.fK)
    a.kind == :mtd && return MTD(base; N_prog, t_max = a.tmax)
    a.kind == :hol && return AdaptiveOnOff(base; N_low = round(Int, a.level * ndet), N_up = ndet,
                                           N_prog, t_max = a.tmax)
    a.kind == :mod && return Stabilising(base; N_stable = round(Int, a.level * ndet), N_prog, t_max = a.tmax)
    a.kind == :met && return Metronomic(base; intervals = [(3.0k, 3.0k + 1) for k in 0:ceil(Int, a.tmax / 3)],
                                        N_prog, t_max = a.tmax)
    a.kind == :single && return FunctionSchedule() do pop, i, t0
        i == 0 && return remake(base; dose = FixedDose(1.0), stop_below = round(Int, 0.1ndet),
                                tmax = t0 + a.tmax, label = "single on")
        i == 1 && pop.log[end].end_reason == :below &&
            return remake(base; dose = FixedDose(0.0), stop_above = ndet, tmax = t0 + a.tmax,
                          label = "single off")
        return nothing
    end
    throw(ArgumentError("unknown arm kind $(a.kind)"))
end

# Make all figure panels from the full data (data/processed/<block>/) into figures/<figure>/.
# Each script runs in its own namespace; a failing script is reported and the others still run.
# Run: julia --project=. code/03_figures/run_figures.jl [regex]     e.g. fig_2 or "fig_[45]"
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FILTER = isempty(ARGS) ? "" : ARGS[1]
get(ENV, "FIG_DATA", "") == "test" && error("FIG_DATA=test is set; use smoke.jl for the test data")

names = sort([d for d in readdir(@__DIR__) if occursin(r"^((fig|supp)_\d+|stabilising)$", d) && occursin(Regex(FILTER), d)])
failed = String[]
for name in names
    script = joinpath(@__DIR__, name, "$name.jl")
    t0 = time()
    try
        m = Module(Symbol(name))                       # each script in its own namespace
        Core.eval(m, :(include(p) = Base.include($m, p)))
        Base.include(m, script)
        dir = joinpath(ROOT, "figures", name)
        new = count(f -> endswith(f, ".png") && mtime(joinpath(dir, f)) >= t0, readdir(dir))
        println("$name: $new panels in $(round(time() - t0; digits = 1)) s")
    catch e
        push!(failed, name)
        println("$name: FAILED after $(round(time() - t0; digits = 1)) s")
        showerror(stdout, e, catch_backtrace()); println()
    end
end
println(isempty(failed) ? "all $(length(names)) figures done" : "failed: $(join(failed, ", "))")
isempty(failed) || exit(1)

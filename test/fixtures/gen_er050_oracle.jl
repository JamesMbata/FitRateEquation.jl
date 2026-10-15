#=
Generate test/fixtures/er050_oracle.toml: a record of what EnzymeRates 0.5.0 (c9ae54a) produced
for every registered FitRateEquation mechanism -- parameter names (fitted, full, dependent),
metabolites, substrates, the elementary-step view, and rate_equation values at fixed points.
test_enzymerates_names.jl replays it against the installed EnzymeRates through FRE's name
layer, so a version bump that renames, reorders or re-derives anything fails loudly.

Must run under an env that resolves EnzymeRates 0.5.0 and FRE commit 1bda0f2 (the last
0.5.0-pinned main), e.g. a checkout of 1bda0f2:

    julia --project=<FRE@1bda0f2> test/fixtures/gen_er050_oracle.jl <out.toml>
=#
using FitRateEquation, EnzymeRates, Random, TOML

pkgversion(EnzymeRates) == v"0.5.0" || error("oracle must be generated on EnzymeRates 0.5.0")
const KEQ = Dict(:G6PD => 13.7, :PGD => 0.079)
out = Dict{String,Any}("enzymerates_version" => "0.5.0", "mechanisms" => Any[])

for enz in (:G6PD, :PGD), v in FitRateEquation.consensus_variants(enz)
    m = v.mech
    fp = collect(EnzymeRates.fitted_params(m))
    dep, _ = EnzymeRates._dependent_param_exprs(typeof(m))
    mets = collect(EnzymeRates.metabolites(m))
    steps = [Dict("lhs" => string.(st.lhs), "rhs" => string.(st.rhs), "is_eq" => st.is_eq,
                  "group" => st.group, "param_index" => st.param_index)
             for st in FitRateEquation._mechanism_steps(m)]
    rng = Xoshiro(hash((enz, v.name)) % UInt32)
    points = Any[]
    for _ in 1:3
        θ = Dict(string(s) => 10.0^(-6 + 3rand(rng)) for s in fp)
        p = NamedTuple{(fp..., :Keq, :E_total)}((map(s -> θ[string(s)], fp)..., KEQ[enz], 1.0))
        for _ in 1:4
            c = Dict(string(s) => 10.0^(-7 + 4rand(rng)) for s in mets)
            concs = NamedTuple{Tuple(mets)}(map(s -> c[string(s)], mets))
            push!(points, Dict("params" => θ, "keq" => KEQ[enz], "concs" => c,
                               "rate" => EnzymeRates.rate_equation(m, concs, p)))
        end
    end
    push!(out["mechanisms"], Dict(
        "enzyme" => string(enz), "variant" => string(v.name),
        "fitted_params" => string.(fp),
        "parameters" => string.(collect(EnzymeRates.parameters(m))),
        "dependent_params" => sort(string.(collect(keys(dep)))),
        "metabolites" => string.(mets),
        "substrates" => string.(collect(EnzymeRates.substrates(m))),
        "steps" => steps, "points" => points))
end

open(ARGS[1], "w") do io; TOML.print(io, out; sorted = true); end
println("wrote ", ARGS[1], " (", length(out["mechanisms"]), " mechanisms)")

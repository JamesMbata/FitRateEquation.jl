using FitRateEquation
using EnzymeRates
using Test

@testset "PGD mechanisms" begin
    vs = FitRateEquation.consensus_variants(:PGD)
    @test [v.name for v in vs] == [:RE_rate_eq, :SS_NADPH_release_rate_eq, Symbol("+NADPH_deadend_rate_eq"), :cha_base]
    v1 = vs[1].mech; v2 = vs[2].mech
    mets1 = Set(EnzymeRates.metabolites(v1))
    for s in (:NADP, :PGA, :Ru5P, :CO2, :NADPH, :ATP)
        @test s in mets1
    end
    # RE: one SS step (catalysis). V2 has THREE SS steps: catalysis gauge anchor +
    # Ru5P-release (rate-limiting, Topham) + NADPH-release (de-conflation). The Ru5P-SS
    # step is required for a working apparent-constant gauge on the 3-product mechanism
    # (see the V2 mechanism comment / spec).
    @test count(st -> st[3] == false, EnzymeRates.reactions(v1)) == 1
    @test count(st -> st[3] == false, EnzymeRates.reactions(v2)) == 3
    # V2 has strictly more free coordinates than V1 (the SS-release reverse DOFs) —
    # the de-conflation is structural, not pin-induced.
    @test length(free_params(v2)) > length(free_params(v1))
    # Both compile a finite forward rate at saturating substrates, zero products.
    p1 = build_params(v1, fill(-3.0, length(free_params(v1))); keq=0.079)
    concs = NamedTuple{Tuple(EnzymeRates.metabolites(v1))}(
        Tuple(s in (:NADP,:PGA) ? 1e-3 : 0.0 for s in EnzymeRates.metabolites(v1)))
    @test isfinite(EnzymeRates.rate_equation(v1, concs, p1))
end

@testset "PGD cha_base: gauge term + SS Ru5P-release + E·PGA NADPH dead-end" begin
    vs = FitRateEquation.consensus_variants(:PGD)
    i = findfirst(v -> Symbol(v.name) === :cha_base, vs)
    @test i !== nothing
    m = vs[i].mech
    # Free-enzyme (constant) denominator term must exist, else no gauge.
    _, den = EnzymeRates._raw_symbolic_rate_polys(typeof(m))
    metset = Set(EnzymeRates.metabolites(m))
    has_const = any(isempty(first(FitRateEquation._split_mono(mono, metset))) for (mono, _) in den)
    @test has_const
    # Exactly one SS release step (Ru5P): is_eq=false, touches :Ru5P.
    steps = FitRateEquation._mechanism_steps(m)
    ss_releases = [st for st in steps if !st.is_eq && (:Ru5P in st.rhs_mets || :Ru5P in st.lhs_mets)]
    @test length(ss_releases) == 1
    # The forward-Ki cross channel [PGA·NADPH] is present (E·PGA dead-end).
    @test any(occursin("PGA", string(mono)) && occursin("NADPH", string(mono)) for (mono, _) in den)
end

using FitRateEquation
using FitRateEquation: v1_mechanism, v2_mechanism
using EnzymeRates
using Test

@testset "mechanisms" begin
    v1 = v1_mechanism(); v2 = v2_mechanism()
    mets1 = Set(EnzymeRates.metabolites(v1))
    # Consensus carries both substrates, both products, and ATP (dead-end regulator)
    for s in (:NADP, :G6P, :NADPH, :PGLn, :ATP)
        @test s in mets1
    end
    # RE_rate_eq = pure RE: every binding+release step is rapid-equilibrium, catalysis is SS.
    res1 = EnzymeRates.reactions(v1)
    @test count(st -> st[3] == false, res1) == 1            # exactly one SS step (catalysis)
    # SS_NADPH_release_rate_eq adds a second SS step (the NADPH release), so SS count is 2.
    res2 = EnzymeRates.reactions(v2_mechanism())
    @test count(st -> st[3] == false, res2) == 2
    # Both compile a finite forward rate.
    p1 = build_params(v1, fill(-3.0, length(free_params(v1))); keq=10.0)
    concs = NamedTuple{Tuple(EnzymeRates.metabolites(v1))}(
        Tuple(s in (:NADP,:G6P) ? 1e-3 : 0.0 for s in EnzymeRates.metabolites(v1)))
    @test isfinite(EnzymeRates.rate_equation(v1, concs, p1))
end

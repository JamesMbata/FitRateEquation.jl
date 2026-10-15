using Test, FitRateEquation, EnzymeRates

@testset "PGD macro_collect reconstruction" begin
    for m in (mech.mech for mech in FitRateEquation.consensus_variants(:PGD))
        free = free_params(m); logθ = -3 .+ 2 .* rand(length(free)); keq = 0.079
        concs = NamedTuple{Tuple(EnzymeRates.metabolites(m))}(
            Tuple(rand()*1e-3 for _ in EnzymeRates.metabolites(m)))
        num, den = EnzymeRates._raw_symbolic_rate_polys(typeof(m))
        vals = FitRateEquation._micro_values(m, logθ; keq=keq)
        vnum = FitRateEquation._eval_poly(num, concs, vals)
        vden = FitRateEquation._eval_poly(den, concs, vals)
        vref = EnzymeRates.rate_equation(m, concs, build_params(m, logθ; keq=keq))
        @test isapprox(vnum / vden, vref; rtol=1e-6)
    end
end

@testset "PGD macro_constants named (forward shape set)" begin
    v1 = FitRateEquation.consensus_variants(:PGD)[1].mech
    free = free_params(v1); logθ = -3 .+ 2 .* rand(length(free))
    mc = macro_constants(v1, logθ; keq=0.079, enzyme=:PGD)
    names = Set(getfield.(mc, :name))
    for n in (:Km_NADP, :Km_PGA)   # forward shape constants (deliverable)
        @test n in names
    end
    # V1 has no [PGA·NADPH] cross term, so it carries NO forward Ki_NADPH coord — the bare
    # [NADPH] term is now a clean reverse-Km lump (lump_NADPH^1), not Ki_NADPH. Forward
    # Ki_NADPH is a V3-only de-conflated ki-ratio coord (2026-06-10).
    @test :Ki_NADPH ∉ names
    @test Symbol("lump_NADPH^1") in names
    km = mc[findfirst(x -> x.name == :Km_PGA, mc)]
    @test km.value > 0 && km.role == :named
end

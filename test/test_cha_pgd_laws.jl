using FitRateEquation
using EnzymeRates
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert

# The PGD :cha_base mechanism (Topham Bi-Ter, base-only: catalysis-SS + Ru5P-release-SS,
# CO2/NADPH release RE). Bound by name, not index.
function _pgd_cha_base_mech()
    vars = FitRateEquation.consensus_variants(:PGD)
    vars[findfirst(v -> Symbol(v.name) === :cha_base, vars)].mech
end

@testset "exactness anchor: cha_rate_PGD == rate_equation (cha_base, rtol 1e-10)" begin
    m = _pgd_cha_base_mech()
    mets = EnzymeRates.metabolites(m)               # (:NADP,:PGA,:CO2,:NADPH,:Ru5P,:ATP)
    # Product grid: forward, +CO2, +Ru5P, +NADPH, all-products, all+ATP.
    grid = [
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=0.0,  NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=5e-5, NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=5e-6, ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=5e-5, NADPH=5e-6, ATP=0.0),
        (; NADP=2e-5, PGA=80e-6, CO2=2e-4, Ru5P=1e-4, NADPH=8e-6, ATP=5e-4),
    ]
    for _ in 1:20
        free = free_params(m)
        logθ = -3 .+ 2 .* rand(length(free))
        keq  = 0.079
        mac  = cha_macro_readoffs_PGD(m, logθ; keq=keq)
        # Characteristic forward rate at these params -> floor for the rel-error denominator
        # (near a Haldane equilibrium null vref can pass exactly through 0).
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
            vcha = cha_rate_PGD(mac; conc...)
            @test isapprox(vcha, vref; rtol=1e-10, atol=1e-10 * vsat)
        end
    end
end

@testset "PGD macro->micro->law identity (no promoted fiber)" begin
    vs = FitRateEquation.consensus_variants(:PGD)
    m  = vs[findfirst(v -> Symbol(v.name) === :cha_base, vs)].mech
    for _ in 1:10
        free = free_params(m); logθ = -3 .+ 2 .* rand(length(free))
        mac  = cha_macro_readoffs_PGD(m, logθ; keq=10.0)
        mic  = cha_micro_from_macro_PGD(mac)
        pt = (; NADP=5e-6, PGA=40e-6, Ru5P=1e-4, CO2=1e-4, NADPH=5e-6, ATP=1e-3)
        @test isapprox(cha_rate_PGD(mic; pt...), cha_rate_PGD(mac; pt...); rtol=1e-10)
    end
end

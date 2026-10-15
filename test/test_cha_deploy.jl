# Closed-form deploy inverse: cha_deploy_micro maps a Cha macro tuple back to the
# mechanism's INDEPENDENT free-params (logθ aligned to free_params(mech)) such that
# EnzymeRates.rate_equation(mech, .) reproduces ChaLaws.cha_rate_* to rtol 1e-9. This is the
# inverse of cha_macro_readoffs_* (cha_invert.jl): we round-trip a real micro point THROUGH
# the macro tuple and back, and check the deployed rate matches the macro law absolutely.
#
# Gauge: the readoff always has kf = k5f = 1 and Et = E_total = 1 (the unit-enzyme gauge),
# and build_params likewise fixes k5f = 1 / E_total = 1. The deploy micro therefore shares
# the macro tuple's gauge automatically, so the 1e-9 ABSOLUTE rate test holds (no shape-only
# fallback needed).
using FitRateEquation
using EnzymeRates
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy

@testset "cha_deploy_micro reproduces the fitted Cha law (G6PD, rtol 1e-9)" begin
    m = FitRateEquation.v2_mechanism(); mets = EnzymeRates.metabolites(m)
    for _ in 1:8
        logθ = -3 .+ 2 .* rand(length(free_params(m))); keq = 10.0
        mac = cha_macro_readoffs_G6PD(m, logθ; keq=keq)
        coords = Dict(s => getfield(mac, s) for s in cha_coords(:G6PD))
        logθd = cha_deploy_micro(:G6PD, m, coords; keq=keq, koffQ=mac.koffQ,
                    release_rate=mac.koffQ, release_eq=mac.Km_NADPH_rev, kr=mac.kr)
        bp = build_params(m, logθd; keq=keq)
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:G6P) ? 1e-2 : 0.0 for s in mets)), bp))
        for conc in [(;NADP=5e-6,G6P=40e-6,NADPH=5e-6,PGLn=1e-4,ATP=1e-3),
                     (;NADP=2e-5,G6P=80e-6,NADPH=8e-6,PGLn=2e-4,ATP=5e-4)]
            cc = NamedTuple{Tuple(mets)}(Tuple(get(conc, s, 0.0) for s in mets))
            vref = ChaLaws.cha_rate_G6PD(mac; conc...)
            vdep = EnzymeRates.rate_equation(m, cc, bp)
            @test isapprox(vdep, vref; rtol=1e-9, atol=1e-9*vsat)
        end
    end
end

@testset "cha_deploy_micro reproduces the fitted Cha law (PGD, rtol 1e-9)" begin
    vs = FitRateEquation.consensus_variants(:PGD)
    m  = vs[findfirst(v -> Symbol(v.name) === :cha_base, vs)].mech
    mets = EnzymeRates.metabolites(m)
    for _ in 1:8
        logθ = -3 .+ 2 .* rand(length(free_params(m))); keq = 10.0
        mac = cha_macro_readoffs_PGD(m, logθ; keq=keq)
        coords = Dict(s => getfield(mac, s) for s in cha_coords(:PGD))
        # PGD promoted release is Ru5P (koff/kon); its release equilibrium KdRu = koff/kon is
        # distinct from Km_NADPH_rev. Pass the actual KdRu so the inverse round-trips exactly.
        logθd = cha_deploy_micro(:PGD, m, coords; keq=keq, koffQ=mac.koff,
                    release_rate=mac.koff, release_eq=mac.koff/mac.kon, kr=mac.kr)
        bp = build_params(m, logθd; keq=keq)
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)), bp))
        for conc in [(;NADP=5e-6,PGA=40e-6,Ru5P=1e-4,CO2=1e-4,NADPH=5e-6,ATP=1e-3),
                     (;NADP=2e-5,PGA=80e-6,Ru5P=2e-4,CO2=2e-4,NADPH=8e-6,ATP=5e-4)]
            cc = NamedTuple{Tuple(mets)}(Tuple(get(conc, s, 0.0) for s in mets))
            vref = ChaLaws.cha_rate_PGD(mac; conc...)
            vdep = EnzymeRates.rate_equation(m, cc, bp)
            @test isapprox(vdep, vref; rtol=1e-9, atol=1e-9*vsat)
        end
    end
end

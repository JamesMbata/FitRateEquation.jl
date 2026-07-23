using FitRateEquation
using EnzymeRates
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert

_p4de_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                vs[findfirst(v -> Symbol(v.name) === :full_re_deadends, vs)].mech)
_p4fr_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                vs[findfirst(v -> Symbol(v.name) === :full_re, vs)].mech)

@testset "exactness gate: cha_rate_PGD_fullRE == rate_equation (:full_re_deadends, rtol 1e-10)" begin
    m    = _p4de_mech()
    mets = EnzymeRates.metabolites(m)               # (:NADP,:PGA,:CO2,:NADPH,:Ru5P) — no ATP
    # Grid that EXERCISES the dead-ends at NADPH=0 (the Weisz 6A–7B regime), plus both-product rows.
    grid = [
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,    Ru5P=0.0,  NADPH=0.0),   # baseline (no products)
        (; NADP=5e-6, PGA=40e-6, CO2=2.2e-3, Ru5P=0.0,  NADPH=0.0),   # CO2 dead-end active, NADPH=0 (6A)
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,    Ru5P=4e-4, NADPH=0.0),   # Ru5P dead-end active, NADPH=0 (7A)
        (; NADP=5e-6, PGA=40e-6, CO2=2.2e-3, Ru5P=4e-4, NADPH=0.0),   # both dead-ends, NADPH=0
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4,   Ru5P=5e-5, NADPH=5e-6),  # all products present
        (; NADP=2e-5, PGA=80e-6, CO2=2e-4,   Ru5P=1e-4, NADPH=8e-6),
    ]
    Random.seed!(41)
    for _ in 1:20
        logθ = -3 .+ 2 .* rand(length(free_params(m)))
        keq  = 0.079
        mac  = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
        @test isfinite(mac.Ki_CO2) && isfinite(mac.Ki_Ru5P)         # dead-ends carried
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc   = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
            vcha = cha_rate_PGD_fullRE(mac; conc...)
            @test isapprox(vcha, vref; rtol=1e-10, atol=1e-10 * vsat)
        end
    end
end

@testset "dead-ends actually inhibit at NADPH=0 (structural: :full_re would be flat)" begin
    m   = _p4de_mech()
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    mac = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=0.079)
    # Tighten the two dead-end constants so the effect is unambiguous.
    mac = merge(mac, (; Ki_CO2=1e-4, Ki_Ru5P=5e-5))
    base = cha_rate_PGD_fullRE(mac; NADP=5e-6, PGA=40e-6)                       # C=R=Q=0
    co2  = cha_rate_PGD_fullRE(mac; NADP=5e-6, PGA=40e-6, CO2=2.2e-3)           # CO2 present, NADPH=0
    ru5p = cha_rate_PGD_fullRE(mac; NADP=5e-6, PGA=40e-6, Ru5P=4e-4)            # Ru5P present, NADPH=0
    @test co2  < base                                # CO2 inhibits at NADPH=0 (the Phase-4 point)
    @test ru5p < base                                # Ru5P inhibits at NADPH=0
end

@testset ":full_re bit-identity: the new terms vanish (Inf) on a :full_re tuple" begin
    # A :full_re readoff has no Ki_CO2/Ki_Ru5P fields -> hasproperty false -> Inf -> terms == 0.
    m   = _p4fr_mech()
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    mac = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=0.079)     # carries Ki_CO2=Ki_Ru5P=Inf
    pt  = (; NADP=5e-6, PGA=8e-5, CO2=1e-3, Ru5P=1e-3, NADPH=2e-5)
    with_inf = cha_rate_PGD_fullRE(mac; pt...)
    # A bare core tuple with NO Ki_* fields at all must give the identical rate.
    bare = (; mac.Kd_NADP, mac.Kd_PGA, mac.alpha, mac.Kd_NADPH, mac.Kd_Ru5P, mac.Kd_CO2,
              mac.kf, mac.kr, mac.Et)
    @test cha_rate_PGD_fullRE(bare; pt...) == with_inf          # bit-identical, not just ≈
end

using FitRateEquation
using FitRateEquation: v2_mechanism
using EnzymeRates
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert

# A random, physically-valid macro tuple for the forward/base-limit shape tests.
# C = 1 + kf/koffQ (corrected two-SS-segment factor; koffQ alone).
function _rand_macro(; koffQ=1.0)
    kf = 100.0 + 200.0 * rand()
    kr = 10.0 + 100.0 * rand()
    konQ = (0.5 + rand()) * koffQ
    (; Kd_NADP=5e-6 + 5e-6*rand(), Kd_G6P=30e-6 + 30e-6*rand(),
       Kd_6PGLn=1e-4 + 2e-4*rand(), alpha=0.5 + rand(),
       Km_NADPH_rev=koffQ/konQ, Ki_NADPH=10e-6 + 10e-6*rand(),
       Ki_ATP=1e-3 + 1e-3*rand(), Ki_ATP_EG=1e-2 + 1e-2*rand(),
       koffQ=koffQ, konQ=konQ, kf=kf, kr=kr, Et=1.0)
end

@testset "forward half-max NADP recovers alpha*Kd_NADP/C (saturating G6P)" begin
    # At P=Q=ATP=0 and saturating G6P, the NADP concentration giving half v_max is the
    # apparent Km_NADP = alpha*Kd_NADP/C, with C = 1 + kf/koffQ.
    for _ in 1:5
        m = _rand_macro()
        C = 1 + m.kf / m.koffQ
        KmA = m.alpha * m.Kd_NADP / C
        Bsat = m.Kd_G6P * 1e6
        vmax  = cha_rate_G6PD(m; NADP=1e3*KmA, G6P=Bsat)
        vhalf = cha_rate_G6PD(m; NADP=KmA, G6P=Bsat)
        @test isapprox(vhalf, vmax/2; rtol=1e-3)
    end
end

@testset "base limit koffQ->oo gives C->1 (Km_NADP -> alpha*Kd_NADP)" begin
    m0 = _rand_macro(; koffQ=1e12)
    C = 1 + m0.kf / m0.koffQ
    @test isapprox(C, 1.0; rtol=1e-6)
    KmA_C1 = m0.alpha * m0.Kd_NADP
    Bsat = m0.Kd_G6P * 1e6
    vmax  = cha_rate_G6PD(m0; NADP=1e3*KmA_C1, G6P=Bsat)
    vhalf = cha_rate_G6PD(m0; NADP=KmA_C1, G6P=Bsat)
    @test isapprox(vhalf, vmax/2; rtol=1e-3)
end

@testset "exactness anchor: cha_rate_G6PD == rate_equation (v2, rtol 1e-10)" begin
    m = v2_mechanism()
    mets = EnzymeRates.metabolites(m)               # (:G6P,:NADP,:NADPH,:PGLn,:ATP)
    # Both-product grid: forward, +product P, +product Q, both, both+ATP (x2).
    grid = [
        (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=0.0,  ATP=0.0),
        (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=1e-4, ATP=0.0),
        (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=0.0,  ATP=0.0),
        (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0),
        (; NADP=1e-6, G6P=10e-6, NADPH=2e-6, PGLn=5e-5, ATP=1e-3),
        (; NADP=2e-5, G6P=80e-6, NADPH=8e-6, PGLn=2e-4, ATP=5e-4),
    ]
    for _ in 1:20
        free = free_params(m)
        logθ = -3 .+ 2 .* rand(length(free))
        keq  = 10.0
        mac  = cha_macro_readoffs_G6PD(m, logθ; keq=keq)
        # Characteristic forward rate at these params -> floor for the rel-error denominator
        # (near a Haldane equilibrium null vref can pass exactly through 0).
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:G6P) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
            vcha = cha_rate_G6PD(mac; conc...)
            @test isapprox(vcha, vref; rtol=1e-10, atol=1e-10 * vsat)
        end
    end
end

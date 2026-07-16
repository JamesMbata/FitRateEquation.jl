using FitRateEquation
using EnzymeRates
using Test
using FitRateEquation.ChaLawsHK1

if FitRateEquation.HK1_AVAILABLE

# Macro tuple read off the literature-pin vectors.
_hk1_macro(alpha) = (; Kd_Glc=50e-6, Kd_ATP=1e-3, Ki_G6P_C=15e-6, Ki_ADP=1.5e-3,
                       Ki_G6P_N=6.9e-6, K_Pi_N=750e-6, alpha=alpha, Keq=2700.0,
                       kf=1.0, k2f=1.0, Et=1.0)

# Map a metabolite NamedTuple to the EnzymeRates param vector for the mechanism, using the
# SAME literal K-names the derivation verified (H1: K15/K21/K30/K39/K48; H3: K14/K20/K26/K34/K43).
function _ka_params(alpha)
    if alpha === :one
        (k1f=1.0, k2f=1.0, K3=50e-6, K15=1e-3, K21=15e-6, K30=1.5e-3, K39=750e-6, K48=6.9e-6,
         Keq=2700.0, E_total=1.0)
    else
        (k1f=1.0, k2f=1.0, K3=50e-6, K14=1e-3, K20=15e-6, K26=1.5e-3, K34=750e-6, K43=6.9e-6,
         Keq=2700.0, E_total=1.0)
    end
end

@testset "cha_rate_HK1 == EnzymeRates.rate_equation (bitwise, both candidates)" begin
    for (alpha_sym, alpha_val) in ((:one, 1.0), (:infinity, Inf))
        m   = build_hk1_mechanism(alpha=alpha_sym, glc_g6p_dead_end=true, glc_adp_dead_end=true)
        mac = _hk1_macro(alpha_val)
        kp  = _ka_params(alpha_sym)
        grid = [
            (Glucose=1e-4, ATP=2e-3, G6P=0.0,  ADP=0.0,  Pi=0.0),
            (Glucose=1e-4, ATP=2e-3, G6P=1e-4, ADP=0.0,  Pi=0.0),
            (Glucose=1e-4, ATP=2e-3, G6P=1e-4, ADP=1e-3, Pi=5e-3),
            (Glucose=5e-5, ATP=1e-3, G6P=3e-5, ADP=2e-3, Pi=9.2e-3),
        ]
        for c in grid
            concs = (Glucose=c.Glucose, ATP=c.ATP, G6P=c.G6P, ADP=c.ADP, Pi=c.Pi)
            v_ka  = EnzymeRates.rate_equation(m, concs, merge(kp, (Keq=2700.0, E_total=1.0)))
            v_cha = cha_rate_HK1(mac; Glucose=c.Glucose, ATP=c.ATP, G6P=c.G6P, ADP=c.ADP, Pi=c.Pi)
            @test isapprox(v_cha, v_ka; rtol=1e-12, atol=1e-300)
        end
    end
end

@testset "Pi-invariance at G6P=0 with k2f=1" begin
    mac = _hk1_macro(1.0)
    v0 = cha_rate_HK1(mac; Glucose=5e-5, ATP=1e-3, G6P=0.0, ADP=0.0, Pi=0.0)
    for pi in (0.0, 1e-3, 9.2e-3)
        v = cha_rate_HK1(mac; Glucose=5e-5, ATP=1e-3, G6P=0.0, ADP=0.0, Pi=pi)
        @test isapprox(v, v0; rtol=1e-12)
    end
end

else
    @testset "HK1 laws (skipped: HK1_AVAILABLE=false)" begin
        @test_skip "build_hk1_mechanism unavailable pending HK1 wiring port"
    end
end

using FitRateEquation
using EnzymeRates
using Test
using FitRateEquation.ChaFit, FitRateEquation.ChaDeploy, FitRateEquation.ChaLawsHK1

if FitRateEquation.HK1_AVAILABLE

const MID = FitRateEquation

# Deploy round-trip: macro coords -> deploy logθ -> build_params -> EnzymeRates.rate_equation must
# reproduce cha_rate_HK1 at the same coords. H1 uses the raw {Ki_G6P_C, Ki_G6P_N}; H4 uses the
# reparameterized {Keff, split_ratio} encoding the SAME {Kc, Kn} — both must deploy to one law
# (H4's deploy back-maps to {Kc,Kn} = H1, so the micro block is identical).
@testset "HK1 deploy micro map reproduces the Cha law (H1 raw & H4 reparam)" begin
    mech = build_hk1_mechanism(alpha=:one, glc_g6p_dead_end=true, glc_adp_dead_end=true)
    c = (Glucose=5e-5, ATP=1e-3, G6P=3e-5, ADP=1e-3, Pi=5e-3)
    Kc, Kn = 15e-6, 6.9e-6                      # C = larger (loose), N = smaller (tight)
    # H1: raw coords
    coords1 = Dict(:Kd_Glc=>50e-6, :Kd_ATP=>1e-3, :Ki_G6P_C=>Kc, :Ki_ADP=>1.5e-3,
                   :Ki_G6P_N=>Kn, :K_Pi_N=>750e-6)
    # H4: reparameterized coords encoding the same {Kc,Kn}
    Keff  = 1 / (1/Kc + 1/Kn); sqrtP = sqrt(Kc * Kn)
    coords4 = Dict(:Kd_Glc=>50e-6, :Kd_ATP=>1e-3, :Keff=>Keff, :Ki_ADP=>1.5e-3,
                   :split_ratio=>sqrtP/Keff, :K_Pi_N=>750e-6)
    for (variant, coords) in ((:H1, coords1), (:H4, coords4))
        logθ = ChaDeploy.cha_deploy_micro(:HK1, mech, coords; keq=2700.0)
        kp   = ChaFit.cha_macro_tuple(:HK1, coords; keq=2700.0, variant=variant)
        p = MID.build_params(mech, logθ; keq=2700.0)
        v_micro = EnzymeRates.rate_equation(mech, c, p)
        v_cha   = cha_rate_HK1(kp; Glucose=c.Glucose, ATP=c.ATP, G6P=c.G6P, ADP=c.ADP, Pi=c.Pi)
        @test isapprox(v_micro, v_cha; rtol=1e-10)
    end
end

else
    @testset "HK1 deploy (skipped: HK1_AVAILABLE=false)" begin
        @test_skip "build_hk1_mechanism unavailable pending HK1 wiring port"
    end
end

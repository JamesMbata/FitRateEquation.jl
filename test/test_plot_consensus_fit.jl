# ##########################################################################################
#            Unit tests for the FitRateEquation Cha fit-vs-data plotter helpers             #
# ##########################################################################################
# The plotter (detect_enzyme, config_for, read_coords, build_cha_adapter, plot_consensus_fit)
# lives in the CairoMakie package extension (FitRateEquationMakieExt), ported in a later task.
# Until that extension exists, every assertion below is a "render" assertion in the sense of
# the port brief (it depends on extension-only helpers) and is gated off; this preserves the
# testset structure/names so a later task can drop the gate without touching the bodies.

using Test, DataFrames, FitRateEquation

const _MAKIE_EXT_LOADED = isdefined(Base, :get_extension) &&
    !isnothing(Base.get_extension(FitRateEquation, :FitRateEquationMakieExt))

@testset "plotter helpers" begin
    if _MAKIE_EXT_LOADED
        ext = Base.get_extension(FitRateEquation, :FitRateEquationMakieExt)
        detect_enzyme     = ext.detect_enzyme
        config_for        = ext.config_for
        read_coords       = ext.read_coords
        build_cha_adapter = ext.build_cha_adapter

        @testset "detect_enzyme" begin
            @test detect_enzyme("fitting/G6PD/rate_eq/consensus_macro/results/2026-06-11_full") == :G6PD
            @test detect_enzyme("/abs/PPP_Experiments/fitting/PGD/rate_eq/consensus_macro/results/2026-06-15") == :PGD
            @test detect_enzyme("fitting/HK1/rate_eq/consensus_macro/results/2026-06-13") == :HK1
            @test_throws ErrorException detect_enzyme("some/unrelated/path")
        end

        @testset "config_for" begin
            @test config_for(:G6PD).name == "G6PD"
            @test config_for(:PGD).name == "PGD"
            @test_throws ErrorException config_for(:NOPE)
        end

        @testset "read_coords" begin
            mc = DataFrame(
                variant = ["SS_NADPH_release_rate_eq", "SS_NADPH_release_rate_eq",
                           "SS_NADPH_release_rate_eq", "SS_NADPH_release_rate_eq",
                           "SS_NADPH_release_rate_eq", "SS_NADPH_release_rate_eq",
                           "SS_NADPH_release_rate_eq", "SS_NADPH_release_rate_eq",
                           "SS_NADPH_release_rate_eq"],
                mode = fill("mode1", 9),
                name = ["Kd_NADP","Kd_G6P","Kd_6PGLn","alpha","Ki_NADPH","Ki_ATP",
                        "Ki_ATP_EG","Km_NADPH_rev","Km_G6P"],           # last is a DERIVED readoff
                value = [5.8e-6, 2.7e-5, 9.9e-4, 3.4, 3.7e-5, 2.2e-3, 6.1e-3, 3.9e-6, 9.3e-5],
            )
            d = read_coords(mc, :G6PD, :SS_NADPH_release_rate_eq, :mode1)
            @test Set(keys(d)) == Set([:Kd_NADP,:Kd_G6P,:Kd_6PGLn,:alpha,:Ki_NADPH,
                                       :Ki_ATP,:Ki_ATP_EG,:Km_NADPH_rev])   # Km_G6P excluded
            @test d[:alpha] == 3.4
            # Missing a required coord -> error naming it
            mc_missing = mc[mc.name .!= "alpha", :]
            @test_throws ErrorException read_coords(mc_missing, :G6PD, :SS_NADPH_release_rate_eq, :mode1)
        end

        @testset "ChaAdapter numerics (G6PD)" begin
            # A realistic G6PD deploy coord set (mode1, 2026-06-11_full).
            coords = Dict(:Kd_NADP=>5.8175e-6, :Kd_G6P=>2.7447e-5, :Kd_6PGLn=>9.896e-4,
                          :alpha=>3.4062, :Ki_NADPH=>3.744e-5, :Ki_ATP=>2.1629e-3,
                          :Ki_ATP_EG=>6.0598e-3, :Km_NADPH_rev=>3.9e-6)
            a = build_cha_adapter(:G6PD, coords, :SS_NADPH_release_rate_eq, 43.743)

            @test a.enzyme === :G6PD
            @test Set(EnzymeRates.metabolites(a)) ==
                  Set([:NADP,:G6P,:NADPH,:PGLn,:ATP])

            # Forward, single-substrate-limited point: positive, finite rate.
            concs = (NADP = 1e-4, G6P = 1e-3, NADPH = 0.0, PGLn = 0.0, ATP = 0.0)
            v = EnzymeRates.rate_equation(a, concs, (Keq = 43.743, E_total = 1.0))
            @test isfinite(v) && v > 0.0

            # FORWARD rate is keq-independent (kr does not enter when products=0), so a
            # different params.Keq gives the SAME forward rate -- even though params.Keq is
            # now honored (see the reverse-arm subtest below).
            v2 = EnzymeRates.rate_equation(a, concs, (Keq = 1.0, E_total = 99.0))
            @test v2 == v

            # Deploy release_rate (1e3) is baked in: forward rate matches a direct cha_rate
            # call at the same macro tuple (guards against drawing at the wrong fiber).
            m_ref = FitRateEquation.ChaFit.cha_macro_tuple(:G6PD, coords;
                        keq = 43.743,
                        release_rate = FitRateEquation.ChaFit.CHA_DEPLOY_RELEASE_RATE,
                        variant = :SS_NADPH_release_rate_eq)
            v_ref = FitRateEquation.ChaLaws.cha_rate_G6PD(m_ref;
                        NADP = 1e-4, G6P = 1e-3, NADPH = 0.0, PGLn = 0.0, ATP = 0.0)
            @test v ≈ v_ref rtol=1e-12

            # PER-FIGURE KEQ IS HONORED on the reverse arm: with products present, kr (Haldane
            # from params.Keq) enters, so different params.Keq -> different rate, and each
            # matches a tuple rebuilt at that keq. This is the whole point of the fix.
            rev = (NADP = 0.0, G6P = 0.0, NADPH = 1e-4, PGLn = 1e-4, ATP = 0.0)
            vA = EnzymeRates.rate_equation(a, rev, (Keq = 43.743, E_total = 1.0))
            vB = EnzymeRates.rate_equation(a, rev, (Keq = 20.0,   E_total = 1.0))
            @test vA < 0.0 && vB < 0.0            # pure-reverse -> negative net rate
            @test vA != vB                        # kr tracks params.Keq (per-figure)
            mB = FitRateEquation.ChaFit.cha_macro_tuple(:G6PD, coords;
                        keq = 20.0,
                        release_rate = FitRateEquation.ChaFit.CHA_DEPLOY_RELEASE_RATE,
                        variant = :SS_NADPH_release_rate_eq)
            vB_ref = FitRateEquation.ChaLaws.cha_rate_G6PD(mB;
                        NADP = 0.0, G6P = 0.0, NADPH = 1e-4, PGLn = 1e-4, ATP = 0.0)
            @test vB ≈ vB_ref rtol=1e-12
        end
    else
        @test_skip "FitRateEquationMakieExt not loaded (Makie plotter ported in a later task)"
    end
end

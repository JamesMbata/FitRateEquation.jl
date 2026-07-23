# ##########################################################################################
#            Unit tests for the FitRateEquation Cha fit-vs-data plotter helpers             #
# ##########################################################################################
# The non-Makie helpers (detect_enzyme, config_for, read_coords, build_cha_adapter) live in
# the main module (src/plot_support.jl) and are reachable without loading CairoMakie, so
# these testsets run in the DEFAULT suite. The rendering method itself
# (FitRateEquation.plot_consensus_fit) is a stub in the main module, implemented only by the
# CairoMakie package extension (ext/FitRateEquationMakieExt.jl) -- see test_plot_render.jl
# (NOT run by runtests.jl) for the actual render check.

using Test, DataFrames, FitRateEquation, EnzymeRates

@testset "plotter helpers" begin
    @testset "detect_enzyme" begin
        @test FitRateEquation.detect_enzyme("fitting/G6PD/rate_eq/consensus_macro/results/2026-06-11_full") == :G6PD
        @test FitRateEquation.detect_enzyme("/abs/PPP_Experiments/fitting/PGD/rate_eq/consensus_macro/results/2026-06-15") == :PGD
        @test FitRateEquation.detect_enzyme("fitting/HK1/rate_eq/consensus_macro/results/2026-06-13") == :HK1
        @test_throws ErrorException FitRateEquation.detect_enzyme("some/unrelated/path")
    end

    @testset "config_for" begin
        @test FitRateEquation.config_for(:G6PD).name == "G6PD"
        @test FitRateEquation.config_for(:PGD).name == "PGD"
        @test_throws ErrorException FitRateEquation.config_for(:NOPE)
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
        d = FitRateEquation.read_coords(mc, :G6PD, :SS_NADPH_release_rate_eq, :mode1)
        @test Set(keys(d)) == Set([:Kd_NADP,:Kd_G6P,:Kd_6PGLn,:alpha,:Ki_NADPH,
                                   :Ki_ATP,:Ki_ATP_EG,:Km_NADPH_rev])   # Km_G6P excluded
        @test d[:alpha] == 3.4
        # Missing a required coord -> error naming it
        mc_missing = mc[mc.name .!= "alpha", :]
        @test_throws ErrorException FitRateEquation.read_coords(mc_missing, :G6PD, :SS_NADPH_release_rate_eq, :mode1)
    end

    @testset "ChaAdapter numerics (G6PD)" begin
        # A realistic G6PD deploy coord set (mode1, 2026-06-11_full).
        coords = Dict(:Kd_NADP=>5.8175e-6, :Kd_G6P=>2.7447e-5, :Kd_6PGLn=>9.896e-4,
                      :alpha=>3.4062, :Ki_NADPH=>3.744e-5, :Ki_ATP=>2.1629e-3,
                      :Ki_ATP_EG=>6.0598e-3, :Km_NADPH_rev=>3.9e-6)
        a = FitRateEquation.build_cha_adapter(:G6PD, coords, :SS_NADPH_release_rate_eq, 43.743)

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

    @testset "ChaAdapter numerics (PGD :full_re dispatches cha_rate_PGD_fullRE)" begin
        # :full_re coords (Run A mode1). The adapter must evaluate through cha_rate_PGD_fullRE,
        # NOT cha_rate_PGD -- the fiber-free full-RE law over its 6 core coords.
        coords = Dict(:Kd_NADP=>2.4638e-6, :Kd_PGA=>1.3504e-5, :alpha=>13.569,
                      :Kd_NADPH=>8.7249e-7, :Kd_Ru5P=>4.6517e-5, :Kd_CO2=>0.011088)
        a = FitRateEquation.build_cha_adapter(:PGD, coords, :full_re, 0.17)
        @test a.enzyme === :PGD && a.variant === :full_re
        @test Set(EnzymeRates.metabolites(a)) == Set([:NADP,:PGA,:Ru5P,:CO2,:NADPH,:ATP])

        concs = (NADP = 3e-5, PGA = 2e-4, Ru5P = 0.0, CO2 = 0.0, NADPH = 0.0, ATP = 0.0)
        v = EnzymeRates.rate_equation(a, concs, (Keq = 0.17, E_total = 1.0))
        @test isfinite(v) && v > 0.0
        # Matches a direct cha_rate_PGD_fullRE call at the DEPLOY-fiber full_re macro tuple.
        m_ref = FitRateEquation.ChaFit.cha_macro_tuple(:PGD, coords;
                    keq = 0.17, release_rate = FitRateEquation.ChaFit.CHA_DEPLOY_RELEASE_RATE,
                    variant = :full_re)
        v_ref = FitRateEquation.ChaLaws.cha_rate_PGD_fullRE(m_ref;
                    NADP = 3e-5, PGA = 2e-4, Ru5P = 0.0, CO2 = 0.0, NADPH = 0.0, ATP = 0.0)
        @test v ≈ v_ref rtol=1e-12

        # ADDITIVE GUARD: a :cha_base PGD adapter still dispatches cha_rate_PGD (unchanged).
        cb = Dict(:Kd_NADP=>3e-6, :Kd_PGA=>1.5e-5, :alpha=>2.0, :Kd_CO2=>1e-4,
                  :Ki_NADPH=>1.7e-5, :Ki_ATP=>1.7e-3, :Ki_ATP_EN=>1e-6, :Km_NADPH_rev=>2e-5)
        acb = FitRateEquation.build_cha_adapter(:PGD, cb, :cha_base, 0.17)
        vcb = EnzymeRates.rate_equation(acb, concs, (Keq = 0.17, E_total = 1.0))
        mcb = FitRateEquation.ChaFit.cha_macro_tuple(:PGD, cb;
                    keq = 0.17, release_rate = FitRateEquation.ChaFit.CHA_DEPLOY_RELEASE_RATE,
                    variant = :cha_base)
        vcb_ref = FitRateEquation.ChaLaws.cha_rate_PGD(mcb;
                    NADP = 3e-5, PGA = 2e-4, Ru5P = 0.0, CO2 = 0.0, NADPH = 0.0, ATP = 0.0)
        @test vcb ≈ vcb_ref rtol=1e-12
    end

    @testset "detect_enzyme maps the :full_re variant to PGD" begin
        # "full_re" is now in the variant->enzyme map, so a run dir need not sit under fitting/PGD/.
        @test FitRateEquation._VARIANT_TO_ENZYME["full_re"] == :PGD
    end
end

@testset "plot stub errors without CairoMakie" begin
    if isnothing(Base.get_extension(FitRateEquation, :FitRateEquationMakieExt))
        @test_throws MethodError plot_consensus_fit(mktempdir())
    else
        @test_skip "CairoMakie loaded; stub not exercised"
    end
end

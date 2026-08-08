# ##########################################################################################
#            Unit tests for the FitRateEquation Cha fit-vs-data plotter helpers             #
# ##########################################################################################
# The non-Makie helpers (detect_enzyme, config_for, read_coords, build_cha_adapter) live in
# the main module (src/plot_support.jl) and are reachable without loading CairoMakie, so
# these testsets run in the DEFAULT suite. The rendering method itself
# (FitRateEquation.plot_consensus_fit) is a stub in the main module, implemented only by the
# CairoMakie package extension (ext/FitRateEquationMakieExt.jl) -- see test_plot_render.jl
# (NOT run by runtests.jl) for the actual render check.

using Test, CSV, DataFrames, FitRateEquation, EnzymeRates

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

@testset "read_fit_corpus" begin
    # A run dir with a snapshot: returned verbatim.
    dir = mktempdir()
    CSV.write(joinpath(dir, "fit_corpus.csv"),
              DataFrame(NADP = [1e-6], G6P = [2e-4], NADPH = [0.0], PGLn = [0.0],
                        ATP = [0.0], Rate = [15.33], source = ["Wang2002|1A"],
                        Apparent_Keq = [43.743], X_axis_label = ["NADP"]))
    df = FitRateEquation.read_fit_corpus(dir)
    @test nrow(df) == 1
    @test df.source[1] == "Wang2002|1A"
    @test df.NADP[1] ≈ 1e-6

    # A run dir with NO snapshot (pre-0.2.0): a clear error, never a silent fallback to
    # the bundled corpus. This is the reported bug.
    empty_dir = mktempdir()
    err = try; FitRateEquation.read_fit_corpus(empty_dir); nothing; catch e; e; end
    @test err isa ErrorException
    @test occursin("fit_corpus.csv", err.msg)
    @test occursin("0.2.0", err.msg)

    # A snapshot with no X_axis_label (HK1-shaped): errors at PLOT time, naming the column.
    hk1_dir = mktempdir()
    CSV.write(joinpath(hk1_dir, "fit_corpus.csv"),
              DataFrame(Glucose = [1e-3], ATP = [1e-3], G6P = [0.0], ADP = [0.0],
                        Pi = [0.0], Rate = [1.0], source = ["Choe|1"],
                        Apparent_Keq = [2700.0]))
    err2 = try; FitRateEquation.read_fit_corpus(hk1_dir); nothing; catch e; e; end
    @test err2 isa ErrorException
    @test occursin("X_axis_label", err2.msg)
end

@testset "read_fit_corpus validates every renderer column" begin
    full = DataFrame(Rate=[1.0], source=["A|1"], Apparent_Keq=[13.7],
                     X_axis_label=["G6P"], G6P=[1e-4])

    # All four present: reads back fine.
    ok = mktempdir()
    CSV.write(joinpath(ok, "fit_corpus.csv"), full)
    @test nrow(FitRateEquation.read_fit_corpus(ok)) == 1

    # Each required column missing in turn: a clear error, NOT a downstream @warn. The
    # message must NAME the missing column — an ErrorException alone would also be satisfied
    # by some unrelated error escaping read_fit_corpus, which is exactly the confusion this
    # validation exists to remove.
    for c in (:Rate, :source, :Apparent_Keq, :X_axis_label)
        d = mktempdir()
        CSV.write(joinpath(d, "fit_corpus.csv"), full[!, Not(c)])
        err = try; FitRateEquation.read_fit_corpus(d); nothing; catch e; e; end
        @test err isa ErrorException
        @test occursin("no `$c` column", err.msg)
    end
end

@testset "read_corpus rejects a metabolite key that collides with a reserved column" begin
    cfg = g6pd_config()
    for reserved in FitRateEquation._RESERVED_CORPUS_COLS
        bad = (; cfg..., metabolites = Dict(reserved => ("[G6P] (uM)", :uM)))
        @test_throws ErrorException FitRateEquation.read_corpus(bad)
    end
    # The bundled config is unaffected.
    @test nrow(FitRateEquation.read_corpus(cfg)) > 0
end

@testset "write_outputs requires the corpus it fit" begin
    # Every argument is valid, so this asserts the corpus guard specifically rather than
    # passing because some earlier argument blew up — and it stays correct wherever in the
    # body the guard sits. Matching the message means the test fails loudly if the guard is
    # ever reworded or removed, instead of being satisfied by an unrelated error.
    d = FitRateEquation.Dataset([(; NADP=1e-5, G6P=1e-4)], [1.0], ["Article|1"], [13.7])
    @test_throws "write_outputs requires `corpus=`" write_outputs(mktempdir(), d, NamedTuple[])
end

@testset "_toml_escape quotes backslashes and quotes" begin
    @test FitRateEquation._toml_escape("C:\\data\\corpus.csv") == "C:\\\\data\\\\corpus.csv"
    @test FitRateEquation._toml_escape("odd\"name.csv")        == "odd\\\"name.csv"
    @test FitRateEquation._toml_escape("/plain/path.csv")      == "/plain/path.csv"
end

@testset "dataset_from_corpus accepts a SubDataFrame view" begin
    cfg = g6pd_config()
    df  = FitRateEquation.read_corpus(cfg)
    @test nrow(df) > 2
    d = FitRateEquation.dataset_from_corpus(view(df, 1:2, :), cfg)
    @test FitRateEquation.nrows(d) == 2
end

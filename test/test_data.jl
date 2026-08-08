using FitRateEquation
using EnzymeRates
using CSV
using DataFrames
using Test


@testset "data" begin
    d = load_dataset(g6pd_config())
    @test nrows(d) > 100
    @test all(occursin("|", g) for g in d.group)   # group key is Article|Fig
    @test isconcretetype(eltype(d.concs))           # Task 11: concs stays concretely typed
end

@testset "read_corpus / dataset_from_corpus" begin
    cfg = g6pd_config()
    df  = FitRateEquation.read_corpus(cfg)

    # Columns: every metabolite symbol, plus the four renderer/loader columns.
    for s in keys(cfg.metabolites)
        @test s in propertynames(df)
    end
    @test :Rate in propertynames(df)
    @test :source in propertynames(df)
    @test :Apparent_Keq in propertynames(df)
    @test :X_axis_label in propertynames(df)      # G6PD corpus has it

    # Bad-rate rows are dropped, exactly as load_dataset drops them.
    @test all(isfinite, df.Rate)
    @test all(!=(0.0), df.Rate)

    # Explicit first-row values, independent of any implementation: the bundled G6PD
    # corpus row 1 is Wang2002|1A, NADP 1 uM, G6P 200 uM, Rate_V 15.33, Keq 43.743.
    @test df.source[1]       == "Wang2002|1A"
    @test df.NADP[1]         ≈ 1e-6            # uM -> M
    @test df.G6P[1]          ≈ 200e-6
    @test df.NADPH[1]        == 0.0
    @test df.Rate[1]         ≈ 15.33
    @test df.Apparent_Keq[1] ≈ 43.743
    @test df.X_axis_label[1] == "NADP"

    # DIFFERENTIAL CHECK: the composition reproduces `load_dataset` field-by-field on
    # every bundled corpus. As of Task 2 `load_dataset` IS this composition, so both
    # sides run the same code path and this is trivially true today. KEEP IT ANYWAY --
    # it is a regression guard, not a tautology: it fires the moment anyone reintroduces
    # a second, independent corpus-loading implementation behind `load_dataset`, which is
    # the exact drift this branch exists to remove. (The in-repo evidence that the loaded
    # VALUES did not move is elsewhere in this file: the literal first-row assertions above
    # -- source/concs/rate/Keq/X_axis_label pinned to the corpus by hand -- and the
    # "read_corpus drop predicate (fixture corpus)" testset below. Not test_byte_identity.jl:
    # by default that compares only variant/mode/name, with value/class/ci behind an opt-in
    # env var, so it does not gate loaded values on a normal run.)
    for c in (g6pd_config(), pgd_config(), hk1_config())
        ref = load_dataset(c)
        got = FitRateEquation.dataset_from_corpus(FitRateEquation.read_corpus(c), c)
        @test nrows(got) == nrows(ref)
        @test got.concs  == ref.concs
        @test got.rate   == ref.rate
        @test got.group  == ref.group
        @test isequal(got.keq, ref.keq)          # isequal: NaN == NaN
        @test eltype(got.concs) === eltype(ref.concs)
        @test isconcretetype(eltype(got.concs))
    end

    # HK1's corpus has no X_axis_label; read_corpus must not invent one.
    @test !(:X_axis_label in propertynames(FitRateEquation.read_corpus(hk1_config())))

    # read_corpus uses load_dataset's `_to_float(x, 0.0)` for concentration cells. The
    # deleted build_plot_df used NaN-then-zero, which differs ONLY for a literal NaN cell.
    # Verified 2026-08-07: no bundled corpus contains one, so the unification is a no-op on
    # real data. This pins that -- if a future corpus introduces a NaN concentration cell,
    # this fires instead of the plot data silently changing.
    for c in (g6pd_config(), pgd_config(), hk1_config())
        raw = CSV.read(c.data_csv, DataFrame)
        for (s, (col, unit)) in c.metabolites
            @test !any(x -> x isa Real && isnan(x), raw[!, col])
        end
    end
end

@testset "read_corpus drop predicate (fixture corpus)" begin
    # The bundled corpora drop zero rows, so they cannot exercise the rate filter.
    # This fixture carries one row per drop path plus two survivors.
    csv = """
    [NADP] (uM),[G6P] (uM),[NADPH] (uM),[PGLn] (uM),[ATP] (uM),Rate_V,Article,Fig,X_axis_label,Apparent_Keq
    1,200,0,0,0,15.33,Wang2002,1A,NADP,43.743
    2,200,0,0,0,0,Wang2002,1A,NADP,43.743
    3,200,0,0,0,,Wang2002,1A,NADP,43.743
    4,200,0,0,0,n/a,Wang2002,1A,NADP,43.743
    5,,0,0,0,-2.5,Wang2002,1B,NADP,43.743
    """
    path = joinpath(mktempdir(), "fixture_corpus.csv")
    write(path, csv)
    cfg = g6pd_config(; data_csv=path)

    df = FitRateEquation.read_corpus(cfg)
    @test nrow(df) == 2                       # zero, blank and non-numeric rates dropped
    @test df.Rate == [15.33, -2.5]            # negative (reverse) rates SURVIVE
    @test df.source == ["Wang2002|1A", "Wang2002|1B"]
    @test df.NADP == [1e-6, 5e-6]             # surviving rows keep their own values
    @test df.G6P[2] == 0.0                    # blank concentration cell -> 0.0

    # The differential check, on a corpus that actually drops rows.
    ref = load_dataset(cfg)
    got = FitRateEquation.dataset_from_corpus(df, cfg)
    @test nrows(got) == nrows(ref) == 2
    @test got.concs == ref.concs
    @test got.rate  == ref.rate
    @test got.group == ref.group
end

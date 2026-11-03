using FitRateEquation
using FitRateEquation: g6pd_config
using EnzymeRates
using CSV, DataFrames
using Test


@testset "fit_consensus_equation outputs (Cha macro-coord)" begin
    outdir = mktempdir()
    res = FitRateEquation._fit_consensus(g6pd_config(); outdir=outdir, n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    for f in ("macro_constants.csv","goodness_of_fit.csv","cv_folds.csv","fit_corpus.csv",
              "identifiable_functions.csv","micro_parameters.jl","report.md")
        @test isfile(joinpath(outdir, f))
    end
    # Deploy-variant-only: G6PD = 1 variant (:SS_NADPH_release_rate_eq) × 2 modes = 2 cells.
    @test length(res) == 2
    @test all(r -> r.variant === :SS_NADPH_release_rate_eq, res)
    # Both CV schemes, clearly labelled, article (headline) first; per-fold losses written.
    gof = CSV.read(joinpath(outdir, "goodness_of_fit.csv"), DataFrame)
    # propertynames, not names: this testset binds a local `names` below, shadowing DataFrames.names.
    @test propertynames(gof) == [:variant, :mode, :in_sample,
                                 :cv_article_mean, :cv_article_se, :cv_article_n_folds,
                                 :cv_figure_mean, :cv_figure_se, :cv_figure_n_folds]
    @test nrow(gof) == 2
    @test all(isfinite, gof.cv_article_mean) && all(isfinite, gof.cv_figure_mean)
    folds = CSV.read(joinpath(outdir, "cv_folds.csv"), DataFrame)
    @test propertynames(folds) == [:variant, :mode, :scheme, :held_out, :n_test_rows, :loss]
    d = load_dataset(g6pd_config())
    for r in res
        sub = folds[string.(folds.mode) .== String(r.mode), :]
        # Article folds first, then figure folds, each in data order.
        @test sub.scheme == vcat(fill("article", r.r.cv.article.n_folds),
                                 fill("figure",  r.r.cv.figure.n_folds))
        fig = sub[sub.scheme .== "figure", :]
        @test fig.held_out == unique(d.group)
        @test fig.n_test_rows == [count(==(g), d.group) for g in unique(d.group)]
        @test fig.loss == [x.loss for x in r.r.cv.figure.per_fold]
        g = gof[string.(gof.mode) .== String(r.mode), :]
        @test g.cv_article_n_folds[1] == r.r.cv.article.n_folds == count(==("article"), sub.scheme)
        @test g.cv_figure_mean[1] == r.r.cv.figure.mean_cv
    end
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    @test :mode in propertynames(mc)
    @test Set(string.(mc.mode)) == Set(["mode1","mode2"])
    # Forward shape constants: the Cha coords (Kd_*/alpha/Ki_*) plus the DERIVED apparent
    # Km_G6P readoff. Km_G6P in-band (~10–200 µM → between 1e-6 and 1e-2 M).
    names = Set(string.(mc.name))
    @test "Km_G6P" in names && "Kd_G6P" in names && "Kd_NADP" in names && "Ki_ATP" in names
    km = mc[(mc.name .== "Km_G6P") .& (string.(mc.mode) .== "mode1"), :]
    @test nrow(km) >= 1
    @test all(1e-6 .< km.value .< 1e-2)
    # Classes are the Cha macro vocabulary (derived apparent Km's carry :derived).
    @test issubset(Set(string.(mc.class)),
        Set(["data_identified","unconstrained","literature_pinned","derived"]))
    # Ki_ATP literature-pinned in Mode 2 (a hard Cha coord pin).
    ki_atp_m2 = mc[(mc.name .== "Ki_ATP") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_atp_m2) >= 1 && all(string.(ki_atp_m2.class) .== "literature_pinned")
    # Forward Ki_NADPH is a Cha coord (E·G6P·NADPH dead-end cross term): :literature_pinned
    # (15 µM) in Mode 2, data-driven in Mode 1; NEVER :conflated_reverse.
    ki_nadph_m2 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_nadph_m2) >= 1 && all(string.(ki_nadph_m2.class) .== "literature_pinned")
    @test all(isapprox.(ki_nadph_m2.value, 15e-6; rtol=1e-6))
    ki_nadph_m1 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode1"), :]
    @test nrow(ki_nadph_m1) >= 1
    @test all(in.(string.(ki_nadph_m1.class), Ref(Set(["data_identified","unconstrained"]))))
    @test !any(string.(mc.class) .== "conflated_reverse")
    # report.md carries the 2-way (mode1<->mode2) agreement + the G6PD koffQ hybrid block.
    rep = read(joinpath(outdir,"report.md"), String)
    @test occursin("agreement across modes", rep)
    @test occursin("leave-one-article-out (headline)", rep)
    @test occursin("leave-one-figure-out", rep)
    @test findfirst("leave-one-article-out", rep).start < findfirst("leave-one-figure-out", rep).start
    @test occursin("mode1 ↔ mode2", rep)
    @test occursin("koffQ hybrid", rep)
    @test occursin("data-identified koffQ", rep)
    # micro_parameters.jl carries the closed-form deploy block for the deploy variant.
    mp = read(joinpath(outdir, "micro_parameters.jl"), String)
    @test occursin("SS_NADPH_release_rate_eq", mp)
    @test occursin("cha_deploy_micro", mp)
end

@testset "anchor_reverse: default off, data-driven NOT DEPLOYABLE banner" begin
    # Default run (anchor_reverse=false since 0.9.0): provenance records it, and each variant ×
    # mode block carries the banner exactly when that fit's forward-NADPH coord classified
    # :unconstrained (_nadph_unidentified) -- never by variant name.
    outdir1 = mktempdir()
    res1 = FitRateEquation._fit_consensus(g6pd_config(); outdir=outdir1,
               n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    @test occursin("anchor_reverse  = false", read(joinpath(outdir1, "provenance.toml"), String))
    mp1  = read(joinpath(outdir1, "micro_parameters.jl"), String)
    rep1 = read(joinpath(outdir1, "report.md"), String)
    for r in res1
        flagged = FitRateEquation._nadph_unidentified(:G6PD, false, r.classed)
        @test occursin("NOT DEPLOYABLE [$(r.variant) $(r.mode)]", mp1) == flagged
    end
    @test occursin("NOT DEPLOYABLE", rep1) ==
          any(r -> FitRateEquation._nadph_unidentified(:G6PD, false, r.classed), res1)
    # Mode 2 pins Ki_NADPH to literature, so its forward-NADPH coord is never unconstrained.
    @test !occursin("NOT DEPLOYABLE [SS_NADPH_release_rate_eq mode2]", mp1)

    # Opt-in anchor: recorded in provenance, Km_NADPH_rev literature-pinned, never a banner.
    outdir2 = mktempdir()
    FitRateEquation._fit_consensus(g6pd_config(); outdir=outdir2, n_restarts=2, maxiter=150,
            maxtime=5.0, seed=1, anchor_reverse=true)
    @test occursin("anchor_reverse  = true", read(joinpath(outdir2, "provenance.toml"), String))
    @test !occursin("NOT DEPLOYABLE", read(joinpath(outdir2, "micro_parameters.jl"), String))
    @test !occursin("NOT DEPLOYABLE", read(joinpath(outdir2, "report.md"), String))
    mc2 = CSV.read(joinpath(outdir2, "macro_constants.csv"), DataFrame)
    kr  = mc2[mc2.name .== "Km_NADPH_rev", :]
    @test nrow(kr) == 2 && all(string.(kr.class) .== "literature_pinned")
end

@testset "provenance records deploy_keq" begin
    d    = FitRateEquation.load_dataset(g6pd_config())
    meta = (n_restarts=1, maxiter=100, maxtime=1.0, seed=1,
            n_rows=FitRateEquation.nrows(d))
    dir  = mktempdir()
    FitRateEquation._write_provenance(dir, d, meta; deploy_keq=13.655)
    txt = read(joinpath(dir, "provenance.toml"), String)
    @test occursin("deploy_keq", txt)
    @test occursin("13.655", txt)
end

@testset "fit_corpus.csv snapshots the corpus actually fit" begin
    bundled = CSV.read(g6pd_config().data_csv, DataFrame)

    # A genuinely different corpus: three of the seven articles. Three (not one) so
    # leave-one-article-out CV still has non-empty train and test sets per fold.
    keep = ["Wang2002", "Ozer2002", "Beutler1986"]
    sub  = bundled[in.(string.(bundled.Article), Ref(Set(keep))), :]
    @test nrow(sub) < nrow(bundled)
    tmpcsv = joinpath(mktempdir(), "subset_corpus.csv")
    CSV.write(tmpcsv, sub)

    outdir = mktempdir()
    FitRateEquation._fit_consensus(g6pd_config(; data_csv=tmpcsv); outdir=outdir,
            n_restarts=2, maxiter=150, maxtime=5.0, seed=1)

    @test isfile(joinpath(outdir, "fit_corpus.csv"))
    fc = CSV.read(joinpath(outdir, "fit_corpus.csv"), DataFrame)

    # THE REGRESSION: the snapshot is the corpus that was passed in, not the default.
    # nrow > 0 first — `all` over an empty collection is true, so an empty snapshot would
    # otherwise satisfy both of the next two assertions.
    @test nrow(fc) > 0
    @test all(a -> a in keep, first.(split.(string.(fc.source), "|")))
    @test nrow(fc) < nrow(bundled)

    # It carries everything the plot renderer needs.
    for c in (:Rate, :source, :Apparent_Keq, :X_axis_label, :NADP, :G6P, :NADPH, :PGLn, :ATP)
        @test c in propertynames(fc)
    end

    # Row count agrees with what the fit actually saw.
    prov = read(joinpath(outdir, "provenance.toml"), String)
    @test occursin("n_rows          = $(nrow(fc))\n", prov)

    # Provenance records where the corpus came from.
    @test occursin("data_csv        = \"$(abspath(tmpcsv))\"", prov)
end

@testset "fit_corpus.csv honors row_filter" begin
    outdir = mktempdir()
    FitRateEquation._fit_consensus(g6pd_config(); outdir=outdir, n_restarts=2, maxiter=150, maxtime=5.0, seed=1,
            variants=[:no_atp], row_filter=FitRateEquation.drop_atp_rows)
    fc = CSV.read(joinpath(outdir, "fit_corpus.csv"), DataFrame)
    @test nrow(fc) > 0
    @test all(fc.ATP .<= 0.0)                 # the filtered-out rows are really gone
    @test occursin("n_rows          = $(nrow(fc))\n",
                   read(joinpath(outdir, "provenance.toml"), String))
end

@testset "fit_corpus.csv column order is deterministic" begin
    cfg  = g6pd_config()
    mets = FitRateEquation.metabolite_syms(cfg)
    @test issorted(mets)
    @test Set(mets) == Set(keys(cfg.metabolites))

    outdir = mktempdir()
    FitRateEquation._fit_consensus(cfg; outdir=outdir, n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    fc = CSV.read(joinpath(outdir, "fit_corpus.csv"), DataFrame)

    # The metabolite columns come first, in metabolite_syms order — not Dict hash order.
    @test propertynames(fc)[1:length(mets)] == mets

    # The fit path's own ordering: `d.concs`'s element type is a NamedTuple whose FIELD ORDER
    # is a type parameter built from the same list. Nothing above observes it — the corpus
    # columns and the adapter can both be sorted while this one silently is not.
    d = FitRateEquation.load_dataset(cfg)
    @test collect(propertynames(d.concs[1])) == mets     # the fit-path type parameter

    # The renderer's adapter sweeps the same order, so the plotter and the fit cannot
    # disagree about which column is which.
    a = FitRateEquation.build_cha_adapter(:G6PD, Dict{Symbol,Float64}(), :SS_NADPH_release_rate_eq, 13.7)
    @test collect(EnzymeRates.metabolites(a)) == mets
end

@testset "absolute outputs record scale + kcat verdict" begin
    out = mktempdir()
    fit_consensus_equation(:g6pd; scale=:absolute, smoke=true, nprocs=1, outdir=out,
        data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"),
        pins=Dict(:Kd_6PGLn=>log10(2.1e-4), :Km_NADPH_rev=>log10(3.9e-6)))
    prov = read(joinpath(out, "provenance.toml"), String)
    # Assert the actual TOML field, not a path substring (the branch dir name contains both
    # "scale" and "absolute", which would make a loose occursin pass spuriously).
    @test occursin(r"scale\s*=\s*\"absolute\"", prov)
    report = read(joinpath(out, "report.md"), String)
    # The kcat VERDICT line (not just the always-present kcatKm_G6P derived row).
    @test occursin(r"kcat.*s⁻¹", report) || occursin("in-band", report) || occursin("OUT OF BAND", report)
end

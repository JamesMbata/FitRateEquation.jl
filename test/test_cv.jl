using FitRateEquation
using FitRateEquation: g6pd_config
using EnzymeRates
using Statistics
using Test


@testset "article CV" begin
    d = load_dataset(g6pd_config())
    arts = FitRateEquation._article.(d.group)
    @test length(unique(arts)) >= 5      # 6 articles in the bundled corpus
    # The fold iterator holds out whole articles: no article in both train and test.
    for fold in FitRateEquation._article_folds(d)
        train_arts = Set(FitRateEquation._article.(d.group[fold.train]))
        test_arts  = Set(FitRateEquation._article.(d.group[fold.test]))
        @test isempty(intersect(train_arts, test_arts))
        @test test_arts == Set([fold.held_out])
    end
end

@testset "_group_folds: one fold per unique group, train/test partition" begin
    cfg = FitRateEquation.g6pd_config(
        data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    folds = FitRateEquation._group_folds(d)
    @test length(folds) == length(unique(d.group))
    for f in folds
        @test sort(vcat(f.train, f.test)) == collect(1:FitRateEquation.nrows(d))
        @test isempty(intersect(f.train, f.test))
        @test all(==(f.held_out), d.group[f.test])
    end
end

@testset "_cv_folds: both schemes, EnzymeRates-compatible figure folds" begin
    d = load_dataset(g6pd_config())
    @test FitRateEquation.CV_SCHEMES == (:article, :figure)   # headline first

    art = FitRateEquation._cv_folds(d, :article)
    fig = FitRateEquation._cv_folds(d, :figure)
    # Article folds: one per article, in first-appearance order.
    @test [f.held_out for f in art] == unique(FitRateEquation._article.(d.group))
    # Figure folds match EnzymeRates' LOOCV exactly: one per `group` in `unique(group)`
    # (data) order, 1-row groups included, held-out rows = that group's rows.
    @test [f.held_out for f in fig] == unique(d.group)
    for f in fig
        @test f.test == findall(==(f.held_out), d.group)
        @test f.train == findall(!=(f.held_out), d.group)
    end
    # More figure folds than article folds whenever some article has >1 figure.
    @test length(fig) > length(art)

    @test_throws ErrorException FitRateEquation._cv_folds(d, :bogus)
end

@testset "_cv_folds: degenerate folds dropped (single-article / single-group corpus)" begin
    full = load_dataset(g6pd_config())
    # One article with several figures (absolute mode's shape): leave-one-article-out has no
    # train set, leave-one-figure-out is a genuine CV.
    d = FitRateEquation._subset(full,
            findall(g -> FitRateEquation._article(g) == "Beutler1986", full.group))
    @test length(unique(d.group)) > 1
    @test isempty(FitRateEquation._cv_folds(d, :article))
    @test length(FitRateEquation._cv_folds(d, :figure)) == length(unique(d.group))

    one = FitRateEquation._subset(d, findall(==(d.group[1]), d.group))
    @test isempty(FitRateEquation._cv_folds(one, :figure))   # sole group => empty train
    @test length(FitRateEquation._group_folds(one)) == 1     # ...which _group_folds keeps
end

@testset "_cv_summary" begin
    per = [(held_out="a", n_test=3, loss=1.0), (held_out="b", n_test=4, loss=3.0)]
    s = FitRateEquation._cv_summary(per)
    @test s.mean_cv == 2.0
    @test s.se == std([1.0, 3.0]) / sqrt(2)
    @test s.n_folds == 2
    @test s.per_fold === per
    e = FitRateEquation._cv_summary(NamedTuple[])
    @test isnan(e.mean_cv) && e.se == 0.0 && e.n_folds == 0
end

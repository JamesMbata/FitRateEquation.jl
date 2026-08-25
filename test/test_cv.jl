using FitRateEquation
using FitRateEquation: g6pd_config
using EnzymeRates
using Test


@testset "article CV" begin
    d = load_dataset(g6pd_config())
    arts = FitRateEquation._article.(d.group)
    @test length(unique(arts)) >= 5      # ~7 articles
    # The fold iterator holds out whole articles: no article in both train and test.
    for fold in FitRateEquation._article_folds(d)
        train_arts = Set(FitRateEquation._article.(d.group[fold.train]))
        test_arts  = Set(FitRateEquation._article.(d.group[fold.test]))
        @test isempty(intersect(train_arts, test_arts))
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
    end
end

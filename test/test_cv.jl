using FitRateEquation
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

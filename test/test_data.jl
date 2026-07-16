using FitRateEquation
using EnzymeRates
using Test


@testset "data" begin
    d = load_dataset(g6pd_config())
    @test nrows(d) > 100
    @test all(occursin("|", g) for g in d.group)   # group key is Article|Fig
end

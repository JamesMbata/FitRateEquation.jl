using Test, FitRateEquation

@testset "run_g6pd smoke writes artifacts" begin
    out = mktempdir()
    res = run_g6pd(smoke=true, nprocs=1, outdir=out)
    @test res !== nothing
    @test isfile(joinpath(out, "macro_constants.csv"))
    @test isfile(joinpath(out, "micro_parameters.jl"))
    @test filesize(joinpath(out, "micro_parameters.jl")) > 0   # deploy block non-empty
end

@testset "run_hk1 errors clearly while guarded" begin
    FitRateEquation.HK1_AVAILABLE || @test_throws ErrorException run_hk1(smoke=true, nprocs=1, outdir=mktempdir())
end

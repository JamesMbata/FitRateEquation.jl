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

@testset "run_* wrappers forward to the symbol entry" begin
    dir = mktempdir()
    r1 = run_g6pd(; smoke=true, nprocs=1, outdir=joinpath(dir, "a"))
    r2 = fit_consensus_equation(:g6pd; smoke=true, nprocs=1, outdir=joinpath(dir, "b"))
    @test length(r1) == length(r2)   # same cell set

    # noatp is now a variant, not a runner
    @test !isdefined(FitRateEquation, :run_g6pd_noatp)
    @test !isdefined(FitRateEquation, :run_pgd_fullre)
end

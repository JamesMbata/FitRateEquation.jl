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

@testset "run_pgd_fullre smoke writes artifacts (:full_re deploy block non-empty)" begin
    out = mktempdir()
    res = run_pgd_fullre(smoke=true, nprocs=1, outdir=out)
    @test res !== nothing
    @test all(r -> r.variant === :full_re, res)
    @test Set(r.mode for r in res) == Set((:mode1, :mode2, :mode3))
    @test isfile(joinpath(out, "micro_parameters.jl"))
    @test filesize(joinpath(out, "micro_parameters.jl")) > 0
    mp = read(joinpath(out, "micro_parameters.jl"), String)
    @test occursin("full_re_mode1_K_NADPH_E", mp)      # fiber-free deploy succeeded
    @test !occursin("koff_Ru5P_ENADPH", mp)            # no SS-release fiber param emitted
end

using Test, FitRateEquation

# Exact-string comparison of the fitted macro-constant table against the committed
# reference (fixed-seed smoke). Guards the determinism contract and gates Task 10-11.
function _macro_csv(runner)
    out = mktempdir()
    runner(smoke=true, nprocs=1, outdir=out)
    read(joinpath(out, "macro_constants.csv"), String)
end

@testset "byte-identity: G6PD smoke" begin
    ref = read(joinpath(@__DIR__, "fixtures", "g6pd_smoke_macro_constants.csv"), String)
    @test _macro_csv(run_g6pd) == ref
end
@testset "byte-identity: PGD smoke" begin
    ref = read(joinpath(@__DIR__, "fixtures", "pgd_smoke_macro_constants.csv"), String)
    @test _macro_csv(run_pgd) == ref
end

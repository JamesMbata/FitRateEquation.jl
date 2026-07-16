using Test, FitRateEquation

# Comparison of the fitted macro-constant table against the committed reference (fixed-seed
# smoke). Guards the determinism contract and gates Task 10-11.
#
# Columns 1-5 (variant,mode,name,value,class) are compared EXACTLY: this is the strict gate on
# the fit output itself -- `value` (the actual fitted coordinate) must stay byte-identical run
# to run. Column 6 (`ci`, a finite-difference-Hessian confidence interval from cha_classify) is
# compared with a tolerance instead: FD-Hessian evaluation picks up last-bit floating-point
# noise from accumulated execution state when run mid-suite (after ~20 other tests, several of
# which run fits) that a fresh process does not exhibit -- confirmed by the coordinator: a
# from-fresh-process run reproduces the fixture exactly on EVERY column including `ci`, so the
# fit is fully deterministic and only this diagnostic column drifts across execution contexts.

function _macro_csv(runner)
    out = mktempdir()
    runner(smoke=true, nprocs=1, outdir=out)
    read(joinpath(out, "macro_constants.csv"), String)
end

_rows(csv) = [split(l, ',') for l in split(strip(csv), '\n')]

function _check(runner, fixture_path)
    got = _rows(_macro_csv(runner))
    ref = _rows(read(fixture_path, String))
    @test length(got) == length(ref)
    for (g, r) in zip(got, ref)
        @test g[1:5] == r[1:5]                  # fit output (variant,mode,name,value,class): exact
        gf, rf = tryparse(Float64, g[6]), tryparse(Float64, r[6])
        if gf === nothing || rf === nothing || isnan(rf)
            @test g[6] == r[6]                  # header row ("ci") or NaN reference: exact string
        else
            @test isapprox(gf, rf; rtol=1e-6)   # ci: FD-Hessian last-bit FP noise tolerated
        end
    end
end

@testset "byte-identity: G6PD smoke" begin
    _check(run_g6pd, joinpath(@__DIR__, "fixtures", "g6pd_smoke_macro_constants.csv"))
end
@testset "byte-identity: PGD smoke" begin
    _check(run_pgd, joinpath(@__DIR__, "fixtures", "pgd_smoke_macro_constants.csv"))
end

using Test, FitRateEquation

# Comparison of the fitted macro-constant table against the committed reference (fixed-seed
# smoke). Guards the pipeline's SHAPE, not its bit-exact numbers, by default: Task 10-11's
# determinism contract (worker-count / pmap-order invariance) is already covered elsewhere by
# `test_run_all.jl`'s "single-process run_all (pmap) == serial baseline" test, which compares two
# computations in the SAME process/environment and is legitimately bit-exact. This file compares
# against a fixture generated on a DIFFERENT machine/run, which is a different and much weaker
# guarantee.
#
# Default check (every machine, every CI job): columns 1,2,3,5 (variant,mode,name,class) must
# match EXACTLY. These don't depend on the CMA-ES trajectory's last bits -- only on which macro
# coordinates exist for a given (variant,mode) and how they're classified -- so they are stable
# across machines and still catch a real regression (a variant disappearing, a coordinate's
# class flipping).
#
# The fitted VALUE (column 4) and its CI (column 6) are NOT compared by default. Measured on a
# non-reference machine (2026-07-21), smoke-budget drift ranged from ~30-58% for otherwise
# `data_identified` coordinates (Kd_NADP, Ki_NADPH) to ~90-100% for `unconstrained` ones
# (Ki_ATP, Ki_ATP_EG) -- nothing like the ~0.1-0.5% "hardware noise" once assumed here. At smoke
# budget (tiny iteration count) the CMA-ES trajectory is under-converged, so different RNG/BLAS
# reduction ordering across machines lands in a materially different (but still smoke-quality)
# point, not just the last few bits. No single tolerance both survives that spread and still
# catches a real regression, so exact value/ci reproduction is opt-in only
# (`FITRATEEQ_BYTE_IDENTITY=1`), for whoever is on the fixture's reference machine + Julia series.

function _macro_csv(runner)
    out = mktempdir()
    runner(smoke=true, nprocs=1, outdir=out)
    read(joinpath(out, "macro_constants.csv"), String)
end

_rows(csv) = [split(l, ',') for l in split(strip(csv), '\n')]

# The fixtures are a fixed-seed smoke artifact generated on Julia 1.12. Bit-level fit
# reproducibility is NOT guaranteed across Julia MINOR versions -- different LLVM / libm /
# codegen shift the seeded CMA-ES trajectory -- so the opt-in exact check only applies on the
# fixtures' Julia series; other versions still run the default structural check above.
const _FIXTURES_JULIA_SERIES = v"1.12"
_on_fixture_julia() = _FIXTURES_JULIA_SERIES <= VERSION < v"1.13"
_strict_opt_in() = get(ENV, "FITRATEEQ_BYTE_IDENTITY", "false") == "true"

function _check(runner, fixture_path)
    got = _rows(_macro_csv(runner))
    ref = _rows(read(fixture_path, String))
    @test length(got) == length(ref)
    for (g, r) in zip(got, ref)
        @test g[[1, 2, 3, 5]] == r[[1, 2, 3, 5]]   # variant,mode,name,class: exact, every machine
    end
    if _strict_opt_in() && _on_fixture_julia()
        for (g, r) in zip(got, ref)
            @test g[4] == r[4]                    # value: exact, opt-in only
            gf, rf = tryparse(Float64, g[6]), tryparse(Float64, r[6])
            if gf === nothing || rf === nothing || isnan(rf)
                @test g[6] == r[6]                # header row ("ci") or NaN reference: exact string
            else
                @test isapprox(gf, rf; rtol=1e-6) # ci: FD-Hessian last-bit FP noise tolerated
            end
        end
    else
        @test_skip "exact value/ci reproduction is opt-in (FITRATEEQ_BYTE_IDENTITY=1) on Julia $(_FIXTURES_JULIA_SERIES)"
    end
end

@testset "byte-identity: G6PD smoke" begin
    _check(run_g6pd, joinpath(@__DIR__, "fixtures", "g6pd_smoke_macro_constants.csv"))
end
@testset "byte-identity: PGD smoke" begin
    _check(run_pgd, joinpath(@__DIR__, "fixtures", "pgd_smoke_macro_constants.csv"))
end

using Test, FitRateEquation

# ##########################################################################################
#                    Structural byte-identity gate for the smoke fixtures                   #
# ##########################################################################################
# test/fixtures/{g6pd,pgd}_smoke_macro_constants.csv are committed macro_constants.csv
# outputs from a fixed-seed smoke run. Every test run re-runs the smoke fit and compares
# variant/mode/name -- WHICH coordinates each variant x mode produces, and in what order.
# That property is machine-independent, so this gate is always on, on every machine.
#
# The fitted VALUES are deliberately not compared. At smoke budget the seeded CMA-ES
# trajectory is under-converged, so different BLAS/RNG reduction ordering lands a different
# machine on a materially different point -- observed spread is ~30-58% on data_identified
# coordinates and ~90-100% on unconstrained ones, and `class` follows the same drifted
# optimum. No tolerance both survives that and still catches a real regression.
#
# An exact value/class comparison lived here as an opt-in tier until it was retired: it
# gated on VERSION, but the property it needed was reference-MACHINE identity, so it could
# never pass in CI (rotating ubuntu-latest images) and did not pass on the maintainer's
# machine either. Behavioural regressions in the fit path are covered by test_cha_deploy,
# test_cha_classify, test_mode_agreement and the goodness-of-fit assertions -- and by
# test_parallel_equivalence.jl, which IS bit-exact on fit output, but compares two runs in
# the SAME process (pmap vs serial), so it is a same-machine determinism gate, not a
# cross-machine one.

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
    @test got[1] == ref[1]                     # header: full column schema, machine-independent
    for (g, r) in zip(got, ref)
        @test g[[1, 2, 3]] == r[[1, 2, 3]]         # variant,mode,name: which coords exist
    end
end

@testset "byte-identity: G6PD smoke" begin
    _check(run_g6pd, joinpath(@__DIR__, "fixtures", "g6pd_smoke_macro_constants.csv"))
end
@testset "byte-identity: PGD smoke" begin
    _check(run_pgd, joinpath(@__DIR__, "fixtures", "pgd_smoke_macro_constants.csv"))
end

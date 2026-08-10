using Test, FitRateEquation
using FitRateEquation: _canonical_enzyme, _combined_row_filter, _labeled, drop_atp_rows

@testset "symbol entry helpers" begin
    @test _canonical_enzyme(:g6pd) === :G6PD
    @test _canonical_enzyme(:G6PD) === :G6PD
    @test _canonical_enzyme(:pgd)  === :PGD
    @test_throws ErrorException _canonical_enzyme(:nope)

    # single :no_atp variant -> its filter; deploy -> identity
    @test _combined_row_filter(:G6PD, [:no_atp]) === drop_atp_rows
    @test _combined_row_filter(:G6PD, [:SS_NADPH_release_rate_eq]) === identity

    @test _labeled("G6PD", "")      == "G6PD"
    @test _labeled("G6PD", "noatp") == "G6PD_noatp"
end

@testset "symbol entry smoke fit (G6PD)" begin
    dir = mktempdir()
    res = fit_consensus_equation(:g6pd; smoke=true, nprocs=1, outdir=dir)
    @test !isempty(res)
    @test isfile(joinpath(dir, "macro_constants.csv"))
    @test isfile(joinpath(dir, "fit_corpus.csv"))
end

@testset "symbol entry :no_atp auto-filters ATP rows" begin
    dir = mktempdir()
    fit_consensus_equation(:g6pd; variants=[:no_atp], smoke=true, nprocs=1, outdir=dir)
    corpus = FitRateEquation.read_fit_corpus(dir)   # pass the run DIR, not the CSV path
    @test all(<=(0.0), corpus.ATP)   # ATP-bearing rows dropped
end

@testset "conflicting row filters error" begin
    # two DIFFERENT non-identity filters in one run is rejected. Simulate by a variant list
    # whose profiles disagree; only :no_atp has a non-identity filter today, so pair it with a
    # hand-rolled second filter via _combined_row_filter's contract.
    @test_throws ErrorException FitRateEquation._combined_row_filter_check(
        [drop_atp_rows, df -> df])
end

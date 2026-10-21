using Test, FitRateEquation
using FitRateEquation: variant_profile, drop_atp_rows

@testset "variant_profile" begin
    # :no_atp carries the ATP row filter and the "noatp" outdir label
    p = variant_profile(:G6PD, :no_atp)
    @test p.row_filter === drop_atp_rows
    @test p.label == "noatp"

    # :full_re labels the outdir but does not filter rows
    p = variant_profile(:PGD, :full_re)
    @test p.row_filter === identity
    @test p.label == "fullre"

    # PGD :no_atp labels the outdir but keeps every row: the ATP-free law is the hypothesis
    # "ATP has no effect", so ATP-bearing rows are the data that test it (unlike G6PD :no_atp).
    # The bundled PGD corpus has none; this matters for a user corpus with ATP (data_csv).
    p = variant_profile(:PGD, :no_atp)
    @test p.row_filter === identity
    @test p.label == "noatp"

    # deploy / any other variant: no filter, no label
    @test variant_profile(:G6PD, :SS_NADPH_release_rate_eq) == (row_filter=identity, label="")
    @test variant_profile(:PGD, :cha_base) == (row_filter=identity, label="")
    @test variant_profile(:G6PD, :no_g6p_atp_deadend) == (row_filter=identity, label="")
end

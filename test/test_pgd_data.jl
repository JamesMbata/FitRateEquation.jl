using FitRateEquation
using Test

@testset "PGD data" begin
    d = load_dataset(pgd_config())
    @test nrows(d) > 100
    @test all(occursin("|", g) for g in d.group)
    # ATP column present and zero-filled (structural regulator, no data yet).
    atp = [getfield(c, :ATP) for c in d.concs]
    @test all(==(0.0), atp)
end

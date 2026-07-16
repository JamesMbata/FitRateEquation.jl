using FitRateEquation
using Test

@testset "config deploy_keq" begin
    g = g6pd_config(); p = pgd_config(); h = hk1_config()
    @test g.deploy_keq == 13.655
    @test p.deploy_keq == 0.17
    @test h.deploy_keq == 2700.0
    # keq_reference is retired on the consensus_macro configs
    @test !hasproperty(g, :keq_reference)
    @test !hasproperty(p, :keq_reference)
    @test !hasproperty(h, :keq_reference)
end

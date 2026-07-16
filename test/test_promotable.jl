using FitRateEquation
using FitRateEquation: Promotable
using Test

@testset "promotable registry is bounded and correct" begin
    @test Promotable.promotable_steps(:G6PD) == [:nadph_release]
    @test Promotable.promotable_steps(:PGD) == Symbol[]
    @test Promotable.fiber_coord(:G6PD, :nadph_release) == :koffQ
    @test_throws KeyError Promotable.promotable_steps(:NOPE)
end

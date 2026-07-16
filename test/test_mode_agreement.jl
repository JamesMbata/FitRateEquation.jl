using FitRateEquation
using Test

@testset "mode agreement" begin
    c1 = [(name=:Km_G6P,  value=4.5e-5, class=:data_identified, ci=0.1),
          (name=:Ki_6PGLn, value=2.0e-4, class=:unconstrained,   ci=NaN)]
    c2 = [(name=:Km_G6P,  value=4.6e-5, class=:data_identified, ci=0.1),
          (name=:Ki_6PGLn, value=2.0e-4, class=:literature_pinned, ci=NaN)]
    ag = mode_agreement(c1, c2)
    @test all(x.name == :Km_G6P for x in ag)          # only both-modes-identified checked
    @test ag[1].agree                                  # 4.5e-5 vs 4.6e-5 within tolerance
    # a forward constant that shifts a full decade between modes fails
    bad = mode_agreement(c1, [(name=:Km_G6P, value=4.5e-4, class=:data_identified, ci=0.1)])
    @test !bad[1].agree
end

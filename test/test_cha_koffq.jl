using FitRateEquation
using FitRateEquation: g6pd_config
using EnzymeRates
using Test

using FitRateEquation.ChaFit
using FitRateEquation.ChaKoffqReport

@testset "koffq_hybrid_report: deploy convention + all-constants-refit profile" begin
    m = FitRateEquation.v2_mechanism(); d = load_dataset(g6pd_config())
    rep = koffq_hybrid_report(m, d)
    @test rep.deploy_value == CHA_DEPLOY_RELEASE_RATE
    @test rep.n_reverse > 0                            # PGLn-bearing rows exist
    @test all(p -> isfinite(p.loss) && p.loss > 0, rep.profile)
    # The corpus does not choose a koffQ: once alpha and Kd_6PGLn can move, koffQ = 1 (the old
    # fit rate) and the deploy convention fit equally well, so the scan's former
    # "data-identified koffQ = 1" was the fit's own rate, not a measurement.
    loss_at(k) = only(p.loss for p in rep.profile if p.koffQ == k)
    @test isapprox(loss_at(1.0), loss_at(CHA_DEPLOY_RELEASE_RATE); rtol=ChaKoffqReport.PROFILE_FLAT_TOL)
    @test rep.spread <= ChaKoffqReport.PROFILE_FLAT_TOL
    # alpha moves with koffQ (alpha/C is what the forward data pin).
    alpha_at(k) = only(p.alpha for p in rep.profile if p.koffQ == k)
    @test isapprox(alpha_at(1.0) / 2, alpha_at(CHA_DEPLOY_RELEASE_RATE) / (1 + 1 / CHA_DEPLOY_RELEASE_RATE); rtol=0.05)
    @test occursin("convention", rep.caveat)
    @test occursin("does not identify", rep.caveat)
end

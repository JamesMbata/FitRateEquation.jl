using FitRateEquation
using EnzymeRates
using Statistics: median
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaKoffqReport

@testset "koffq_hybrid_report: report-only hybrid (deploy swept + diagnostic + CI + gap)" begin
    m = FitRateEquation.v2_mechanism(); d = load_dataset(g6pd_config())
    rep = koffq_hybrid_report(m, d; koffQ_deploy=1.0e3)
    @test rep.deploy_value == 1.0e3
    @test rep.data_identified_value > 0
    @test rep.ci[1] < rep.ci[2]                       # a real (wide) interval
    @test isapprox(rep.gap_dex, log10(rep.data_identified_value) - log10(rep.deploy_value); atol=1e-9)
    @test rep.n_reverse > 0                            # PGLn-bearing rows exist
    @test occursin("forward gate", rep.caveat)        # honest framing present
end

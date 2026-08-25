using FitRateEquation
using FitRateEquation: g6pd_config
using FitRateEquation.ChaFit
using FitRateEquation.ChaLaws
using Statistics: median
using Test

@testset "absolute mode recovers planted kcat" begin
    # Build synthetic single-scale FORWARD data from a known kcat + the deployed variant, then
    # assert absolute mode recovers kcat within 10% — the end-to-end acceptance gate.
    d0 = load_dataset(g6pd_config())
    m  = FitRateEquation.v2_mechanism()
    # Forward-only rows (no products): PGLn = NADPH = 0. This is the regime that makes the
    # fiber-free CHA_ABS_RELEASE_RATE limit exact (the design's forward-only assumption).
    fwd = [i for i in 1:FitRateEquation.nrows(d0)
           if d0.concs[i].PGLn == 0.0 && d0.concs[i].NADPH == 0.0]
    concs = d0.concs[fwd]; group = d0.group[fwd]; keq = d0.keq[fwd]

    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    kcat_true = 178.0; Et = fill(5e-9, length(concs))
    rates = Float64[]
    for (i, cc) in enumerate(concs)
        mtup = ChaFit.cha_macro_tuple(:G6PD, coords; keq=13.655, kf=kcat_true, Et=Et[i],
                    release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
        push!(rates, ChaLaws.cha_rate_G6PD(mtup; ChaFit._cha_row_kwargs(:G6PD, cc)...))
    end
    dsyn = Dataset(concs, rates, group, keq, Et)

    fit = cha_fit_candidate(:G6PD, m, dsyn; n_restarts=8, maxiter=400, maxtime=90.0,
                            seed=1, keq=13.655, scale=:absolute,
                            pins=Dict(:Kd_6PGLn=>log10(2e-4), :Km_NADPH_rev=>log10(3.9e-6)))
    @test isapprox(fit.coords[:kcat], kcat_true; rtol=0.1)   # within 10%
end

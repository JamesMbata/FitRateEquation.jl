using FitRateEquation
using EnzymeRates
using Statistics: median
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaClassify

# Small helpers (copied from test_cha_fit.jl; file-local names to avoid a method-overwrite
# warning when both files' top-level helpers land in the same runtests.jl Main scope).
_classify_f(concs, sym) = hasproperty(concs, sym) ? getfield(concs, sym) : 0.0
# Copy a Dataset but swap in new rates (positional Dataset(concs, rate, group, keq)).
_classify_with_rates(d, r) = Dataset(d.concs, collect(float.(r)), d.group, d.keq)

@testset "classify_cha labels identified vs pinned vs flat (G6PD)" begin
    # SEEDED planted point (same scheme as the Task 8 recovery test): under the per-(Article,
    # Fig) mean-centered (Vmax-gauged) log-ratio loss the corpus robustly identifies the two
    # substrate Kd's and the ATP dead-end Ki, while Km_NADPH_rev (the bare-[NADPH] reverse
    # release term) is flat on this forward + sparse-product-inhibition coverage. Seed 21 lands
    # the planted draw on a corpus-identifiable minimum (verified in Task 8).
    m = FitRateEquation.v2_mechanism()
    d0 = load_dataset(g6pd_config()); keq = median(d0.keq)
    # Plant a FIXED, physiologically identifiable macro tuple (basis-robust: cha_deploy_micro
    # handles the free-param basis, so this no longer depends on free_params order or a tuned
    # random seed — both of which drifted when upstream EnzymeRates renamed the micro params).
    # Substrate Kd's sit in the corpus-covered range (identifiable); Km_NADPH_rev is the
    # bare-[NADPH] reverse term (flat on this forward + sparse-product coverage).
    target = Dict(:Kd_NADP=>5e-6, :Kd_G6P=>4e-5, :Kd_6PGLn=>1e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>1e-2, :Km_NADPH_rev=>5e-2)
    logθ = FitRateEquation.ChaDeploy.cha_deploy_micro(:G6PD, m, target; keq=keq, koffQ=1e3,
                            release_rate=1e3, release_eq=target[:Km_NADPH_rev])
    planted = cha_macro_readoffs_G6PD(m, logθ; keq=keq)
    syn = [ChaLaws.cha_rate_G6PD(planted;
               NADP=_classify_f(d0.concs[i],:NADP), G6P=_classify_f(d0.concs[i],:G6P),
               NADPH=_classify_f(d0.concs[i],:NADPH), PGLn=_classify_f(d0.concs[i],:PGLn),
               ATP=_classify_f(d0.concs[i],:ATP)) for i in 1:nrows(d0)]
    d = _classify_with_rates(d0, syn)

    pins = Dict(:Ki_ATP => log10(1.5e-3))
    fit = cha_fit_candidate(:G6PD, m, d; n_restarts=8, maxiter=400, maxtime=60.0,
                            seed=1, keq=keq, pins=pins)
    idf = cha_identifiable_functions(:G6PD, m, d, fit.coords; keq=keq, pins=pins)
    sigma2 = fit.loss / max(nrows(d) - idf.rank, 1)
    classed = classify_cha(:G6PD, m, d, fit.coords, pins, idf; keq=keq, sigma2=sigma2)
    cl(name) = classed[findfirst(c -> c.name === name, classed)]

    # Robustly-identified substrate constants: data_identified with finite positive CI.
    @test cl(:Kd_NADP).class === :data_identified
    @test isfinite(cl(:Kd_NADP).ci) && cl(:Kd_NADP).ci > 0
    @test cl(:Kd_G6P).class === :data_identified
    @test isfinite(cl(:Kd_G6P).ci) && cl(:Kd_G6P).ci > 0
    # Pinned coord -> literature_pinned, at the pin value.
    @test cl(:Ki_ATP).class === :literature_pinned
    @test isapprox(cl(:Ki_ATP).value, 1.5e-3; rtol=1e-8)
    # Known-flat coord under all-free fitting: NOT data_identified.
    @test cl(:Km_NADPH_rev).class !== :data_identified
end

using FitRateEquation
using EnzymeRates
using Statistics: median
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit

_p4de_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                vs[findfirst(v -> Symbol(v.name) === :full_re_deadends, vs)].mech)
_p4de_f(cc, s) = hasproperty(cc, s) ? getfield(cc, s) : 0.0
_p4de_with_rates(d, r) = Dataset(d.concs, collect(float.(r)), d.group, d.keq)

@testset "cha_coords(:PGD, :full_re_deadends) = 6 core + Ki_CO2/Ki_Ru5P (8)" begin
    cs = cha_coords(:PGD, :full_re_deadends)
    @test cs == [:Kd_NADP, :Kd_PGA, :alpha, :Kd_NADPH, :Kd_Ru5P, :Kd_CO2, :Ki_CO2, :Ki_Ru5P]
    lo, hi = cha_coord_bounds(:PGD, :full_re_deadends)
    @test length(lo) == length(hi) == 8
    for (s, l, h) in zip(cs, lo, hi)
        if s === :alpha
            @test (l, h) == (-2.0, 2.0)
        else
            @test (l, h) == (-9.0, 0.0)                # Ki_CO2/Ki_Ru5P get the Kd default band
        end
    end
end

@testset "haldane/macro-tuple/apparent-Km reuse the :full_re fiber-free path (C=1)" begin
    coords = Dict(:Kd_NADP=>1e-5, :Kd_PGA=>4e-5, :alpha=>1.4, :Kd_NADPH=>1e-6,
                  :Kd_Ru5P=>5e-5, :Kd_CO2=>1e-4, :Ki_CO2=>2e-4, :Ki_Ru5P=>8e-5)
    kr = cha_haldane_kr(:PGD, coords; keq=0.079, release_rate=1.0, kf=1.0,
                        variant=:full_re_deadends)
    @test isapprox(kr, 1.0*coords[:Kd_NADPH]*coords[:Kd_Ru5P]*coords[:Kd_CO2] /
                       (coords[:Kd_NADP]*coords[:Kd_PGA]*coords[:alpha]*0.079); rtol=1e-12)
    tup = cha_macro_tuple(:PGD, coords; keq=0.079, variant=:full_re_deadends)
    @test isfinite(tup.Ki_CO2) && tup.Ki_CO2 == coords[:Ki_CO2]
    @test isfinite(tup.Ki_Ru5P) && tup.Ki_Ru5P == coords[:Ki_Ru5P]
    @test cha_apparent_km(:PGD, coords, :Km_PGA; variant=:full_re_deadends) ≈
          coords[:alpha]*coords[:Kd_PGA]                                       # C=1
end

@testset "loss dispatch: :full_re_deadends -> cha_rate_PGD_fullRE (planted L≈0)" begin
    d0 = load_dataset(pgd_config()); keq = median(d0.keq)
    m  = _p4de_mech()
    Random.seed!(42)
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    planted = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
    syn = [cha_rate_PGD_fullRE(planted; NADP=_p4de_f(d0.concs[i],:NADP),
             PGA=_p4de_f(d0.concs[i],:PGA), Ru5P=_p4de_f(d0.concs[i],:Ru5P),
             CO2=_p4de_f(d0.concs[i],:CO2), NADPH=_p4de_f(d0.concs[i],:NADPH))
           for i in 1:nrows(d0)]
    d = _p4de_with_rates(d0, syn)
    coords = Dict(s => getfield(planted, s) for s in cha_coords(:PGD, :full_re_deadends))
    L = cha_centered_logratio_loss(:PGD, m, d, coords; keq=keq, variant=:full_re_deadends)
    @test L < 1e-8                                    # the fullRE law reproduces its own planted data
end

@testset "end-to-end cha_fit_candidate(:PGD, …; variant=:full_re_deadends) recovers Ki_CO2/Ki_Ru5P" begin
    d0 = load_dataset(pgd_config()); keq = median(d0.keq)
    m  = _p4de_mech()
    Random.seed!(7)
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    planted = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
    # Plant the two dead-end constants at physiologically sensible values the corpus can see
    # (Weisz 6A–7B carry CO2 up to 2.2 mM and Ru5P up to ~1 mM), so they are identifiable here.
    planted = merge(planted, (; Ki_CO2=3e-4, Ki_Ru5P=2e-4))
    syn = [cha_rate_PGD_fullRE(planted; NADP=_p4de_f(d0.concs[i],:NADP),
             PGA=_p4de_f(d0.concs[i],:PGA), Ru5P=_p4de_f(d0.concs[i],:Ru5P),
             CO2=_p4de_f(d0.concs[i],:CO2), NADPH=_p4de_f(d0.concs[i],:NADPH))
           for i in 1:nrows(d0)]
    d = _p4de_with_rates(d0, syn)
    fit = cha_fit_candidate(:PGD, m, d; n_restarts=12, maxiter=600, maxtime=90.0, seed=1,
                            keq=keq, variant=:full_re_deadends)
    @test fit.loss < 1e-3
    # apparent Km's recovered (fiber-free C=1); the dead-end constants recovered within ~0.5 log.
    for which in (:Km_NADP, :Km_PGA)
        got = cha_apparent_km(:PGD, fit.coords, which; variant=:full_re_deadends)
        exp = cha_apparent_km(:PGD, Dict(s=>getfield(planted,s) for s in
                    cha_coords(:PGD,:full_re_deadends)), which; variant=:full_re_deadends)
        @test abs(log10(got) - log10(exp)) < 0.5
    end
    @test abs(log10(fit.coords[:Ki_CO2])  - log10(3e-4)) < 0.6
    @test abs(log10(fit.coords[:Ki_Ru5P]) - log10(2e-4)) < 0.6
end

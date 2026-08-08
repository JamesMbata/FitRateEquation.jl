using FitRateEquation
using EnzymeRates
using Statistics: median
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy
using FitRateEquation.ChaClassify

# File-local helpers (prefixed _p3_ to avoid Main-scope method-overwrite warnings).
_p3_re_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                 vs[findfirst(v -> Symbol(v.name) === :RE_rate_eq, vs)].mech)
_p3_base_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                   vs[findfirst(v -> Symbol(v.name) === :cha_base, vs)].mech)
_p3_f(concs, sym) = hasproperty(concs, sym) ? getfield(concs, sym) : 0.0
_p3_with_rates(d, r) = Dataset(d.concs, collect(float.(r)), d.group, d.keq)

@testset "cha_deploy_micro reproduces the fitted Cha law (PGD :full_re vs V1, rtol 1e-9)" begin
    m = _p3_re_mech()                       # V1 = fully-RE topology (:RE_rate_eq), ATP dead-ends ON
    mets = EnzymeRates.metabolites(m)
    Random.seed!(3)
    for _ in 1:8
        logθ = -3 .+ 2 .* rand(length(free_params(m))); keq = 0.079
        mac = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
        # V1 carries BOTH ATP dead-ends as free-params, so the deploy coords must include the
        # readoff-recovered effector constants alongside the 6 core Kd's (no NADPH dead-end in V1).
        coord_syms = vcat(cha_coords(:PGD, :full_re), [:Ki_ATP, :Ki_ATP_EN])
        coords = Dict(s => getfield(mac, s) for s in coord_syms)
        logθd = cha_deploy_micro(:PGD, m, coords; keq=keq, variant=:full_re)
        bp = build_params(m, logθd; keq=keq)
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)), bp))
        for conc in [(;NADP=5e-6,PGA=40e-6,Ru5P=1e-4,CO2=1e-4,NADPH=5e-6,ATP=1e-3),
                     (;NADP=2e-5,PGA=80e-6,Ru5P=2e-4,CO2=2e-4,NADPH=8e-6,ATP=5e-4)]
            cc = NamedTuple{Tuple(mets)}(Tuple(get(conc, s, 0.0) for s in mets))
            vref = cha_rate_PGD_fullRE(mac; conc...)
            vdep = EnzymeRates.rate_equation(m, cc, bp)
            @test isapprox(vdep, vref; rtol=1e-9, atol=1e-9*vsat)
        end
    end
end

@testset ":full_re deploy map is fiber-free and covers V1's free_params" begin
    m = _p3_re_mech()
    coords = Dict(:Kd_NADP=>1e-5, :Kd_PGA=>4e-5, :alpha=>1.4, :Kd_NADPH=>1e-6,
                  :Kd_Ru5P=>5e-5, :Kd_CO2=>1e-4, :Ki_ATP=>1e-3, :Ki_ATP_EN=>2e-3)
    micro = FitRateEquation.ChaDeploy._deploy_micro_map(:PGD, coords;
                release_rate=1.0, release_eq=1.0, mech=m, variant=:full_re)
    @test !haskey(micro, :koff_Ru5P_ENADPH) && !haskey(micro, :kon_Ru5P_ENADPH)   # fiber-free
    @test micro[:K_NADP_EPGA] ≈ coords[:alpha]*coords[:Kd_NADP]                    # detailed balance
    @test micro[:K_NADPH_E] == coords[:Kd_NADPH]
    @test micro[:K_Ru5P_ENADPH] == coords[:Kd_Ru5P]
    @test Set(keys(micro)) == Set(free_params(m))                                  # covers V1 exactly
    @test length(cha_deploy_micro(:PGD, m, coords; keq=0.079, variant=:full_re)) ==
          length(free_params(m))
    # ADDITIVE GUARD: cha_base PGD deploy map still targets the SS-release fiber, unchanged.
    vb = _p3_base_mech()
    base = Dict(:Kd_NADP=>1e-5,:Kd_PGA=>4e-5,:alpha=>1.4,:Kd_CO2=>1e-4,:Ki_NADPH=>1e-5,
                :Ki_ATP=>1e-3,:Ki_ATP_EN=>2e-3,:Km_NADPH_rev=>5e-2)
    mb = FitRateEquation.ChaDeploy._deploy_micro_map(:PGD, base;
                release_rate=1e3, release_eq=1.0, mech=vb)
    @test haskey(mb, :koff_Ru5P_ENADPH) && haskey(mb, :kon_Ru5P_ENADPH)
end

_p3_syn(d, planted) = [cha_rate_PGD_fullRE(planted;
        NADP=_p3_f(d.concs[i],:NADP), PGA=_p3_f(d.concs[i],:PGA), Ru5P=_p3_f(d.concs[i],:Ru5P),
        CO2=_p3_f(d.concs[i],:CO2), NADPH=_p3_f(d.concs[i],:NADPH), ATP=_p3_f(d.concs[i],:ATP))
    for i in 1:nrows(d)]

@testset "classify_cha(:PGD, :full_re) yields a well-formed identifiability spectrum" begin
    m = _p3_re_mech()
    d0 = load_dataset(pgd_config()); keq = median(d0.keq)
    Random.seed!(23)
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    planted = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
    d = _p3_with_rates(d0, _p3_syn(d0, planted))
    fit = cha_fit_candidate(:PGD, m, d; n_restarts=8, maxiter=400, maxtime=60.0, seed=1,
                            keq=keq, variant=:full_re)
    idf = cha_identifiable_functions(:PGD, m, d, fit.coords; keq=keq, variant=:full_re)
    cs = cha_coords(:PGD, :full_re)
    @test length(idf.eigvals) == length(cs) == 6
    @test all(isfinite, idf.eigvals)
    @test idf.idx == collect(1:6)             # mode-1: nothing pinned
    @test 0 <= idf.rank <= 6
    sigma2 = fit.loss / max(nrows(d) - idf.rank, 1)
    classed = classify_cha(:PGD, m, d, fit.coords, Dict{Symbol,Float64}(), idf;
                           keq=keq, sigma2=sigma2, variant=:full_re)
    @test Set(getfield.(classed, :name)) == Set(cs)
    for c in classed
        @test c.class in (:data_identified, :unconstrained)
        @test (c.class === :data_identified) ? (isfinite(c.ci) && c.ci > 0) : isnan(c.ci)
    end
    classed2 = classify_cha(:PGD, m, d, fit.coords, Dict{Symbol,Float64}(), idf;
                            keq=keq, sigma2=sigma2, variant=:full_re)
    @test all(classed[i].class === classed2[i].class for i in eachindex(classed))
end

@testset "classify_cha(:PGD, :full_re) honors a coord pin (literature_pinned path)" begin
    m = _p3_re_mech()
    d0 = load_dataset(pgd_config()); keq = median(d0.keq)
    Random.seed!(23)
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    planted = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
    d = _p3_with_rates(d0, _p3_syn(d0, planted))
    pins = Dict(:Kd_CO2 => log10(1e-4))
    fit = cha_fit_candidate(:PGD, m, d; n_restarts=8, maxiter=400, maxtime=60.0, seed=1,
                            keq=keq, pins=pins, variant=:full_re)
    idf = cha_identifiable_functions(:PGD, m, d, fit.coords; keq=keq, pins=pins, variant=:full_re)
    @test idf.idx == [j for (j, s) in enumerate(cha_coords(:PGD, :full_re)) if s !== :Kd_CO2]
    sigma2 = fit.loss / max(nrows(d) - idf.rank, 1)
    classed = classify_cha(:PGD, m, d, fit.coords, pins, idf; keq=keq, sigma2=sigma2, variant=:full_re)
    cco2 = classed[findfirst(c -> c.name === :Kd_CO2, classed)]
    @test cco2.class === :literature_pinned
    @test isapprox(cco2.value, 1e-4; rtol=1e-8)
end

_p3_fullre_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                     vs[findfirst(v -> Symbol(v.name) === :full_re, vs)].mech)

@testset "cha_deploy_micro round-trips the registered :full_re core mech (rtol 1e-9)" begin
    m = _p3_fullre_mech()                    # effectors-off core: 6 RE free_params
    mets = EnzymeRates.metabolites(m)
    Random.seed!(5)
    for _ in 1:8
        logθ = -3 .+ 2 .* rand(length(free_params(m))); keq = 0.079
        mac = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
        coords = Dict(s => getfield(mac, s) for s in cha_coords(:PGD, :full_re))   # 6 core only
        logθd = cha_deploy_micro(:PGD, m, coords; keq=keq, variant=:full_re)
        bp = build_params(m, logθd; keq=keq)
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)), bp))
        for conc in [(;NADP=5e-6,PGA=40e-6,Ru5P=1e-4,CO2=1e-4,NADPH=5e-6),
                     (;NADP=2e-5,PGA=80e-6,Ru5P=2e-4,CO2=2e-4,NADPH=8e-6)]
            cc = NamedTuple{Tuple(mets)}(Tuple(get(conc, s, 0.0) for s in mets))
            @test isapprox(EnzymeRates.rate_equation(m, cc, bp),
                           cha_rate_PGD_fullRE(mac; conc...); rtol=1e-9, atol=1e-9*vsat)
        end
    end
end

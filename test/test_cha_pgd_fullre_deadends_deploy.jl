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

_p4de_mech() = (vs = FitRateEquation.consensus_variants(:PGD);
                vs[findfirst(v -> Symbol(v.name) === :full_re_deadends, vs)].mech)
_p4de_f(cc, s) = hasproperty(cc, s) ? getfield(cc, s) : 0.0
_p4de_with_rates(d, r) = Dataset(d.concs, collect(float.(r)), d.group, d.keq)

@testset "cha_deploy_micro round-trips :full_re_deadends (rtol 1e-9)" begin
    m    = _p4de_mech()
    mets = EnzymeRates.metabolites(m)
    Random.seed!(9)
    for _ in 1:8
        logθ = -3 .+ 2 .* rand(length(free_params(m))); keq = 0.079
        mac  = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
        coords = Dict(s => getfield(mac, s) for s in cha_coords(:PGD, :full_re_deadends))  # 8 coords
        logθd  = cha_deploy_micro(:PGD, m, coords; keq=keq, variant=:full_re_deadends)
        bp = build_params(m, logθd; keq=keq)
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)), bp))
        for conc in [(; NADP=5e-6, PGA=40e-6, CO2=2.2e-3, Ru5P=0.0,  NADPH=0.0),   # CO2 dead-end
                     (; NADP=5e-6, PGA=40e-6, CO2=0.0,    Ru5P=4e-4, NADPH=0.0),   # Ru5P dead-end
                     (; NADP=2e-5, PGA=80e-6, CO2=2e-4,   Ru5P=1e-4, NADPH=8e-6)]  # all present
            cc = NamedTuple{Tuple(mets)}(Tuple(get(conc, s, 0.0) for s in mets))
            @test isapprox(EnzymeRates.rate_equation(m, cc, bp),
                           cha_rate_PGD_fullRE(mac; conc...); rtol=1e-9, atol=1e-9*vsat)
        end
    end
end

@testset ":full_re_deadends deploy map fills K_CO2_E/K_Ru5P_E and stays fiber-free" begin
    m = _p4de_mech()
    coords = Dict(:Kd_NADP=>1e-5, :Kd_PGA=>4e-5, :alpha=>1.4, :Kd_NADPH=>1e-6, :Kd_Ru5P=>5e-5,
                  :Kd_CO2=>1e-4, :Ki_CO2=>2e-4, :Ki_Ru5P=>8e-5)
    micro = FitRateEquation.ChaDeploy._deploy_micro_map(:PGD, coords;
                release_rate=1.0, release_eq=1.0, mech=m, variant=:full_re_deadends)
    @test micro[:K_CO2_E]  == coords[:Ki_CO2]
    @test micro[:K_Ru5P_E] == coords[:Ki_Ru5P]
    @test !haskey(micro, :koff_Ru5P_ENADPH) && !haskey(micro, :kon_Ru5P_ENADPH)   # fiber-free
    @test Set(keys(micro)) == Set(free_params(m))                                  # covers the mech
end

@testset "classify_cha(:PGD, :full_re_deadends) yields a well-formed 8-coord spectrum" begin
    m  = _p4de_mech()
    d0 = load_dataset(pgd_config()); keq = median(d0.keq)
    Random.seed!(23)
    logθ = -3 .+ 2 .* rand(length(free_params(m)))
    planted = merge(cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq), (; Ki_CO2=3e-4, Ki_Ru5P=2e-4))
    syn = [cha_rate_PGD_fullRE(planted; NADP=_p4de_f(d0.concs[i],:NADP), PGA=_p4de_f(d0.concs[i],:PGA),
             Ru5P=_p4de_f(d0.concs[i],:Ru5P), CO2=_p4de_f(d0.concs[i],:CO2),
             NADPH=_p4de_f(d0.concs[i],:NADPH)) for i in 1:nrows(d0)]
    d = _p4de_with_rates(d0, syn)
    fit = cha_fit_candidate(:PGD, m, d; n_restarts=8, maxiter=400, maxtime=60.0, seed=1,
                            keq=keq, variant=:full_re_deadends)
    idf = cha_identifiable_functions(:PGD, m, d, fit.coords; keq=keq, variant=:full_re_deadends)
    cs  = cha_coords(:PGD, :full_re_deadends)
    @test length(idf.eigvals) == length(cs) == 8
    @test all(isfinite, idf.eigvals)
    @test idf.idx == collect(1:8)                     # mode-1: nothing pinned
    @test 0 <= idf.rank <= 8
    sigma2 = fit.loss / max(nrows(d) - idf.rank, 1)
    classed = classify_cha(:PGD, m, d, fit.coords, Dict{Symbol,Float64}(), idf;
                           keq=keq, sigma2=sigma2, variant=:full_re_deadends)
    @test Set(getfield.(classed, :name)) == Set(cs)
    for c in classed
        @test c.class in (:data_identified, :unconstrained)
    end
end

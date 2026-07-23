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

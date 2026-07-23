using FitRateEquation
using EnzymeRates
using Statistics: median
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit

# The registered fully-RE PGD mechanism (:RE_rate_eq = V1) — for reading self-consistent
# planted points off EnzymeRates in Task 3.
function _pgd_re_mech()
    vs = FitRateEquation.consensus_variants(:PGD)
    vs[findfirst(v -> Symbol(v.name) === :RE_rate_eq, vs)].mech
end

# A representative fully-RE coord dict (µM-scale substrate Kd's, competitive NADPH Kq).
_fr_coords() = Dict(:Kd_NADP=>1e-5, :Kd_PGA=>4e-5, :alpha=>1.4,
                    :Kd_NADPH=>1e-6, :Kd_Ru5P=>5e-5, :Kd_CO2=>1e-4)

@testset "cha_coords(:PGD, :full_re) is the 6 core fiber-free coords" begin
    cs = cha_coords(:PGD, :full_re)
    @test cs == [:Kd_NADP, :Kd_PGA, :alpha, :Kd_NADPH, :Kd_Ru5P, :Kd_CO2]
    for bad in (:kf, :kr, :Et, :koff, :kon, :Km_NADPH_rev, :Ki_NADPH, :Ki_ATP)
        @test !(bad in cs)
    end
    # cha_base PGD coords unchanged (additive guard).
    @test cha_coords(:PGD) == [:Kd_NADP, :Kd_PGA, :alpha, :Kd_CO2, :Ki_NADPH, :Ki_ATP,
                               :Ki_ATP_EN, :Km_NADPH_rev]
end

@testset "cha_coord_bounds(:PGD, :full_re): alpha bounded, Kd's in 1nM..1M" begin
    cs = cha_coords(:PGD, :full_re)
    lo, hi = cha_coord_bounds(:PGD, :full_re)
    @test length(lo) == length(hi) == length(cs)
    for (i, s) in enumerate(cs)
        if s === :alpha
            @test (lo[i], hi[i]) == (-2.0, 2.0)
        else
            @test (lo[i], hi[i]) == (-9.0, 0.0)
        end
    end
end

@testset "cha_haldane_kr(:PGD, :full_re) is the fully-RE Haldane (no release fiber)" begin
    coords = _fr_coords(); keq, kf = 0.079, 1.0
    kr = cha_haldane_kr(:PGD, coords; keq=keq, release_rate=1e3, kf=kf, variant=:full_re)
    expected = coords[:Kd_NADPH]*coords[:Kd_Ru5P]*coords[:Kd_CO2]*kf /
               (coords[:Kd_NADP]*coords[:Kd_PGA]*coords[:alpha]*keq)
    @test isapprox(kr, expected; rtol=1e-12)
    # release_rate is inert in the fully-RE Haldane (fiber-free).
    @test cha_haldane_kr(:PGD, coords; keq=keq, release_rate=1.0, kf=kf, variant=:full_re) == kr
end

@testset "cha_macro_tuple(:PGD, :full_re): fiber-free tuple, Haldane kr, no koff/kon" begin
    coords = _fr_coords(); keq = 0.079
    m = cha_macro_tuple(:PGD, coords; keq=keq, variant=:full_re)
    # Exact fields the fully-RE law consumes; NO koff/kon/Km_NADPH_rev fiber symbols.
    @test Set(keys(m)) == Set((:Kd_NADP,:Kd_PGA,:alpha,:Kd_NADPH,:Kd_Ru5P,:Kd_CO2,
                               :kf,:kr,:Et,:Keq))
    @test m.kf == 1.0 && m.Et == 1.0 && m.Keq == keq
    krH = cha_haldane_kr(:PGD, coords; keq=keq, release_rate=1e3, kf=1.0, variant=:full_re)
    @test isapprox(m.kr, krH; rtol=1e-12)
    # The tuple drives cha_rate_PGD_fullRE and is Haldane-balanced (v≈0 at equilibrium).
    A,B = 1e-5,4e-5; Q,R = 1e-6,1e-3; Cc = keq*A*B/(Q*R)
    v  = cha_rate_PGD_fullRE(m; NADP=A, PGA=B, NADPH=Q, Ru5P=R, CO2=Cc)
    vf = cha_rate_PGD_fullRE(m; NADP=A, PGA=B)
    @test abs(v) < 1e-10*abs(vf)
    # Config-gated dead-ends appended ONLY when present as coords (default absent).
    @test !haskey(m, :Ki_ATP)
    m2 = cha_macro_tuple(:PGD, merge(coords, Dict(:Ki_ATP=>1e-3)); keq=keq, variant=:full_re)
    @test haskey(m2, :Ki_ATP) && m2.Ki_ATP == 1e-3
end

@testset "cha_apparent_km(:PGD, :full_re): C=1 so Km == alpha*Kd" begin
    coords = _fr_coords()
    @test cha_apparent_km(:PGD, coords, :Km_NADP; variant=:full_re) ≈ coords[:alpha]*coords[:Kd_NADP]
    @test cha_apparent_km(:PGD, coords, :Km_PGA;  variant=:full_re) ≈ coords[:alpha]*coords[:Kd_PGA]
    # Fiber-free: independent of release_rate (C=1 exactly, unlike the cha_base 1+kf/r).
    @test cha_apparent_km(:PGD, coords, :Km_NADP; variant=:full_re, release_rate=1.0) ==
          cha_apparent_km(:PGD, coords, :Km_NADP; variant=:full_re, release_rate=1e9)
    # cha_base default path unchanged: C = 1 + kf/CHA_DEPLOY_RELEASE_RATE ≠ 1.
    base = Dict(:Kd_NADP=>1e-5, :Kd_PGA=>4e-5, :alpha=>1.4, :Kd_CO2=>1e-4)
    @test cha_apparent_km(:PGD, base, :Km_NADP) ≈ 1.4*1e-5 / (1 + 1/CHA_DEPLOY_RELEASE_RATE)
end

@testset "cha_specificity(:PGD): kcat/Km = kf/(alpha*Kd), fiber-invariant" begin
    coords = _fr_coords()
    @test cha_specificity(:PGD, coords, :Km_NADP) ≈ 1.0/(coords[:alpha]*coords[:Kd_NADP])
    @test cha_specificity(:PGD, coords, :Km_PGA)  ≈ 1.0/(coords[:alpha]*coords[:Kd_PGA])
end

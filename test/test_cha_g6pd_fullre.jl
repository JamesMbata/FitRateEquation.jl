using FitRateEquation
using FitRateEquation: variant_profile, drop_atp_rows, modes_for
using EnzymeRates
using DataFrames, CSV
using Random
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy

# G6PD :full_re -- random RE substrate binding, catalysis the only SS step, RE product release
# (PGLn off first, NADPH last), no dead ends, no ATP term. File-local helpers prefixed _gfr_.

_gfr_mech() = (vs = FitRateEquation.consensus_variants(:G6PD);
               vs[findfirst(v -> Symbol(v.name) === :full_re, vs)].mech)

# A representative coord dict (µM-scale Kd's, alpha off 1).
_gfr_coords() = Dict(:Kd_NADP => 1e-5, :Kd_G6P => 4e-5, :alpha => 1.5,
                     :Kd_NADPH => 2e-6, :Kd_6PGLn => 1e-3)

const _GFR_KEQ = 13.655

@testset "cha_coords(:G6PD, :full_re) is the 5 fiber-free coords" begin
    @test cha_coords(:G6PD, :full_re) == [:Kd_NADP, :Kd_G6P, :alpha, :Kd_NADPH, :Kd_6PGLn]
    # Existing G6PD variants unchanged (additive guard).
    @test cha_coords(:G6PD) == [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Ki_ATP,
                                :Ki_ATP_EG, :Km_NADPH_rev]
    @test cha_coords(:G6PD, :no_atp) ==
          [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Km_NADPH_rev]
end

@testset ":full_re is registered under consensus_variants(:G6PD): no dead ends, no ATP" begin
    m = _gfr_mech()
    @test Set(EnzymeRates.metabolites(m)) == Set([:NADP, :G6P, :NADPH, :PGLn])
    fp = Set(free_params(m))
    @test fp == Set([:K_NADP_E, :K_G6P_E, :K_G6P_ENADP, :K_NADPH_E, :K_PGLn_ENADPH])
end

@testset "closed-form law == er_rate_equation of the registered mechanism (200 pts, rtol 1e-10)" begin
    m    = _gfr_mech()
    mets = EnzymeRates.metabolites(m)
    rng  = MersenneTwister(20261112)
    lu() = 10.0^(-7 + 5 * rand(rng))                 # log-uniform 1e-7 .. 1e-2 M
    for _ in 1:10                                     # 10 random coord sets x 20 points = 200
        coords = Dict(:Kd_NADP => 10.0^(-7 + 5rand(rng)), :Kd_G6P => 10.0^(-7 + 5rand(rng)),
                      :alpha => 10.0^(-2 + 4rand(rng)), :Kd_NADPH => 10.0^(-7 + 5rand(rng)),
                      :Kd_6PGLn => 10.0^(-7 + 5rand(rng)))
        keq = 10.0^(-1 + 3rand(rng))
        logθ = cha_deploy_micro(:G6PD, m, coords; keq=keq, variant=:full_re)
        bp   = build_params(m, logθ; keq=keq)
        mac  = cha_macro_tuple(:G6PD, coords; keq=keq, variant=:full_re)
        # Readoff round-trip: the micro point reads back to the same macro tuple.
        rd = cha_macro_readoffs_G6PD_fullRE(m, logθ; keq=keq)
        for s in cha_coords(:G6PD, :full_re)
            @test isapprox(getfield(rd, s), coords[s]; rtol=1e-12)
        end
        @test isapprox(rd.kr, mac.kr; rtol=1e-12)
        for _ in 1:20
            conc = (; NADP=lu(), G6P=lu(), NADPH=lu(), PGLn=lu())
            cc   = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = er_rate_equation(m, cc, bp)
            vcha = cha_rate_G6PD_fullRE(mac; conc...)
            # atol guards only the cancellation floor near equilibrium: 1e-12 of the larger
            # one-way flux at this point.
            fwd  = abs(cha_rate_G6PD_fullRE(mac; NADP=conc.NADP, G6P=conc.G6P,
                                            NADPH=conc.NADPH, PGLn=0.0))
            @test isapprox(vcha, vref; rtol=1e-10, atol=1e-12 * fwd)
        end
    end
end

@testset "Haldane: v == 0 when NADPH*PGLn/(NADP*G6P) == Keq" begin
    mac = cha_macro_tuple(:G6PD, _gfr_coords(); keq=_GFR_KEQ, variant=:full_re)
    A, B, Q = 2e-5, 5e-5, 3e-6
    P = _GFR_KEQ * A * B / Q
    v = cha_rate_G6PD_fullRE(mac; NADP=A, G6P=B, NADPH=Q, PGLn=P)
    # Forward term kf*Et*gAB/D at the same point (the characteristic scale).
    c = _gfr_coords()
    gAB = A * B / (c[:alpha] * c[:Kd_NADP] * c[:Kd_G6P])
    gQP = Q * P / (c[:Kd_NADPH] * c[:Kd_6PGLn])
    D = 1 + A / c[:Kd_NADP] + B / c[:Kd_G6P] + gAB + Q / c[:Kd_NADPH] + gQP
    @test abs(v) < 1e-12 * (mac.kf * mac.Et * gAB / D)
end

@testset "finite rate with NADPH = 0, PGLn = 0, or both" begin
    mac = cha_macro_tuple(:G6PD, _gfr_coords(); keq=_GFR_KEQ, variant=:full_re)
    for (Q, P) in ((0.0, 1e-4), (1e-5, 0.0), (0.0, 0.0))
        v = cha_rate_G6PD_fullRE(mac; NADP=2e-5, G6P=5e-5, NADPH=Q, PGLn=P)
        @test isfinite(v) && v > 0
    end
    # All four metabolites at zero: no flux, still finite.
    @test cha_rate_G6PD_fullRE(mac) == 0.0
    # Pure reverse (no substrates): finite and negative.
    vr = cha_rate_G6PD_fullRE(mac; NADPH=1e-5, PGLn=1e-4)
    @test isfinite(vr) && vr < 0
end

@testset "saturating forward rate (NADP = G6P = 1 M, no products) -> kf*Et" begin
    mac = cha_macro_tuple(:G6PD, _gfr_coords(); keq=_GFR_KEQ, variant=:full_re)
    v = cha_rate_G6PD_fullRE(mac; NADP=1.0, G6P=1.0)
    @test abs(v - mac.kf * mac.Et) <= 1e-4 * mac.kf * mac.Et
end

@testset "fiber-free: apparent Km = alpha*Kd (C = 1), release rate inert" begin
    c = _gfr_coords()
    @test cha_apparent_km(:G6PD, c, :Km_NADP; variant=:full_re) == c[:alpha] * c[:Kd_NADP]
    @test cha_apparent_km(:G6PD, c, :Km_G6P;  variant=:full_re) == c[:alpha] * c[:Kd_G6P]
    @test cha_apparent_km(:G6PD, c, :Km_G6P; variant=:full_re, release_rate=1.0) ==
          cha_apparent_km(:G6PD, c, :Km_G6P; variant=:full_re, release_rate=1e9)
    mac = cha_macro_tuple(:G6PD, c; keq=_GFR_KEQ, variant=:full_re)
    @test !hasproperty(mac, :koffQ) && !hasproperty(mac, :konQ)
    @test mac.kr ≈ mac.kf * c[:Kd_NADPH] * c[:Kd_6PGLn] /
                   (c[:alpha] * c[:Kd_NADP] * c[:Kd_G6P] * _GFR_KEQ)
end

@testset "deploy micro names (fiber-free, K_G6P_ENADP = alpha*Kd_G6P)" begin
    c = _gfr_coords()
    micro = FitRateEquation.ChaDeploy._deploy_micro_map(:G6PD, c; release_rate=1e3,
                release_eq=1.0, mech=_gfr_mech(), variant=:full_re)
    @test micro == Dict(:K_NADP_E => c[:Kd_NADP], :K_G6P_E => c[:Kd_G6P],
                        :K_G6P_ENADP => c[:alpha] * c[:Kd_G6P],
                        :K_NADPH_E => c[:Kd_NADPH], :K_PGLn_ENADPH => c[:Kd_6PGLn])
end

@testset "variant_profile / modes_for" begin
    p = variant_profile(:G6PD, :full_re)
    @test p.label == "fullre"
    @test p.row_filter === drop_atp_rows
    df = DataFrame(ATP=[0.0, 1e-3, 0.0], Rate=[1.0, 2.0, 3.0])
    @test p.row_filter(df).Rate == [1.0, 3.0]
    @test modes_for(:G6PD, :full_re) == (:mode1,)
end

@testset "smoke fit writes the 5 coords + derived Km_NADP, Km_G6P" begin
    out = mktempdir()
    fit_consensus_equation(:g6pd; variants=[:full_re], smoke=true, nprocs=1, outdir=out)
    mc = CSV.read(joinpath(out, "macro_constants.csv"), DataFrame)
    names = Set(String.(mc.name))
    for s in ("Kd_NADP", "Kd_G6P", "alpha", "Kd_NADPH", "Kd_6PGLn", "Km_NADP", "Km_G6P")
        @test s in names
    end
    @test all(==("full_re"), String.(mc.variant))
    @test all(==("mode1"), String.(mc.mode))
    fc = CSV.read(joinpath(out, "fit_corpus.csv"), DataFrame)
    @test !hasproperty(fc, :ATP) || all(<=(0.0), fc.ATP)
    # The deploy block names the five fiber-free micro constants, no koff/kon.
    mp = read(joinpath(out, "micro_parameters.jl"), String)
    for s in ("K_NADP_E", "K_G6P_E", "K_G6P_ENADP", "K_NADPH_E", "K_PGLn_ENADPH")
        @test occursin("full_re_mode1_$(s) = ", mp)
    end
    @test !occursin("koff", mp) && !occursin("deploy unavailable", mp)
    # Fiber-free report: no koffQ block, the :full_re NADPH note instead of the Ki_NADPH one.
    rep = read(joinpath(out, "report.md"), String)
    @test !occursin("koffQ hybrid", rep)
    @test occursin("NADPH product inhibition (G6PD `:full_re`)", rep)
    @test !occursin("NOT DEPLOYABLE (", rep)
end

# =========================================================================================
#                    ATP-free G6PD Cha variant (:no_atp) — unit tests
# =========================================================================================
# Proves the ATP-free law is the exact Ki_ATP,Ki_ATP_EG -> Inf limit of the deployed Cha law,
# and that the reduced-coordinate fit / deploy / pin wiring is consistent.
using FitRateEquation
using EnzymeRates
using Random
using Test
using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy

@testset "no_atp mechanism registration" begin
    vs = FitRateEquation.consensus_variants(:G6PD)
    names = [Symbol(v.name) for v in vs]
    @test :no_atp in names
    # Deploy variant indices are unchanged (v1/v2 stay at 1/2).
    @test Symbol(vs[1].name) === :RE_rate_eq
    @test Symbol(vs[2].name) === :SS_NADPH_release_rate_eq
    # The ATP-free mechanism has no ATP metabolite at all.
    mna = vs[findfirst(v -> Symbol(v.name) === :no_atp, vs)].mech
    @test !(:ATP in EnzymeRates.metabolites(mna))
    # The NADPH product-inhibition dead-end is retained (NADPH still a metabolite).
    @test :NADPH in EnzymeRates.metabolites(mna)
end

@testset "no_atp exactness gate: cha_rate ≡ rate_equation(no_atp mechanism)" begin
    Random.seed!(20260712)
    vs  = FitRateEquation.consensus_variants(:G6PD)
    m   = vs[findfirst(v -> Symbol(v.name) === :no_atp, vs)].mech
    mets = EnzymeRates.metabolites(m)          # no :ATP
    grid = [
        (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=0.0),
        (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=1e-4),
        (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=0.0),
        (; NADP=1e-6, G6P=10e-6, NADPH=2e-6, PGLn=5e-5),
        (; NADP=2e-5, G6P=80e-6, NADPH=8e-6, PGLn=2e-4),
    ]
    maxrel = 0.0
    for _ in 1:25
        free = free_params(m)
        logθ = -3 .+ 2 .* rand(length(free))
        keq  = 10.0
        mac  = ChaInvert.cha_macro_readoffs_G6PD(m, logθ; keq=keq)
        @test isinf(mac.Ki_ATP) && isinf(mac.Ki_ATP_EG)   # ATP terms structurally absent
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:G6P) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
            vcha = ChaLaws.cha_rate_G6PD(mac; conc...)
            scale = max(abs(vref), abs(vcha), 1e-6 * vsat)
            maxrel = max(maxrel, abs(vcha - vref) / scale)
        end
    end
    @test maxrel <= 1e-12
end

@testset "no_atp reduced coords + macro tuple" begin
    c = ChaFit.cha_coords(:G6PD, :no_atp)
    @test :Ki_ATP ∉ c && :Ki_ATP_EG ∉ c
    @test Set(c) == Set([:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Km_NADPH_rev])
    # A reduced coords dict (no ATP keys) assembles a tuple with Ki_ATP/Ki_ATP_EG = Inf,
    # and the resulting law equals the full law evaluated at ATP = 0.
    coords = Dict(:Kd_NADP=>5e-6, :Kd_G6P=>25e-6, :Kd_6PGLn=>1e-3, :alpha=>4.0,
                  :Ki_NADPH=>20e-6, :Km_NADPH_rev=>3.9e-6)
    m = ChaFit.cha_macro_tuple(:G6PD, coords; keq=13.7, variant=:no_atp)
    @test isinf(m.Ki_ATP) && isinf(m.Ki_ATP_EG)
    r_noatp = ChaLaws.cha_rate_G6PD(m; NADP=5e-6, G6P=40e-6, NADPH=3e-6, ATP=2e-3)
    r_atp0  = ChaLaws.cha_rate_G6PD(m; NADP=5e-6, G6P=40e-6, NADPH=3e-6, ATP=0.0)
    @test r_noatp == r_atp0   # ATP-blind: rate independent of ATP
end

@testset "no_atp pins + modes" begin
    @test FitRateEquation.modes_for(:G6PD, :no_atp) == (:mode1, :mode2)
    p1 = ChaFit.resolve_cha_pins(:G6PD, :no_atp, :mode1)
    @test Set(keys(p1)) == Set([:Km_NADPH_rev])          # all-modes reverse anchor only
    p2 = ChaFit.resolve_cha_pins(:G6PD, :no_atp, :mode2)
    @test :Ki_NADPH in keys(p2)                          # lit cross-check
    @test :Km_NADPH_rev in keys(p2)
    @test :Ki_ATP ∉ keys(p2)                             # not a coord -> auto-skipped
end

@testset "no_atp deploy inverse round-trips" begin
    Random.seed!(20260713)
    vs = FitRateEquation.consensus_variants(:G6PD)
    m  = vs[findfirst(v -> Symbol(v.name) === :no_atp, vs)].mech
    mets = EnzymeRates.metabolites(m)
    coords = Dict(:Kd_NADP=>5e-6, :Kd_G6P=>25e-6, :Kd_6PGLn=>1e-3, :alpha=>4.0,
                  :Ki_NADPH=>20e-6, :Km_NADPH_rev=>3.9e-6)
    keq = 13.7
    rr  = ChaFit.CHA_DEPLOY_RELEASE_RATE
    logθ = ChaDeploy.cha_deploy_micro(:G6PD, m, coords; keq=keq, koffQ=rr, release_rate=rr)
    @test length(logθ) == length(free_params(m))
    mac = ChaFit.cha_macro_tuple(:G6PD, coords; keq=keq, variant=:no_atp,
                                 release_rate=rr, release_eq=coords[:Km_NADPH_rev])
    p = build_params(m, logθ; keq=keq)
    for conc in [(; NADP=5e-6, G6P=40e-6, NADPH=0.0, PGLn=0.0),
                 (; NADP=1e-5, G6P=80e-6, NADPH=4e-6, PGLn=1e-4)]
        cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
        vref = EnzymeRates.rate_equation(m, cc, p)
        vcha = ChaLaws.cha_rate_G6PD(mac; conc...)
        @test isapprox(vref, vcha; rtol=1e-9, atol=1e-30)
    end
end

@testset "runner: no_atp cell list" begin
    cells = FitRateEquation._cells(:G6PD; variants=[:no_atp])
    @test length(cells) == 2
    @test all(c -> c[1] === :no_atp, cells)
    @test [c[3] for c in cells] == [:mode1, :mode2]
    # Default behavior unchanged: the deploy variant still drives the normal cell list.
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:G6PD),
              FitRateEquation._cells(:G6PD))
end

@testset "drop_atp_rows filter" begin
    concs = [(; NADP=5e-6, G6P=4e-5, NADPH=0.0, PGLn=0.0, ATP=0.0),
             (; NADP=5e-6, G6P=4e-5, NADPH=0.0, PGLn=0.0, ATP=2e-3),
             (; NADP=1e-5, G6P=8e-5, NADPH=0.0, PGLn=0.0, ATP=0.0)]
    d = Dataset(collect(concs), [1.0, 0.8, 1.2],
                            ["A|1", "A|1", "A|1"], [13.7, 13.7, 13.7])
    d2 = FitRateEquation.drop_atp_rows(d)
    @test nrows(d2) == 2
    @test all(nt -> get(nt, :ATP, 0.0) <= 0.0, d2.concs)
    @test d2.rate == [1.0, 1.2]
end

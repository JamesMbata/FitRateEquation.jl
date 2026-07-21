# =========================================================================================
#     Dead-end-dropped G6PD Cha variants: no_g6p_nadph_deadend / no_g6p_atp_deadend /
#                              no_g6p_both_deadends — unit tests
# =========================================================================================
# Each variant is an exact Ki -> Inf limit of the deployed Cha law (SS_NADPH_release_rate_eq):
#   no_g6p_nadph_deadend  drops E·G6P·NADPH  (Ki_NADPH -> Inf),   keeps both ATP dead-ends
#   no_g6p_atp_deadend    drops E·G6P·ATP    (Ki_ATP_EG -> Inf),  keeps E·ATP + NADPH dead-end
#   no_g6p_both_deadends  drops both of the above,                keeps E·ATP only
# Mirrors test_cha_noatp.jl's testset structure, table-driven over the 3 variants since the
# 5 testset shapes are near-identical (differing only in which forms/coords are present).
using FitRateEquation
using EnzymeRates
using Random
using Test
using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy

# Per-variant expectations: which of {Ki_NADPH, Ki_ATP, Ki_ATP_EG} survive as coords, and the
# anchor_reverse convention each is fit with (informational here; enforced at run_all call time).
const _DEADEND_VARIANTS = [
    (name=:no_g6p_nadph_deadend, has_nadph=false, has_atp=true,  has_atp_eg=true),
    (name=:no_g6p_atp_deadend,   has_nadph=true,  has_atp=true,  has_atp_eg=false),
    (name=:no_g6p_both_deadends, has_nadph=false, has_atp=true,  has_atp_eg=false),
]

_expected_coords(v) = Set(vcat(
    [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Km_NADPH_rev],
    v.has_nadph  ? [:Ki_NADPH]  : Symbol[],
    v.has_atp    ? [:Ki_ATP]    : Symbol[],
    v.has_atp_eg ? [:Ki_ATP_EG] : Symbol[],
))

# Both-product grid including ATP (all 3 variants keep ATP as a metabolite, unlike :no_atp) —
# same shape as test_cha_laws.jl's v2 exactness grid.
const _GRID = [
    (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=0.0,  ATP=0.0),
    (; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=1e-4, ATP=0.0),
    (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=0.0,  ATP=0.0),
    (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0),
    (; NADP=1e-6, G6P=10e-6, NADPH=2e-6, PGLn=5e-5, ATP=1e-3),
    (; NADP=2e-5, G6P=80e-6, NADPH=8e-6, PGLn=2e-4, ATP=5e-4),
]

for v in _DEADEND_VARIANTS
    vname = v.name

    @testset "$vname mechanism registration" begin
        vs = FitRateEquation.consensus_variants(:G6PD)
        names = [Symbol(x.name) for x in vs]
        @test vname in names
        # Pre-existing deploy variant indices are unchanged.
        @test Symbol(vs[1].name) === :RE_rate_eq
        @test Symbol(vs[2].name) === :SS_NADPH_release_rate_eq
        m = vs[findfirst(x -> Symbol(x.name) === vname, vs)].mech
        @test :ATP in EnzymeRates.metabolites(m)     # unlike :no_atp, ATP stays a metabolite
        @test :NADPH in EnzymeRates.metabolites(m)
    end

    @testset "$vname exactness gate: cha_rate ≡ rate_equation" begin
        Random.seed!(20260720)
        vs = FitRateEquation.consensus_variants(:G6PD)
        m  = vs[findfirst(x -> Symbol(x.name) === vname, vs)].mech
        mets = EnzymeRates.metabolites(m)
        maxrel = 0.0
        for _ in 1:25
            free = free_params(m)
            logθ = -3 .+ 2 .* rand(length(free))
            keq  = 10.0
            mac  = ChaInvert.cha_macro_readoffs_G6PD(m, logθ; keq=keq)
            v.has_nadph  || @test isinf(mac.Ki_NADPH)
            v.has_atp_eg || @test isinf(mac.Ki_ATP_EG)
            vsat = abs(EnzymeRates.rate_equation(m,
                NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:G6P) ? 1e-2 : 0.0 for s in mets)),
                build_params(m, logθ; keq=keq)))
            for conc in _GRID
                cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
                vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
                vcha = ChaLaws.cha_rate_G6PD(mac; conc...)
                scale = max(abs(vref), abs(vcha), 1e-6 * vsat)
                maxrel = max(maxrel, abs(vcha - vref) / scale)
            end
        end
        @test maxrel <= 1e-12
    end

    @testset "$vname reduced coords + macro tuple" begin
        c = ChaFit.cha_coords(:G6PD, vname)
        @test Set(c) == _expected_coords(v)
        v.has_nadph  || @test :Ki_NADPH  ∉ c
        v.has_atp_eg || @test :Ki_ATP_EG ∉ c

        coords = Dict{Symbol,Float64}(:Kd_NADP=>5e-6, :Kd_G6P=>25e-6, :Kd_6PGLn=>1e-3,
                                       :alpha=>4.0, :Km_NADPH_rev=>3.9e-6)
        v.has_nadph  && (coords[:Ki_NADPH]  = 20e-6)
        v.has_atp    && (coords[:Ki_ATP]    = 1.5e-3)
        v.has_atp_eg && (coords[:Ki_ATP_EG] = 3.0e-2)

        m = ChaFit.cha_macro_tuple(:G6PD, coords; keq=13.7, variant=vname)
        v.has_nadph  || @test isinf(m.Ki_NADPH)
        v.has_atp_eg || @test isinf(m.Ki_ATP_EG)
        v.has_atp    && @test m.Ki_ATP == 1.5e-3
        # Neither ATP nor NADPH remains a pure regulator once any release/dead-end channel
        # stays live, so (unlike :no_atp's clean rate-invariance check) the exactness gate
        # above -- which anchors against the TRUE micro mechanism lacking the dropped step --
        # is the load-bearing correctness proof here; this testset only checks the readoffs.
    end

    @testset "$vname pins + modes" begin
        @test FitRateEquation.modes_for(:G6PD, vname) == (:mode1, :mode2)
        p1 = ChaFit.resolve_cha_pins(:G6PD, vname, :mode1)
        @test Set(keys(p1)) == Set([:Km_NADPH_rev])   # all-modes reverse anchor only
        p1_free = ChaFit.resolve_cha_pins(:G6PD, vname, :mode1; anchor_reverse=false)
        @test isempty(p1_free)
        p2 = ChaFit.resolve_cha_pins(:G6PD, vname, :mode2)
        @test (:Ki_NADPH in keys(p2)) == v.has_nadph
        @test (:Ki_ATP in keys(p2))   == v.has_atp
        @test :Ki_ATP_EG ∉ keys(p2)   # never a literature pin target (K9 has no lit anchor)
        @test :Km_NADPH_rev in keys(p2)
    end

    @testset "$vname deploy inverse round-trips" begin
        Random.seed!(20260721)
        vs = FitRateEquation.consensus_variants(:G6PD)
        m  = vs[findfirst(x -> Symbol(x.name) === vname, vs)].mech
        mets = EnzymeRates.metabolites(m)
        coords = Dict{Symbol,Float64}(:Kd_NADP=>5e-6, :Kd_G6P=>25e-6, :Kd_6PGLn=>1e-3,
                                       :alpha=>4.0, :Km_NADPH_rev=>3.9e-6)
        v.has_nadph  && (coords[:Ki_NADPH]  = 20e-6)
        v.has_atp    && (coords[:Ki_ATP]    = 1.5e-3)
        v.has_atp_eg && (coords[:Ki_ATP_EG] = 3.0e-2)
        keq = 13.7
        rr  = ChaFit.CHA_DEPLOY_RELEASE_RATE
        logθ = ChaDeploy.cha_deploy_micro(:G6PD, m, coords; keq=keq, koffQ=rr, release_rate=rr)
        @test length(logθ) == length(free_params(m))
        mac = ChaFit.cha_macro_tuple(:G6PD, coords; keq=keq, variant=vname,
                                     release_rate=rr, release_eq=coords[:Km_NADPH_rev])
        p = build_params(m, logθ; keq=keq)
        for conc in [(; NADP=5e-6, G6P=40e-6, NADPH=0.0, PGLn=0.0, ATP=0.0),
                     (; NADP=1e-5, G6P=80e-6, NADPH=4e-6, PGLn=1e-4, ATP=5e-4)]
            cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, p)
            vcha = ChaLaws.cha_rate_G6PD(mac; conc...)
            @test isapprox(vref, vcha; rtol=1e-9, atol=1e-30)
        end
    end

    @testset "runner: $vname cell list" begin
        cells = FitRateEquation._cells(:G6PD; variants=[vname])
        @test length(cells) == 2
        @test all(c -> c[1] === vname, cells)
        @test [c[3] for c in cells] == [:mode1, :mode2]
    end
end

@testset "deploy variant cell list unaffected by new registrations" begin
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:G6PD),
              FitRateEquation._cells(:G6PD))
end

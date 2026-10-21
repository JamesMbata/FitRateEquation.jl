# =========================================================================================
#                    ATP-free PGD Cha variant (:no_atp) — unit tests
# =========================================================================================
# Proves the ATP-free PGD law is the exact Ki_ATP,Ki_ATP_EN -> Inf limit of the consensus
# :cha_base law, and that the reduced-coordinate fit / deploy / pin / mode wiring is consistent.
using FitRateEquation
using FitRateEquation: pgd_config
using EnzymeRates
using CSV, DataFrames
using Random
using Test
using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert
using FitRateEquation.ChaFit
using FitRateEquation.ChaDeploy

_pgd_mech(name::Symbol) = (vs = FitRateEquation.consensus_variants(:PGD);
                           vs[findfirst(v -> Symbol(v.name) === name, vs)].mech)

const _PGD_NOATP_COORDS = Dict(:Kd_NADP=>3e-6, :Kd_PGA=>2e-5, :alpha=>2.0, :Kd_CO2=>1e-2,
                               :Ki_NADPH=>17e-6, :Km_NADPH_rev=>5e-6)

@testset "PGD no_atp mechanism registration" begin
    names = [Symbol(v.name) for v in FitRateEquation.consensus_variants(:PGD)]
    @test :no_atp in names
    @test :cha_base in names && :full_re in names
    m = _pgd_mech(:no_atp)
    @test !(:ATP in EnzymeRates.metabolites(m))            # ATP is not a metabolite at all
    @test :NADPH in EnzymeRates.metabolites(m)             # NADPH dead-end on E·PGA kept
end

@testset "PGD no_atp exactness gate: cha_rate_PGD ≡ rate_equation(no_atp mechanism)" begin
    Random.seed!(20261021)
    m    = _pgd_mech(:no_atp)
    mets = EnzymeRates.metabolites(m)                      # no :ATP
    grid = [
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=0.0,  NADPH=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=5e-5, NADPH=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=5e-6),
        (; NADP=2e-5, PGA=80e-6, CO2=2e-4, Ru5P=1e-4, NADPH=8e-6),
    ]
    for _ in 1:20
        logθ = -3 .+ 2 .* rand(length(free_params(m)))
        keq  = 0.1606
        mac  = cha_macro_readoffs_PGD(m, logθ; keq=keq)
        @test isinf(mac.Ki_ATP) && isinf(mac.Ki_ATP_EN)    # ATP terms structurally absent
        vsat = abs(er_rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc   = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = er_rate_equation(m, cc, build_params(m, logθ; keq=keq))
            @test isapprox(cha_rate_PGD(mac; conc...), vref; rtol=1e-10, atol=1e-10 * vsat)
        end
    end
end

@testset "PGD no_atp reduced coords + macro tuple" begin
    c = cha_coords(:PGD, :no_atp)
    @test Set(c) == Set([:Kd_NADP, :Kd_PGA, :alpha, :Kd_CO2, :Ki_NADPH, :Km_NADPH_rev])
    lo, hi = cha_coord_bounds(:PGD, :no_atp)
    @test length(lo) == length(hi) == 6
    mt = cha_macro_tuple(:PGD, _PGD_NOATP_COORDS; keq=0.17, variant=:no_atp)
    @test isinf(mt.Ki_ATP) && isinf(mt.Ki_ATP_EN)
    pt = (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=5e-5, NADPH=3e-6)
    @test cha_rate_PGD(mt; pt..., ATP=2e-3) == cha_rate_PGD(mt; pt..., ATP=0.0)   # ATP-blind
end

@testset "PGD no_atp loss == cha_base loss on the ATP-free corpus (any ATP constants)" begin
    # Every PGD corpus row has [ATP] = 0, so the two laws share one loss surface over the six
    # shared coords: the cha_base ATP constants are inert there, whatever their value.
    d = load_dataset(pgd_config())
    @test all(c -> getfield(c, :ATP) == 0.0, d.concs)
    L0 = cha_centered_logratio_loss(:PGD, _pgd_mech(:no_atp), d, _PGD_NOATP_COORDS;
                                    variant=:no_atp)
    for (ka, ken) in ((1e-3, 1e-6), (1e-9, 1.0))
        cb = merge(_PGD_NOATP_COORDS, Dict(:Ki_ATP=>ka, :Ki_ATP_EN=>ken))
        @test cha_centered_logratio_loss(:PGD, _pgd_mech(:cha_base), d, cb) == L0
    end
end

@testset "PGD no_atp pins + modes" begin
    @test FitRateEquation.modes_for(:PGD, :no_atp) == (:mode1, :mode2)
    @test isempty(resolve_cha_pins(:PGD, :no_atp, :mode1))
    p2 = resolve_cha_pins(:PGD, :no_atp, :mode2)
    @test Set(keys(p2)) == Set([:Ki_NADPH])
    @test isapprox(p2[:Ki_NADPH], log10(17e-6); atol=1e-9)
    cells = FitRateEquation._cells(:PGD; variants=[:no_atp])
    @test [c[1] for c in cells] == [:no_atp, :no_atp]
    @test [c[3] for c in cells] == [:mode1, :mode2]
    @test FitRateEquation.variant_profile(:PGD, :no_atp) == (row_filter=identity, label="noatp")
end

@testset "PGD no_atp deploy inverse round-trips" begin
    m    = _pgd_mech(:no_atp)
    mets = EnzymeRates.metabolites(m)
    keq  = 0.17
    rr   = CHA_DEPLOY_RELEASE_RATE
    logθ = cha_deploy_micro(:PGD, m, _PGD_NOATP_COORDS; keq=keq, koffQ=rr, release_rate=rr,
                            variant=:no_atp)
    @test length(logθ) == length(free_params(m))
    mac = cha_macro_tuple(:PGD, _PGD_NOATP_COORDS; keq=keq, variant=:no_atp, release_rate=rr)
    p   = build_params(m, logθ; keq=keq)
    for conc in [(; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=0.0),
                 (; NADP=1e-5, PGA=80e-6, CO2=1e-4, Ru5P=5e-5, NADPH=4e-6)]
        cc = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
        @test isapprox(er_rate_equation(m, cc, p), cha_rate_PGD(mac; conc...);
                       rtol=1e-9, atol=1e-30)
    end
end

@testset "PGD no_atp plot adapter evaluates the ATP-free law" begin
    a = FitRateEquation.build_cha_adapter(:PGD, _PGD_NOATP_COORDS, :no_atp, 0.17)
    pt = (NADP=5e-6, PGA=40e-6, CO2=0.0, Ru5P=0.0, NADPH=3e-6, ATP=0.0)
    mt = cha_macro_tuple(:PGD, _PGD_NOATP_COORDS; keq=0.17, variant=:no_atp,
                         release_rate=CHA_DEPLOY_RELEASE_RATE)
    @test er_rate_equation(a, pt, (Keq=0.17, E_total=1.0)) ≈ cha_rate_PGD(mt; pt...) rtol=1e-12
end

@testset "PGD no_atp end-to-end run writes a two-mode, ATP-free, deployable output" begin
    outdir = mktempdir()
    res = FitRateEquation._fit_consensus(pgd_config(); outdir=outdir, variants=[:no_atp],
                                         n_restarts=1, maxiter=40, maxtime=30.0, seed=1)
    @test [r.mode for r in res] == [:mode1, :mode2]
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    @test all(==("no_atp"), string.(mc.variant))
    @test !any(in(("Ki_ATP", "Ki_ATP_EN")), string.(mc.name))
    ki2 = mc[(string.(mc.name) .== "Ki_NADPH") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki2) == 1 && string(ki2.class[1]) == "literature_pinned"
    mp = read(joinpath(outdir, "micro_parameters.jl"), String)
    @test occursin("no_atp_mode1_", mp) && occursin("no_atp_mode2_", mp)
    @test !occursin("deploy unavailable", mp)
end

@testset "PGD full_re end-to-end run: one mode, no agreement section, C=1 apparent Km" begin
    outdir = mktempdir()
    res = FitRateEquation._fit_consensus(pgd_config(); outdir=outdir, variants=[:full_re],
                                         n_restarts=1, maxiter=40, maxtime=30.0, seed=1)
    @test [r.mode for r in res] == [:mode1]
    rep = read(joinpath(outdir, "report.md"), String)
    @test !occursin("Forward-constant agreement", rep)    # nothing to compare with one mode
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    v(n) = mc.value[findfirst(==(n), string.(mc.name))]
    @test v("Km_PGA") == v("alpha") * v("Kd_PGA")        # fully-RE law is fiber-free: C = 1
    @test !occursin("deploy unavailable", read(joinpath(outdir, "micro_parameters.jl"), String))
end

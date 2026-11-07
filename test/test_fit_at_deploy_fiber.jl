# ##########################################################################################
#                  The fitted law IS the deployed law (same promoted release rate)          #
# ##########################################################################################
# Regression for the fit/deploy fiber split: the fit ran at koffQ = 1 (C = 2) while
# micro_parameters.jl and the plots used the same alpha at CHA_DEPLOY_RELEASE_RATE = 1e3
# (C ~ 1). In the forward law only alpha/C matters, so the deployed law was not the fitted
# law: on the G6PD corpus its loss was ~30% above the reported in-sample loss and its apparent
# Km's were twice the fit's.
#
# The check goes through the artifacts a downstream user consumes: parse the emitted
# micro_parameters.jl, build the law with EnzymeRates' own rate_equation (not the cha_rate_*
# path the fit used), and compare its centered loss with goodness_of_fit.csv's in_sample.
using FitRateEquation
using FitRateEquation: g6pd_config, pgd_config, ChaFit, _SIGN_PENALTY
using EnzymeRates
using CSV, DataFrames
using Test

# One (variant, mode) block of an emitted micro_parameters.jl: free-param symbol => value.
function _emitted_block(outdir, variant::Symbol, mode::Symbol)
    prefix = "$(variant)_$(mode)_"
    block = Dict{Symbol,Float64}()
    for line in readlines(joinpath(outdir, "micro_parameters.jl"))
        startswith(line, prefix) || continue
        name, val = strip.(split(line, " = "))
        block[Symbol(name[length(prefix)+1:end])] = parse(Float64, val)
    end
    block
end

# Centered log-ratio loss of any rate function `v(concs, keq)`, written out independently of
# ChaFit: per-figure mean-centering, sign/zero penalty, total / n.
function _independent_centered_loss(v, d)
    total = 0.0
    for g in unique(d.group)
        idx = findall(==(g), d.group)
        keq = only(unique(d.keq[idx]))
        lr = Float64[]
        for i in idx
            p = v(d.concs[i], keq)
            if !isfinite(p) || p == 0 || sign(p) != sign(d.rate[i])
                total += _SIGN_PENALTY
            else
                push!(lr, log(abs(p)) - log(abs(d.rate[i])))
            end
        end
        isempty(lr) || (total += sum(abs2, lr .- sum(lr) / length(lr)))
    end
    total / length(d.rate)
end

# The main fits only (no CV folds: this test is about the fit -> artifacts path) at a small
# budget, written out by the real write_outputs.
function _fit_and_write(cfg, enzyme, variants, outdir)
    d     = FitRateEquation.load_dataset(cfg)
    cells = FitRateEquation._cells(enzyme; variants=variants)
    mechs = Dict(v => FitRateEquation._mech_for(enzyme, v) for v in variants)
    tasks = filter(t -> t.kind === :main, FitRateEquation._build_tasks(cells, d; enzyme=enzyme))
    raw   = [FitRateEquation._run_fit_task(t, d, mechs; n_restarts=2, maxiter=150, maxtime=5.0,
                                           enzyme=enzyme) for t in tasks]
    res   = FitRateEquation._reduce_cells(raw, cells, d, mechs; enzyme=enzyme)
    FitRateEquation.write_outputs(outdir, d, res; name=String(cfg.name), enzyme=enzyme,
                                  deploy_keq=cfg.deploy_keq, corpus=FitRateEquation.read_corpus(cfg),
                                  data_csv=cfg.data_csv)
    d, mechs, res
end

@testset "one release rate: fit fiber == deploy fiber" begin
    for enzyme in (:G6PD, :PGD)
        @test ChaFit._default_release_rate(enzyme) === ChaFit.CHA_DEPLOY_RELEASE_RATE
    end
    # Absolute mode keeps its own fiber-free limit (C = 1 at any kcat in bounds).
    @test ChaFit.CHA_ABS_RELEASE_RATE == 1.0e8
end

@testset "emitted constants reproduce the reported in-sample loss ($enzyme)" for (enzyme, cfg, variants) in (
        (:G6PD, g6pd_config(), [:no_g6p_both_deadends, :no_g6p_atp_deadend]),
        (:PGD,  pgd_config(),  [:cha_base, :no_atp]))
    outdir = mktempdir()
    d, mechs, res = _fit_and_write(cfg, enzyme, variants, outdir)
    gof = CSV.read(joinpath(outdir, "goodness_of_fit.csv"), DataFrame)
    mc  = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    for r in res
        reported = only(gof[(gof.variant .== String(r.variant)) .& (gof.mode .== String(r.mode)), :in_sample])
        @test reported == r.r.fit.loss

        # (1) micro_parameters.jl through EnzymeRates' own rate equation.
        mech  = mechs[r.variant]
        block = _emitted_block(outdir, r.variant, r.mode)
        fp    = free_params(mech)
        @test Set(keys(block)) == Set(fp)
        logθ  = [log10(block[s]) for s in fp]
        mets  = Tuple(EnzymeRates.metabolites(mech))
        emitted = (c, keq) -> er_rate_equation(mech,
            NamedTuple{mets}(Tuple(hasproperty(c, s) ? getproperty(c, s) : 0.0 for s in mets)),
            build_params(mech, logθ; keq=keq))
        @test isapprox(_independent_centered_loss(emitted, d), reported; rtol=1e-6)

        # (2) macro_constants.csv through the Cha law at the deploy release rate.
        coords = FitRateEquation.read_coords(mc, enzyme, r.variant, r.mode)
        @test isapprox(ChaFit.cha_centered_logratio_loss(enzyme, mech, d, coords; variant=r.variant,
                                                         release_rate=ChaFit.CHA_DEPLOY_RELEASE_RATE),
                       reported; rtol=1e-6)

        # (3) the plotted curve: the adapter the renderer draws through.
        adapter = FitRateEquation.build_cha_adapter(enzyme, coords, r.variant, cfg.deploy_keq)
        plotted = (c, keq) -> EnzymeRates.rate_equation(adapter, c, (; Keq=keq))
        @test isapprox(_independent_centered_loss(plotted, d), reported; rtol=1e-6)
    end
end

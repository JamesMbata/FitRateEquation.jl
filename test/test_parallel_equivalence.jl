using FitRateEquation
using FitRateEquation: g6pd_config
using EnzymeRates
using Distributed
using Statistics
using Test


# Replicates the serial fit_consensus_equation cell loop on the CHA PATH (using the now-Cha `_fit_and_cv`
# + ChaClassify) so the new pmap path can be checked against it. With a maxiter-bound,
# seeded CMA-ES the fit is a pure function of (data, seed, maxiter), so the pmap path and
# this serial baseline are BIT-IDENTICAL — the integration test below asserts exact equality
# on fits, losses, and macro constants (not merely tolerance). The worker-COUNT invariance of
# the reduction is additionally pinned by the deterministic task-generation test. Mirrors
# `_reduce_cells`'s classify block exactly (ChaClassify.cha_identifiable_functions +
# classify_cha, same sigma2).
function _serial_baseline(cfg; n_restarts, maxiter, maxtime, seed)
    d = load_dataset(cfg)
    enzyme = Symbol(cfg.name)
    keq = median(d.keq)
    cells = FitRateEquation._cells(enzyme)
    results = NamedTuple[]
    for (ci, (variant, mech, mode)) in enumerate(cells)
        r = FitRateEquation._fit_and_cv(variant, mech, d; mode=mode, enzyme=enzyme,
                                       n_restarts=n_restarts, maxiter=maxiter,
                                       maxtime=maxtime, seed=seed + ci)
        pins = r.pins
        idf = FitRateEquation.ChaClassify.cha_identifiable_functions(enzyme, mech, d, r.fit.coords;
                                                                    keq=keq, pins=pins)
        sigma2 = r.fit.loss / max(nrows(d) - idf.rank, 1)
        classed = FitRateEquation.ChaClassify.classify_cha(enzyme, mech, d, r.fit.coords, pins, idf;
                                                          keq=keq, variant=variant, mode=mode,
                                                          sigma2=sigma2)
        push!(results, (variant=variant, mode=mode, r=r, idf=idf, classed=classed))
    end
    results
end

@testset "task generation (deterministic, worker-invariant)" begin
    d = load_dataset(g6pd_config())
    cells = FitRateEquation._cells()
    # Every cell is cross-validated under BOTH schemes: its fold tasks are the article folds,
    # then the figure folds (CV_SCHEMES order), each in `_cv_folds` order.
    folds = [(scheme=sc, fold=f) for sc in FitRateEquation.CV_SCHEMES
                                 for f in FitRateEquation._cv_folds(d, sc)]
    n_folds = length(folds)
    @test n_folds == length(FitRateEquation._cv_folds(d, :article)) +
                     length(FitRateEquation._cv_folds(d, :figure))
    tasks = FitRateEquation._build_tasks(cells, d; seed=1)

    # Deploy-variant-only: G6PD = 1 variant (:SS_NADPH_release_rate_eq) × 2 modes = 2 cells;
    # n_cells × (1 + n_article_folds + n_figure_folds) tasks.
    @test length(cells) == 2
    @test all(v === :SS_NADPH_release_rate_eq for (v, _, _) in cells)
    @test length(tasks) == length(cells) * (1 + n_folds)

    allrows = collect(1:nrows(d))
    for (ci, (variant, mech, mode)) in enumerate(cells)
        cell_tasks = [t for t in tasks if t.ci == ci]
        @test length(cell_tasks) == 1 + n_folds

        # Exactly one :main per cell: full train rows, empty test, seed = base + ci.
        mains = [t for t in cell_tasks if t.kind === :main]
        @test length(mains) == 1
        m = mains[1]
        @test m.variant === variant && m.mode === mode
        @test m.seed == 1 + ci
        @test m.train_idx == allrows
        @test isempty(m.test_idx)

        # Folds follow CV_SCHEMES × _cv_folds order, share the cell seed, carry the right subsets.
        cell_folds = [t for t in cell_tasks if t.kind === :fold]
        @test length(cell_folds) == n_folds
        for (f, sf) in zip(cell_folds, folds)
            @test f.kind === :fold
            @test f.seed == 1 + ci
            @test f.scheme === sf.scheme
            @test f.held_out == sf.fold.held_out
            @test f.train_idx == sf.fold.train
            @test f.test_idx == sf.fold.test
        end

        # Pins are the resolved Cha values, shared per cell. Default (anchor_reverse=false):
        # Mode 1 pins nothing, Mode 2 pins Ki_ATP + Ki_NADPH. The pmap path threads exactly
        # these, so the reduction is worker-invariant.
        @test m.pins == FitRateEquation.ChaFit.resolve_cha_pins(:G6PD, variant, mode)
        @test !haskey(m.pins, :Km_NADPH_rev)                          # reverse anchor is opt-in
        if mode === :mode2
            @test Set(keys(m.pins)) == Set([:Ki_ATP, :Ki_NADPH])
        else
            @test isempty(m.pins)
        end
        @test all(t.pins == m.pins for t in cell_tasks)
    end
end

@testset "single-process fit_consensus_equation (pmap) == serial baseline (bit-identical, smoke)" begin
    # Single-process pool == the serial path; never spawns workers, so it does not
    # contend with any other worktree's workers. The budget is maxiter-bound (seeded
    # CMA-ES), so fit_consensus_equation's pmap path is bit-identical to the serial baseline.
    @test nworkers() == 1
    budget = (n_restarts=2, maxiter=150, maxtime=120.0)

    par = FitRateEquation._fit_consensus(g6pd_config(); outdir=mktempdir(), n_restarts=budget.n_restarts,
                  maxiter=budget.maxiter, maxtime=budget.maxtime, seed=1)
    ser = _serial_baseline(g6pd_config(); n_restarts=budget.n_restarts,
                           maxiter=budget.maxiter, maxtime=budget.maxtime, seed=1)

    @test length(par) == length(ser) == 2
    for (p, s) in zip(par, ser)
        @test (p.variant, p.mode) == (s.variant, s.mode)
        # Fit is a pure function of (data, seed, maxiter): bit-identical coords and loss.
        @test keys(p.r.fit.coords) == keys(s.r.fit.coords)
        @test all(isequal(p.r.fit.coords[k], s.r.fit.coords[k]) for k in keys(p.r.fit.coords))
        @test isequal(p.r.fit.loss, s.r.fit.loss)
        # CV reduction (fold order, mean, se, per-fold losses) reproduces _cha_loocv exactly,
        # for both schemes.
        @test keys(p.r.cv) == keys(s.r.cv) == FitRateEquation.CV_SCHEMES
        for sc in FitRateEquation.CV_SCHEMES
            pc, sc_ = p.r.cv[sc], s.r.cv[sc]
            @test pc.n_folds == sc_.n_folds > 0
            @test isequal(pc.mean_cv, sc_.mean_cv)
            @test isequal(pc.se, sc_.se)
            @test [a.held_out for a in pc.per_fold] == [a.held_out for a in sc_.per_fold]
            @test [a.n_test for a in pc.per_fold] == [a.n_test for a in sc_.per_fold]
            @test isequal([a.loss for a in pc.per_fold], [a.loss for a in sc_.per_fold])
        end
        # Macro constants: same names, classes, and (bit-identical) values.
        @test [m.name for m in p.classed] == [m.name for m in s.classed]
        @test [m.class for m in p.classed] == [m.class for m in s.classed]
        @test isequal([m.value for m in p.classed], [m.value for m in s.classed])
    end
end

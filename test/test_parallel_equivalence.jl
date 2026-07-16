using FitRateEquation
using EnzymeRates
using Distributed
using Statistics
using Test


# Replicates the serial run_all cell loop on the CHA PATH (using the now-Cha `_fit_and_cv`
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
    folds = FitRateEquation._article_folds(d)
    n_articles = length(folds)
    tasks = FitRateEquation._build_tasks(cells, d; seed=1)

    # Deploy-variant-only: G6PD = 1 variant (:SS_NADPH_release_rate_eq) × 2 modes = 2 cells;
    # n_cells × (1 + n_articles) = 2 × (1 + 7) = 16 for the G6PD corpus.
    @test length(cells) == 2
    @test all(v === :SS_NADPH_release_rate_eq for (v, _, _) in cells)
    @test length(tasks) == length(cells) * (1 + n_articles)

    allrows = collect(1:nrows(d))
    for (ci, (variant, mech, mode)) in enumerate(cells)
        cell_tasks = [t for t in tasks if t.ci == ci]
        @test length(cell_tasks) == 1 + n_articles

        # Exactly one :main per cell: full train rows, empty test, seed = base + ci.
        mains = [t for t in cell_tasks if t.kind === :main]
        @test length(mains) == 1
        m = mains[1]
        @test m.variant === variant && m.mode === mode
        @test m.seed == 1 + ci
        @test m.train_idx == allrows
        @test isempty(m.test_idx)

        # Folds follow _article_folds order, share the cell seed, carry the right subsets.
        cell_folds = [t for t in cell_tasks if t.kind === :fold]
        @test length(cell_folds) == n_articles
        for (f, fold) in zip(cell_folds, folds)
            @test f.kind === :fold
            @test f.seed == 1 + ci
            @test f.article == fold.article
            @test f.train_idx == fold.train
            @test f.test_idx == fold.test
        end

        # Pins + anchors are the resolved Cha values, shared per cell. G6PD always-anchors
        # the conflating reverse channel Km_NADPH_rev (ALL modes); Mode-2 additionally pins
        # Ki_ATP + Ki_NADPH. G6PD carries NO Km_PGA anchor (that is PGD-only). The pmap path
        # threads exactly these, so the reduction is worker-invariant.
        @test m.pins == FitRateEquation.ChaFit.resolve_cha_pins(:G6PD, variant, mode)
        @test m.anchors === FitRateEquation.cha_anchors(:G6PD, mode)   # nothing for G6PD
        @test haskey(m.pins, :Km_NADPH_rev)                           # always-anchored (all modes)
        if mode === :mode2
            @test haskey(m.pins, :Ki_ATP) && haskey(m.pins, :Ki_NADPH)
        else
            @test Set(keys(m.pins)) == Set([:Km_NADPH_rev])           # Mode 1: only the reverse anchor
        end
        @test all(t.pins == m.pins for t in cell_tasks)
        @test all(t.anchors === m.anchors for t in cell_tasks)
    end
end

@testset "single-process run_all (pmap) == serial baseline (bit-identical, smoke)" begin
    # Single-process pool == the serial path; never spawns workers, so it does not
    # contend with any other worktree's workers. The budget is maxiter-bound (seeded
    # CMA-ES), so run_all's pmap path is bit-identical to the serial baseline.
    @test nworkers() == 1
    budget = (n_restarts=2, maxiter=150, maxtime=120.0)

    par = run_all(g6pd_config(); outdir=mktempdir(), n_restarts=budget.n_restarts,
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
        # CV reduction (fold order, mean, se) reproduces _cha_loocv exactly.
        @test isequal(p.r.cv.mean_cv, s.r.cv.mean_cv)
        @test isequal(p.r.cv.se, s.r.cv.se)
        @test [a.article for a in p.r.cv.per_article] == [a.article for a in s.r.cv.per_article]
        @test isequal([a.loss for a in p.r.cv.per_article], [a.loss for a in s.r.cv.per_article])
        # Macro constants: same names, classes, and (bit-identical) values.
        @test [m.name for m in p.classed] == [m.name for m in s.classed]
        @test [m.class for m in p.classed] == [m.class for m in s.classed]
        @test isequal([m.value for m in p.classed], [m.value for m in s.classed])
    end
end

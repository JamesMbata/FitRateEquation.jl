# =========================================================================================
#         Pipeline: fit the MACRO constants DIRECTLY via the Cha path (Phase B)
# =========================================================================================
#
# The Cha law is the PROMOTED law; the pure-RE variant is just its koffQ->infinity limit, so
# the pipeline runs ONLY the deploy variant per enzyme. Each cell is a (deploy_variant, mech,
# mode); fitting optimizes the NAMED forward shape constants (cha_coords) directly through
# `cha_rate_*` (ChaFit), classifies them in macro-coordinate space (ChaClassify), and deploys a
# closed-form micro representative (ChaDeploy). Mode 2/3 anchor the product side to literature
# via hard coord-pins (resolve_cha_pins) plus, for PGD Km_PGA, a soft/high-weight apparent-Km
# anchor (cha_anchors -> ChaFit anchors). The pmap/CachingPool/per-cell-seed determinism
# scaffolding (`seed+ci`, order-preserving reduce) is preserved byte-identical in STRUCTURE;
# only the fit/loss/classify/deploy CALLS changed (coeff_* retired on this path; Task 16 deletes).

# The deploy variant per enzyme — the single PROMOTED Cha law the pipeline fits. The pure-RE
# variant is its release_rate->infinity limit, so it is NOT run separately.
deploy_variant(enzyme::Symbol) =
    enzyme === :G6PD ? :SS_NADPH_release_rate_eq :
    enzyme === :PGD  ? :cha_base :
    enzyme === :HK1  ? :H1 :
    error("deploy_variant: unknown enzyme $enzyme (expected :G6PD, :PGD, or :HK1)")

# The deploy mechanism: select it from consensus_variants(enzyme) by the deploy-variant name.
function _deploy_mech(enzyme::Symbol)
    v = deploy_variant(enzyme)
    for x in consensus_variants(enzyme)
        Symbol(x.name) === v && return x.mech
    end
    error("_deploy_mech: deploy variant $v not found in consensus_variants($enzyme)")
end

# Variants the run sweeps. G6PD/PGD run only their single deploy variant; HK1 runs H1 (alpha=:one,
# the raw {Ki_G6P_C, Ki_G6P_N} parameterization, 3 modes) and H4 (the SAME alpha=1 law but
# reparameterized in the data-identifiable {Keff, split_ratio}, mode1-only). H3 (alpha=:infinity)
# was REMOVED: the reverse-rate turnover requires the [G6P]² term H3 deletes
# (notes/2026-06-13_hk1_g6p_ridge_resolution_report.md).
run_variants(enzyme::Symbol) =
    enzyme === :HK1 ? [:H1, :H4] : [deploy_variant(enzyme)]

# G6PD dead-end-dropped ablation variants (src/enzymes/g6pd.jl) for which
# `anchor_reverse=false` is now a supported default, not a diagnostic: fit on the smaller,
# more-physiological "mydata" corpus (134 rows), `no_g6p_atp_deadend` jointly data-identifies
# BOTH Km_NADPH_rev and Ki_NADPH without the anchor (2026-07-21) -- the deployed law
# (`SS_NADPH_release_rate_eq`) and the raw conflating `RE_rate_eq` do NOT show this on the
# full 565-row historical corpus (Ki_NADPH stays unconstrained), so the anchor requirement is
# unchanged for those. This has not been confirmed on the full corpus, so the exemption stays
# scoped to exactly these three variants -- see AGENTS.md for the corpus caveat.
const _G6PD_ANCHOR_OPTIONAL_VARIANTS = (:no_g6p_nadph_deadend, :no_g6p_atp_deadend,
                                        :no_g6p_both_deadends)

# Whether `variant` still REQUIRES the reverse anchor for deployability (the deployed law and
# the raw conflating RE law, for G6PD; always true for PGD/HK1, which have no such anchor).
_requires_reverse_anchor(enzyme::Symbol, variant::Symbol) =
    !(enzyme === :G6PD && variant in _G6PD_ANCHOR_OPTIONAL_VARIANTS)

# Default `anchor_reverse` for a `run_all` call: false only when EVERY variant being fit is one
# of the anchor-optional ablation variants; true otherwise (unchanged default everywhere else,
# including any call that mixes an ablation variant with the deploy variant).
_default_anchor_reverse(enzyme::Symbol, variants::AbstractVector{Symbol}) =
    !(enzyme === :G6PD && !isempty(variants) && all(v -> v in _G6PD_ANCHOR_OPTIONAL_VARIANTS, variants))

function _mech_for(enzyme::Symbol, variant::Symbol)
    for x in consensus_variants(enzyme)
        Symbol(x.name) === variant && return x.mech
    end
    error("_mech_for: variant $variant not found in consensus_variants($enzyme)")
end

# Per-(enzyme, mode) SOFT/OVERRIDE anchors on the DERIVED apparent Michaelis constants (the
# Km_PGA gap; PGD only). Km_PGA is the apparent constant `alpha·Kd_PGA/C`, NOT a cha_coord, so
# it can never be a hard pin — it is realized through the ChaFit `anchors` penalty:
#   G6PD: nothing (all modes).
#   PGD :mode1: nothing (unanchored data baseline).
#   PGD :mode2: soft pull to the 38–80µM band midpoint (59µM, weight 1.0 — PROVISIONAL per
#               spec §11 / open item T1-a; resolved against the Phase-C flux gate).
#   PGD :mode3: the explicit Km_PGA "hard override" realized as a HIGH-weight anchor (38µM,
#               weight 100.0) — because Km_PGA cannot be a hard coord-pin.
function cha_anchors(enzyme::Symbol, mode::Symbol)
    enzyme === :PGD || return nothing
    mode === :mode2 && return Dict(:Km_PGA => (target = log10(59e-6), weight = 1.0))
    mode === :mode3 && return Dict(:Km_PGA => (target = log10(38e-6), weight = 100.0))
    return nothing                                   # :mode1 (and any other) is unanchored
end

# One (variant, mode) cell fit + leave-one-article-out CV, via the Cha path. Mode 2/3 hard-pin
# the literature coords (resolve_cha_pins) and, for PGD, soft/override-anchor Km_PGA (the
# DERIVED apparent constant, via `anchors`). Mode 1 fits free over the always-anchored reverse
# coords only. Fit the macro coords, then CV with the same pins/anchors. (Used by tests; the
# pmap path threads the same pins/anchors per cell.)
function _fit_and_cv(variant::Symbol, mech, d::Dataset;
                     mode::Symbol=:mode2, n_restarts::Int=8, maxiter::Int=1_000_000,
                     maxtime::Real=20.0, seed::Int=1, enzyme::Symbol=_enzyme_of(mech),
                     anchor_reverse::Bool=true)
    keq     = enzyme === :HK1 ? median(d.keq) : nothing
    pins    = ChaFit.resolve_cha_pins(enzyme, variant, mode; anchor_reverse=anchor_reverse)
    anchors = cha_anchors(enzyme, mode)
    fit = ChaFit.cha_fit_candidate(enzyme, mech, d; n_restarts=n_restarts, maxiter=maxiter,
                                   maxtime=maxtime, seed=seed, keq=keq, pins=pins,
                                   anchors=anchors, variant=variant)
    cv  = _cha_loocv(enzyme, mech, d; n_restarts=n_restarts, maxiter=maxiter, maxtime=maxtime,
                     seed=seed, keq=keq, pins=pins, anchors=anchors, variant=variant)
    (variant=variant, mode=mode, mech=mech, pins=pins, anchors=anchors, fit=fit, cv=cv)
end

# Leave-one-article-out CV in Cha macro-coordinate space (the twin of loocv_by_article): each
# fold refits the cha_coords on the train rows (same pins/anchors), then scores the held-out
# article with cha_centered_logratio_loss at the train-fit coords.
function _cha_loocv(enzyme::Symbol, mech, d::Dataset; n_restarts::Int,
                    maxiter::Int, maxtime::Real, seed::Int, keq::Union{Nothing,Real},
                    pins::Dict{Symbol,Float64}, anchors, variant::Symbol=:_deploy)
    per = NamedTuple[]
    for fold in _article_folds(d)
        dtr = _subset(d, fold.train); dte = _subset(d, fold.test)
        fit = ChaFit.cha_fit_candidate(enzyme, mech, dtr; n_restarts=n_restarts, maxiter=maxiter,
                                       maxtime=maxtime, seed=seed, keq=keq, pins=pins,
                                       anchors=anchors, variant=variant)
        push!(per, (article=fold.article,
                    loss=ChaFit.cha_centered_logratio_loss(enzyme, mech, dte, fit.coords; keq=keq, variant=variant)))
    end
    losses = getfield.(per, :loss)
    (per_article=per,
     mean_cv = isempty(losses) ? NaN : mean(losses),
     se      = isempty(losses) ? 0.0 : std(losses)/sqrt(length(losses)))
end

"Compare the macro constants that are `data_identified` in BOTH modes; flag any whose
 Mode-1 and Mode-2 point estimates differ by more than `tol_dex` decades (log10). A
 disagreement means a Mode-2 'pinned' constant was not actually flat (it coupled into the
 forward fit) — surfaced in report.md, never hidden."
function mode_agreement(classed_mode1, classed_mode2; tol_dex::Real=0.5)
    by1 = Dict(c.name => c for c in classed_mode1)
    out = NamedTuple[]
    for c2 in classed_mode2
        # The agreement check validates the FORWARD shape constants (named Michaelis/
        # inhibition constants). Anonymous denominator-monomial lumps are higher-order
        # products, not pinnable, and carry no information about whether a pin was flat —
        # exclude them so noisy lumps cannot fire a spurious distortion WARNING.
        startswith(String(c2.name), "lump_") && continue
        c1 = get(by1, c2.name, nothing)
        c1 === nothing && continue
        (c1.class === :data_identified && c2.class === :data_identified) || continue
        ddex = abs(log10(c2.value) - log10(c1.value))
        push!(out, (name=c2.name, mode1=c1.value, mode2=c2.value,
                    delta_dex=ddex, agree=(ddex <= tol_dex)))
    end
    out
end

# The (deploy_variant × per-enzyme-modes) cell list. The Cha law is the PROMOTED law, so the
# pipeline runs ONLY the deploy variant per enzyme (the pure-RE variant is its release-rate->inf
# limit). G6PD = 1 variant × 2 modes = 2 cells; PGD = 1 variant × 3 modes = 3 cells. Order is
# the canonical reduction order, so any code that groups task results by `ci` reproduces the
# serial cell sequence.
_cells(enzyme::Symbol=:G6PD; variants::Vector{Symbol}=run_variants(enzyme)) =
    [(v, _mech_for(enzyme, v), mode)
     for v in variants for mode in modes_for(enzyme, v)]

# ---------------------------------------------------------------------------------------
#   Row filter for the ATP-free fit: drop any row carrying ATP (ATP > 0). The :no_atp law
#   is ATP-blind, so ATP-inhibited rows would become forced-misfit residuals biasing the
#   core constants. Used by run_g6pd_noatp.jl as run_all's `row_filter`.
# ---------------------------------------------------------------------------------------
function drop_atp_rows(d::Dataset)
    keep = [i for i in 1:nrows(d) if get(d.concs[i], :ATP, 0.0) <= 0.0]
    Dataset(d.concs[keep], d.rate[keep], d.group[keep], d.keq[keep])
end

# Flat list of independent fit-tasks across all cells: one `:main` per cell, then one
# `:fold` per article fold of that cell, in `_article_folds(d)` order. Each task is a
# small descriptor carrying only identity + row-index subsets + the resolved pins/anchors; the
# Dataset and mechanisms travel once via the CachingPool, not per task. Seed is a pure
# function of cell identity (`seed + ci`) — exactly the `_fit_and_cv` scheme, where the main
# fit and every fold of a cell share that one seed — so it is independent of worker count and
# dispatch order.
function _build_tasks(cells, d::Dataset; seed::Int=1, enzyme::Symbol=:G6PD, anchor_reverse::Bool=true)
    folds = _article_folds(d)
    allrows = collect(1:nrows(d))
    tasks = NamedTuple[]
    for (ci, (variant, mech, mode)) in enumerate(cells)
        pins    = ChaFit.resolve_cha_pins(enzyme, variant, mode; anchor_reverse=anchor_reverse)
        anchors = cha_anchors(enzyme, mode)
        s = seed + ci
        push!(tasks, (ci=ci, variant=variant, mode=mode, kind=:main, article="",
                      seed=s, train_idx=allrows, test_idx=Int[], pins=pins, anchors=anchors))
        for fold in folds
            push!(tasks, (ci=ci, variant=variant, mode=mode, kind=:fold, article=fold.article,
                          seed=s, train_idx=fold.train, test_idx=fold.test, pins=pins, anchors=anchors))
        end
    end
    tasks
end

# Worker-pure: run one fit-task and return only what the master needs to reduce. Selects
# the mechanism by variant from the captured `mechs` (EnzymeMechanism serializes fine, so
# it is shipped via the pool, not rebuilt). Fits the macro COORDINATES over the task's pins
# + anchors via ChaFit. `:main` returns the cell's coord-fit; `:fold` scores the held-out
# article with `cha_centered_logratio_loss` at the train-fit coords — identical semantics to
# `_cha_loocv` (macro-coord space).
function _run_fit_task(t, d::Dataset, mechs; n_restarts::Int, maxiter::Int, maxtime::Real,
                       enzyme::Symbol=:G6PD)
    mech = mechs[t.variant]
    dtr  = _subset(d, t.train_idx)
    keq  = enzyme === :HK1 ? median(d.keq) : nothing
    fit  = ChaFit.cha_fit_candidate(enzyme, mech, dtr; n_restarts=n_restarts, maxiter=maxiter,
                                    maxtime=maxtime, seed=t.seed, keq=keq, pins=t.pins,
                                    anchors=t.anchors, variant=t.variant)
    if t.kind === :main
        (ci=t.ci, kind=:main, fit=fit)
    else
        (ci=t.ci, kind=:fold, article=t.article,
         loss=ChaFit.cha_centered_logratio_loss(enzyme, mech, _subset(d, t.test_idx), fit.coords; keq=keq, variant=t.variant))
    end
end

# Master-side reduction: group the order-preserved `raw` results by cell and rebuild each
# cell's `r` NamedTuple (`per_article`, `mean_cv`, `se` reproduce `_cha_loocv` exactly, same
# isempty NaN/0.0 guards), then run the cheap macro-coord identifiability + classification
# serially. Cells are emitted in `ci` order so `write_outputs` is untouched.
function _reduce_cells(raw, cells, d::Dataset, mechs; seed::Int=1, enzyme::Symbol=:G6PD,
                       anchor_reverse::Bool=true)
    keq = median(d.keq)
    results = NamedTuple[]
    for (ci, (variant, mech, mode)) in enumerate(cells)
        pins    = ChaFit.resolve_cha_pins(enzyme, variant, mode; anchor_reverse=anchor_reverse)
        anchors = cha_anchors(enzyme, mode)
        mres = raw[findfirst(x -> x.ci == ci && x.kind === :main, raw)]
        fit  = mres.fit
        per  = [(article=x.article, loss=x.loss)
                for x in raw if x.ci == ci && x.kind === :fold]
        losses = getfield.(per, :loss)
        cv = (per_article=per,
              mean_cv = isempty(losses) ? NaN : mean(losses),
              se      = isempty(losses) ? 0.0 : std(losses)/sqrt(length(losses)))
        r = (variant=variant, mode=mode, mech=mech, pins=pins, anchors=anchors, fit=fit, cv=cv)
        idf = ChaClassify.cha_identifiable_functions(enzyme, mech, d, r.fit.coords;
                    keq=(enzyme === :HK1 ? median(d.keq) : nothing), pins=pins, variant=variant)
        # Residual variance σ̂² = in-sample loss / dof, for the calibrated macro-constant CIs.
        sigma2 = r.fit.loss / max(nrows(d) - idf.rank, 1)
        classed = ChaClassify.classify_cha(enzyme, mech, d, r.fit.coords, pins, idf;
                                            keq=keq, sigma2=sigma2, variant=variant, mode=mode)
        push!(results, (variant=variant, mode=mode, r=r, idf=idf, classed=classed))
    end
    results
end

"Run the deploy variant × per-enzyme modes end-to-end and write all artifacts. The
 `n_cells × (1 + n_articles)` independent fits are dispatched via `pmap` over a
 `CachingPool` (a 1-process pool is the serial path); `pmap` order preservation + the
 identity-derived per-cell seed make the reduced output byte-identical at any worker count.
 `anchor_reverse` is a G6PD-only switch — see `run_g6pd`; its default is variant-aware
 (`_default_anchor_reverse`): `false` when every variant in `variants` is one of the
 anchor-optional ablations (`_G6PD_ANCHOR_OPTIONAL_VARIANTS`), `true` otherwise. Explicitly
 passing `false` for a variant that still requires the anchor (the deploy variant,
 `:RE_rate_eq`) reproduces the conflating fit and marks that variant's output NOT DEPLOYABLE.
 No-op for PGD/HK1 (no always-on reverse anchor)."
function run_all(cfg; outdir::AbstractString, n_restarts::Int=8, maxiter::Int=1_000_000,
                 maxtime::Real=20.0, seed::Int=1,
                 variants::Vector{Symbol}=run_variants(Symbol(cfg.name)),
                 row_filter=identity,
                 anchor_reverse::Bool=_default_anchor_reverse(Symbol(cfg.name), variants))
    enzyme = Symbol(cfg.name)
    d = row_filter(load_dataset(cfg))
    deploy_keq = cfg.deploy_keq
    mkpath(outdir)
    cells = _cells(enzyme; variants=variants)
    mechs = Dict(v => _mech_for(enzyme, v) for v in variants)

    tasks = _build_tasks(cells, d; seed=seed, enzyme=enzyme, anchor_reverse=anchor_reverse)
    pool  = CachingPool(workers())
    raw   = pmap(pool, tasks) do t           # captures d, mechs, enzyme (cached per worker)
        _run_fit_task(t, d, mechs; n_restarts=n_restarts, maxiter=maxiter, maxtime=maxtime,
                      enzyme=enzyme)
    end

    results = _reduce_cells(raw, cells, d, mechs; seed=seed, enzyme=enzyme,
                            anchor_reverse=anchor_reverse)
    meta = (n_restarts=n_restarts, maxiter=maxiter, maxtime=maxtime, seed=seed,
            n_rows=nrows(d), anchor_reverse=anchor_reverse, variants=variants)
    write_outputs(outdir, d, results; meta=meta, name=String(cfg.name), enzyme=enzyme,
                 deploy_keq=deploy_keq, anchor_reverse=anchor_reverse)
    results
end

# Short git SHA of the repo containing `dir`, or "unknown" off a checkout that has no
# `.git` (an installed package depot, `Pkg.test`'s sandbox, etc.). Stderr is piped to
# `devnull` so a missing-repo `git` failure never leaks "fatal: not a git repository"
# onto the caller's stderr/test output — the try/catch already handles the outcome, the
# subprocess just should not narrate the miss.
function _git_sha(dir)
    try
        return readchomp(pipeline(`git -C $dir rev-parse --short HEAD`; stderr=devnull))
    catch
        return "unknown"
    end
end

# Self-certifying provenance: the package's own version + the EnzymeRates dependency SHA
# + the exact fit budget/seed + corpus size, so a run dir is interpretable (full vs smoke,
# which code) without inferring from dir name.
#
# The package's own identity is recorded via `pkgversion`, NOT a self-repo git SHA: once
# FitRateEquation is an installed/tested PACKAGE (not the standalone consensus_macro
# script this was ported from), `@__DIR__` under `Pkg.test`'s sandbox — and any installed
# depot — has no `.git` next to the source, so shelling out to `git` for our own version
# always failed (silently returning "unknown") while still leaking git's stderr. The
# package version is the correct, always-available identity for this field. EnzymeRates
# is a `[sources]`-pinned dependency that DOES live in a real git checkout in typical dev
# setups, so its SHA is still meaningful there; `_git_sha` is kept for it, with the same
# stderr suppression applied.
function _write_provenance(outdir, d, meta; deploy_keq::Union{Nothing,Real}=nothing)
    open(joinpath(outdir, "provenance.toml"), "w") do io
        println(io, "# FitRateEquation run provenance (auto-written)")
        println(io, "timestamp       = \"$(Libc.strftime("%Y-%m-%dT%H:%M:%S%z", time()))\"")
        println(io, "julia_version   = \"$(VERSION)\"")
        println(io, "package_version = \"$(pkgversion(@__MODULE__))\"")
        println(io, "enzymerates_sha = \"$(_git_sha(pkgdir(EnzymeRates)))\"")
        println(io, "n_restarts      = $(meta.n_restarts)")
        println(io, "maxiter         = $(meta.maxiter)")
        println(io, "maxtime         = $(meta.maxtime)")
        println(io, "seed            = $(meta.seed)")
        println(io, "n_rows          = $(meta.n_rows)")
        deploy_keq === nothing || println(io, "deploy_keq      = $(deploy_keq)")
        # Self-describing run: record the reverse-anchor state and the fitted variant(s) so a
        # conflated (anchor_reverse=false) diagnostic run dir is never mistaken for a deploy run.
        hasproperty(meta, :anchor_reverse) &&
            println(io, "anchor_reverse  = $(meta.anchor_reverse)")
        hasproperty(meta, :variants) &&
            println(io, "variants        = $(collect(String.(meta.variants)))")
        println(io, "smoke           = $(meta.maxiter <= 150)")
    end
end

# The DERIVED apparent Michaelis constants to surface per enzyme (NOT cha_coords — they are
# alpha·Kd/C readoffs via ChaFit.cha_apparent_km). Reported alongside the classed coords so
# downstream readers still see Km_G6P (G6PD) / Km_PGA (PGD). Class is :derived (a readoff, not
# a fit/pinned coord). The apparent Km is read at ChaFit.CHA_DEPLOY_RELEASE_RATE so it matches
# the DEPLOYED law (micro_parameters.jl), and each is paired with the fiber-INVARIANT
# specificity constant kcat/Km (`kcatKm_*`, class :derived) — the koffQ-robust readoff.
_apparent_kms(enzyme::Symbol) =
    enzyme === :G6PD ? (:Km_G6P,) :
    enzyme === :PGD  ? (:Km_PGA, :Km_NADP) :
    enzyme === :HK1  ? (:Km_Glc, :Km_ATP) :
    ()

function _apparent_km_rows(enzyme::Symbol, coords::AbstractDict)
    rows = NamedTuple[]
    for which in _apparent_kms(enzyme)
        km = ChaFit.cha_apparent_km(enzyme, coords, which)   # at deploy koffQ -> deployed law
        push!(rows, (name=which, value=km, class=:derived, ci=NaN))
        # HK1 apparent Km == Kd (no fiber); specificity is undefined there, so skip it.
        if enzyme !== :HK1
            spec = ChaFit.cha_specificity(enzyme, coords, which)   # kcat/Km, koffQ-invariant
            push!(rows, (name=Symbol("kcatKm_", String(which)[4:end]),
                         value=spec, class=:derived, ci=NaN))
        end
    end
    rows
end

# H4-only DERIVED rows: the reparameterized {Keff, split_ratio} coords back-map to the physical
# G6P dissociation constants {Ki_G6P_C (loose, C-half), Ki_G6P_N (tight, N-half)} and the product
# √(Kc·Kn). Surfaced as :derived so downstream readers see the physical constants alongside the
# data-identified Keff/split_ratio. (H1 reports Ki_G6P_C/Ki_G6P_N directly as its classed coords.)
function _hk1_h4_derived_rows(enzyme::Symbol, variant::Symbol, coords::AbstractDict)
    (enzyme === :HK1 && variant === :H4 && haskey(coords, :Keff)) || return NamedTuple[]
    Keff = coords[:Keff]; ratio = coords[:split_ratio]
    sqrtP = Keff * ratio; P = sqrtP^2; sumK = P / Keff
    sq = sqrt(max(sumK^2 - 4P, 0.0)); KiC = (sumK + sq) / 2; KiN = (sumK - sq) / 2
    [(name=:Ki_G6P_C, value=KiC,    class=:derived, ci=NaN),
     (name=:Ki_G6P_N, value=KiN,    class=:derived, ci=NaN),
     (name=:sqrt_KcKn, value=sqrtP, class=:derived, ci=NaN)]
end

function write_outputs(outdir, d, results; meta=nothing, name::AbstractString="G6PD",
                       enzyme::Symbol=:G6PD, deploy_keq::Real=median(d.keq),
                       anchor_reverse::Bool=true)
    meta === nothing || _write_provenance(outdir, d, meta; deploy_keq=deploy_keq)
    keq = deploy_keq
    # macro_constants.csv (keyed by variant × mode): the classed cha_coords PLUS the derived
    # apparent Michaelis constants (Km_G6P / Km_PGA), so downstream readers see the named
    # forward Km's even though they are readoffs, not fit coords.
    rows = NamedTuple[]
    for res in results
        for m in res.classed
            push!(rows, (variant=res.variant, mode=res.mode, name=m.name,
                         value=m.value, class=m.class, ci=m.ci))
        end
        for m in _apparent_km_rows(enzyme, res.r.fit.coords)
            push!(rows, (variant=res.variant, mode=res.mode, name=m.name,
                         value=m.value, class=m.class, ci=m.ci))
        end
        for m in _hk1_h4_derived_rows(enzyme, res.variant, res.r.fit.coords)
            push!(rows, (variant=res.variant, mode=res.mode, name=m.name,
                         value=m.value, class=m.class, ci=m.ci))
        end
    end
    CSV.write(joinpath(outdir, "macro_constants.csv"), DataFrame(rows))
    # goodness_of_fit.csv
    CSV.write(joinpath(outdir, "goodness_of_fit.csv"),
        DataFrame([(variant=res.variant, mode=res.mode, in_sample=res.r.fit.loss,
                    cv_mean=res.r.cv.mean_cv, cv_se=res.r.cv.se) for res in results]))
    # identifiable_functions.csv (eigen-spectrum per variant × mode; `index` = eigen rank)
    spec = NamedTuple[]
    for res in results, (k, λ) in enumerate(res.idf.eigvals)
        push!(spec, (variant=res.variant, mode=res.mode, index=k, eigval=λ,
                     stiff=(k <= res.idf.rank)))
    end
    CSV.write(joinpath(outdir, "identifiable_functions.csv"), DataFrame(spec))
    # micro_parameters.jl — the DEPLOY block. The fitted macro coords are mapped to a positive
    # micro-k representative by the CLOSED-FORM ChaDeploy.cha_deploy_micro (the inverse of the
    # Cha readoff — no LM/root-finding, so the de-conflation "blocker E" cannot recur). Emitted
    # for the deploy variant (the only variant the pipeline now runs), one `{variant}_{mode}_{sym}`
    # line per free param. Guarded so a deploy failure never aborts write_outputs.
    open(joinpath(outdir, "micro_parameters.jl"), "w") do io
        println(io, "# Consensus $(name) micro parameters (one point on a flat manifold;")
        println(io, "# only the macro constants labeled data_identified are determined).")
        println(io, "# DEPLOY block: closed-form ChaDeploy.cha_deploy_micro of the fitted Cha")
        println(io, "# macro coords (the promoted Cha law's micro representative). The pure-RE")
        println(io, "# variant is this law's release-rate->infinity limit (not emitted separately).")
        for res in results
            # Per-variant guardrail: an identifiability failure, not a parsimony trade-off --
            # dropping the anchor doesn't simplify THIS variant's mechanism, it leaves Ki_NADPH
            # railed at the same parameter count. Only fires for variants that still require the
            # anchor (_requires_reverse_anchor) -- never for the anchor-optional ablations.
            if !anchor_reverse && _requires_reverse_anchor(enzyme, res.variant)
                println(io, "# ##########################################################################")
                println(io, "# NOT DEPLOYABLE [$(res.variant) $(res.mode)] — Ki_NADPH is structurally")
                println(io, "# undetermined without this anchor (railed). Do NOT wire into the ODE.")
                println(io, "# ##########################################################################")
            end
            logθ = nothing
            try
                logθ = ChaDeploy.cha_deploy_micro(enzyme, res.r.mech, res.r.fit.coords; keq=keq,
                                                  koffQ=ChaFit.CHA_DEPLOY_RELEASE_RATE,
                                                  release_rate=ChaFit.CHA_DEPLOY_RELEASE_RATE)
            catch err
                println(io, "# --- $(res.variant) [$(res.mode)] : deploy unavailable ($(err)) ---")
            end
            if logθ !== nothing
                free = free_params(res.r.mech)
                println(io, "# --- $(res.variant) [$(res.mode)] (deploy: cha_deploy_micro) ---")
                for (s, lv) in zip(free, logθ)
                    println(io, "$(res.variant)_$(res.mode)_$(s) = $(10.0^lv)")
                end
            end
        end
    end
    # report.md (variant × mode tables + Mode1<->Mode2 agreement + enzyme-specific blocks)
    open(joinpath(outdir, "report.md"), "w") do io
        println(io, "# Consensus macro-constant extraction ($(name))\n")
        for res in results
            println(io, "## $(res.variant) — $(res.mode)\n")
            # Per-variant: only the variants that still require the anchor get flagged (the
            # anchor-optional ablations run anchor-off as a supported default, not a diagnostic).
            if !anchor_reverse && _requires_reverse_anchor(enzyme, res.variant)
                println(io, "> **⚠ NOT DEPLOYABLE (`anchor_reverse=false`).** `Km_NADPH_rev` was left ",
                            "unanchored, so the forward `Ki_NADPH` is non-identifiable (railed) here ",
                            "— an identifiability failure, not a parsimony trade-off (same parameter ",
                            "count, just an undetermined value). This variant's `micro_parameters.jl` ",
                            "block is **not deployable**. Compare fit quality only.\n")
            end
            println(io, "CV (leave-one-article-out): $(res.r.cv.mean_cv) ± $(res.r.cv.se)\n")
            println(io, "| macro | value | class | ci |\n|---|---|---|---|")
            for m in res.classed
                println(io, "| $(m.name) | $(m.value) | $(m.class) | $(m.ci) |")
            end
            for m in _apparent_km_rows(enzyme, res.r.fit.coords)
                println(io, "| $(m.name) | $(m.value) | $(m.class) | $(m.ci) |")
            end
            for m in _hk1_h4_derived_rows(enzyme, res.variant, res.r.fit.coords)
                println(io, "| $(m.name) | $(m.value) | $(m.class) | $(m.ci) |")
            end
            println(io)
        end
        # Forward-constant agreement across modes over the classed cha_coords. PGD has 3 modes
        # (mode1/2/3) → pairwise (1↔2, 1↔3, 2↔3); G6PD has 2 (mode1/2) → 1↔2. Only both-modes
        # `:data_identified` coords are compared; pinned/derived rows never enter mode_agreement
        # (class ≠ :data_identified), so an intended data-vs-physiology contrast is NOT mislabeled
        # a flatness WARNING.
        modes = unique(r.mode for r in results)
        mode_pairs = [(modes[i], modes[j]) for i in 1:length(modes) for j in (i+1):length(modes)]
        println(io, "## Forward-constant agreement across modes\n")
        for (ma, mb) in mode_pairs
            println(io, "### $(ma) ↔ $(mb)\n")
            for variant in unique(r.variant for r in results)
                ia = findfirst(r -> r.variant==variant && r.mode==ma, results)
                ib = findfirst(r -> r.variant==variant && r.mode==mb, results)
                (ia === nothing || ib === nothing) && continue
                ag = mode_agreement(results[ia].classed, results[ib].classed)
                println(io, "#### $(variant)\n")
                println(io, "| macro | $(ma) | $(mb) | Δdex | agree |\n|---|---|---|---|---|")
                for a in ag
                    println(io, "| $(a.name) | $(a.mode1) | $(a.mode2) | $(round(a.delta_dex, digits=3)) | $(a.agree) |")
                end
                any(!a.agree for a in ag) &&
                    println(io, "\n**WARNING:** a forward constant moved beyond tolerance between $(ma) and $(mb) — a mode pin was not flat.")
                println(io)
            end
        end
        # Enzyme-specific resolved-decision blocks.
        enzyme === :G6PD && _write_g6pd_koffq_block(io, results, d, keq; anchor_reverse=anchor_reverse)
        enzyme === :PGD  && _write_pgd_km_pga_block(io, results)
        print(io, _report_note(enzyme, results))
    end
end

# G6PD koffQ HYBRID block: deploy the promoted release at a healthy swept default and ALSO
# report the (weak, wide-CI) data-identified koffQ from the reverse-weighted diagnostic, with
# the deploy<->data gap and the honest caveat. Guarded so a report-only refit failure never
# aborts write_outputs.
function _write_g6pd_koffq_block(io, results, d, keq; anchor_reverse::Bool=true)
    println(io, "## G6PD koffQ hybrid (promoted NADPH-release fiber)\n")
    idx = findfirst(r -> r.variant === deploy_variant(:G6PD), results)
    res = results[idx === nothing ? 1 : idx]
    rvariant = res.variant
    try
        hr = ChaKoffqReport.koffq_hybrid_report(res.r.mech, d; keq=keq,
                                                variant=rvariant, koffQ_deploy=1.0e3,
                                                anchor_reverse=anchor_reverse)
        println(io, "| quantity | value |\n|---|---|")
        println(io, "| deploy koffQ (swept) | $(hr.deploy_value) |")
        println(io, "| data-identified koffQ | $(hr.data_identified_value) |")
        println(io, "| log10 CI | [$(hr.ci[1]), $(hr.ci[2])] |")
        println(io, "| gap (dex, data − deploy) | $(hr.gap_dex) |")
        println(io, "| reverse rows used | $(hr.n_reverse) |\n")
        println(io, "> $(hr.caveat)\n")
    catch err
        println(io, "> koffQ hybrid report unavailable ($(err))\n")
    end
end

# PGD Km_PGA gap warning: the apparent Km_PGA in Mode 1 (unanchored — the corpus pulls it high)
# vs the anchored Mode 2 readoff. Prints a LOUD warning quantifying the Mode-1-vs-anchored gap
# (spec §6). Guarded.
function _write_pgd_km_pga_block(io, results)
    println(io, "## PGD Km_PGA gap (unanchored corpus vs literature anchor)\n")
    try
        i1 = findfirst(r -> r.mode === :mode1, results)
        i2 = findfirst(r -> r.mode === :mode2, results)
        (i1 === nothing || i2 === nothing) && error("missing Mode 1 or Mode 2 cell")
        km1 = ChaFit.cha_apparent_km(:PGD, results[i1].r.fit.coords, :Km_PGA)
        km2 = ChaFit.cha_apparent_km(:PGD, results[i2].r.fit.coords, :Km_PGA)
        gap = log10(km1) - log10(km2)
        println(io, "| quantity | value |\n|---|---|")
        println(io, "| Km_PGA Mode 1 (unanchored) | $(km1) M |")
        println(io, "| Km_PGA Mode 2 (anchored 59µM) | $(km2) M |")
        println(io, "| gap (dex, Mode1 − Mode2) | $(gap) |\n")
        println(io, "> **WARNING — Km_PGA is NOT data-identified in the forward corpus.** The ",
                    "unanchored Mode-1 apparent Km_PGA (", round(km1*1e6; sigdigits=3), " µM) ",
                    "differs from the literature-anchored Mode-2 value (",
                    round(km2*1e6; sigdigits=3), " µM) by ", round(gap; sigdigits=3),
                    " dex. The corpus pulls Km_PGA high (it lacks sub-Km [6PG] coverage); the ",
                    "38–80 µM band is imposed via the anchor (Mode 2/3), not recovered from data.\n")
    catch err
        println(io, "> Km_PGA gap report unavailable ($(err))\n")
    end
end

# Enzyme-specific closing note for report.md. PGD/G6PD: the forward-Ki de-conflation caveat
# (the G6PD note reports the ACTUAL fitted Mode-1 cross-term Ki_NADPH — fit directly, not
# hardcoded — at full budget the corpus reads it ABOVE the 9–24 µM literature band). HK1: the
# H1/H3 candidate + pinned-feedback note. Any other enzyme returns "" (no note).
function _report_note(enzyme::Symbol, results)
    enzyme === :PGD && return string(
        "> **Forward `Ki_NADPH` de-conflation (PGD):** on the Cha law the forward product-\n",
        "> inhibition `Ki_NADPH` is read as the E·PGA dead-end cross term, decoupled from the\n",
        "> bare-[NADPH] reverse-release Km (`Km_NADPH_rev`). It is literature-pinned (Cottreau\n",
        "> 17 µM) in Mode 2/3; Mode 1 reports it diagnostic/unconstrained (the forward-only\n",
        "> cross-term de-conflation was refuted on the PGD corpus — pin-only 17 µM).\n",
        "> FROZEN: do not add SS steps to chase this (see test_rec4_topology_freeze.jl).\n")
    enzyme === :HK1 && return string(
        "> **HK1 consensus (H1 α=1 two G6P sites / H3 α=∞ single net site).** The two\n",
        "> data-identified substrate constants (`Km_Glc`, `Km_ATP`) are fit from the corpus; the\n",
        "> product-feedback constants (`Ki_G6P_C`, `Ki_G6P_N`, `K_Pi_N`, `Ki_ADP`) are literature-\n",
        "> pinned in mode2/mode3. The two G6P sites are symmetric quadratic roots (non-identifiable\n",
        "> together); mode1 (no pins) is expected to show that C/N ridge — the honest baseline.\n")
    enzyme === :G6PD && return string(
        "> **Forward `Ki_NADPH` de-conflation (G6PD):** the Cha law reads the forward product-\n",
        "> inhibition `Ki_NADPH` from the E·G6P·NADPH dead-end cross term, with `Km_NADPH_rev`\n",
        "> anchored (3.9 µM) to break the trade-off with the bare-[NADPH] productive-release\n",
        "> reverse channel. Cross-term data-id lands Ki_NADPH ", _g6pd_ki_nadph_desc(results),
        " in Mode 1 (the full-budget corpus reads ABOVE the 9–24 µM literature band); Mode 2\n",
        "> pins it to the literature 15 µM.\n")
    return ""
end

# =========================================================================================
#                    Exported runners: run_g6pd / run_pgd / run_hk1
# =========================================================================================
#
# Thin wrappers around `run_all` that pick the budget (smoke vs full), the default outdir,
# and spin up workers via `setup_workers`. These are the library replacement for the old
# per-enzyme launcher scripts (run_g6pd.jl / run_pgd.jl / run_hk1.jl).

function _budget(smoke::Bool)
    smoke ? (n_restarts=2, maxiter=150, maxtime=120.0) :
            (n_restarts=48, maxiter=1000, maxtime=300.0)
end

function _default_outdir(enzyme::AbstractString, smoke::Bool)
    joinpath(pwd(), "results", "$(enzyme)_" * Dates.format(Dates.today(), "yyyy-mm-dd") * (smoke ? "_smoke" : ""))
end

function _run_enzyme(cfg, enzyme::AbstractString; outdir=nothing, smoke::Bool=false, nprocs=nothing,
                     variants=nothing, row_filter=nothing, anchor_reverse::Bool=true)
    b = _budget(smoke)
    od = isnothing(outdir) ? _default_outdir(enzyme, smoke) : outdir
    setup_workers(nprocs)
    @info "FitRateEquation run starting" enzyme nworkers=nworkers() smoke outdir=od anchor_reverse
    # `variants`/`row_filter` are forwarded to `run_all` only when the caller supplies them
    # (e.g. `run_g6pd_noatp`'s :no_atp variant + ATP-row filter), so the plain per-enzyme
    # runners (run_g6pd/run_pgd/run_hk1) keep run_all's own defaults untouched.
    extra = NamedTuple()
    variants   === nothing || (extra = merge(extra, (variants=variants,)))
    row_filter === nothing || (extra = merge(extra, (row_filter=row_filter,)))
    run_all(cfg; outdir=od, n_restarts=b.n_restarts, maxiter=b.maxiter, maxtime=b.maxtime,
            anchor_reverse=anchor_reverse, extra...)
end

"""
    run_g6pd(; outdir=nothing, smoke=false, nprocs=nothing, anchor_reverse=true)

Run the deploy-variant × mode consensus macro-constant extraction for G6PD end-to-end and
write the six artifacts (macro_constants.csv, goodness_of_fit.csv,
identifiable_functions.csv, micro_parameters.jl, report.md, provenance.toml) to `outdir`
(default: `./results/G6PD_<date>[_smoke]`). `smoke=true` uses a tiny fit budget for a fast
sanity check. `nprocs` overrides the local worker-count default (see `setup_workers`); a
SLURM allocation always overrides `nprocs`. Returns the `run_all` results.

`anchor_reverse` (default `true`) controls the G6PD reverse-channel anchor. **The deployed
law REQUIRES `anchor_reverse=true`** — it anchors `Km_NADPH_rev` (3.9 µM) in every mode to
de-conflate the forward `Ki_NADPH` from the reverse-release Km. `anchor_reverse=false` leaves
`Km_NADPH_rev` free, deliberately reintroducing that conflation (forward `Ki_NADPH` becomes
non-identifiable). It is a **conflation/identifiability DIAGNOSTIC only**: the run is tagged
`NOT DEPLOYABLE` in `micro_parameters.jl` and `report.md`, and the anchor state is recorded in
`provenance.toml`. Use it with `variants=[:RE_rate_eq]` via `run_all` to reproduce the
original full-RE conflating fit.
"""
run_g6pd(; outdir=nothing, smoke=false, nprocs=nothing, anchor_reverse=true) =
    _run_enzyme(g6pd_config(), "G6PD"; outdir, smoke, nprocs, anchor_reverse)

"""
    run_g6pd_noatp(; outdir=nothing, smoke=false, nprocs=nothing, data_csv=nothing)

The ATP-free (`:no_atp`) G6PD variant: fits `run_all` with `variants=[:no_atp]` and
`row_filter=drop_atp_rows`, so ATP-bearing rows (ATP > 0) are dropped from the corpus before
the ATP-blind `:no_atp` mechanism is fit (the ATP-tolerant deploy variant is `run_g6pd`; this
is the library replacement for the standalone `run_g6pd_noatp.jl` launcher). `data_csv`, if
given, overrides the bundled G6PD corpus (see `g6pd_config`); otherwise the default corpus is
used. Default outdir is labeled `G6PD_noatp_<date>[_smoke]` (distinct from plain `run_g6pd`'s
`G6PD_<date>[_smoke]`) so the two runs never collide in `./results/`.
"""
function run_g6pd_noatp(; outdir=nothing, smoke=false, nprocs=nothing, data_csv=nothing)
    cfg = data_csv === nothing ? g6pd_config() : g6pd_config(; data_csv=data_csv)
    _run_enzyme(cfg, "G6PD_noatp"; outdir, smoke, nprocs, variants=[:no_atp], row_filter=drop_atp_rows)
end

"""
    run_pgd(; outdir=nothing, smoke=false, nprocs=nothing)

As `run_g6pd`, for PGD.
"""
run_pgd(;  outdir=nothing, smoke=false, nprocs=nothing) = _run_enzyme(pgd_config(),  "PGD";  outdir, smoke, nprocs)

"""
    run_hk1(; outdir=nothing, smoke=false, nprocs=nothing)

As `run_g6pd`, for HK1. Errors clearly if HK1 wiring is unavailable on this EnzymeRates
build (`FitRateEquation.HK1_AVAILABLE == false`; a deferred port — see AGENTS.md) rather
than crashing deeper in the pipeline.
"""
function run_hk1(; outdir=nothing, smoke=false, nprocs=nothing)
    HK1_AVAILABLE || error("HK1 is not available on this EnzymeRates build (deferred port). See AGENTS.md.")
    _run_enzyme(hk1_config(), "HK1"; outdir, smoke, nprocs)
end

# Format the actual Mode-1 G6PD cross-term Ki_NADPH (value + class + CI) for the report note.
function _g6pd_ki_nadph_desc(results)
    i1 = findfirst(r -> r.mode === :mode1, results)
    i1 === nothing && return "in the forward (tens-of-µM) band"
    j = findfirst(x -> x.name === :Ki_NADPH, results[i1].classed)
    j === nothing && return "in the forward (tens-of-µM) band"
    c = results[i1].classed[j]
    c.class === :data_identified && isfinite(c.ci) ?
        string("~", round(c.value * 1e6; sigdigits = 3), " µM (`data_identified`, ±",
               round(c.ci * 1e6; sigdigits = 2), " µM)") :
        string("~", round(c.value * 1e6; sigdigits = 3), " µM (`", c.class, "`)")
end

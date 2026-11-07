# =========================================================================================
#         G6PD koffQ report: the deploy convention + an all-constants-refit profile
# =========================================================================================
#
# REPORT-ONLY. This module refits the corpus purely to REPORT the koffQ situation; it does NOT
# alter any deployed law or fitted coords.
#
#   - The fit and the deployed law share one koffQ, CHA_DEPLOY_RELEASE_RATE (kf = 1 gauge
#     units). That value is a CONVENTION (release fast relative to catalysis, C = 1 + kf/koffQ
#     ~ 1), not a measurement.
#
#   - koffQ is free only if the other constants move with it. In the forward law (P = 0) the
#     shape depends on alpha/C, so alpha moves with C. In the reverse rows koffQ enters through
#     the kr*r/koffQ denominator term, which Kd_6PGLn can absorb when the rows sit at a single
#     [NADPH], as the 22 Beutler 1986 Fig 2a delta-lactone rows do (all at 50 uM NADPH). So the
#     profile below refits EVERY free constant at each koffQ. Scanning koffQ with the constants
#     held fixed does not measure koffQ: it just returns whatever rate those constants were fit
#     at. (The pre-0.10.0 report did that and called the result "data-identified koffQ = 1",
#     which was the old fit rate.)
#
#   - The PGLn-bearing rows are up-weighted so any reverse-data handle has its best chance to
#     show. On the current corpus the profile is flat: the data do not choose a koffQ.
# =========================================================================================

module ChaKoffqReport

using ..ChaFit
using ..ChaLaws
using ..FitRateEquation
using CMAEvolutionStrategy: minimize, xbest, fbest

export koffq_hybrid_report

# Number of extra copies of each PGLn-bearing row in the profile's dataset. K is a DIAGNOSTIC
# weight, not a fit hyperparameter: it does not touch the deployed law or the fitted coords.
const REVERSE_UPWEIGHT_K = 20

# The koffQ grid (log10, kf = 1 gauge units) and the relative loss rise that counts as flat.
const PROFILE_LOG10_KOFFQ = (-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0)
const PROFILE_FLAT_TOL = 0.01

# Build a reverse-up-weighted copy of `d`: every PGLn>0 row is duplicated K extra times. Rows
# keep their group key, so the duplicates stay inside their own (Article,Fig) groups (the
# centering is per group; replicating rows within a group scales its SSE contribution by ~K).
function _reverse_upweighted(d::Dataset, rev_idx::Vector{Int}, K::Int)
    concs = copy(d.concs); rate = copy(d.rate); group = copy(d.group); keq = copy(d.keq)
    for _ in 1:K, i in rev_idx
        push!(concs, d.concs[i]); push!(rate, d.rate[i])
        push!(group, d.group[i]); push!(keq, d.keq[i])
    end
    Dataset(concs, rate, group, keq)
end

# PGLn>0 row indices (the reverse rows). Guards for the field's presence.
_pgln_rows(d::Dataset) =
    [i for i in 1:nrows(d)
     if hasproperty(d.concs[i], :PGLn) && d.concs[i].PGLn > 0]

# Refit every free (non-pinned) coord at one koffQ, warm-started from the deploy-rate fit with
# alpha moved to hold alpha/C. Same per-figure keq as the main fit.
function _refit_at(mech, d_w::Dataset, base_coords::AbstractDict, pins::AbstractDict, koffQ::Real;
                   variant::Symbol)
    syms = cha_coords(:G6PD, variant)
    lo, hi = cha_coord_bounds(:G6PD, variant)
    C_ratio = (1 + 1 / koffQ) / (1 + 1 / CHA_DEPLOY_RELEASE_RATE)
    u0 = [log10(s === :alpha ? base_coords[s] * C_ratio : base_coords[s]) for s in syms]
    pinned(u) = [haskey(pins, s) ? pins[s] : x for (s, x) in zip(syms, u)]
    objective = ChaFit._cha_objective(:G6PD, d_w, syms, pins; variant=variant, release_rate=koffQ)
    o = minimize(objective, clamp.(u0, lo, hi), 0.1; lower=lo, upper=hi, maxiter=300, maxtime=30.0,
                 seed=UInt64(1), verbosity=0)
    coords = Dict(syms .=> 10 .^ pinned(collect(xbest(o))))
    (koffQ=float(koffQ), loss=fbest(o), alpha=coords[:alpha], Kd_6PGLn=coords[:Kd_6PGLn])
end

# -----------------------------------------------------------------------------------------
#   koffq_hybrid_report(mech, d; variant, anchor_reverse) -> NamedTuple
#
#   deploy_value : the koffQ the law is fit and deployed at (CHA_DEPLOY_RELEASE_RATE).
#   profile      : per koffQ, (koffQ, loss, alpha, Kd_6PGLn) with every free coord refit on the
#                  reverse-up-weighted corpus.
#   flat_range   : (lo, hi) koffQ values whose loss is within PROFILE_FLAT_TOL of the best.
#   spread       : (max loss - min loss) / min loss across the profile.
#   caveat       : plain-language reading of the profile (with the actual numbers interpolated).
#   n_reverse    : count of PGLn>0 rows up-weighted.
# -----------------------------------------------------------------------------------------
function koffq_hybrid_report(mech, d::Dataset;
                             variant::Symbol = :SS_NADPH_release_rate_eq,
                             anchor_reverse::Bool = false)
    rev_idx = _pgln_rows(d)
    n_reverse = length(rev_idx)
    pins = resolve_cha_pins(:G6PD, variant, :mode1; anchor_reverse=anchor_reverse)
    base = cha_fit_candidate(:G6PD, mech, d; n_restarts=6, maxiter=300, maxtime=60.0, seed=1,
                             pins=pins, variant=variant)
    d_w = _reverse_upweighted(d, rev_idx, REVERSE_UPWEIGHT_K)
    profile = [_refit_at(mech, d_w, base.coords, pins, 10.0^lk; variant=variant)
               for lk in PROFILE_LOG10_KOFFQ]

    losses = [p.loss for p in profile]
    best = minimum(losses)
    spread = (maximum(losses) - best) / best
    flat = [p.koffQ for p in profile if p.loss <= best * (1 + PROFILE_FLAT_TOL)]
    flat_range = (minimum(flat), maximum(flat))

    r3(x) = round(x; sigdigits=3)
    convention = string("koffQ = ", r3(CHA_DEPLOY_RELEASE_RATE), " is a convention (release fast ",
        "relative to catalysis, C = 1 + kf/koffQ ~ ", round(1 + 1 / CHA_DEPLOY_RELEASE_RATE; sigdigits=4),
        " on the kf = 1 gauge), not a measurement. ")
    caveat = spread <= PROFILE_FLAT_TOL ?
        string(convention, "The corpus does not identify it: with every free constant refit at ",
            "each koffQ from ", r3(first(profile).koffQ), " to ", r3(last(profile).koffQ), " and the ",
            n_reverse, " PGLn-bearing rows up-weighted, the loss moves by only ", r3(100 * spread),
            "%. alpha and Kd_6PGLn move with koffQ, so they are reported at this convention.") :
        string(convention, "The all-constants-refit profile is not flat (loss spread ",
            r3(100 * spread), "% with the ", n_reverse, " PGLn-bearing rows up-weighted); koffQ values ",
            "within ", r3(100 * PROFILE_FLAT_TOL), "% of the best loss span ", r3(flat_range[1]), " to ",
            r3(flat_range[2]), ". Read the table before relying on the convention.")

    (deploy_value = float(CHA_DEPLOY_RELEASE_RATE),
     profile = profile,
     flat_range = flat_range,
     spread = spread,
     caveat = caveat,
     n_reverse = n_reverse)
end

end # module ChaKoffqReport

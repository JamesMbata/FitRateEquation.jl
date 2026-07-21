# =========================================================================================
#         G6PD koffQ HYBRID report: swept-deploy default + reverse-weighted diagnostic
# =========================================================================================
#
# REPORT-ONLY. This module reads/refits the corpus purely to REPORT the koffQ situation; it
# does NOT alter any deployed law or fitted forward coords. The resolved decision is the
# HYBRID treatment of G6PD's promoted fiber `koffQ` (the catalytic NADPH-release rate):
#
#   - DEPLOY koffQ at a healthy SWEPT default that is flux-neutral under the forward+product-
#     inhibition deploy gate. On the centered log-ratio loss the overall Vmax is gauged out
#     and the observable law depends on koffQ only through the Cha C-factor C = 1 + kf/koffQ
#     (cha_invert.jl). With the fit gauge kf = 1, a release rate koffQ >> kf drives C -> 1, so
#     release is NOT rate-limiting and the forward shape is koffQ-invariant on the fiber
#     (cha_micro_from_macro_G6PD co-adjusts alpha so apparent Km's stay fixed). The default
#     `koffQ_deploy = 1.0e3` (in the kf=1 gauge units) gives C = 1.001 -- release is ~1000x
#     faster than catalysis, so deploy is effectively at the fast-release limit and flux-
#     neutral. Any koffQ >= ~1e2 is equally valid; 1e3 is a round, safely-fast default.
#
#   - ADDITIONALLY run a reverse-weighted DIAGNOSTIC refit and REPORT the data-identified
#     koffQ with its (wide) CI, flagging the deploy<->data gap. koffQ is silent to the forward
#     gate (the gate's rows are P=0 / forward + product-inhibition; koffQ only enters the
#     bare-[NADPH] productive-release reverse channel), but is WEAKLY data-identified by the
#     18 PGLn-bearing reverse rows the gate discards (Beutler1986|3, Gordon1994|1A,
#     Gordon1994|1B). The refit up-weights those rows; the handle is single-regime, single-pH
#     and thin, so the CI is WIDE on purpose -- that width IS the finding.
# =========================================================================================

module ChaKoffqReport

using ..ChaFit
using ..ChaLaws
using ..FitRateEquation
using Statistics: median

export koffq_hybrid_report

# Number of times the PGLn-bearing groups' rows are replicated in the reverse-weighted
# dataset. Replicating the 18 reverse rows K=20 times makes those 3 groups dominate the
# centered log-ratio sum (the forward groups still mean-center their own residuals to ~0 at
# the data-identified forward coords, so the koffQ-sensitive reverse SSE is what the 1-D scan
# actually minimizes). K is a DIAGNOSTIC weight, not a fit hyperparameter -- it does not touch
# the deployed law or the forward coords (those come from the unweighted Mode-1 fit).
const REVERSE_UPWEIGHT_K = 20

# Build a reverse-up-weighted copy of `d`: every PGLn>0 row is duplicated K extra times so the
# PGLn-bearing groups dominate the per-group centered loss. Forward rows are kept once. We
# duplicate at the ROW level but preserve each row's group key, so the duplicated reverse rows
# stay inside their own (Article,Fig) groups (the centering is per group; replicating rows
# within a group scales that group's SSE contribution by ~K).
function _reverse_upweighted(d::Dataset, rev_idx::Vector{Int}, K::Int)
    concs = copy(d.concs); rate = copy(d.rate); group = copy(d.group); keq = copy(d.keq)
    for _ in 1:K, i in rev_idx
        push!(concs, d.concs[i]); push!(rate, d.rate[i])
        push!(group, d.group[i]); push!(keq, d.keq[i])
    end
    Dataset(concs, rate, group, keq)
end

# PGLn>0 row indices (the koffQ handle). Guards for the field's presence.
_pgln_rows(d::Dataset) =
    [i for i in 1:nrows(d)
     if hasproperty(d.concs[i], :PGLn) && d.concs[i].PGLn > 0]

# Reverse-weighted 1-D loss as a function of lk = log10(koffQ). Evaluates the Cha centered
# log-ratio loss on the reverse-up-weighted dataset, holding the data-identified forward
# `coords` fixed and sweeping ONLY the promoted release rate. The release equilibrium
# `release_eq` is held at coords[:Km_NADPH_rev] (the NADPH-release RE equilibrium, a free
# coord on the G6PD fiber); kr defaults to the Haldane value inside cha_macro_tuple.
function _reverse_loss(mech, d_w::Dataset, coords::AbstractDict, keq::Real, lk::Real)
    cha_centered_logratio_loss(:G6PD, mech, d_w, coords; keq=keq,
                               release_rate=10.0^lk,
                               release_eq=coords[:Km_NADPH_rev])
end

# Hand-rolled 1-D minimizer over lk in [lo, hi]: coarse grid then parabolic refine around the
# best grid node (no heavy deps). Returns (lk_hat, f_hat, grid_lk, grid_f).
function _minimize_1d(f, lo::Real, hi::Real; ngrid::Int=71)
    grid = collect(range(lo, hi; length=ngrid))
    fv = [f(x) for x in grid]
    j = argmin(fv)
    lk_hat = grid[j]; f_hat = fv[j]
    # Parabolic refinement using the two neighbors (only when an interior bracket exists).
    if 1 < j < length(grid)
        x0, x1, x2 = grid[j-1], grid[j], grid[j+1]
        y0, y1, y2 = fv[j-1], fv[j], fv[j+1]
        denom = (y0 - 2y1 + y2)
        if denom > 0                                   # convex bracket -> vertex
            xv = x1 - 0.5 * (x2 - x0) * (y2 - y0) / (2 * denom)
            if lo <= xv <= hi
                fxv = f(xv)
                if fxv < f_hat
                    lk_hat = xv; f_hat = fxv
                end
            end
        end
    end
    (lk_hat, f_hat, grid, fv)
end

# Profile-likelihood-style CI on lk = log10(koffQ): the interval where the reverse-weighted
# loss rises by `rise` above the minimum. We walk the grid out from the minimizer to either
# side until the loss exceeds f_hat + rise (linear-interpolating the crossing), railing the CI
# bound to the scan edge if the loss never rises enough on that side (an HONEST "unbounded on
# this side" outcome -- the handle is too weak to close the interval). The threshold `rise`
# uses a chi-square-1 deviance scale on the per-row centered loss; since the loss here is an
# averaged SSE (not a true 2*negloglik), `rise` is a small fixed deviance proxy -- documented
# as a profile WIDTH, not a calibrated likelihood interval.
function _profile_ci(grid::Vector{Float64}, fv::Vector{Float64}, lk_hat::Real, f_hat::Real;
                     rise::Real)
    thr = f_hat + rise
    n = length(grid)
    jstar = argmin(abs.(grid .- lk_hat))
    # Walk left.
    lo = grid[1]
    for j in jstar:-1:2
        if fv[j-1] >= thr
            # interpolate crossing between grid[j-1] (>=thr) and grid[j] (<thr)
            t = (thr - fv[j]) / (fv[j-1] - fv[j])
            lo = grid[j] + t * (grid[j-1] - grid[j])
            break
        end
        j == 2 && (lo = grid[1])                       # never crossed -> rail to edge
    end
    # Walk right.
    hi = grid[end]
    for j in jstar:(n-1)
        if fv[j+1] >= thr
            t = (thr - fv[j]) / (fv[j+1] - fv[j])
            hi = grid[j] + t * (grid[j+1] - grid[j])
            break
        end
        j == n-1 && (hi = grid[end])                   # never crossed -> rail to edge
    end
    (lo, hi)
end

# -----------------------------------------------------------------------------------------
#   koffq_hybrid_report(mech, d; keq, variant, koffQ_deploy) -> NamedTuple
#
#   deploy_value          : the healthy SWEPT koffQ used in deploy (= koffQ_deploy).
#   data_identified_value : koffQ from the reverse-weighted diagnostic refit.
#   ci                    : (lo, hi) CI on log10(data_identified_value) from the 1-D profile.
#   gap_dex               : log10(data_identified_value) - log10(deploy_value).
#   caveat                : honest framing string (with the actual numbers interpolated).
#   n_reverse             : count of PGLn>0 rows used.
# -----------------------------------------------------------------------------------------
function koffq_hybrid_report(mech, d::Dataset;
                             keq::Real = median(d.keq),
                             variant::Symbol = :SS_NADPH_release_rate_eq,
                             koffQ_deploy::Real = 1.0e3,
                             anchor_reverse::Bool = true)
    rev_idx = _pgln_rows(d)
    n_reverse = length(rev_idx)

    # (1) Standard Mode-1 G6PD fit -> data-identified FORWARD coords. The forward coords are
    #     fit on the FULL unweighted corpus (the deploy fit); koffQ does not enter them (it is
    #     gauged into C, and the fit uses the default healthy release_rate).
    base = cha_fit_candidate(:G6PD, mech, d; n_restarts=6, maxiter=300, maxtime=60.0, seed=1,
                             keq=keq,
                             pins=resolve_cha_pins(:G6PD, variant, :mode1; anchor_reverse=anchor_reverse))

    # (2) Reverse-weighted diagnostic: up-weight the PGLn>0 groups, then 1-D minimize the
    #     reverse-weighted loss over lk = log10(koffQ) holding the forward coords fixed.
    d_w = _reverse_upweighted(d, rev_idx, REVERSE_UPWEIGHT_K)
    f1d = lk -> _reverse_loss(mech, d_w, base.coords, keq, lk)
    lo_scan, hi_scan = -1.0, 6.0
    lk_hat, f_hat, grid, fv = _minimize_1d(f1d, lo_scan, hi_scan; ngrid=71)
    data_identified_value = 10.0^lk_hat

    # (3) CI from the 1-D profile: the lk interval where the reverse-weighted loss rises by a
    #     small deviance proxy above the minimum (profile-likelihood-style). Wide by design.
    ci = _profile_ci(grid, fv, lk_hat, f_hat; rise=0.10 * abs(f_hat) + 1e-3)

    gap_dex = log10(data_identified_value) - log10(koffQ_deploy)

    # Honest framing carried verbatim (numbers interpolated). "forward gate" appears.
    caveat = string(
        "koffQ is effectively free under the forward gate (deploy swept to ",
        round(koffQ_deploy; sigdigits=3), " in kf=1 gauge units, C=1+kf/koffQ~",
        round(1 + 1 / koffQ_deploy; sigdigits=4), ", flux-neutral); a reverse-weighted refit ",
        "exposes a weak (~0.4 log-unit within-group), single-regime, single-pH handle from the ",
        n_reverse, " reverse (PGLn-bearing) rows -- data-identified koffQ ~ ",
        round(data_identified_value; sigdigits=3), " (log10 CI [",
        round(ci[1]; sigdigits=3), ", ", round(ci[2]; sigdigits=3), "], gap ",
        round(gap_dex; sigdigits=3), " dex vs deploy) -- too thin to overturn the silent-",
        "variant null-space.")

    (deploy_value = float(koffQ_deploy),
     data_identified_value = data_identified_value,
     ci = ci,
     gap_dex = gap_dex,
     caveat = caveat,
     n_reverse = n_reverse)
end

end # module ChaKoffqReport

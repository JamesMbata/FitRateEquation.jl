# =========================================================================================
#         Direct macro-constant fit coordinates for the Cha-form rate laws
# =========================================================================================
#
# The new Cha fit optimizes the NAMED forward shape constants DIRECTLY through `cha_rate_*`
# (retiring the coeff-space coordinate set). Under the per-(Article,Fig) mean-centered
# log-ratio loss the overall Vmax is GAUGED OUT, so the kinetic-scale parameters are FIXED,
# not fit:
#   - kf = 1.0, Et = 1.0          (kcat=1 / unit-enzyme gauge)
#   - the promoted SS-release rate (`release_rate`: koffQ for G6PD, the Ru5P-release koff for
#     PGD) is a SWEPT FIBER held at a healthy default -- NOT a fit coord,
#   - the promoted-release on-rate (konQ / kon) is DERIVED from `release_rate / release_eq`,
#   - kr (reverse catalysis) is HALDANE-determined from `keq` -- NOT a free fit coord.
#
# So the FREE fit coordinates (`cha_coords`) are exactly the data-identifiable forward shape
# constants. `cha_macro_tuple` assembles the FULL named tuple `cha_rate_*` consumes from
# (free coords + gauge/fiber/Haldane inputs); `cha_haldane_kr` returns the kr making the
# law's apparent Keq equal `keq`.
#
# RELEASE-EQUILIBRIUM (`release_eq`) -- the promoted step's koff/kon -- is enzyme-specific:
#   - G6PD: the single promoted release IS the NADPH release, so the release equilibrium
#     EQUALS Km_NADPH_rev (a free coord). `release_eq` defaults to coords[:Km_NADPH_rev]
#     and konQ = release_rate / Km_NADPH_rev. One constant serves both Haldane and konQ.
#   - PGD: the promoted release is Ru5P, whose equilibrium KdRu = koff/kon is DISTINCT from
#     Km_NADPH_rev (the NADPH-release RE equilibrium, K8). KdRu is a NUISANCE DOF:
#     unidentifiable on forward+product-inhibition data (the parity gate), reported
#     unconstrained, never pinned (mirrors the enzymes/pgd.jl reverse-Km Dalziel framing).
#     It is therefore a FIXED healthy fiber `release_eq` (default `CHA_KDRU_DEFAULT`), and
#     Km_NADPH_rev (a free coord) enters the Haldane formula AND the bare-[NADPH] reverse
#     term INDEPENDENTLY of KdRu.
#
# Haldane closed forms (verified vs the readoff to rtol 1e-8; derivations/cha_derive_*.py
# Property 4, on the physical fiber kon = release_rate / release_eq):
#   G6PD: kr = Kd_6PGLn * Km_NADPH_rev * kf / (Kd_NADP * Kd_G6P * alpha * keq)
#   PGD:  kr = Kd_CO2  * Km_NADPH_rev * KdRu * kf / (Kd_NADP * Kd_PGA * alpha * keq)
# (PGD: BOTH Km_NADPH_rev [coord] and KdRu [release_eq] appear.)
# =========================================================================================

module ChaFit

export cha_coords, cha_macro_tuple, cha_haldane_kr, CHA_KDRU_DEFAULT
export cha_centered_logratio_loss
export cha_fit_candidate, cha_coord_bounds
export cha_apparent_km, cha_specificity, CHA_DEPLOY_RELEASE_RATE
export resolve_cha_pins

# Parent-relative imports so ChaFit composes BOTH as a Main-level module (tests, where the
# parents are Main.ChaLaws / Main.FitRateEquation) AND as a FitRateEquation submodule
# (pipeline). The test header includes cha_laws.jl before cha_fit.jl, so `..ChaLaws` exists.
using ..ChaLaws
using ..ChaLawsHK1
using ..FitRateEquation

# CMA-ES multi-start machinery (same package symbols coeff_fit.jl uses). ChaFit is its own
# module so these are imported here directly rather than inherited from FitRateEquation.
using CMAEvolutionStrategy: minimize, xbest, fbest
using Random, LinearAlgebra
using Statistics: median

# Healthy Ru5P-release equilibrium nuisance default (KdRu = koff/kon) for PGD. Ru5P release is
# the Topham 1986 rate-limiting promoted step; on the forward+product-inhibition parity gate
# the reverse arm is ~0, so KdRu is inert there and the exact value affects only the
# diagnostic both-product shape. Chosen as a neutral mid-band dimensionless release
# equilibrium (observed readoff KdRu spans ~1e-2..1e0 in these rate-constant-ratio units).
const CHA_KDRU_DEFAULT = 1.0

# HK1 candidate alpha (the N/C-half G6P negative-cooperativity axis) is a FIXED structural
# value, NOT a fit coord: H1 = two independent G6P sites (alpha=1), H3 = mutual exclusion
# (alpha=Inf, drops the [G6P]^2 coupling). Threaded by variant through cha_macro_tuple.
_hk1_variant_alpha(variant::Symbol) =
    variant === :H1 ? 1.0 :
    variant === :H4 ? 1.0 :    # H4 = H1 (alpha=1) reparameterized in {Keff, split_ratio}
    variant === :H3 ? Inf  :
    error("_hk1_variant_alpha: HK1 variant must be :H1, :H4, or :H3 (got $variant)")

# -----------------------------------------------------------------------------------------
#   Free fit coordinates: the data-identifiable forward shape constants. EXCLUDES the gauge
#   (kf, Et), the swept fiber (release_rate = koffQ / koff and its on-rate konQ / kon), and
#   the Haldane-determined reverse catalysis kr.
# -----------------------------------------------------------------------------------------
function cha_coords(enzyme::Symbol, variant::Symbol=:_deploy)
    if enzyme === :G6PD
        return variant === :no_atp ?
            [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Km_NADPH_rev] :
        variant === :no_g6p_nadph_deadend ?
            [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_ATP, :Ki_ATP_EG, :Km_NADPH_rev] :
        variant === :no_g6p_atp_deadend ?
            [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Ki_ATP, :Km_NADPH_rev] :
        variant === :no_g6p_both_deadends ?
            [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_ATP, :Km_NADPH_rev] :
            [:Kd_NADP, :Kd_G6P, :Kd_6PGLn, :alpha, :Ki_NADPH, :Ki_ATP, :Ki_ATP_EG,
             :Km_NADPH_rev]
    elseif enzyme === :PGD
        return variant === :full_re ?
            # Fully-RE (fiber-free): NO promoted-release fiber DOF, NO separate forward Ki_NADPH.
            # Kd_NADPH is the single competitive NADPH constant (Km_NADPH_rev ≡ Kd_NADPH); Kd_Ru5P
            # is a real RE coord. Effector coords (:Ki_ATP/:Ki_ATP_EN/:Ki_NADPH) are appended by
            # the config only when the dead-ends are enabled (default OFF).
            [:Kd_NADP, :Kd_PGA, :alpha, :Kd_NADPH, :Kd_Ru5P, :Kd_CO2] :
            [:Kd_NADP, :Kd_PGA, :alpha, :Kd_CO2, :Ki_NADPH, :Ki_ATP, :Ki_ATP_EN,
             :Km_NADPH_rev]
    elseif enzyme === :HK1
        # H4 reparameterizes the two G6P dissociation constants {Ki_G6P_C, Ki_G6P_N} into the
        # data-identifiable pair {Keff, split_ratio}, where Keff = 1/(1/Kc+1/Kn) (effective G6P
        # feedback, forward-identified) and split_ratio = √(Kc·Kn)/Keff (the reverse-turnover-
        # identified split). The back-map (cha_macro_tuple, variant=:H4) reconstructs Kc,Kn.
        # H1/H3 keep the raw {Ki_G6P_C, Ki_G6P_N} coords.
        # See notes/2026-06-13_hk1_g6p_ridge_resolution_report.md.
        return variant === :H4 ?
            [:Kd_Glc, :Kd_ATP, :Keff, :Ki_ADP, :split_ratio, :K_Pi_N] :
            [:Kd_Glc, :Kd_ATP, :Ki_G6P_C, :Ki_ADP, :Ki_G6P_N, :K_Pi_N]
    else
        error("cha_coords: unknown enzyme $enzyme (expected :G6PD, :PGD, or :HK1)")
    end
end

# -----------------------------------------------------------------------------------------
#   Haldane kr: the reverse catalysis rate making the law's apparent Keq equal `keq`.
#   `release_eq` is the promoted-step release equilibrium (koff/kon):
#     G6PD -> Km_NADPH_rev (single NADPH-release promoted step),
#     PGD  -> KdRu (Ru5P-release equilibrium, distinct from Km_NADPH_rev).
# -----------------------------------------------------------------------------------------
function cha_haldane_kr(enzyme::Symbol, coords::AbstractDict; keq::Real, release_rate::Real,
                        kf::Real, release_eq::Real = _default_release_eq(enzyme, coords),
                        variant::Symbol = :_deploy)
    # Fully-RE PGD: no promoted-release fiber. Haldane comes straight from the RE dissociation
    # constants: kr = kf·Kd_NADPH·Kd_Ru5P·Kd_CO2 / (Keq·α·Kd_NADP·Kd_PGA). release_rate/release_eq
    # are inert here (fiber-free); they are accepted for signature parity with the cha_base call.
    if enzyme === :PGD && variant === :full_re
        return coords[:Kd_NADPH] * coords[:Kd_Ru5P] * coords[:Kd_CO2] * kf /
               (coords[:Kd_NADP] * coords[:Kd_PGA] * coords[:alpha] * keq)
    end
    if enzyme === :G6PD
        # Keq_app = Kd_6PGLn * Km_NADPH_rev * kf / (Kd_NADP * Kd_G6P * alpha * kr).
        # release_eq == Km_NADPH_rev on the G6PD fiber.
        return coords[:Kd_6PGLn] * release_eq * kf /
               (coords[:Kd_NADP] * coords[:Kd_G6P] * coords[:alpha] * keq)
    elseif enzyme === :PGD
        # Keq_app = Kd_CO2 * Km_NADPH_rev * KdRu * kf / (Kd_NADP * Kd_PGA * alpha * kr).
        # Km_NADPH_rev (coord) and KdRu (release_eq) enter independently.
        return coords[:Kd_CO2] * coords[:Km_NADPH_rev] * release_eq * kf /
               (coords[:Kd_NADP] * coords[:Kd_PGA] * coords[:alpha] * keq)
    else
        error("cha_haldane_kr: unknown enzyme $enzyme (expected :G6PD or :PGD)")
    end
end

# Default release equilibrium per enzyme (see module header).
function _default_release_eq(enzyme::Symbol, coords::AbstractDict)
    enzyme === :G6PD && return coords[:Km_NADPH_rev]
    enzyme === :PGD  && return CHA_KDRU_DEFAULT
    enzyme === :HK1  && return 1.0
    error("_default_release_eq: unknown enzyme $enzyme")
end

# -----------------------------------------------------------------------------------------
#   Assemble the FULL named tuple `cha_rate_*` consumes from (free coords + gauge/fiber/
#   Haldane inputs). Defaults give the FIT path (kf=Et=1, healthy release_rate/release_eq,
#   Haldane kr); the round-trip test overrides them with the readoff's actual values.
# -----------------------------------------------------------------------------------------
function cha_macro_tuple(enzyme::Symbol, coords::AbstractDict; keq::Real,
                         kf::Real = 1.0, Et::Real = 1.0,
                         release_rate::Real = _default_release_rate(enzyme),
                         release_eq::Real = _default_release_eq(enzyme, coords),
                         kr::Union{Nothing,Real} = nothing,
                         variant::Symbol = :_deploy)
    # HK1 has no Haldane kr (reverse arm is internal via Keq), and its alpha is a FIXED
    # per-variant structural value, NOT a coord — so return early, never calling cha_haldane_kr.
    if enzyme === :HK1
        if variant === :H4
            # Back-map the reparameterized coords {Keff, split_ratio} → {Ki_G6P_C, Ki_G6P_N}.
            # √P = Keff·split_ratio; sumK = Kc+Kn = P/Keff; Kc,Kn = roots of x²−sumK·x+P.
            # split_ratio ≥ 2 (bound) guarantees disc ≥ 0; convention Kc = larger (loose) root.
            Keff  = coords[:Keff]
            ratio = coords[:split_ratio]
            sqrtP = Keff * ratio
            P     = sqrtP^2
            sumK  = P / Keff
            disc  = sumK^2 - 4P
            sq    = sqrt(max(disc, 0.0))
            KiC   = (sumK + sq) / 2
            KiN   = (sumK - sq) / 2
            return (; Kd_Glc   = coords[:Kd_Glc],
                      Kd_ATP   = coords[:Kd_ATP],
                      Ki_G6P_C = KiC,
                      Ki_ADP   = coords[:Ki_ADP],
                      Ki_G6P_N = KiN,
                      K_Pi_N   = coords[:K_Pi_N],
                      alpha    = 1.0,
                      Keq = keq, kf = kf, k2f = kf, Et = Et)
        end
        return (; Kd_Glc   = coords[:Kd_Glc],
                  Kd_ATP   = coords[:Kd_ATP],
                  Ki_G6P_C = coords[:Ki_G6P_C],
                  Ki_ADP   = coords[:Ki_ADP],
                  Ki_G6P_N = coords[:Ki_G6P_N],
                  K_Pi_N   = coords[:K_Pi_N],
                  alpha    = _hk1_variant_alpha(variant),
                  Keq = keq, kf = kf, k2f = kf, Et = Et)   # k2f == kf == 1: Pi competitor-only
    end
    # PGD fully-RE variant: NO promoted SS-release fiber (no koff/kon; C=1). The reverse arm is
    # carried entirely by the Haldane kr; the product Kd's (Kd_NADPH/Kd_Ru5P/Kd_CO2) are real
    # coords. Mirrors the HK1 early-return so the generic release-fiber path below is untouched.
    # Effector dead-ends are appended ONLY when present as coords (default OFF → law uses Inf).
    if enzyme === :PGD && variant === :full_re
        krv = kr === nothing ?
            cha_haldane_kr(enzyme, coords; keq=keq, release_rate=release_rate, kf=kf,
                           release_eq=release_eq, variant=variant) : kr
        tup = (; Kd_NADP  = coords[:Kd_NADP],
                 Kd_PGA   = coords[:Kd_PGA],
                 alpha    = coords[:alpha],
                 Kd_NADPH = coords[:Kd_NADPH],
                 Kd_Ru5P  = coords[:Kd_Ru5P],
                 Kd_CO2   = coords[:Kd_CO2],
                 kf = kf, kr = krv, Et = Et, Keq = keq)
        for s in (:Ki_ATP, :Ki_ATP_EN, :Ki_NADPH)
            haskey(coords, s) && (tup = merge(tup, NamedTuple{(s,)}((coords[s],))))
        end
        return tup
    end
    krv = kr === nothing ?
        cha_haldane_kr(enzyme, coords; keq=keq, release_rate=release_rate, kf=kf,
                       release_eq=release_eq) : kr
    if enzyme === :G6PD
        koffQ = release_rate
        konQ  = koffQ / release_eq           # release_eq == Km_NADPH_rev on the G6PD fiber
        return (; Kd_NADP   = coords[:Kd_NADP],
                  Kd_G6P    = coords[:Kd_G6P],
                  Kd_6PGLn  = coords[:Kd_6PGLn],
                  alpha     = coords[:alpha],
                  Km_NADPH_rev = coords[:Km_NADPH_rev],
                  Ki_NADPH  = get(coords, :Ki_NADPH,  Inf),
                  Ki_ATP    = get(coords, :Ki_ATP,    Inf),
                  Ki_ATP_EG = get(coords, :Ki_ATP_EG, Inf),
                  koffQ = koffQ, konQ = konQ, kf = kf, kr = krv, Et = Et,
                  Keq = keq)
    elseif enzyme === :PGD
        koff = release_rate
        kon  = koff / release_eq             # release_eq == KdRu (Ru5P-release equilibrium)
        return (; Kd_NADP   = coords[:Kd_NADP],
                  Kd_PGA    = coords[:Kd_PGA],
                  alpha     = coords[:alpha],
                  Kd_CO2    = coords[:Kd_CO2],
                  Km_NADPH_rev = coords[:Km_NADPH_rev],
                  Ki_NADPH  = coords[:Ki_NADPH],
                  Ki_ATP    = coords[:Ki_ATP],
                  Ki_ATP_EN = coords[:Ki_ATP_EN],
                  koff = koff, kon = kon, kf = kf, kr = krv, Et = Et,
                  Keq = keq)
    else
        error("cha_macro_tuple: unknown enzyme $enzyme (expected :G6PD, :PGD, or :HK1)")
    end
end

# -----------------------------------------------------------------------------------------
#   Cha-form centered log-ratio loss. The Cha twin of coeff_centered_logratio_loss
#   (coeff_fit.jl): IDENTICAL per-(Article,Fig) centered-loss arithmetic (sign/finite
#   penalty, per-group mean-centering, total/n), only the forward-rate evaluation differs --
#   the Cha law cha_rate_* through the macro tuple assembled from the named coordinates, the
#   gauge (kf, Et), the swept fiber (release_rate), its release equilibrium (release_eq, used
#   for the on-rate and -- on the G6PD fiber -- the Haldane), and the Haldane kr.
#
#   Defaults match cha_macro_tuple so the loss and the tuple agree: kf=Et=1, healthy
#   release_rate, per-enzyme release_eq, Haldane kr.
#
#   `keq`: `nothing` (default) -> PER-FIGURE keq, built once per (Article,Fig) group from
#   that figure's own (asserted-uniform) `d.keq`. A scalar `keq` reproduces the old
#   single-macro-tuple-for-the-whole-dataset behavior EXACTLY (every figure builds the same
#   tuple from the same scalar). The within-(Article,Fig)-uniform keq is cancelled by the
#   per-group centering (same argument as coeff_centered_logratio_loss -- see coeff_fit.jl
#   header) ONLY when keq truly is uniform within the dataset; per-figure keq is exactly the
#   generalization that keeps that per-group-uniform property while letting keq vary ACROSS
#   figures.
# -----------------------------------------------------------------------------------------
function cha_centered_logratio_loss(enzyme::Symbol, mech, d::Dataset,
                                    coords::AbstractDict; keq::Union{Nothing,Real}=nothing,
                                    kf::Real = 1.0, Et::Real = 1.0,
                                    release_rate::Real = _default_release_rate(enzyme),
                                    release_eq::Real = _default_release_eq(enzyme, coords),
                                    kr::Union{Nothing,Real} = nothing,
                                    variant::Symbol = :_deploy)
    cha_rate_enz = enzyme === :G6PD ? ChaLaws.cha_rate_G6PD :
                   enzyme === :PGD  ? ChaLaws.cha_rate_PGD :
                   enzyme === :HK1  ? ChaLawsHK1.cha_rate_HK1 :
                   error("cha_centered_logratio_loss: unknown enzyme $enzyme")
    n = nrows(d)
    logratio = fill(NaN, n)
    penalty = 0.0
    groups = unique(d.group)
    # Single pass per group: gather idx ONCE (was findall'd twice -- once here, once again in
    # a second centering loop). Prediction and per-group centering are computed together, but
    # each group's variance contribution is STASHED in group_variances rather than added to
    # `total` inline, so the final fold order below is byte-identical to the old two-loop
    # structure (`total = penalty` first, THEN group variances added in group order) -- adding
    # inline here would interleave penalty and variance terms into `total` in a different order
    # and risk perturbing the last bits (floating-point addition is not associative).
    group_variances = Vector{Float64}(undef, length(groups))
    for (gi, g) in enumerate(groups)
        idx = findall(==(g), d.group)
        # keq for this figure: scalar override if given, else the figure's own (uniform) d.keq.
        keq_g = if keq === nothing
            ks = unique(d.keq[idx])
            length(ks) == 1 ||
                error("cha_centered_logratio_loss: keq not uniform within figure $g: $ks")
            ks[1]
        else
            keq
        end
        m = cha_macro_tuple(enzyme, coords; keq=keq_g, kf=kf, Et=Et,
                            release_rate=release_rate, release_eq=release_eq, kr=kr,
                            variant=variant)
        for i in idx
            v = cha_rate_enz(m; _cha_row_kwargs(enzyme, d.concs[i])...)
            o = d.rate[i]
            if !isfinite(v) || v == 0 || sign(v) != sign(o)
                penalty += _SIGN_PENALTY
                logratio[i] = NaN
            else
                logratio[i] = log(abs(v)) - log(abs(o))
            end
        end
        vals = filter(isfinite, logratio[idx])
        if isempty(vals)
            group_variances[gi] = 0.0
        else
            μ = sum(vals) / length(vals)
            group_variances[gi] = sum(x -> (x - μ)^2, vals)
        end
    end
    total = penalty
    for v in group_variances
        total += v
    end
    total / n
end

# Map a per-row concentration NamedTuple to the keyword args of the enzyme's Cha law,
# pulling each metabolite by field name (absent field => 0.0, the law's "species absent").
_cha_field(cc, s) = hasproperty(cc, s) ? getproperty(cc, s) : 0.0
function _cha_row_kwargs(enzyme::Symbol, cc)
    if enzyme === :G6PD
        return (NADP = _cha_field(cc, :NADP), G6P = _cha_field(cc, :G6P),
                NADPH = _cha_field(cc, :NADPH), PGLn = _cha_field(cc, :PGLn),
                ATP = _cha_field(cc, :ATP))
    elseif enzyme === :PGD
        return (NADP = _cha_field(cc, :NADP), PGA = _cha_field(cc, :PGA),
                Ru5P = _cha_field(cc, :Ru5P), CO2 = _cha_field(cc, :CO2),
                NADPH = _cha_field(cc, :NADPH), ATP = _cha_field(cc, :ATP))
    elseif enzyme === :HK1
        return (Glucose = _cha_field(cc, :Glucose), ATP = _cha_field(cc, :ATP),
                G6P = _cha_field(cc, :G6P), ADP = _cha_field(cc, :ADP),
                Pi = _cha_field(cc, :Pi))
    else
        error("_cha_row_kwargs: unknown enzyme $enzyme")
    end
end

# FIT-fiber promoted-release rate default (the swept fiber used while FITTING). Under the
# per-(Article,Fig) mean-centered log-ratio loss the overall scale is gauged out, and the only
# fiber-INVARIANT shape observable the forward+product-inhibition corpus constrains is the
# SPECIFICITY constant kcat/Km = kf/(alpha*Kd) (`cha_specificity`). The apparent kcat and the
# apparent Km individually are NOT release-rate-invariant: along the fiber they trade off as
# kcat = kf*r/(kf+r) and Km = alpha*Kd*r/(kf+r) = alpha*Kd/C (C = 1+kf/r), holding kcat/Km
# fixed. So the FIT can sit at any r (koffQ is unidentifiable) — but a Km/kcat READOFF must use
# the SAME r the law is DEPLOYED at (`CHA_DEPLOY_RELEASE_RATE`), or it reports a constant the
# deployed model does not exhibit. (At r=1 -> C=2 -> Km=alpha*Kd/2; at the deploy r=1e3 ->
# C~1.001 -> Km~alpha*Kd, the textbook catalysis-limited value.)
_default_release_rate(enzyme::Symbol) =
    enzyme === :G6PD ? 1.0 :
    enzyme === :PGD  ? 1.0 :
    enzyme === :HK1  ? 1.0 :
    error("_default_release_rate: unknown enzyme $enzyme")

# The promoted-release rate the pipeline DEPLOYS at (see the cha_deploy_micro call in run.jl,
# koffQ = release_rate = this value). The apparent-Km readoff defaults to THIS so the reported
# Km describes the law actually written to model_parameters.jl, not the fit-fiber default.
# Single source of truth: the deploy call and the readoff both read it, so they cannot drift.
const CHA_DEPLOY_RELEASE_RATE = 1.0e3

# -----------------------------------------------------------------------------------------
#   Biophysical log10 bounds (lo, hi) aligned to cha_coords(enzyme). The Kd_*/Ki_*/Km_*_rev
#   shape constants get -9..0 (1 nM .. 1 M), mirroring coord_bounds (coeff_fit.jl). :alpha is
#   a dimensionless interaction factor, NOT a dissociation constant, so it gets a bounded
#   interaction range -2..2 (0.01 .. 100) instead.
# -----------------------------------------------------------------------------------------
function cha_coord_bounds(enzyme::Symbol, variant::Symbol=:_deploy)
    coords = cha_coords(enzyme, variant)
    lo = Float64[]; hi = Float64[]
    for s in coords
        if s === :alpha
            push!(lo, -2.0); push!(hi, 2.0)
        elseif s === :split_ratio
            # √(Kc·Kn)/Keff ∈ [2, 1000] in log10. The floor log10(2) is the real-roots
            # constraint (split_ratio = 2 ⇒ Kc = Kn, the single-site point; > 2 ⇒ two sites).
            push!(lo, log10(2.0)); push!(hi, 3.0)
        else
            push!(lo, -9.0); push!(hi, 0.0)
        end
    end
    lo, hi
end

# -----------------------------------------------------------------------------------------
#   APPARENT Michaelis constant readoff. The substrate Km's are NOT fit coords -- they are
#   DERIVED from the Cha shape constants and the gauge:
#       Km = alpha * Kd / C,   C = 1 + kf/release_rate
#   (release_rate is the promoted SS-release rate -- koffQ for G6PD, the Ru5P-release koff for
#   PGD; on the gauge kf=1.0). `which` selects:
#       :Km_PGA  -> alpha*Kd_PGA/C   (PGD apparent 6PG Michaelis constant)
#       :Km_NADP -> alpha*Kd_NADP/C  (PGD or G6PD apparent NADP Michaelis constant)
#       :Km_G6P  -> alpha*Kd_G6P/C   (G6PD apparent G6P Michaelis constant)
#   `release_rate` DEFAULTS to `CHA_DEPLOY_RELEASE_RATE` (the koffQ the law is deployed at), so
#   the readoff describes the DEPLOYED law. The apparent Km is fiber-DEPENDENT (it carries the
#   1+kf/release_rate factor), so reporting it at any other release_rate -- e.g. the fit-fiber
#   _default_release_rate(enzyme) -- yields a Km the deployed model does NOT exhibit (a factor
#   ~2 at r=1 vs the catalysis-limited deploy). For the fiber-INVARIANT constant use
#   `cha_specificity` (kcat/Km). The Mode-2/3 ANCHOR penalty also reads at this deploy fiber
#   (see `_cha_anchor_penalty`), so it pulls the DEPLOYED apparent Km onto the literature target.
# -----------------------------------------------------------------------------------------
function cha_apparent_km(enzyme::Symbol, coords::AbstractDict, which::Symbol;
                         kf::Real = 1.0,
                         release_rate::Real = CHA_DEPLOY_RELEASE_RATE,
                         variant::Symbol = :_deploy)
    # HK1: C = 1 (no SS-release fiber) and gamma = 1, so apparent Km == binary Kd directly.
    if enzyme === :HK1
        which === :Km_Glc && return coords[:Kd_Glc]
        which === :Km_ATP && return coords[:Kd_ATP]
        error("cha_apparent_km(:HK1): expected :Km_Glc or :Km_ATP (got $which)")
    end
    # Fully-RE PGD is fiber-free (HK1 precedent): C = 1, so apparent Km == alpha*Kd exactly.
    C = (enzyme === :PGD && variant === :full_re) ? 1.0 : 1 + kf / release_rate
    kd = which === :Km_PGA  ? coords[:Kd_PGA] :
         which === :Km_NADP ? coords[:Kd_NADP] :
         which === :Km_G6P  ? coords[:Kd_G6P] :
         error("cha_apparent_km: unknown apparent Km $which (expected :Km_PGA/:Km_NADP/:Km_G6P)")
    coords[:alpha] * kd / C
end

# SPECIFICITY constant kcat/Km = kf/(alpha*Kd) for a substrate -- the koffQ-fiber INVARIANT.
# Apparent kcat and Km each slide along the release-rate fiber (kcat=kf*r/(kf+r),
# Km=alpha*Kd*r/(kf+r)) but their ratio kf/(alpha*Kd) does not, so this is the constant the
# forward corpus actually pins regardless of where the unidentifiable koffQ is parked. `which`
# selects the substrate exactly as `cha_apparent_km`. Units: 1/M on the kf=1 gauge.
function cha_specificity(enzyme::Symbol, coords::AbstractDict, which::Symbol; kf::Real = 1.0)
    enzyme === :HK1 &&
        error("cha_specificity: not defined for HK1 (no SS-release fiber; apparent Km == Kd)")
    kd = which === :Km_PGA  ? coords[:Kd_PGA] :
         which === :Km_NADP ? coords[:Kd_NADP] :
         which === :Km_G6P  ? coords[:Kd_G6P] :
         error("cha_specificity: unknown substrate $which (expected :Km_PGA/:Km_NADP/:Km_G6P)")
    kf / (coords[:alpha] * kd)
end

# Pins: override certain coordinates at fixed log10 values during the fit (Mode-2 twin of
# _coeff_loss_with_pins, coeff_fit.jl). `anchors` carries soft-anchor targets on DERIVED
# apparent Km's (Task 11): a Dict like `Dict(:Km_PGA => (target=<log10 M>, weight=<w>))`. The
# penalty is STRICTLY ADDITIVE on top of the unchanged centered log-ratio base loss -- only
# applied when `anchors` carries entries -- so anchors=nothing (or empty) is a perfect no-op.
function _cha_loss_with_pins(enzyme, mech, d, u, coords_syms, pins, anchors;
                             keq::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy)
    if !isempty(pins)
        u = collect(u)
        for (idx, k) in enumerate(coords_syms)
            haskey(pins, k) && (u[idx] = pins[k])
        end
    end
    coords_dict = Dict(coords_syms .=> 10 .^ u)
    L = cha_centered_logratio_loss(enzyme, mech, d, coords_dict; keq=keq, variant=variant)
    L += _cha_anchor_penalty(enzyme, coords_dict, anchors)
    L
end

# Additive soft-anchor penalty on the apparent Michaelis constants. For each `which =>
# (target, weight)` in `anchors`, adds weight*(log10(apparent_Km) - target)^2. anchors=nothing
# (or empty) contributes EXACTLY 0.0 -- the base centered-loss arithmetic is untouched.
function _cha_anchor_penalty(enzyme, coords_dict, anchors)
    anchors === nothing && return 0.0
    pen = 0.0
    for (which, spec) in anchors
        # Anchor on the DEPLOY fiber (CHA_DEPLOY_RELEASE_RATE -- the release rate the law is
        # actually deployed at), so the penalty pulls the DEPLOYED apparent Km to `target`.
        # (Through 2026-06-11 this read the FIT fiber `_default_release_rate`, which left the
        # deployed apparent Km ~(C_fit/C_deploy)x the target -- ~2x for PGD Km_PGA. Re-anchored
        # on the deploy fiber 2026-06-15 so a PGD consensus deploy lands Km_PGA on the literature
        # band, not 2x it; this is a re-fit -- Mode-2/3 PGD results change, Mode-1 is unaffected.)
        km = cha_apparent_km(enzyme, coords_dict, which)   # release_rate defaults to CHA_DEPLOY_RELEASE_RATE
        pen += spec.weight * (log10(km) - spec.target)^2
    end
    pen
end

# -----------------------------------------------------------------------------------------
#   CMA-ES multi-start basin-hopping over cha_coords(enzyme) in LOG10 space. Mirrors
#   coeff_fit_candidate (coeff_fit.jl) EXACTLY: hermetic per-restart seed hash((seed,r)),
#   MersenneTwister-seeded uniform x0, maxiter/maxtime bounds, verbosity=0, best-of-restarts
#   pick, restarts NamedTuple vector, and the post-fit pin-overwrite. The determinism scheme
#   is preserved so the pipeline's pmap/seed contract is unchanged. Only the objective differs
#   -- it threads through the Cha law (cha_centered_logratio_loss) over the named coords.
#   `anchors` (Task 11) carries soft-anchor targets on the DERIVED apparent Michaelis
#   constants -- e.g. `Dict(:Km_PGA => (target=log10(59e-6), weight=1.0))` pulls the PGD
#   apparent Km_PGA toward the 38-80µM literature band midpoint (Mode 2). The penalty is
#   STRICTLY ADDITIVE on the centered base loss and applied only when anchors is non-empty, so
#   anchors=nothing reproduces the Mode-1 fit bit-for-bit. The default anchor weight (1.0 at
#   the call site) is PROVISIONAL -- to be resolved against the Phase-C flux gate (spec §11 /
#   open item T1-a). Returns (coords::Dict, loss, restarts).
# -----------------------------------------------------------------------------------------
function cha_fit_candidate(enzyme::Symbol, mech, d::Dataset; n_restarts::Int=8,
                           maxiter::Int=1_000_000, maxtime::Real=20.0, seed::Int=1,
                           keq::Union{Nothing,Real}=nothing,
                           pins::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
                           anchors=nothing, variant::Symbol=:_deploy)
    coords_syms = cha_coords(enzyme, variant)
    lo, hi = cha_coord_bounds(enzyme, variant)
    objective = u -> _cha_loss_with_pins(enzyme, mech, d, u, coords_syms, pins, anchors;
                                         keq=keq, variant=variant)
    best = (u=fill(NaN, length(coords_syms)), loss=Inf)
    endpoints = Vector{Float64}[]
    losses = Float64[]
    for r in 1:n_restarts
        s   = hash((seed, r)) % UInt64
        rng = Random.MersenneTwister(s)
        x0  = lo .+ (hi .- lo) .* rand(rng, length(coords_syms))
        o   = minimize(objective, x0, 0.1; lower=lo, upper=hi,
                       maxiter=maxiter, maxtime=maxtime, seed=s, verbosity=0)
        u   = collect(xbest(o))
        push!(endpoints, u)
        push!(losses, fbest(o))
        if fbest(o) < best.loss
            best = (u=u, loss=fbest(o))
        end
    end
    restarts = NamedTuple[(loss=losses[i], dist=norm(endpoints[i] .- best.u))
                          for i in eachindex(endpoints)]
    # Overwrite pinned coordinates with their pinned values (CMA-ES searched them but the
    # loss ignored them -- the raw optimizer value is meaningless).
    final = copy(best.u)
    for (idx, k) in enumerate(coords_syms)
        haskey(pins, k) && (final[idx] = pins[k])
    end
    coords = Dict(coords_syms .=> 10 .^ final)
    (coords=coords, loss=best.loss, restarts=restarts)
end

# -----------------------------------------------------------------------------------------
#   ERROR-on-no-op GUARD (carried over from pins.jl::resolve_pins). Every name a Cha pin
#   intends to clamp MUST be a member of cha_coords(enzyme); a pin that silently no-ops would
#   let the report label a constant :literature_pinned while the fit ignored it (report and
#   fit decouple). All names resolve_cha_pins emits ARE coords on the happy path, so this only
#   fires on a future coord-set change.
# -----------------------------------------------------------------------------------------
function _assert_pin_is_coord(enzyme::Symbol, name::Symbol, variant::Symbol=:_deploy)
    name in cha_coords(enzyme, variant) && return nothing
    error("resolve_cha_pins: intended pin :$name (enzyme=$enzyme) is NOT a member of " *
          "cha_coords($enzyme) — the pin would be a silent no-op while the report still labels " *
          "it :literature_pinned at the anchor (report and fit disagree). Add :$name to " *
          "cha_coords, or remove it from the Cha pin set.")
end

# -----------------------------------------------------------------------------------------
#   Cha-path MODE-PIN resolution (the hard coord-pins). Returns Dict{cha_coord => log10(M)}.
#   Literature values come from FitRateEquation._lit_values(enzyme) (macro-name => log10(M)).
#
#   Pin rules (only names that ARE members of cha_coords(enzyme) are emitted; the guard
#   _assert_pin_is_coord protects against a future coord-set change):
#
#   ALL MODES (incl. Mode 1) — anchor the conflating reverse channel:
#     - G6PD: pin :Km_NADPH_rev to lit[:Km_NADPH_rev] (3.9µM). This protects the Ki_NADPH
#       cross-term identifiability (Task 8: leaving Km_NADPH_rev free makes the [G6P·NADPH]
#       dead-end Ki_NADPH trade off with the bare-[NADPH] productive-release reverse channel,
#       railing Ki_NADPH; anchoring Km_NADPH_rev de-conflates it to ~24µM). Kd_6PGLn is left
#       FREE (a harmless 6PGL-binding nuisance → classifies :unconstrained on forward data).
#       The `anchor_reverse` kwarg (default true) can turn this off: `anchor_reverse=false`
#       leaves Km_NADPH_rev FREE, deliberately reintroducing the forward/reverse Ki_NADPH
#       conflation. This is a DIAGNOSTIC ONLY — the deployed law REQUIRES the anchor; a fit
#       with anchor_reverse=false is not deployable (Ki_NADPH goes non-identifiable).
#     - PGD: nothing always-anchored. PGD's reverse constants have no literature values
#       (Dalziel nuisance) so they are left free/unconstrained.
#
#   MODE 2 AND MODE 3 — additionally pin the forward inhibition constants to literature:
#     - both enzymes: :Ki_ATP → lit[:Ki_ATP].
#     - both enzymes: :Ki_NADPH → lit[:Ki_NADPH] (G6PD 15µM; PGD 17µM). In Mode 1 Ki_NADPH is
#       left FREE: G6PD Mode-1 data-identifies it (with Km_NADPH_rev anchored); PGD Mode-1
#       reports it diagnostic/unconstrained (the cross-term de-conflation was refuted there).
#
#   NOT pinned here: :Km_PGA — an APPARENT constant (alpha·Kd_PGA/C), not a coord. Its Mode-3
#   PGD hard override and Mode-2 PGD soft anchor go through the `anchors` mechanism in run.jl,
#   NEVER as a hard coord-pin (mirrors pins.jl::resolve_coord_pins which keeps Km_PGA on the
#   coord side; here Km_PGA is not even a coord).
# -----------------------------------------------------------------------------------------
function resolve_cha_pins(enzyme::Symbol, variant::Symbol, mode::Symbol; anchor_reverse::Bool=true)
    lit    = FitRateEquation._lit_values(enzyme)
    coords = cha_coords(enzyme, variant)
    pins   = Dict{Symbol,Float64}()

    # Emit a pin only after asserting the name is a real coord (ERROR-on-no-op). The literature
    # value MUST exist for any name we intend to pin; a missing lit entry is also a no-op risk.
    function _pin!(name::Symbol)
        _assert_pin_is_coord(enzyme, name, variant)
        haskey(lit, name) || error("resolve_cha_pins: intended pin :$name (enzyme=$enzyme, " *
            "mode=$mode) has no literature value in _lit_values($enzyme) — cannot anchor it.")
        pins[name] = lit[name]
        return nothing
    end

    # HK1 per-mode pin sets. Mode 1: nothing pinned (all forward shape constants free).
    # Mode 2: pin the N-half regulatory constants (Ki_G6P_N, K_Pi_N) to literature. Mode 3:
    # additionally pin the C-half product-inhibition constants (Ki_G6P_C, Ki_ADP).
    if enzyme === :HK1
        # H4 is the data-driven reparameterized variant {Keff, split_ratio}: NO pins in any mode
        # (the literature pin names Ki_G6P_N/Ki_G6P_C are not H4 coords by construction).
        variant === :H4 && return pins
        if mode === :mode2 || mode === :mode3
            _pin!(:Ki_G6P_N); _pin!(:K_Pi_N)
        end
        if mode === :mode3
            _pin!(:Ki_G6P_C); _pin!(:Ki_ADP)
        end
        if mode ∉ (:mode1, :mode2, :mode3)
            error("resolve_cha_pins: unknown mode :$mode (expected :mode1/:mode2/:mode3)")
        end
        return pins
    end

    # ALL MODES: anchor the conflating reverse channel where it is a coord with a lit value.
    # `anchor_reverse=false` leaves it FREE (diagnostic — reintroduces the Ki_NADPH conflation).
    anchor_reverse && enzyme === :G6PD && (:Km_NADPH_rev in coords) && _pin!(:Km_NADPH_rev)

    # MODE 2 / MODE 3: additionally pin the forward inhibition constants to literature.
    if mode === :mode2 || mode === :mode3
        (:Ki_ATP   in coords) && _pin!(:Ki_ATP)
        (:Ki_NADPH in coords) && _pin!(:Ki_NADPH)
    elseif mode !== :mode1
        error("resolve_cha_pins: unknown mode :$mode (expected :mode1, :mode2, or :mode3)")
    end

    pins
end

end # module ChaFit

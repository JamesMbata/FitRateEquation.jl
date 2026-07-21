# =========================================================================================
#         Closed-form macro <-> micro map for the Cha-form G6PD law (no LM bridge)
# =========================================================================================
#
# `cha_macro_readoffs_G6PD(m, logθ; keq)` reads a complete Cha macro tuple OFF a real micro
# mechanism `m` (v2_mechanism()) at parameters logθ. The resulting tuple, fed to
# `ChaLaws.cha_rate_G6PD` (the corrected two-SS-segment super-node law), reproduces
# `EnzymeRates.rate_equation` exactly (A2 exactness anchor, rtol 1e-10).
#
# `cha_micro_from_macro_G6PD(mac; koffQ)` is the closed-form fiber map: it slides the
# promoted release rate koffQ to a new value, co-adjusting (alpha, konQ) so that every
# DATA-IDENTIFIABLE macro constant (Kd's, Ki's, apparent Km's which depend only on alpha/C,
# and Km_NADPH_rev) is invariant. cha_rate_G6PD is therefore koffQ-invariant on the fiber.
#
# Symbol/constant bindings (introspected from v2_mechanism(), 2026-06-10):
#   catalysis E_NB<-->E_C  : kf = k5f (gauge=1), kr = k5r (Haldane-dependent)
#   NADPH release E_H<-->E+NADPH (SS) : koffQ = k7f (release/off), konQ = k7r (rebind/on)
#       => Km_NADPH_rev = koffQ/konQ = k7f/k7r
#   PGLn release E_C<->E_H+PGLn (RE)  : Kd_6PGLn = K6 (step-6 binding constant)
#       NOTE: this is the law's KdP, NOT macro_constants' :Ki_6PGLn (= KdP*koff/kr).
#   CORRECTED Cha C-factor: C = 1 + kf/koffQ  (koffQ alone, NOT kr+koffQ).
# =========================================================================================

module ChaInvert

using ..FitRateEquation
using EnzymeRates

export cha_macro_readoffs_G6PD, cha_micro_from_macro_G6PD

# Bind the Cha rate-constant symbols to their UPSTREAM-EnzymeRates semantic names by SHAPE
# (the imposed topology makes each unique). Upstream renamed rate constants from step-index
# (k5f, k7f, K6) to composition-semantic (k_EG6PNADP_to_ENADPHPGLn, koff_NADPH_E, K_PGLn_ENADPH),
# so we identify each by pattern rather than constructing `k<idx>f`. fitted_params holds the
# independents; the Haldane-reverse chemistry rate is a DEPENDENT param (via _dependent_param_exprs).
_cha_all_ksyms(m) = vcat(collect(EnzymeRates.fitted_params(m)),
                         [s for (s, _) in EnzymeRates._dependent_param_exprs(typeof(m))[1]])
_is_kto(s) = (t = String(s); startswith(t, "k_") && occursin("_to_", t))   # chemistry k_<from>_to_<to>
_ktoks(s)  = split(String(s), "_")                                          # K_<met>_<form> tokens
function _only_sym(pool, pred, what)
    h = filter(pred, pool)
    length(h) == 1 || error("Cha symbol bind: expected exactly 1 $what, got $(h)")
    first(h)
end
_opt_sym(pool, pred) = (h = filter(pred, pool); isempty(h) ? nothing : first(h))

# Bind (a) the catalysis kf (unique forward chemistry `k_*_to_*` = gauge) and its Haldane
# reverse kr (the other `k_*_to_*`, dependent), (b) the SS NADPH-release koffQ/konQ (unique
# `koff_*`/`kon_*`), (c) the PGLn RE-binding K (`K_PGLn_*`), and (d) the ATP dead-end on E·G6P
# (`K_ATP..._*G6P*`; nothing for the ATP-free mechanism).
function _g6pd_cha_symbols(m)
    fp   = collect(EnzymeRates.fitted_params(m))
    allk = _cha_all_ksyms(m)
    kf_sym    = gauge_param(m)
    kr_sym    = _only_sym(allk, s -> _is_kto(s) && s != kf_sym, "reverse chemistry rate")
    koffQ_sym = _only_sym(fp, s -> startswith(String(s), "koff"), "SS release off-rate")
    konQ_sym  = _only_sym(fp, s -> startswith(String(s), "kon"),  "SS release on-rate")
    pgln_sym  = _only_sym(fp, s -> startswith(String(s), "K") && occursin("PGLn", String(s)),
                          "PGLn binding K")
    atpEG_sym = _opt_sym(fp, s -> (t = _ktoks(s); length(t) == 3 && t[1] == "K" &&
                                   occursin("ATP", t[2]) && occursin("G6P", t[3])))
    (kf_sym, kr_sym, koffQ_sym, konQ_sym, pgln_sym, atpEG_sym)
end

# Read a complete Cha macro tuple off the micro mechanism `m` at `logθ`.
function cha_macro_readoffs_G6PD(m, logθ; keq::Real)
    mc    = FitRateEquation.macro_constants(m, logθ; keq=keq, enzyme=:G6PD)
    named = Dict(x.name => x.value for x in mc)
    vals  = FitRateEquation._micro_values(m, logθ; keq=keq)

    kf_sym, kr_sym, koffQ_sym, konQ_sym, pgln_sym, atpEG_sym = _g6pd_cha_symbols(m)
    kf    = vals[kf_sym]
    kr    = vals[kr_sym]
    koffQ = vals[koffQ_sym]
    konQ  = vals[konQ_sym]

    # Corrected Cha C-factor (verified: the two-SS-segment ternary coeff is gAB*C with
    # C = 1 + kf/koffQ, koffQ alone -- see cha_derive_g6pd.py Property 1).
    C = 1 + kf / koffQ

    Kd_NADP   = named[:Kd_NADP]     # = 1/g[NADP]
    Kd_G6P    = named[:Kd_G6P]      # = 1/g[G6P]
    Ki_NADPH  = get(named, :Ki_NADPH, Inf)  # = g[G6P]/g[G6P·NADPH]; Inf when no E·G6P·NADPH dead-end
    Ki_ATP    = get(named, :Ki_ATP, Inf)              # Inf when no free-E ATP dead-end
    Ki_ATP_EG = atpEG_sym === nothing ? Inf : vals[atpEG_sym]  # Inf when no E·G6P ATP dead-end
    Kd_6PGLn  = vals[pgln_sym]      # = K6 (step-6 PGLn binding), NOT macro :Ki_6PGLn

    # Km_NADP = alpha*Kd_NADP/C  =>  alpha = (Km_NADP/Kd_NADP)*C  (same from Km_G6P/Kd_G6P).
    alpha = (named[:Km_NADP] / Kd_NADP) * C

    Km_NADPH_rev = koffQ / konQ     # productive NADPH release equilibrium of v2
    Et = vals[:E_total]

    (; Kd_NADP, Kd_G6P, Kd_6PGLn, alpha, Km_NADPH_rev, Ki_NADPH, Ki_ATP, Ki_ATP_EG,
       koffQ, konQ, kf, kr, Et,
       Keq = get(vals, :Keq, NaN),
       Km_NADP_apparent = named[:Km_NADP],
       Km_G6P_apparent  = named[:Km_G6P])
end

# Closed-form fiber map: re-express the same observable macro law with a new promoted
# release rate `koffQ`. alpha is re-scaled so alpha/C is invariant (apparent Km's depend
# only on alpha/C), and konQ is set to hold Km_NADPH_rev = koffQ/konQ fixed.
function cha_micro_from_macro_G6PD(mac; koffQ::Real)
    kf = mac.kf
    kr = mac.kr
    C_macro = 1 + mac.kf / mac.koffQ
    C_new   = 1 + kf / koffQ
    alpha   = (mac.alpha / C_macro) * C_new
    konQ    = koffQ / mac.Km_NADPH_rev

    (; mac.Kd_NADP, mac.Kd_G6P, mac.Kd_6PGLn, alpha, mac.Km_NADPH_rev,
       mac.Ki_NADPH, mac.Ki_ATP, mac.Ki_ATP_EG, koffQ, konQ, kf, kr, mac.Et,
       Keq = get(mac, :Keq, NaN),
       Km_NADP_apparent = get(mac, :Km_NADP_apparent, NaN),
       Km_G6P_apparent  = get(mac, :Km_G6P_apparent, NaN))
end

export cha_macro_readoffs_PGD, cha_micro_from_macro_PGD

# Bind the PGD :cha_base Cha symbols to their upstream semantic names by SHAPE (as for G6PD):
# (a) catalysis kf (unique `k_*_to_*` = gauge) + Haldane reverse kr (the other `k_*_to_*`);
# (b) the SS Ru5P-release koff/kon (unique `koff_*`/`kon_*`); (c) CO2-release binding (`K_CO2_*`);
# (d) the bare-[NADPH] reverse-release Km (`K_NADPH_*` on the FREE enzyme, form token "E",
# distinct from the E·PGA NADPH dead-end); (e) the ATP dead-end on E·NADP (`K_ATP..._*NADP*`).
function _pgd_cha_symbols(m)
    fp   = collect(EnzymeRates.fitted_params(m))
    allk = _cha_all_ksyms(m)
    kf_sym   = gauge_param(m)
    kr_sym   = _only_sym(allk, s -> _is_kto(s) && s != kf_sym, "reverse chemistry rate")
    koff_sym = _only_sym(fp, s -> startswith(String(s), "koff"), "SS release off-rate")
    kon_sym  = _only_sym(fp, s -> startswith(String(s), "kon"),  "SS release on-rate")
    co2_sym  = _only_sym(fp, s -> startswith(String(s), "K") && occursin("CO2", String(s)),
                         "CO2 binding K")
    nadphRev_sym = _only_sym(fp, s -> (t = _ktoks(s); length(t) == 3 && t[1] == "K" &&
                             occursin("NADPH", t[2]) && t[3] == "E"), "NADPH reverse-release Km")
    atpEN_sym = _only_sym(fp, s -> (t = _ktoks(s); length(t) == 3 && t[1] == "K" &&
                          occursin("ATP", t[2]) && occursin("NADP", t[3])), "ATP-EN dead-end K")
    (kf_sym, kr_sym, koff_sym, kon_sym, co2_sym, nadphRev_sym, atpEN_sym)
end

# Read a complete Cha macro tuple off the PGD micro mechanism `m` at `logθ`.
function cha_macro_readoffs_PGD(m, logθ; keq::Real)
    mc    = FitRateEquation.macro_constants(m, logθ; keq=keq, enzyme=:PGD)
    named = Dict(x.name => x.value for x in mc)
    vals  = FitRateEquation._micro_values(m, logθ; keq=keq)

    kf_sym, kr_sym, koff_sym, kon_sym, co2_sym, nadphRev_sym, atpEN_sym = _pgd_cha_symbols(m)
    kf   = vals[kf_sym]
    kr   = vals[kr_sym]
    koff = vals[koff_sym]
    kon  = vals[kon_sym]

    # Cha C-factor (verified by the exactness gate: two-SS-segment ternary coeff is gAB*C
    # with C = 1 + kf/koff, koff alone -- see cha_derive_pgd.py Property 1).
    C = 1 + kf / koff

    Kd_NADP      = named[:Kd_NADP]      # = 1/g[NADP]
    Kd_PGA       = named[:Kd_PGA]       # = 1/g[PGA]
    Ki_NADPH     = named[:Ki_NADPH]     # = g[PGA]/g[PGA·NADPH] (forward dead-end cross term)
    Ki_ATP       = named[:Ki_ATP]       # = 1/g[ATP] = K9 (ATP dead-end on free E)
    Ki_ATP_EN    = vals[atpEN_sym]      # = K10 (ATP dead-end on E.NADP; distinct from Ki_ATP)
    Kd_CO2       = vals[co2_sym]        # = K6 (CO2-release RE binding constant)
    Km_NADPH_rev = vals[nadphRev_sym]   # = K8 (bare-[NADPH] reverse-release Km)

    # Km_NADP = alpha*Kd_NADP/C  =>  alpha = (Km_NADP/Kd_NADP)*C  (same from Km_PGA/Kd_PGA).
    alpha = (named[:Km_NADP] / Kd_NADP) * C

    Et = vals[:E_total]

    (; Kd_NADP, Kd_PGA, alpha, Kd_CO2, Km_NADPH_rev, Ki_NADPH, Ki_ATP, Ki_ATP_EN,
       koff, kon, kf, kr, Et,
       Keq = get(vals, :Keq, NaN),
       Km_NADP_apparent = named[:Km_NADP],
       Km_PGA_apparent  = named[:Km_PGA])
end

# Closed-form macro->micro map for PGD. Unlike G6PD, PGD ships BASE-ONLY: there is no
# promoted silent fiber to sweep, so this is the IDENTITY-class map. The promoted SS step is
# Ru5P release (koff/kon). By default `k_off_Ru5P = mac.koff`, which round-trips `mac`
# bit-identically. A non-default value re-scales alpha to hold alpha/C invariant (apparent
# Km's depend only on alpha/C; C = 1 + kf/koff) and the Ru5P-rebind on-rate kon to hold the
# Ru5P-release equilibrium koff/kon fixed -- EXACTLY analogous to the G6PD koffQ fiber map.
# PGD never sweeps it in practice (no silent fiber); the kwarg exists only for symmetry.
function cha_micro_from_macro_PGD(mac; k_off_Ru5P::Real = mac.koff)
    kf   = mac.kf
    kr   = mac.kr
    koff = k_off_Ru5P
    C_macro = 1 + mac.kf / mac.koff
    C_new   = 1 + kf / koff
    alpha   = (mac.alpha / C_macro) * C_new
    # Hold the promoted-release equilibrium koff/kon fixed (analogous to Km_NADPH_rev).
    Km_Ru5P_rev_macro = mac.koff / mac.kon
    kon = koff / Km_Ru5P_rev_macro

    (; mac.Kd_NADP, mac.Kd_PGA, alpha, mac.Kd_CO2, mac.Km_NADPH_rev,
       mac.Ki_NADPH, mac.Ki_ATP, mac.Ki_ATP_EN, koff, kon, kf, kr, mac.Et,
       Keq = get(mac, :Keq, NaN),
       Km_NADP_apparent = get(mac, :Km_NADP_apparent, NaN),
       Km_PGA_apparent  = get(mac, :Km_PGA_apparent, NaN))
end

end # module ChaInvert

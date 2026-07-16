# HK1 consensus wiring — Cha-form two-half hexokinase. The C-half is the catalytic
# random Bi-Bi site (Glc, ATP in; G6P, ADP out) with the two product dead-ends
# (E·Glc·G6P, E·Glc·ADP); the N-half is the G6P/Pi regulatory site. The two G6P-site
# candidates differ ONLY by the alpha structural axis: H1 = two independent G6P sites
# (alpha=1), H3 = mutual exclusion / Liu (alpha=infinity, drops the [G6P]^2 coupling).
# alpha is threaded per-variant through cha_fit.jl::_hk1_variant_alpha — it is NOT a fit
# coord, so the two variants share the SAME six forward shape constants (cha_coords(:HK1)).
#
# Unlike G6PD/PGD there is no King-Altman macro readoff for HK1 (Cha-form macro/micro
# hybrid; the macro constants are fit DIRECTLY via cha_rate_HK1, Q3 resolved). So the
# ki_map / kd_map / km_ratio alias maps are EMPTY (unused), and the reverse arm is
# internal via Keq (no Haldane kr). lit_values store log10(M), matching enzymes/g6pd.jl
# and enzymes/pgd.jl, so resolve_cha_pins' _pin! reads them as log10 anchors directly.

include(joinpath(@__DIR__, "..", "..", "HK1", "king_altman", "hk1_candidate_equations.jl"))

# This file is included BOTH inside the FitRateEquation module (FitRateEquation_impl.jl, for
# the pipeline runner) AND at Main level (test_hk1_fit.jl header). `EnzymeWiring` and
# `register_enzyme!` are exported (reachable bare in both contexts), but `MonoKey` is NOT
# exported — so resolve it against whichever module owns the definition.
const _HK1_CM = isdefined(@__MODULE__, :MonoKey) ? (@__MODULE__) : FitRateEquation
const _HK1_MonoKey = _HK1_CM.MonoKey

const _HK1_MECH_H1 = build_hk1_mechanism(alpha=:one, glc_g6p_dead_end=true, glc_adp_dead_end=true)
# H4 = the SAME alpha=1 King-Altman law as H1; it differs ONLY in the fit parameterization
# (cha_coords(:HK1,:H4) = {Keff, split_ratio} reparameterization of {Ki_G6P_C, Ki_G6P_N}).
# H3 (alpha=:infinity) was removed — the reverse-rate turnover requires the [G6P]² term H3
# deletes (mechanistically refuted; see notes/2026-06-13_hk1_g6p_ridge_resolution_report.md).
const _HK1_MECH_H4 = _HK1_MECH_H1

register_enzyme!(EnzymeWiring(
    :HK1, :Glucose,
    [(name=:H1, mech=_HK1_MECH_H1),
     (name=:H4, mech=_HK1_MECH_H4)],
    Dict{Symbol,Float64}(
        :Ki_G6P_N => log10(6.9e-6),
        :K_Pi_N   => log10(750e-6),
        :Ki_G6P_C => log10(15e-6),
        :Ki_ADP   => log10(1.5e-3),
    ),
    Dict{Symbol,Vector{Symbol}}(
        :H1 => [:Ki_G6P_N, :K_Pi_N, :Ki_G6P_C, :Ki_ADP],
        :H4 => Symbol[],   # H4 is data-driven (no pins)
    ),
    Dict{_HK1_MonoKey,Symbol}(),                  # ki_map  — unused (no King-Altman readoff)
    Dict{_HK1_MonoKey,Symbol}(),                  # kd_map  — unused
    Dict{_HK1_MonoKey,Tuple{Symbol,_HK1_MonoKey}}(),  # km_ratio — unused
    _HK1_MonoKey([:ATP => 1, :Glucose => 1]),     # substrate_pair (nominal, sorted)
))

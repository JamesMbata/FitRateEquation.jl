# G6PD consensus wiring (random Bi-Bi RE; SS NADPH-release de-conflation variant).
# ATP competitive vs NADP -> binds open-dinucleotide forms [E, E_G].
# NADPH dead-end: on [:E_G] ONLY (option B, 2026-06-09; E_N + free-E E·NADPH dropped).
# Rationale: the free-E E·NADPH dead-end is the REVERSE of the productive NADPH release
# (double-counts NADPH binding to free E) and is non-identifiable (the decoupled forward
# DOF never reaches the forward band — the PGD-V3 result). Keeping it makes the bare-[NADPH]
# competitive coefficient a 2-knob lump (release K7 + dead-end K10) that cannot be cleanly
# ODE-pinned. Dropping it leaves the release K7 as the SOLE bare-[NADPH] term => forward
# Ki_NADPH is a clean single knob (= K7), pinnable in the ODE the PGD k7r way. The E_G
# dead-end (E·G6P·NADPH, noncompetitive vs G6P) is the physiologically load-bearing one
# (Wang 2002 Fig 4; §15.1) and is KEPT. See docs/G6PD_session_context.md (2026-06-09).

const _G6PD_ATP_FORMS   = [:E, :E_G]
const _G6PD_NADPH_FORMS = [:E_G]

function _g6pd_consensus(nadph_release_op::Symbol)
    extra = _deadends([(_G6PD_ATP_FORMS, :ATP), (_G6PD_NADPH_FORMS, :NADPH)])
    _mech([:NADP, :G6P], [:NADPH, :PGLn], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :G6P],    [:E_G],  :(⇌)),
        ([:E_N, :G6P],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_H, :PGLn],  :(⇌)),               # PGLn release (RE, both variants)
        ([:E_H], [:E, :NADPH],   nadph_release_op),   # NADPH release: ⇌ / <-->
    ], extra); regs=[:ATP])
end

# ATP-free consensus: the random-Bi-Bi SS-NADPH-release skeleton with the ATP dead-ends
# DROPPED (the NADPH dead-end on E·G6P is kept). This is the mechanism behind the :no_atp
# fit variant; the ATP terms are additive free-E-pool dead-ends, so this is the exact
# Ki_ATP,Ki_ATP_EG -> ∞ limit of _g6pd_consensus(:(<-->)). No ATP regulator.
function _g6pd_consensus_noatp()
    extra = _deadends([(_G6PD_NADPH_FORMS, :NADPH)])
    _mech([:NADP, :G6P], [:NADPH, :PGLn], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :G6P],    [:E_G],  :(⇌)),
        ([:E_N, :G6P],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_H, :PGLn],  :(⇌)),               # PGLn release (RE)
        ([:E_H], [:E, :NADPH],   :(<-->)),            # NADPH release (SS, promoted)
    ], extra))
end

# Three dead-end-dropped ablation variants, each an exact ∞-limit of _g6pd_consensus(:(<-->))
# (the deployed SS_NADPH_release_rate_eq law): unlike :no_atp/:no_g6p_nadph_deadend below,
# these do NOT remove ATP as a metabolite (E·ATP, K8=Ki_ATP, is kept in all three), so
# `regs=[:ATP]` is retained throughout. See docs/G6PD_kinetics_literature.md §2.3/§15.4/§17.

# Drop ONLY the E·G6P·NADPH dead-end (K10=Ki_NADPH -> ∞); keep both ATP dead-end forms
# (K8=Ki_ATP, K9=Ki_ATP_EG) unchanged. Tests whether the bare-[NADPH] reverse-release channel
# alone (anchor_reverse=false at fit time) can carry the forward NADPH product inhibition.
function _g6pd_consensus_no_g6p_nadph_deadend()
    extra = _deadends([(_G6PD_ATP_FORMS, :ATP)])
    _mech([:NADP, :G6P], [:NADPH, :PGLn], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :G6P],    [:E_G],  :(⇌)),
        ([:E_N, :G6P],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_H, :PGLn],  :(⇌)),               # PGLn release (RE)
        ([:E_H], [:E, :NADPH],   :(<-->)),            # NADPH release (SS, promoted)
    ], extra); regs=[:ATP])
end

# Drop ONLY the E·G6P·ATP dead-end (K9=Ki_ATP_EG -> ∞): ATP dead-end restricted to free-E form
# only ([:E], i.e. K8=Ki_ATP kept). NADPH dead-end (K10=Ki_NADPH) unchanged. Literature Ki_ATP_EG
# (Ninfali 1996, ATP·Mg uncompetitive vs G6P) is ≈30 mM, superphysiological — expect droppable.
function _g6pd_consensus_no_g6p_atp_deadend()
    extra = _deadends([([:E], :ATP), (_G6PD_NADPH_FORMS, :NADPH)])
    _mech([:NADP, :G6P], [:NADPH, :PGLn], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :G6P],    [:E_G],  :(⇌)),
        ([:E_N, :G6P],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_H, :PGLn],  :(⇌)),               # PGLn release (RE)
        ([:E_H], [:E, :NADPH],   :(<-->)),            # NADPH release (SS, promoted)
    ], extra); regs=[:ATP])
end

# Drop BOTH: E·G6P·NADPH (K10 -> ∞) AND E·G6P·ATP (K9 -> ∞). Keeps only the free-E ATP
# dead-end (K8=Ki_ATP). Combines both ablations above in one mechanism.
function _g6pd_consensus_no_g6p_both_deadends()
    extra = _deadends([([:E], :ATP)])
    _mech([:NADP, :G6P], [:NADPH, :PGLn], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :G6P],    [:E_G],  :(⇌)),
        ([:E_N, :G6P],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_H, :PGLn],  :(⇌)),               # PGLn release (RE)
        ([:E_H], [:E, :NADPH],   :(<-->)),            # NADPH release (SS, promoted)
    ], extra); regs=[:ATP])
end

const _G6PD_KI_MAP = Dict{MonoKey,Symbol}(
    [:ATP   => 1] => :Ki_ATP,
    [:PGLn  => 1] => :Ki_6PGLn,
)
# Forward Ki_NADPH is the noncompetitive dead-end on E·G6P: read as the cross-term
# ratio g[G6P]/g[G6P·NADPH], NOT the bare [NADPH] term (= reverse release Km). The cross
# monomial is sorted by symbol (:G6P < :NADPH) to match `_split_mono`.
const _G6PD_KI_RATIO = Dict{Symbol,Tuple{MonoKey,MonoKey}}(
    :Ki_NADPH => ([:G6P => 1], [:G6P => 1, :NADPH => 1]),
)
const _G6PD_KD_MAP = Dict{MonoKey,Symbol}(
    [:NADP => 1] => :Kd_NADP,
    [:G6P  => 1] => :Kd_G6P,
)
const _G6PD_KM_RATIO = Dict{MonoKey,Tuple{Symbol,MonoKey}}(
    [:NADP => 1] => (:Km_NADP, [:G6P => 1]),
    [:G6P  => 1] => (:Km_G6P,  [:NADP => 1]),
)
const _G6PD_SUBSTRATE_PAIR = MonoKey([:G6P => 1, :NADP => 1])   # sorted (G6P < NADP)

register_enzyme!(EnzymeWiring(
    :G6PD, :G6P,
    [(name=:RE_rate_eq,                mech=_g6pd_consensus(:(⇌))),
     (name=:SS_NADPH_release_rate_eq,  mech=_g6pd_consensus(:(<-->))),
     (name=:no_atp,                    mech=_g6pd_consensus_noatp()),
     (name=:no_g6p_nadph_deadend,      mech=_g6pd_consensus_no_g6p_nadph_deadend()),
     (name=:no_g6p_atp_deadend,        mech=_g6pd_consensus_no_g6p_atp_deadend()),
     (name=:no_g6p_both_deadends,      mech=_g6pd_consensus_no_g6p_both_deadends())],
    Dict{Symbol,Float64}(
        :Ki_6PGLn     => log10(2.0e-4),
        :Ki_ATP       => log10(1.5e-3),
        :Km_NADPH_rev => log10(3.9e-6),
        # Forward product-inhibition Ki_NADPH for the literature-anchored deploy
        # (~15 µM; kinetics-literature band 9–24 µM). Distinct from the conflated-reverse
        # value the FIT reads (~2.3 µM); see `conflated` below and 3.4b classification.
        :Ki_NADPH     => log10(15e-6),
    ),
    Dict{Symbol,Vector{Symbol}}(
        :SS_NADPH_release_rate_eq => [:Ki_6PGLn, :Ki_ATP, :Km_NADPH_rev],
        :RE_rate_eq               => [:Ki_6PGLn, :Ki_ATP],
    ),
    _G6PD_KI_MAP, _G6PD_KD_MAP, _G6PD_KM_RATIO, _G6PD_SUBSTRATE_PAIR;
    # Forward Ki_NADPH is now read from the noncompetitive dead-end CROSS term on E·G6P
    # (g[G6P]/g[G6P·NADPH]), distinct from the bare [NADPH] reverse-release Km. This
    # de-conflates the forward Ki from the reverse Km, so it is no longer reported as
    # :conflated_reverse. See docs/G6PD_session_context.md and the forward-Ki cross-term plan.
    ki_ratio = _G6PD_KI_RATIO,
))

# Back-compat: tests reference `FitRateEquation.LIT_VALUES[:Km_NADPH_rev]` directly.
const LIT_VALUES = _wiring(:G6PD).lit_values

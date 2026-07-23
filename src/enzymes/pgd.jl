# PGD consensus wiring — random-binding / ordered-release Bi-Ter.
# Substrate binding random (NADP, PGA); product release ordered CO2 -> Ru5P -> NADPH
# (NADPH last: bound first in reverse on free E => competitive vs NADP, and forward
# Ki_NADPH ≡ reverse Km_NADPH in V1, separated in V2). ATP competitive vs PGA =>
# binds the PGA-free forms [E, E_N]. No product dead-ends (Problem 5 / Sessions 41-42).

const _PGD_ATP_FORMS = [:E, :E_N]

# `ru5p_op`/`nadph_op` toggle the Ru5P and NADPH release steps RE(⇌) vs SS(<-->).
# V1 = both RE. V2 = both SS: NADPH-release SS de-conflates forward Ki_NADPH from reverse
# Km_NADPH, and Ru5P-release SS is REQUIRED for a working apparent-constant gauge — the
# 3-product ordered release (CO2->Ru5P->NADPH) otherwise leaves two RE steps between the
# catalysis-SS gauge anchor and the NADPH-SS step (SS->RE->RE->SS), which has NO free-enzyme
# denominator term (gauged_denominator_coeffs needs one). Making Ru5P-release SS too
# (SS->RE->SS->SS) restores the free-E term + a finite rate, and is biochemically correct
# (Topham 1986: E·Ru5P isomerization/release is rate-limiting). The Ru5P-release reverse
# constant is a nuisance DOF (unidentifiable on forward data; reported unconstrained, not
# pinned). See the spec's "Why V2 also makes Ru5P-release SS (3-product gauge constraint)".
function _pgd_consensus(ru5p_op::Symbol, nadph_op::Symbol)
    extra = _deadends([(_PGD_ATP_FORMS, :ATP)])
    _mech([:NADP, :PGA], [:Ru5P, :CO2, :NADPH], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),               # NADP binds free E   ┐ random
        ([:E, :PGA],    [:E_G],  :(⇌)),               # PGA binds free E    │ substrate
        ([:E_N, :PGA],  [:E_NB], :(⇌)),               # PGA binds E·NADP    │ binding
        ([:E_G, :NADP], [:E_NB], :(⇌)),               # NADP binds E·PGA    ┘
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_1, :CO2],   :(⇌)),               # CO2 release  (RE, first out)
        ([:E_1], [:E_2, :Ru5P],  ru5p_op),            # Ru5P release: ⇌ (V1) / <--> (V2)
        ([:E_2], [:E, :NADPH],   nadph_op),           # NADPH release: ⇌ (V1) / <--> (V2), last
    ], extra); regs=[:ATP])
end

# V3 forward de-conflation: NADPH binds the E·PGA form [:E_G] as a rapid-equilibrium
# DEAD-END (its own micro-K, pendant node). On E·PGA ⇒ noncompetitive vs 6PG
# (Cottreau/Weisz). The free-E dead-end (E·NADPH) is DROPPED (2026-06-10, mirroring the
# G6PD option-B change): it is the REVERSE of the productive NADPH release and lumps the
# bare-[NADPH] competitive coeff into a 2-knob term (release + dead-end). With only the
# E·PGA dead-end, the bare [NADPH] term is a CLEAN reverse Km and the [PGA·NADPH] cross
# term is the CLEAN forward-Ki slot (read as the ki-ratio g[PGA]/g[PGA·NADPH], see
# _PGD_KI_RATIO). Built on the RE base (no SS steps) so the free-E gauge term survives
# and the REC-4 freeze holds.
const _PGD_NADPH_FWD_FORMS = [:E_G]

function _pgd_v3()
    extra = _deadends([(_PGD_ATP_FORMS, :ATP), (_PGD_NADPH_FWD_FORMS, :NADPH)])
    _mech([:NADP, :PGA], [:Ru5P, :CO2, :NADPH], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :PGA],    [:E_G],  :(⇌)),
        ([:E_N, :PGA],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),
        ([:E_C], [:E_1, :CO2],   :(⇌)),
        ([:E_1], [:E_2, :Ru5P],  :(⇌)),
        ([:E_2], [:E, :NADPH],   :(⇌)),
    ], extra); regs=[:ATP])
end

# Cha base: Topham asymmetric-random Bi-Ter with the ONE gauge-mandated SS release
# (Ru5P, Topham rate-limiting) + the V3 E·PGA NADPH dead-end (clean fwd-Ki/rev-Km channel).
# CO2 + NADPH release stay RE (PGD NADPH release is fast >800/s — base-only, no promotion).
function _pgd_cha_base()
    extra = _deadends([(_PGD_ATP_FORMS, :ATP), (_PGD_NADPH_FWD_FORMS, :NADPH)])
    _mech([:NADP, :PGA], [:Ru5P, :CO2, :NADPH], vcat([
        ([:E, :NADP],   [:E_N],  :(⇌)),
        ([:E, :PGA],    [:E_G],  :(⇌)),
        ([:E_N, :PGA],  [:E_NB], :(⇌)),
        ([:E_G, :NADP], [:E_NB], :(⇌)),
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_1, :CO2],   :(⇌)),               # CO2 release  (RE, first out)
        ([:E_1], [:E_2, :Ru5P],  :(<-->)),            # Ru5P release (SS, Topham rate-limiting)
        ([:E_2], [:E, :NADPH],   :(⇌)),               # NADPH release (RE, last)
    ], extra); regs=[:ATP])
end

# Fully-RE CORE (:full_re): V1's topology (random RE substrate binding, SS catalysis gauge,
# ordered RE product release CO2->Ru5P->NADPH) with the ATP dead-end effectors OFF (the
# config-gated default of the fully-RE law). free_params = the 6 RE dissociation constants only
# (no K_ATPinh_*), so the :full_re fit's 6 core coords map bijectively (cha_deploy_micro branch).
function _pgd_fullre_core()
    _mech([:NADP, :PGA], [:Ru5P, :CO2, :NADPH], [
        ([:E, :NADP],   [:E_N],  :(⇌)),               # NADP binds free E   ┐ random
        ([:E, :PGA],    [:E_G],  :(⇌)),               # PGA binds free E    │ substrate
        ([:E_N, :PGA],  [:E_NB], :(⇌)),               # PGA binds E·NADP    │ binding
        ([:E_G, :NADP], [:E_NB], :(⇌)),               # NADP binds E·PGA    ┘
        ([:E_NB],       [:E_C],  :(<-->)),            # catalysis (forced SS, gauge anchor)
        ([:E_C], [:E_1, :CO2],   :(⇌)),               # CO2 release  (RE, first out)
        ([:E_1], [:E_2, :Ru5P],  :(⇌)),               # Ru5P release (RE)
        ([:E_2], [:E, :NADPH],   :(⇌)),               # NADPH release (RE, last)
    ]; regs=Symbol[])
end

const _PGD_KI_MAP = Dict{MonoKey,Symbol}(
    [:Ru5P  => 1] => :Ki_Ru5P,
    [:CO2   => 1] => :Ki_CO2,
    [:ATP   => 1] => :Ki_ATP,
)
# Forward Ki_NADPH is the noncompetitive dead-end on E·PGA: read as the cross-term ratio
# g[PGA]/g[PGA·NADPH], NOT the bare [NADPH] term (= reverse release Km). The cross
# monomial is sorted by symbol (:NADPH < :PGA) to match `_split_mono`.
const _PGD_KI_RATIO = Dict{Symbol,Tuple{MonoKey,MonoKey}}(
    :Ki_NADPH => ([:PGA => 1], [:NADPH => 1, :PGA => 1]),
)
const _PGD_KD_MAP = Dict{MonoKey,Symbol}(
    [:NADP => 1] => :Kd_NADP,
    [:PGA  => 1] => :Kd_PGA,
)
const _PGD_KM_RATIO = Dict{MonoKey,Tuple{Symbol,MonoKey}}(
    [:NADP => 1] => (:Km_NADP, [:PGA  => 1]),
    [:PGA  => 1] => (:Km_PGA,  [:NADP => 1]),
)
const _PGD_SUBSTRATE_PAIR = MonoKey([:NADP => 1, :PGA => 1])   # sorted (NADP < PGA)

##########################################################################################
# FROZEN DECISION (panel review 2026-06-07; notes/2026-06-07_pgd-ki-nadph-deconflation-
# panel-review.md): do NOT add steady-state release steps to de-conflate forward Ki_NADPH.
# Proven inert — the Villet-downweight probe showed Ki_NADPH stays ~0.66 µM with all reverse
# rows removed. The blocker is the shared E·NADPH symbol + corpus geometry, not SS-step count.
# The only de-conflation route is a DECOUPLED forward symbol (V3, +NADPH_deadend_rate_eq) or
# bench data. See test_rec4_topology_freeze.jl.
##########################################################################################
register_enzyme!(EnzymeWiring(
    :PGD, :PGA,
    [(name=:RE_rate_eq,               mech=_pgd_consensus(:(⇌), :(⇌))),
     (name=:SS_NADPH_release_rate_eq, mech=_pgd_consensus(:(<-->), :(<-->))),
     (name=Symbol("+NADPH_deadend_rate_eq"), mech=_pgd_v3()),
     (name=:cha_base, mech=_pgd_cha_base()),
     (name=:full_re,  mech=_pgd_fullre_core())],
    Dict{Symbol,Float64}(
        :Ki_ATP => log10(1.7e-3),
        :Km_PGA => log10(38e-6),   # lit low-end of the 38–80 µM band (Mode-3 override; provisional)
        # Forward product-inhibition Ki_NADPH for the literature-anchored deploy (Cottreau
        # 17 µM). On V3 the forward Ki is now read from the noncompetitive E·PGA dead-end
        # cross term (ki-ratio g[PGA]/g[PGA·NADPH]), de-conflated from the reverse Km.
        :Ki_NADPH => log10(17e-6),
    ),
    Dict{Symbol,Vector{Symbol}}(
        :RE_rate_eq               => [:Ki_ATP],
        :SS_NADPH_release_rate_eq => [:Ki_ATP],
        Symbol("+NADPH_deadend_rate_eq") => [:Ki_ATP],
        :cha_base                 => [:Ki_ATP],
    ),
    _PGD_KI_MAP, _PGD_KD_MAP, _PGD_KM_RATIO, _PGD_SUBSTRATE_PAIR;
    # Forward Ki_NADPH read from the noncompetitive E·PGA dead-end CROSS term
    # (g[PGA]/g[PGA·NADPH]) on V3, distinct from the bare [NADPH] reverse-release Km. V1/V2
    # have no [PGA·NADPH] monomial, so they carry no forward Ki_NADPH coord (handled by the
    # coord gating). This replaces the obsolete ki_micro_direct free-E dead-end read.
    ki_ratio = _PGD_KI_RATIO,
))

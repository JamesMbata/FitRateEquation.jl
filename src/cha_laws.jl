# =========================================================================================
#         Cha-form partial-equilibrium G6PD rate law (closed-form, macro-parameterized)
# =========================================================================================
#
# Ported from the A1 derivation (derivations/cha_derive_g6pd.py), corrected to the v2
# TWO-STEADY-STATE topology. The v2 consensus mechanism promotes BOTH catalysis (step 5,
# kf/kr) and NADPH-release (step 7, koffQ/konQ) to steady state, with PGLn release
# (step 6) RAPID-EQUILIBRIUM between them. E_C and E_H therefore form ONE SS "super-node"
# S = {E_C <-> E_H}, internally RE-distributed (r = PGLn/Kd_6PGLn) but separated from the
# free-E RE pool by the SS NADPH-release barrier. The derivation's numeric gate
# (Property V2) matches EnzymeRates.rate_equation for v2_mechanism() to <=1e-12.
#
# A = NADP, B = G6P, P = PGLn (6-phosphogluconolactone), Q = NADPH.
#   kf = k5f (catalysis fwd, gauge=1), kr = k5r (Haldane-dependent),
#   koffQ = k7f (NADPH release/off), konQ = k7r (rebind/on), Km_NADPH_rev = koffQ/konQ.
#
# The reverse arm is FULLY determined by (kr, koffQ, konQ, Kd_6PGLn) -- do NOT inject a
# separate Keq factor (that double-specifies Haldane). `Keq` on the macro tuple is for
# reporting only; it is unused by the rate.
#
# Both NADPH channels coexist (spec section 5.1): the B*Q dead-end (forward Ki_NADPH on
# E.G6P) AND the productive bare-Q release (Km_NADPH_rev = koffQ/konQ, from konQ*Q feeding
# the super-node). The A1 single-node collapse wrongly dropped the bare-Q term.
# =========================================================================================

module ChaLaws

export cha_rate_G6PD

# Macro tuple `m` fields consumed here:
#   Kd_NADP, Kd_G6P, Kd_6PGLn, alpha, Ki_NADPH, Ki_ATP, Ki_ATP_EG, koffQ, konQ, kf, kr, Et
# (`Km_NADPH_rev` optional; if `konQ` is absent it is derived as koffQ/Km_NADPH_rev.
#  `Keq` optional, ignored.)  Concentrations are keyword args (default 0 => species absent).
# Ki_ATP is the ATP dead-end on FREE E (K8); Ki_ATP_EG is the DISTINCT ATP dead-end on
# E.G6P (K9) -- the two are independent binding constants in v2 (g6pd.jl forms [:E,:E_G]).
function cha_rate_G6PD(m; NADP=0.0, G6P=0.0, NADPH=0.0, PGLn=0.0, ATP=0.0, BPG=0.0, PGA=0.0)
    A = NADP; B = G6P; Q = NADPH; P = PGLn

    Kd_NADP   = m.Kd_NADP
    Kd_G6P    = m.Kd_G6P
    Kd_6PGLn  = m.Kd_6PGLn
    alpha     = m.alpha
    Ki_NADPH  = m.Ki_NADPH
    Ki_ATP    = m.Ki_ATP
    Ki_ATP_EG = m.Ki_ATP_EG
    koffQ     = m.koffQ
    kf        = m.kf
    kr        = m.kr
    Et        = m.Et
    konQ      = hasproperty(m, :konQ) ? m.konQ : koffQ / m.Km_NADPH_rev
    # Optional 2,3-BPG effector (silent-variants Probe B; inert unless macro carries Ki_BPG
    # AND BPG>0). 2,3-BPG is COMPETITIVE vs G6P (Özer 2001) -> free-E dead-end E.BPG, so the
    # term is +BPG/Ki_BPG on the free-E pool (NOT the noncompetitive (B/Kd_G6P)(BPG/Ki) form;
    # competitive-vs-G6P raises Km_G6P, not a ternary dead-end). Conservative lower bound on
    # the effect vs a noncompetitive reading.
    Ki_BPG    = hasproperty(m, :Ki_BPG) ? m.Ki_BPG : Inf
    # Optional 6-phosphogluconate (PGA) effector (Probe C; inert unless macro carries Ki_PGA
    # AND PGA>0). PGA is COMPETITIVE vs NADP⁺ (Ulusu 2012) -> binds the NADP⁺ subsite on free E
    # (the "1") and on E·G6P (the "B/Kd_G6P") -> term (PGA/Ki_PGA)(1 + B/Kd_G6P); raises Km_NADP.
    Ki_PGA    = hasproperty(m, :Ki_PGA) ? m.Ki_PGA : Inf

    # Forward-competent ternary fraction (relative to free E).
    gAB = A * B / (alpha * Kd_NADP * Kd_G6P)

    # Free-enzyme RE pool (constant free-E term = 1). NO bare-P term (6PGL is on the
    # super-node side, not free E).
    D_pool = 1.0 +
             A / Kd_NADP +
             B / Kd_G6P +
             gAB +
             (B / Kd_G6P) * (Q / Ki_NADPH) +    # NADPH dead-end on E.G6P (forward Ki)
             ATP / Ki_ATP +                     # ATP dead-end on free E (K8)
             (B / Kd_G6P) * (ATP / Ki_ATP_EG) + # ATP dead-end on E.G6P (K9, distinct)
             BPG / Ki_BPG +                     # 2,3-BPG dead-end on free E (comp. vs G6P)
             (PGA / Ki_PGA) * (1.0 + B / Kd_G6P) # 6PG dead-end on NADP⁺ subsite (comp. vs NADP⁺)

    # Super-node S = {E_C <-> E_H}, r = P/Kd_6PGLn; SS occupancy and net flux (per free E).
    r     = P / Kd_6PGLn
    drain = kr * r + koffQ
    S_over_E = (kf * gAB + konQ * Q) * (1 + r) / drain
    v_over_E = (kf * koffQ * gAB - kr * konQ * r * Q) / drain

    # rate = Et * (v/E) / (Et/E);  Et/E = D_pool + S/E.
    return Et * v_over_E / (D_pool + S_over_E)
end

export cha_rate_PGD
export cha_rate_PGD_fullRE

# =========================================================================================
#         Cha-form partial-equilibrium PGD rate law (Topham Bi-Ter, base-only)
# =========================================================================================
#
# A = NADP, B = PGA (6-phosphogluconate); products released ORDERED CO2 -> Ru5P -> NADPH.
# Random RE substrate binding + ordered product release. The Cha BASE promotes exactly ONE
# release to steady state besides catalysis: Ru5P release (Topham 1986, rate-limiting). CO2
# and NADPH release stay RAPID-EQUILIBRIUM. There is NO promoted silent fiber (base-only).
#
# Topology (mirrors the G6PD two-SS-segment super-node, derivations/cha_derive_pgd.py):
#   catalysis E_NB --> E_C   SS  (kf = k5f gauge, kr = k5r Haldane-dep),
#   E_C <-> E_1 + CO2        RE  (Kd_CO2 = K6, the CO2-release equilibrium),
#   E_1 --> E_2 + Ru5P       SS  (koff = k7f release, kon = k7r rebind),  the promoted step,
#   E_2 <-> E + NADPH        RE  (Km_NADPH_rev = K8, the NADPH-release equilibrium).
# The SS super-node is S = {E_C <-> E_1}, RE-linked by CO2 (s = CO2/Kd_CO2). E_2 (the NADPH-
# bound free form) is in the FREE-E RE pool, so the bare-[NADPH] denominator term is the
# reverse release Km (= K8) and the Ru5P-rebind reverse arm enters via the SS kon*Ru5P source.
#
# Both NADPH channels are distinct symbols: the [PGA*NADPH] cross-term dead-end (forward
# Ki_NADPH on E.PGA, K11) and the bare-[NADPH] productive reverse release (Km_NADPH_rev = K8).
# The reverse arm is FULLY determined by (kr, koff, kon, Kd_CO2, Km_NADPH_rev); do NOT inject
# a separate Keq factor. `Keq` on the macro tuple is for reporting only; it is unused.
#
# Macro tuple `m` fields consumed here:
#   Kd_NADP, Kd_PGA, alpha, Kd_CO2, Km_NADPH_rev, Ki_NADPH, Ki_ATP, Ki_ATP_EN,
#   koff, kon, kf, kr, Et.  Concentrations are keyword args (default 0 => species absent).
function cha_rate_PGD(m; NADP=0.0, PGA=0.0, Ru5P=0.0, CO2=0.0, NADPH=0.0, ATP=0.0)
    A = NADP; B = PGA; Q = NADPH

    Kd_NADP      = m.Kd_NADP
    Kd_PGA       = m.Kd_PGA
    alpha        = m.alpha
    Kd_CO2       = m.Kd_CO2
    Km_NADPH_rev = m.Km_NADPH_rev    # = K8 (bare-[NADPH] reverse release Km), the law's KdQ
    Ki_NADPH     = m.Ki_NADPH        # = K11 ([PGA*NADPH] dead-end cross term, forward Ki)
    Ki_ATP       = m.Ki_ATP          # = K9  (ATP dead-end on free E)
    Ki_ATP_EN    = m.Ki_ATP_EN       # = K10 (ATP dead-end on E.NADP, distinct from Ki_ATP)
    koff         = m.koff            # = k7f (Ru5P release/off)
    kon          = m.kon             # = k7r (Ru5P rebind/on)
    kf           = m.kf              # = k5f (catalysis fwd, gauge)
    kr           = m.kr              # = k5r (catalysis rev, Haldane-dependent)
    Et           = m.Et

    # Forward-competent ternary fraction (relative to free E).
    gAB = A * B / (alpha * Kd_NADP * Kd_PGA)

    # Free-enzyme RE pool (constant free-E term = 1). E_2 (NADPH-bound free form) carries the
    # bare-[NADPH] term (1/Km_NADPH_rev); the [PGA*NADPH] dead-end carries the forward Ki.
    D_pool = 1.0 +
             A / Kd_NADP +
             B / Kd_PGA +
             gAB +
             Q / Km_NADPH_rev +                 # E_2 (bare [NADPH], reverse release Km = K8)
             (B / Kd_PGA) * (Q / Ki_NADPH) +    # NADPH dead-end on E.PGA (forward Ki, K11)
             ATP / Ki_ATP +                     # ATP dead-end on free E (K9)
             (A / Kd_NADP) * (ATP / Ki_ATP_EN)  # ATP dead-end on E.NADP (K10, distinct)

    # Super-node S = {E_C <-> E_1}, RE-linked by CO2 (s = CO2/Kd_CO2). SS occupancy and net
    # flux per free E. Sources: catalysis kf*gAB into E_C, Ru5P-rebind kon*Ru5P*(Q/Km_NADPH_rev)
    # into E_1 (E_2 = Q/Km_NADPH_rev is the NADPH-bound free form). Drains: kr from E_C, koff
    # from E_1.
    s     = CO2 / Kd_CO2
    drain = kr * s + koff
    E1    = (kf * gAB + kon * Ru5P * (Q / Km_NADPH_rev)) / drain
    S_over_E = E1 * (1.0 + s)
    v_over_E = kf * gAB - kr * s * E1

    # rate = Et * (v/E) / (Et/E);  Et/E = D_pool + S/E.
    return Et * v_over_E / (D_pool + S_over_E)
end

# =========================================================================================
#     Cha-form FULLY-RAPID-EQUILIBRIUM PGD rate law (only catalysis is steady-state)
# =========================================================================================
#
# Random RE substrate binding + ordered product release CO2 -> Ru5P -> NADPH, with the
# SINGLE steady-state step being catalysis (the gauge). Unlike `cha_rate_PGD` (cha_base),
# there is NO promoted SS-release fiber (no koff/kon, no super-node): every product release
# is a rapid-equilibrium dissociation. NADPH is released LAST, so on the reverse/rebind side
# it binds FREE E first => the bare-[NADPH] term is a COMPETITIVE product-inhibition term
# (Kd_NADPH), which also serves as the reverse-release Km (one constant, both directions).
# Ru5P binds only E·NADPH (Kd_Ru5P); CO2 binds only E·NADPH·Ru5P (Kd_CO2). The product
# central complex E·CO2·Ru5P·NADPH carries the reverse arm.
#
#   A = NADP, B = PGA (6-phosphogluconate), Q = NADPH, R = Ru5P, C = CO2.
#   kf = catalysis fwd (gauge=1), kr = catalysis rev (Haldane-determined by the caller).
#   Keq is carried for provenance only; the reverse arm is fully determined by kr.
#
# Optional config-gated dead-ends (default absent => Inf => term vanishes; the fit/readoff
# supply them only when the mechanism carries them): Ki_ATP (ATP on free E), Ki_ATP_EN (ATP
# on E·NADP), Ki_NADPH (NADPH on E·PGA, the forward cross-term). These are identical in form
# to cha_rate_PGD's dead-ends and independent of the RE-vs-SS release change.
function cha_rate_PGD_fullRE(m; NADP=0.0, PGA=0.0, Ru5P=0.0, CO2=0.0, NADPH=0.0, ATP=0.0)
    A = NADP; B = PGA; Q = NADPH; R = Ru5P; Cc = CO2

    Kd_NADP  = m.Kd_NADP
    Kd_PGA   = m.Kd_PGA
    alpha    = m.alpha
    Kd_NADPH = m.Kd_NADPH       # competitive free-E NADPH constant (= reverse-release Km)
    Kd_Ru5P  = m.Kd_Ru5P
    Kd_CO2   = m.Kd_CO2
    kf       = m.kf
    kr       = m.kr
    Et       = m.Et

    Ki_ATP    = hasproperty(m, :Ki_ATP)    ? m.Ki_ATP    : Inf   # ATP dead-end on free E
    Ki_ATP_EN = hasproperty(m, :Ki_ATP_EN) ? m.Ki_ATP_EN : Inf   # ATP dead-end on E·NADP
    Ki_NADPH  = hasproperty(m, :Ki_NADPH)  ? m.Ki_NADPH  : Inf   # NADPH dead-end on E·PGA

    # Forward-competent ternary fraction (relative to free E).
    gAB = A * B / (alpha * Kd_NADP * Kd_PGA)

    # Partition function (relative to free E). Substrate side + ordered-RE product side (nested,
    # NADPH last) + config-gated dead-ends.
    D = 1.0 +
        A / Kd_NADP +
        B / Kd_PGA +
        gAB +
        Q / Kd_NADPH +                                  # E·NADPH (competitive)
        Q * R / (Kd_NADPH * Kd_Ru5P) +                  # E·NADPH·Ru5P
        Q * R * Cc / (Kd_NADPH * Kd_Ru5P * Kd_CO2) +    # E·NADPH·Ru5P·CO2 (product central)
        ATP / Ki_ATP +                                  # ATP dead-end on free E
        (A / Kd_NADP) * (ATP / Ki_ATP_EN) +             # ATP dead-end on E·NADP
        (B / Kd_PGA) * (Q / Ki_NADPH)                   # NADPH dead-end on E·PGA (fwd Ki)

    # Net flux through the single SS catalytic step: kf*[E·NADP·PGA] - kr*[E·CO2·Ru5P·NADPH].
    num = kf * gAB - kr * Q * R * Cc / (Kd_NADPH * Kd_Ru5P * Kd_CO2)
    return Et * num / D
end

end # module ChaLaws

module ChaLawsHK1

export cha_rate_HK1

# Macro tuple `m` fields consumed here:
#   Kd_Glc, Kd_ATP, Ki_G6P_C, Ki_ADP, Ki_G6P_N, K_Pi_N, alpha, Keq
#   (optional) gamma (default 1.0), kf (default 1.0), k2f (default kf), Et (default 1.0).
# `alpha` may be a real (e.g. 1.0) or Inf / :infinity (mutual exclusion, drops the N.G6P
#  state on C-half-G6P-bearing forms). Concentrations are keyword args (default 0 => absent).
function cha_rate_HK1(m; Glucose=0.0, ATP=0.0, G6P=0.0, ADP=0.0, Pi=0.0)
    KG = m.Kd_Glc;    KA = m.Kd_ATP
    KC = m.Ki_G6P_C;  KD = m.Ki_ADP
    KN = m.Ki_G6P_N;  KP = m.K_Pi_N
    Keq = m.Keq
    γ   = hasproperty(m, :gamma) ? m.gamma : 1.0
    kf  = hasproperty(m, :kf)    ? m.kf    : 1.0
    k2f = hasproperty(m, :k2f)   ? m.k2f   : kf
    Et  = hasproperty(m, :Et)    ? m.Et    : 1.0

    # alpha -> the N.G6P-on-C-half-G6P weight 1/(alpha*KN); :infinity / Inf => 0 (exclusion).
    α = m.alpha
    invαKN = (α === :infinity || (α isa Real && isinf(α))) ? 0.0 :
             1.0 / ((α isa Real ? α : float(α)) * KN)

    g = Glucose; a = ATP; p = G6P; d = ADP; π = Pi

    # N-half 3-state factors (free / C-bearing-G6P).
    F3_noG6P   = 1.0 + p / KN     + π / KP
    F3_withG6P = 1.0 + p * invαKN + π / KP

    # C-half RE pools relative to free E (beta=1 product ternary).
    Z_noG6P   = 1.0 + g / KG + a / KA + (g * a) / (γ * KG * KA) +
                d / KD + (g * d) / (KG * KD)
    Z_withG6P = p / KC + (g * p) / (KG * KC) + (p * d) / (KC * KD)

    DEN = Z_noG6P * F3_noG6P + Z_withG6P * F3_withG6P
    NUM = (kf + k2f * π / KP) * (g * a - p * d / Keq) / (KG * KA)
    return Et * NUM / DEN
end

end # module ChaLawsHK1

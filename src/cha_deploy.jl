# =========================================================================================
#         Closed-form DEPLOY inverse: macro Cha tuple -> mechanism micro free-params
# =========================================================================================
#
# `cha_deploy_micro(enzyme, mech, coords; keq, koffQ, release_rate, release_eq, kr)` is the
# exact INVERSE of `cha_macro_readoffs_*` (cha_invert.jl). Given a Cha macro tuple (the
# data-identifiable forward shape constants `coords` plus the gauge/fiber/Haldane inputs), it
# assigns each INDEPENDENT free-param of `mech` (the symbols in `free_params(mech)`) its
# macro-derived value and returns the aligned `logθ` (log10). `build_params(mech, logθ; keq)`
# then fills the Wegscheider/Haldane-DEPENDENT params (G6PD/PGD: K3 Wegscheider, k5r Haldane)
# from those independents + keq. No optimization / LM / root-finding: the map is closed form,
# field by field. This is the deploy block `write_outputs` emits (Task 15), retiring the
# multi-start LM `macro_bridge` (spec §2/§5.3 — "blocker E cannot recur" because it is closed
# form, not a fitted inverse).
#
# Micro <-> macro correspondence (introspected, the inverse of cha_macro_readoffs_*):
#   G6PD (free_params K1,K2,K4,K6,k7f,k7r,K8,K9,K10; dependent K3,k5r):
#       K1 = Kd_NADP                  (NADP binding, g[NADP] = 1/K1)
#       K2 = Kd_G6P                   (G6P  binding, g[G6P]  = 1/K2)
#       K4 = alpha * Kd_NADP          (ternary interaction: alpha = K4/K1, verified field-exact)
#       K6 = Kd_6PGLn                 (PGLn RE-binding, the law's KdP)
#       k7f = koffQ = release_rate    (promoted NADPH release/off rate)
#       k7r = konQ  = release_rate / release_eq   (NADPH rebind/on; release_eq = Km_NADPH_rev)
#       K8  = Ki_ATP                  (ATP dead-end on free E)
#       K9  = Ki_ATP_EG               (ATP dead-end on E.G6P, distinct from Ki_ATP)
#       K10 = Ki_NADPH                (NADPH dead-end on E.G6P, the forward cross-term Ki)
#       [:no_atp variant] the two ATP dead-ends above are ABSENT (no K8/K9 ATP slot) and the
#         sole remaining NADPH dead-end occupies K8 instead -- see the inline comment in the
#         G6PD branch of `_deploy_micro_map` below.
#   PGD (free_params K1,K2,K4,K6,k7f,k7r,K8,K9,K10,K11; dependent K3,k5r):
#       K1 = Kd_NADP                  K2  = Kd_PGA
#       K4 = alpha * Kd_NADP          K6  = Kd_CO2
#       k7f = koff = release_rate     k7r = kon = release_rate / release_eq  (KdRu = koff/kon)
#       K8  = Km_NADPH_rev            (bare-[NADPH] reverse-release Km)
#       K9  = Ki_ATP    K10 = Ki_ATP_EN    K11 = Ki_NADPH
#
# Gauge: the readoff (and the macro tuple) always carries kf = k5f = 1 and Et = E_total = 1
# (unit-enzyme gauge); build_params likewise fixes k5f = 1 / E_total = 1, so the deploy micro
# shares the macro tuple's gauge and the ABSOLUTE rate matches (no shape-only rescale). The
# Haldane reverse catalysis k5r is NOT set here -- build_params derives it from the
# independents + keq. `kr` is accepted for interface symmetry with cha_macro_tuple but unused
# (build_params owns the Haldane); `koffQ` defaults to a healthy value but is normally the
# actual release-off rate passed by the caller (== release_rate).
# =========================================================================================

module ChaDeploy

using ..ChaFit
using ..ChaLaws
using ..FitRateEquation
using EnzymeRates

export cha_deploy_micro

# Build the symbol => value map for the mechanism's INDEPENDENT free-params from the macro
# tuple. `release_rate` is the promoted release off-rate (k7f / koff), `release_eq` its
# release equilibrium (koff/kon): G6PD -> Km_NADPH_rev, PGD -> KdRu. The on-rate kon = k7r is
# release_rate / release_eq.
function _deploy_micro_map(enzyme::Symbol, coords::AbstractDict; release_rate::Real,
                           release_eq::Real, mech=nothing)
    if enzyme === :HK1
        # free_params role order (Task 0): [k2f, Kd_Glc, Kd_ATP, Ki_G6P_C, Ki_ADP, K_Pi_N, Ki_G6P_N].
        # k1f is the gauge (dropped from free_params); k2f is hard-gauged to 1.0 (Pi competitor-only).
        fp = free_params(mech)
        @assert length(fp) == 7 "HK1 deploy expects 7 free params (k1f gauge-dropped), got $(length(fp)): $fp"
        @assert fp[1] === :k2f "HK1 free_params[1] must be :k2f, got $(fp[1])"
        # H4 carries the reparameterized G6P coords {Keff, split_ratio}; back-map to {Kc, Kn}
        # (identical to cha_macro_tuple's H4 branch) so H4 deploys the same micro block as H1.
        if haskey(coords, :Keff)
            Keff = coords[:Keff]; ratio = coords[:split_ratio]
            sqrtP = Keff * ratio; P = sqrtP^2; sumK = P / Keff
            sq = sqrt(max(sumK^2 - 4P, 0.0))
            KiC = (sumK + sq) / 2; KiN = (sumK - sq) / 2
        else
            KiC = coords[:Ki_G6P_C]; KiN = coords[:Ki_G6P_N]
        end
        vals = (1.0, coords[:Kd_Glc], coords[:Kd_ATP], KiC,
                coords[:Ki_ADP], coords[:K_Pi_N], KiN)
        return Dict{Symbol,Float64}(fp[i] => vals[i] for i in 1:7)
    end
    konv = release_rate / release_eq
    if enzyme === :G6PD
        # Upstream semantic names (see cha_invert). Dead-ends are named by composition, so
        # there is no position shuffling: the NADPH dead-end on E·G6P (K_NADPH_EG6P) is always
        # present; the two ATP dead-ends are absent for the :no_atp variant. BASIS: upstream's
        # free binding param is G6P→E·NADP (K_G6P_ENADP), not the old NADP→E·G6P; by the
        # detailed-balance box K1·K3=K2·K4 it equals alpha·Kd_G6P (alpha = K_G6P_ENADP/Kd_G6P).
        micro = Dict{Symbol,Float64}(
            :K_NADP_E      => coords[:Kd_NADP],
            :K_G6P_E       => coords[:Kd_G6P],
            :K_G6P_ENADP   => coords[:alpha] * coords[:Kd_G6P],
            :K_PGLn_ENADPH => coords[:Kd_6PGLn],
            :koff_NADPH_E  => release_rate,
            :kon_NADPH_E   => konv,
        )
        # Each dead-end is guarded independently: :no_atp drops both ATP entries; the
        # no_g6p_*_deadend variants drop exactly one of {Ki_NADPH, Ki_ATP_EG} while keeping
        # the other ATP entry (Ki_ATP) -- these must not be coupled behind one shared guard.
        haskey(coords, :Ki_NADPH)  && (micro[:K_NADPH_EG6P]  = coords[:Ki_NADPH])
        haskey(coords, :Ki_ATP)    && (micro[:K_ATPinh_E]    = coords[:Ki_ATP])
        haskey(coords, :Ki_ATP_EG) && (micro[:K_ATPinh_EG6P] = coords[:Ki_ATP_EG])
        return micro
    elseif enzyme === :PGD
        # PGD's free basis keeps NADP→E·PGA (K_NADP_EPGA = old K4 = alpha·Kd_NADP), so unlike
        # G6PD there is no basis change.
        return Dict{Symbol,Float64}(
            :K_NADP_E         => coords[:Kd_NADP],
            :K_PGA_E          => coords[:Kd_PGA],
            :K_NADP_EPGA      => coords[:alpha] * coords[:Kd_NADP],
            :K_CO2_ENADPHRu5P => coords[:Kd_CO2],
            :koff_Ru5P_ENADPH => release_rate,
            :kon_Ru5P_ENADPH  => konv,
            :K_NADPH_E        => coords[:Km_NADPH_rev],   # bare-[NADPH] reverse-release Km (free E)
            :K_ATPinh_E       => coords[:Ki_ATP],
            :K_ATPinh_ENADP   => coords[:Ki_ATP_EN],
            :K_NADPH_EPGA     => coords[:Ki_NADPH],        # forward NADPH dead-end on E·PGA
        )
    else
        error("cha_deploy_micro: unknown enzyme $enzyme (expected :G6PD, :PGD, or :HK1)")
    end
end

# Closed-form deploy: return logθ (log10) aligned to free_params(mech), assigning each
# INDEPENDENT free-param its macro-derived value. build_params(mech, logθ; keq) fills the
# Wegscheider/Haldane-dependent params (K3, k5r). `koffQ` / `kr` are accepted for interface
# symmetry with cha_macro_tuple; `koffQ` defaults to a healthy release-off rate and is the
# promoted release rate used for k7f when `release_rate` is omitted, while `kr` is owned by
# build_params (Haldane) and unused here.
function cha_deploy_micro(enzyme::Symbol, mech, coords::AbstractDict; keq::Real,
                          koffQ::Real = 1.0,
                          release_rate::Real = koffQ,
                          release_eq::Real = ChaFit._default_release_eq(enzyme, coords),
                          kr::Union{Nothing,Real} = nothing)
    micro = _deploy_micro_map(enzyme, coords; release_rate=release_rate, release_eq=release_eq, mech=mech)
    fp = free_params(mech)
    logθ = Vector{Float64}(undef, length(fp))
    for (i, s) in enumerate(fp)
        haskey(micro, s) ||
            error("cha_deploy_micro: no macro->micro mapping for free param $s of $enzyme")
        logθ[i] = log10(micro[s])
    end
    logθ
end

end # module ChaDeploy

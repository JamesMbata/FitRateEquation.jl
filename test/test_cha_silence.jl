# =========================================================================================
#       In-Julia koffQ silence + identification guard for the Cha-form G6PD law
# =========================================================================================
#
# VERIFIED FINDING (A1-A3, derivations/cha_derive_g6pd.py Property 5;
# notes/2026-06-10_g6pd_cha_supernode_derivation.md): on the corrected two-SS-segment
# super-node topology of v2_mechanism(), the promoted NADPH-release rate koffQ is
#
#   * SILENT on the pure-Vmax / uniform-rescale fiber (v -> lambda*v cancels in any
#     gauge-invariant rate ratio), and SILENT on P=0 data (forward-only NADPH=0 AND
#     NADPH-only PGLn=0) IN THE GAUGE-INVARIANT SHAPE: sliding koffQ on the
#     macro-constant-preserving fiber `cha_micro_from_macro_G6PD` leaves the
#     gauge-invariant ratio v(point)/v(reference) invariant (relvar <= 1e-12).
#
#   * DATA-IDENTIFIED on P>0 data (PGLn-product-inhibition PGLn>0, NADPH=0; and both-product
#     PGLn>0, NADPH>0): the SAME macro-preserving fiber MOVES even the gauge-invariant SHAPE
#     (empirically ratio relvar ~0.8 to ~8). Mechanism: catalysis-reverse kr re-populates
#     E_C through the SS NADPH-release barrier, generating a PGLn-cross family + reverse P*Q
#     term that pin koffQ. THE PLAN'S ORIGINAL A4 ("koffQ silent on a both-product grid") IS
#     FALSE and is deliberately NOT asserted here.
#
# SUBTLETY (investigated, A4): the `cha_micro_from_macro_G6PD` fiber is NOT a pure-Vmax
# gauge -- it holds apparent macro constants (Kd's, alpha/C, Km_NADPH_rev) fixed but slides
# the ABSOLUTE Vmax (kf is held while C=1+kf/koffQ changes). So on P=0 data the ABSOLUTE
# rate moves (~1.36 relvar) while the gauge-invariant SHAPE v(pt)/v(ref) is silent to 1e-16.
# This matches test_cha_invert.jl's framing (fiber preserves apparent constants, not the
# Vmax). The loss quotients out the Vmax scale, so the SHAPE ratio is the observable the
# silence/identification statement is about -- hence testsets 1 and 2 assert on the ratio.
#
# Testset 2 is the regression tripwire: if a future law re-collapses to the wrong
# single-node topology, koffQ would go spuriously silent (in shape) on P>0 data and testset
# 2 would fail loudly.
# =========================================================================================

using FitRateEquation
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert

# Representative macro tuple (same style as test_cha_invert.jl read-offs). cha_rate_G6PD
# consumes: Kd_NADP, Kd_G6P, Kd_6PGLn, alpha, Ki_NADPH, Ki_ATP, Ki_ATP_EG, koffQ, konQ
# (or Km_NADPH_rev), kf, kr, Et. cha_micro_from_macro_G6PD needs Km_NADPH_rev to re-derive
# konQ on the fiber, so we supply both. Ki_ATP_EG is the distinct ATP dead-end on E.G6P.
const MAC = (; Kd_NADP=7.5e-6, Kd_G6P=45e-6, Kd_6PGLn=2e-4, alpha=1.0,
             Km_NADPH_rev=3.9e-6, Ki_NADPH=15e-6, Ki_ATP=1.5e-3, Ki_ATP_EG=2.0e-3,
             koffQ=1e3, kf=200.0, kr=50.0, Keq=9.9, Et=1.0)

# konQ on MAC itself so it is a self-consistent law-ready tuple (= koffQ/Km_NADPH_rev).
const MAC_LAW = merge(MAC, (; konQ = MAC.koffQ / MAC.Km_NADPH_rev))

const KOFFQ_SWEEP = (1e1, 1e2, 1e3, 1e4, 1e5, 1e6)

relvar(vs) = (maximum(vs) - minimum(vs)) / abs(sum(vs) / length(vs))

@testset "koffQ silent on P=0 data (forward-only + NADPH-only)" begin
    # All points have PGLn=0.0 (P=0): mix forward-only (NADPH=0) and NADPH-only (NADPH>0),
    # varying NADP / G6P / ATP. The fiber slides the absolute Vmax, so the OBSERVABLE we
    # assert silence on is the gauge-invariant shape v(point)/v(reference) (the Vmax scale
    # is quotiented by the loss). It must be koffQ-invariant on the macro-preserving fiber.
    pref = (; NADP=5e-6, G6P=40e-6, NADPH=0.0, PGLn=0.0, ATP=0.0)       # P=0 reference
    pts = [(; NADP=2e-5,  G6P=80e-6, NADPH=0.0,  PGLn=0.0, ATP=1e-3),  # forward-only + ATP
           (; NADP=1e-6,  G6P=10e-6, NADPH=0.0,  PGLn=0.0, ATP=0.0),   # forward-only low
           (; NADP=5e-6,  G6P=40e-6, NADPH=5e-6, PGLn=0.0, ATP=0.0),   # NADPH-only
           (; NADP=2e-5,  G6P=80e-6, NADPH=2e-5, PGLn=0.0, ATP=1e-3),  # NADPH-only + ATP
           (; NADP=1e-6,  G6P=10e-6, NADPH=1e-5, PGLn=0.0, ATP=0.0)]   # NADPH-only low
    maxrv = 0.0
    for p in pts
        rs = [cha_rate_G6PD(cha_micro_from_macro_G6PD(MAC_LAW; koffQ=kq); p...) /
              cha_rate_G6PD(cha_micro_from_macro_G6PD(MAC_LAW; koffQ=kq); pref...)
              for kq in KOFFQ_SWEEP]
        rv = relvar(rs)
        maxrv = max(maxrv, rv)
        @test rv <= 1e-12
    end
    @info "P=0 shape silence: max ratio relvar over grid = $maxrv (expect <= 1e-12)"
end

@testset "koffQ data-identified on P>0 data (PGLn-inhibition + both-product)" begin
    # All points have PGLn>0 (P>0): mix PGLn-only (NADPH=0) and both-product (NADPH>0).
    # The SAME fiber MUST move even the gauge-invariant SHAPE v(point)/v(reference) (where
    # reference is a P=0 forward-only point). 1e-2 is a safe floor (empirical ~0.8-8) that
    # still fails loudly if the topology ever regresses to spurious silence.
    pref = (; NADP=5e-6, G6P=40e-6, NADPH=0.0, PGLn=0.0, ATP=0.0)       # P=0 reference
    pts = [(; NADP=5e-6,  G6P=40e-6, NADPH=0.0,  PGLn=1e-4, ATP=0.0),   # PGLn-only
           (; NADP=2e-5,  G6P=80e-6, NADPH=0.0,  PGLn=2e-4, ATP=1e-3),  # PGLn-only + ATP
           (; NADP=5e-6,  G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0),   # both products
           (; NADP=2e-5,  G6P=80e-6, NADPH=2e-5, PGLn=2e-4, ATP=1e-3)]  # both products + ATP
    minrv = Inf
    for p in pts
        rs = [cha_rate_G6PD(cha_micro_from_macro_G6PD(MAC_LAW; koffQ=kq); p...) /
              cha_rate_G6PD(cha_micro_from_macro_G6PD(MAC_LAW; koffQ=kq); pref...)
              for kq in KOFFQ_SWEEP]
        rv = relvar(rs)
        minrv = min(minrv, rv)
        @test rv >= 1e-2
    end
    @info "P>0 shape identification: min ratio relvar over grid = $minrv (expect >= 1e-2)"
end

@testset "koffQ silent on the Vmax/uniform-rescale fiber" begin
    # Pure-Vmax rescale: scale ALL rate constants (kf, kr, koffQ, konQ) by lambda, so
    # v -> lambda*v. The gauge-invariant ratio v(point)/v(reference) cancels lambda and
    # must be lambda-invariant to <= 1e-12 on a both-product grid. cha_rate_G6PD reads
    # konQ directly, so we scale it explicitly in the rescaled tuple.
    base = cha_micro_from_macro_G6PD(MAC_LAW; koffQ=1e3)
    pref = (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0)   # reference point
    pts  = [(; NADP=5e-6,  G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0),
            (; NADP=2e-5,  G6P=80e-6, NADPH=2e-5, PGLn=2e-4, ATP=1e-3),
            (; NADP=1e-6,  G6P=10e-6, NADPH=1e-5, PGLn=5e-5, ATP=0.0)]
    # reference ratios at lambda = 1.
    ratios_ref = [cha_rate_G6PD(base; p...) / cha_rate_G6PD(base; pref...) for p in pts]
    maxdev = 0.0
    for lam in (0.1, 1.0, 10.0, 100.0)
        scaled = merge(base, (; kf = lam*base.kf, kr = lam*base.kr,
                                koffQ = lam*base.koffQ, konQ = lam*base.konQ))
        vp_ref = cha_rate_G6PD(scaled; pref...)
        for (i, p) in enumerate(pts)
            ratio = cha_rate_G6PD(scaled; p...) / vp_ref
            dev = abs(ratio - ratios_ref[i]) / abs(ratios_ref[i])
            maxdev = max(maxdev, dev)
            @test isapprox(ratio, ratios_ref[i]; rtol=1e-12)
        end
    end
    @info "Vmax-fiber silence: max gauge-ratio deviation = $maxdev (expect <= 1e-12)"
end

using FitRateEquation
using Test

@testset "mode set + deploy-variant cells are per-enzyme and per-variant" begin
    @test FitRateEquation.modes_for(:G6PD) == (:mode1, :mode2)
    # PGD: Mode 3 was retired (2026-10-21). A variant gets Mode 2 only when it carries a
    # forward Ki_NADPH coord to pin; :full_re has none, so its Mode 2 would duplicate Mode 1.
    @test FitRateEquation.modes_for(:PGD)            == (:mode1, :mode2)
    @test FitRateEquation.modes_for(:PGD, :cha_base) == (:mode1, :mode2)
    @test FitRateEquation.modes_for(:PGD, :no_atp)   == (:mode1, :mode2)
    @test FitRateEquation.modes_for(:PGD, :full_re)  == (:mode1,)
    g = FitRateEquation._cells(:G6PD); p = FitRateEquation._cells(:PGD)
    # Deploy-variant-only: one variant per enzyme, one cell per mode.
    @test length(g) == 2 && length(p) == 2
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:G6PD), g)
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:PGD), p)
    @test [c[3] for c in g] == [:mode1, :mode2]
    @test [c[3] for c in p] == [:mode1, :mode2]
    @test [c[3] for c in FitRateEquation._cells(:PGD; variants=[:full_re])] == [:mode1]
end

@testset "PGD Mode 2 pins ONLY the forward Ki_NADPH (no Ki_ATP pin, no Km_PGA pull)" begin
    rp = FitRateEquation.ChaFit.resolve_cha_pins
    @test isempty(rp(:PGD, :cha_base, :mode1))
    p2 = rp(:PGD, :cha_base, :mode2)
    @test Set(keys(p2)) == Set([:Ki_NADPH])
    @test isapprox(p2[:Ki_NADPH], log10(17e-6); atol=1e-9)
    # G6PD Mode 2 is unchanged: reverse anchor + Ki_ATP + Ki_NADPH.
    @test Set(keys(rp(:G6PD, :SS_NADPH_release_rate_eq, :mode2))) ==
          Set([:Km_NADPH_rev, :Ki_ATP, :Ki_NADPH])
    # Mode 3 no longer exists for either enzyme.
    @test_throws ErrorException rp(:PGD, :cha_base, :mode3)
    @test_throws ErrorException rp(:G6PD, :SS_NADPH_release_rate_eq, :mode3)
end

@testset "PGD Mode-2 fit == a fit with only Ki_NADPH pinned (no Km_PGA penalty)" begin
    d    = load_dataset(FitRateEquation.pgd_config())
    mech = FitRateEquation._deploy_mech(:PGD)
    r2   = FitRateEquation._fit_and_cv(:cha_base, mech, d; mode=:mode2,
                                       n_restarts=1, maxiter=40, maxtime=10.0, seed=3)
    ref  = FitRateEquation.ChaFit.cha_fit_candidate(:PGD, mech, d; n_restarts=1, maxiter=40,
               maxtime=10.0, seed=3, pins=Dict(:Ki_NADPH => log10(17e-6)), variant=:cha_base)
    @test r2.fit.loss == ref.loss
    @test r2.fit.coords == ref.coords
end

@testset "anchor_reverse=false frees the G6PD Km_NADPH_rev anchor (all modes)" begin
    rp = FitRateEquation.ChaFit.resolve_cha_pins
    # DEFAULT (anchor_reverse=true): the reverse channel is anchored in EVERY mode, and Mode 1
    # pins ONLY that anchor. This is the deployed-law de-conflation contract — must be unchanged.
    for variant in (:RE_rate_eq, :SS_NADPH_release_rate_eq)
        @test Set(keys(rp(:G6PD, variant, :mode1))) == Set([:Km_NADPH_rev])
        @test Set(keys(rp(:G6PD, variant, :mode2))) == Set([:Km_NADPH_rev, :Ki_ATP, :Ki_NADPH])
        # anchor_reverse=false: Km_NADPH_rev is left FREE → the conflation returns.
        @test isempty(rp(:G6PD, variant, :mode1; anchor_reverse=false))
        @test Set(keys(rp(:G6PD, variant, :mode2; anchor_reverse=false))) == Set([:Ki_ATP, :Ki_NADPH])
    end
    # The ATP-free variant carries the same reverse anchor and responds to the flag identically.
    @test Set(keys(rp(:G6PD, :no_atp, :mode1))) == Set([:Km_NADPH_rev])
    @test isempty(rp(:G6PD, :no_atp, :mode1; anchor_reverse=false))
    @test !haskey(rp(:G6PD, :no_atp, :mode2; anchor_reverse=false), :Km_NADPH_rev)

    # The 3 dead-end-dropped variants carry the same all-modes reverse anchor and respond to
    # the flag identically; each also auto-skips pinning whichever Ki it structurally lacks
    # (:no_g6p_nadph_deadend has no :Ki_NADPH coord, :no_g6p_atp_deadend/:no_g6p_both_deadends
    # have no :Ki_ATP_EG coord — Ki_ATP_EG is never a literature pin target regardless).
    for (variant, has_nadph, has_atp) in [
        (:no_g6p_nadph_deadend, false, true),
        (:no_g6p_atp_deadend,   true,  true),
        (:no_g6p_both_deadends, false, true),
    ]
        @test Set(keys(rp(:G6PD, variant, :mode1))) == Set([:Km_NADPH_rev])
        @test isempty(rp(:G6PD, variant, :mode1; anchor_reverse=false))
        m2 = rp(:G6PD, variant, :mode2)
        @test (:Ki_NADPH in keys(m2)) == has_nadph
        @test (:Ki_ATP in keys(m2))   == has_atp
        @test :Ki_ATP_EG ∉ keys(m2)
        @test !haskey(rp(:G6PD, variant, :mode2; anchor_reverse=false), :Km_NADPH_rev)
    end

    # PGD has no always-on reverse anchor, so the flag is a no-op for it.
    pairs = [(:PGD, :cha_base, :mode1), (:PGD, :cha_base, :mode2)]
    for (enz, variant, mode) in pairs
        @test rp(enz, variant, mode; anchor_reverse=false) == rp(enz, variant, mode)
    end
end

@testset "_default_anchor_reverse: ablation variants default false, everything else true" begin
    dar = FitRateEquation._default_anchor_reverse
    # Anchor-optional ablation variants (alone or together): default false.
    @test dar(:G6PD, [:no_g6p_atp_deadend]) == false
    @test dar(:G6PD, [:no_g6p_nadph_deadend]) == false
    @test dar(:G6PD, [:no_g6p_both_deadends]) == false
    @test dar(:G6PD, [:no_g6p_nadph_deadend, :no_g6p_both_deadends]) == false
    # The deploy variant, the raw conflating RE law, and :no_atp are unaffected (still true) --
    # :no_atp is untested for anchor-off joint identifiability, so it keeps the safe default.
    @test dar(:G6PD, [:SS_NADPH_release_rate_eq]) == true
    @test dar(:G6PD, [:RE_rate_eq]) == true
    @test dar(:G6PD, [:no_atp]) == true
    # Mixing an ablation variant with one that still requires the anchor: conservative true.
    @test dar(:G6PD, [:no_g6p_atp_deadend, :SS_NADPH_release_rate_eq]) == true
    # Empty variant list: conservative true.
    @test dar(:G6PD, Symbol[]) == true
    # PGD has no anchor-optional variants at all: always true.
    @test dar(:PGD, [:cha_base]) == true
end

@testset "_requires_reverse_anchor: NOT DEPLOYABLE banner scope matches the default's scope" begin
    rra = FitRateEquation._requires_reverse_anchor
    @test rra(:G6PD, :SS_NADPH_release_rate_eq) == true
    @test rra(:G6PD, :RE_rate_eq) == true
    @test rra(:G6PD, :no_atp) == true
    @test rra(:G6PD, :no_g6p_nadph_deadend) == false
    @test rra(:G6PD, :no_g6p_atp_deadend) == false
    @test rra(:G6PD, :no_g6p_both_deadends) == false
    @test rra(:PGD, :cha_base) == true
end

using FitRateEquation
using Test

@testset "mode set + deploy-variant cells are per-enzyme" begin
    @test FitRateEquation.modes_for(:G6PD) == (:mode1, :mode2)
    @test FitRateEquation.modes_for(:PGD)  == (:mode1, :mode2, :mode3)
    g = FitRateEquation._cells(:G6PD); p = FitRateEquation._cells(:PGD)
    # Deploy-variant-only: one variant per enzyme, one cell per mode.
    @test length(g) == 2 && length(p) == 3
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:G6PD), g)
    @test all(c -> c[1] === FitRateEquation.deploy_variant(:PGD), p)
    @test all(c -> c[3] in (:mode1, :mode2), g)
    @test any(c -> c[3] === :mode3, p)
    @test !any(c -> c[3] === :mode3, g)   # G6PD has no Mode 3
end

@testset "PGD Mode 3 overrides Km_PGA via a HIGH-weight anchor (not a coord pin)" begin
    # Km_PGA is the DERIVED apparent constant alpha·Kd_PGA/C, NOT a cha_coord, so the Mode-3
    # override is realized as a high-weight ChaFit anchor — never a hard coord pin.
    pins3 = FitRateEquation.ChaFit.resolve_cha_pins(:PGD, :cha_base, :mode3)
    @test !haskey(pins3, :Km_PGA)
    a3 = FitRateEquation.cha_anchors(:PGD, :mode3)
    @test a3 !== nothing && haskey(a3, :Km_PGA)
    @test a3[:Km_PGA].weight == 100.0
    @test isapprox(a3[:Km_PGA].target, log10(38e-6); atol=1e-9)
    # Mode 2 has a SOFT (weight-1) Km_PGA anchor toward the band midpoint (59 µM).
    a2 = FitRateEquation.cha_anchors(:PGD, :mode2)
    @test a2 !== nothing && a2[:Km_PGA].weight == 1.0
    @test isapprox(a2[:Km_PGA].target, log10(59e-6); atol=1e-9)
    # Mode 1 is unanchored; G6PD is unanchored in every mode.
    @test FitRateEquation.cha_anchors(:PGD, :mode1) === nothing
    @test FitRateEquation.cha_anchors(:G6PD, :mode2) === nothing
    # Ki_ATP is a hard Cha coord pin in Mode 2 (and 3); Mode 1 pins only the reverse channel.
    @test haskey(FitRateEquation.ChaFit.resolve_cha_pins(:PGD, :cha_base, :mode2), :Ki_ATP)
    @test !haskey(FitRateEquation.ChaFit.resolve_cha_pins(:PGD, :cha_base, :mode1), :Ki_ATP)
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

    # PGD/HK1 have no always-on reverse anchor, so the flag is a no-op for them.
    pairs = [(:PGD, :cha_base, :mode1), (:PGD, :cha_base, :mode2)]
    FitRateEquation.HK1_AVAILABLE && append!(pairs, [(:HK1, :H1, :mode2), (:HK1, :H1, :mode3)])
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
    # PGD/HK1 have no anchor-optional variants at all: always true.
    @test dar(:PGD, [:cha_base]) == true
    @test dar(:HK1, [:H1]) == true
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

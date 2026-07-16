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

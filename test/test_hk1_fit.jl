using FitRateEquation
using EnzymeRates
using Test
using FitRateEquation.ChaFit

if FitRateEquation.HK1_AVAILABLE

@testset "cha_coords(:HK1) is the six forward shape constants (H1 default)" begin
    @test ChaFit.cha_coords(:HK1) ==
          [:Kd_Glc, :Kd_ATP, :Ki_G6P_C, :Ki_ADP, :Ki_G6P_N, :K_Pi_N]
    @test ChaFit.cha_coords(:HK1, :H1) == ChaFit.cha_coords(:HK1)
end

@testset "cha_coords(:HK1, :H4) is the reparameterized {Keff, split_ratio} set" begin
    @test ChaFit.cha_coords(:HK1, :H4) ==
          [:Kd_Glc, :Kd_ATP, :Keff, :Ki_ADP, :split_ratio, :K_Pi_N]
end

@testset "cha_coord_bounds aligns to coords; H4 split_ratio has the real-roots floor" begin
    lo, hi = ChaFit.cha_coord_bounds(:HK1)
    @test length(lo) == length(ChaFit.cha_coords(:HK1))
    @test all(lo .== -9.0) && all(hi .== 0.0)        # H1: no :alpha/:split_ratio coord
    lo4, hi4 = ChaFit.cha_coord_bounds(:HK1, :H4)
    i = findfirst(==(:split_ratio), ChaFit.cha_coords(:HK1, :H4))
    @test lo4[i] == log10(2.0) && hi4[i] == 3.0      # √P/Keff ∈ [2, 1000]
end

@testset "per-variant alpha: H1/H4 two sites (α=1), H3 exclusion (α=∞)" begin
    @test ChaFit._hk1_variant_alpha(:H1) == 1.0
    @test ChaFit._hk1_variant_alpha(:H4) == 1.0
    @test ChaFit._hk1_variant_alpha(:H3) == Inf
    coords = Dict(:Kd_Glc=>50e-6, :Kd_ATP=>1e-3, :Ki_G6P_C=>15e-6, :Ki_ADP=>1.5e-3,
                  :Ki_G6P_N=>6.9e-6, :K_Pi_N=>750e-6)
    m1 = ChaFit.cha_macro_tuple(:HK1, coords; keq=2700.0, variant=:H1)
    m3 = ChaFit.cha_macro_tuple(:HK1, coords; keq=2700.0, variant=:H3)
    @test m1.alpha == 1.0
    @test isinf(m3.alpha)
    using FitRateEquation.ChaLawsHK1
    f(m; g6p) = cha_rate_HK1(m; Glucose=4.4e-3, ATP=8.3e-3, G6P=g6p, ADP=1.68e-3, Pi=9.2e-3)
    @test f(m1; g6p=1e-4) < f(m3; g6p=1e-4)
end

@testset "H4 back-map {Keff, split_ratio} → {Ki_G6P_C, Ki_G6P_N} is the exact inverse" begin
    # Pick a {Kc (loose), Kn (tight)} pair; reparameterize forward then back via cha_macro_tuple.
    Kc, Kn = 1432e-6, 24.3e-6
    Keff  = 1 / (1/Kc + 1/Kn)
    sqrtP = sqrt(Kc * Kn)
    coords4 = Dict(:Kd_Glc=>20e-6, :Kd_ATP=>1.1e-3, :Keff=>Keff, :Ki_ADP=>0.9e-3,
                   :split_ratio=>sqrtP/Keff, :K_Pi_N=>1.8e-3)
    m4 = ChaFit.cha_macro_tuple(:HK1, coords4; keq=2700.0, variant=:H4)
    @test m4.alpha == 1.0
    @test isapprox(m4.Ki_G6P_C, Kc; rtol=1e-9)        # convention: C = larger (loose) root
    @test isapprox(m4.Ki_G6P_N, Kn; rtol=1e-9)
    # ⇒ H4 evaluates the SAME law as H1 at {Kc,Kn} (loss-preserving relabel).
    coords1 = Dict(:Kd_Glc=>20e-6, :Kd_ATP=>1.1e-3, :Ki_G6P_C=>Kc, :Ki_ADP=>0.9e-3,
                   :Ki_G6P_N=>Kn, :K_Pi_N=>1.8e-3)
    m1 = ChaFit.cha_macro_tuple(:HK1, coords1; keq=2700.0, variant=:H1)
    using FitRateEquation.ChaLawsHK1
    g(m) = cha_rate_HK1(m; Glucose=4.4e-3, ATP=8.3e-3, G6P=1e-4, ADP=1.68e-3, Pi=9.2e-3)
    @test isapprox(g(m4), g(m1); rtol=1e-12)
end

@testset "apparent Km for HK1 equals the binary Kd (C=1, gamma=1)" begin
    coords = Dict(:Kd_Glc=>50e-6, :Kd_ATP=>1e-3, :Ki_G6P_C=>15e-6, :Ki_ADP=>1.5e-3,
                  :Ki_G6P_N=>6.9e-6, :K_Pi_N=>750e-6)
    @test ChaFit.cha_apparent_km(:HK1, coords, :Km_Glc) == 50e-6
    @test ChaFit.cha_apparent_km(:HK1, coords, :Km_ATP) == 1e-3
end

@testset "resolve_cha_pins(:HK1, …) per-mode pin sets; H4 is pin-free" begin
    @test isempty(ChaFit.resolve_cha_pins(:HK1, :H1, :mode1))
    @test Set(keys(ChaFit.resolve_cha_pins(:HK1, :H1, :mode2))) == Set([:Ki_G6P_N, :K_Pi_N])
    @test Set(keys(ChaFit.resolve_cha_pins(:HK1, :H1, :mode3))) ==
          Set([:Ki_G6P_N, :K_Pi_N, :Ki_G6P_C, :Ki_ADP])
    @test isempty(ChaFit.resolve_cha_pins(:HK1, :H4, :mode1))   # data-driven, no pins
end

@testset "run plumbing: HK1 runs H1 (3 modes) + H4 (mode1); H3 removed" begin
    cells = FitRateEquation._cells(:HK1)
    @test length(cells) == 4                          # H1×{mode1,2,3} + H4×{mode1}
    @test Set(c[1] for c in cells) == Set([:H1, :H4])
    @test Set(c[3] for c in cells) == Set([:mode1, :mode2, :mode3])
    @test count(c -> c[1] === :H4, cells) == 1        # H4 is mode1-only
    @test length(FitRateEquation._cells(:G6PD)) == 2   # G6PD unchanged (1 variant × 2 modes)
end

else
    @testset "HK1 fit (skipped: HK1_AVAILABLE=false)" begin
        @test_skip "HK1 wiring unavailable pending HK1 port"
    end
end

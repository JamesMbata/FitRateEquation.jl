using FitRateEquation
using EnzymeRates
using Statistics
using Test


@testset "_fit_and_cv Cha macro-coord fit (both modes)" begin
    d = load_dataset(g6pd_config())
    variant = FitRateEquation.deploy_variant(:G6PD)        # :SS_NADPH_release_rate_eq
    mech = FitRateEquation._deploy_mech(:G6PD)
    coords_syms = FitRateEquation.ChaFit.cha_coords(:G6PD)
    # Mode 2: the Ki_ATP COORD is literature-pinned; the fit holds it at its lit value.
    r2 = FitRateEquation._fit_and_cv(variant, mech, d;
                                    mode=:mode2, n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    @test r2.mode === :mode2
    @test isfinite(r2.fit.loss) && isfinite(r2.cv.mean_cv)
    # The fit result carries coords (Dict{Symbol,Float64}) + loss + restarts.
    @test r2.fit.coords isa AbstractDict
    @test Set(keys(r2.fit.coords)) == Set(coords_syms)
    @test !isempty(r2.pins)
    @test haskey(r2.pins, :Ki_ATP)               # hard Cha coord pin, log10 value
    @test haskey(r2.pins, :Km_NADPH_rev)         # always-anchored reverse channel (all modes)
    for (k, v) in r2.pins
        @test isapprox(log10(r2.fit.coords[k]), v; atol=1e-6)
    end
    # G6PD carries no Km_PGA anchor in any mode.
    @test FitRateEquation.cha_anchors(:G6PD, :mode2) === nothing
    # Mode 1: only the always-anchored reverse channel (Km_NADPH_rev) is pinned; Ki_ATP free.
    r1 = FitRateEquation._fit_and_cv(variant, mech, d;
                                    mode=:mode1, n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    @test r1.mode === :mode1
    @test haskey(r1.pins, :Km_NADPH_rev) && !haskey(r1.pins, :Ki_ATP)
    @test r1.fit.coords isa AbstractDict
    @test isfinite(r1.fit.loss) && isfinite(r1.cv.mean_cv)
end

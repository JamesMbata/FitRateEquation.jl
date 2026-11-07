using FitRateEquation
using FitRateEquation: g6pd_config, pgd_config, load_dataset, nrows, Dataset
using FitRateEquation.ChaFit
using FitRateEquation.ChaFit: cha_coords, cha_coord_bounds, resolve_cha_pins, _cha_objective
using CMAEvolutionStrategy: minimize, xbest, fbest
using LinearAlgebra: norm
using Random
using Test

include("frozen_cha_loss.jl")
using .FrozenChaLoss

# ##########################################################################################
#                          What the fit objective must keep and must gain                    #
# ##########################################################################################
# The rewrite is only allowed to make the loss faster. So every value is compared with `==`
# against the frozen pre-rewrite code, in this process, and a seeded fit must follow the same
# CMA-ES path. `mech` is unused by the loss arithmetic; the frozen code takes it only for
# signature parity.

const _NOPINS = Dict{Symbol,Float64}()

_abs_data() = load_dataset(g6pd_config(data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv")))

# One figure wholly wrong-signed (every row penalized, so its group is skipped by the variance
# fold) and one row of another figure wrong-signed (penalty mixed with finite rows).
function _with_wrong_signs(d)
    rate = copy(d.rate)
    groups = unique(d.group)
    rate[d.group .== groups[1]] .*= -1
    rate[findfirst(==(groups[2]), d.group)] *= -1
    Dataset(d.concs, rate, d.group, d.keq, d.Et)
end

_points(lo, hi, n; seed=1) = (rng = MersenneTwister(seed);
                              [lo .+ (hi .- lo) .* rand(rng, length(lo)) for _ in 1:n])

const _RELATIVE_CASES = [
    (:G6PD, :SS_NADPH_release_rate_eq), (:G6PD, :RE_rate_eq), (:G6PD, :no_atp),
    (:G6PD, :no_g6p_nadph_deadend), (:G6PD, :no_g6p_atp_deadend), (:G6PD, :no_g6p_both_deadends),
    (:PGD, :cha_base), (:PGD, :no_atp), (:PGD, :full_re)]

_corpus(enzyme) = load_dataset(enzyme === :G6PD ? g6pd_config() : pgd_config())

@testset "objective equals the frozen loss bit for bit" begin
    @testset "$enzyme $variant" for (enzyme, variant) in _RELATIVE_CASES
        d = _corpus(enzyme)
        syms = cha_coords(enzyme, variant)
        lo, hi = cha_coord_bounds(enzyme, variant)
        o = _cha_objective(enzyme, d, syms, _NOPINS; variant=variant)
        for u in _points(lo, hi, 20)
            @test o(u) == FrozenChaLoss._cha_loss_with_pins(enzyme, nothing, d, u, syms, _NOPINS;
                                                            variant=variant)
        end
    end
    @testset "pinned coordinates are ignored at their pinned values" begin
        for enzyme in (:G6PD, :PGD)
            variant = enzyme === :G6PD ? :SS_NADPH_release_rate_eq : :cha_base
            d = _corpus(enzyme)
            syms = cha_coords(enzyme, variant)
            lo, hi = cha_coord_bounds(enzyme, variant)
            pins = resolve_cha_pins(enzyme, variant, :mode2)
            @test !isempty(pins)
            o = _cha_objective(enzyme, d, syms, pins; variant=variant)
            for u in _points(lo, hi, 10)
                v = copy(u)
                @test o(u) == FrozenChaLoss._cha_loss_with_pins(enzyme, nothing, d, u, syms, pins;
                                                                variant=variant)
                @test u == v     # the search vector CMA-ES owns is not written to
            end
        end
    end
    @testset "scalar keq override" begin
        d = _corpus(:G6PD); syms = cha_coords(:G6PD)
        lo, hi = cha_coord_bounds(:G6PD)
        o = _cha_objective(:G6PD, d, syms, _NOPINS; keq=13.7)
        for u in _points(lo, hi, 5)
            @test o(u) == FrozenChaLoss._cha_loss_with_pins(:G6PD, nothing, d, u, syms, _NOPINS; keq=13.7)
        end
    end
    @testset "wrong-signed rows take the same penalty path" begin
        for (enzyme, variant) in ((:G6PD, :SS_NADPH_release_rate_eq), (:PGD, :cha_base))
            d = _with_wrong_signs(_corpus(enzyme))
            syms = cha_coords(enzyme, variant)
            lo, hi = cha_coord_bounds(enzyme, variant)
            o = _cha_objective(enzyme, d, syms, _NOPINS; variant=variant)
            for u in _points(lo, hi, 5)
                l = o(u)
                @test l == FrozenChaLoss._cha_loss_with_pins(enzyme, nothing, d, u, syms, _NOPINS;
                                                             variant=variant)
                @test l >= ChaFit._SIGN_PENALTY / nrows(d)   # the penalty path really ran
            end
        end
    end
    @testset "absolute scale" begin
        d = _abs_data()
        syms = cha_coords(:G6PD, :_deploy; scale=:absolute)
        lo, hi = cha_coord_bounds(:G6PD, :_deploy; scale=:absolute)
        for pins in (_NOPINS, Dict(:Kd_6PGLn => log10(2.1e-4)))
            o = _cha_objective(:G6PD, d, syms, pins; scale=:absolute)
            for u in _points(lo, hi, 20)
                @test o(u) == FrozenChaLoss._cha_loss_with_pins(:G6PD, nothing, d, u, syms, pins;
                                                                scale=:absolute)
            end
        end
    end
    @testset "explicit release rate, as the koffQ profile uses it" begin
        d = _corpus(:G6PD); syms = cha_coords(:G6PD)
        lo, hi = cha_coord_bounds(:G6PD)
        for r in (0.1, 1.0, 1e6)
            o = _cha_objective(:G6PD, d, syms, _NOPINS; release_rate=r)
            for u in _points(lo, hi, 5)
                @test o(u) == FrozenChaLoss._cha_loss_with_pins(:G6PD, nothing, d, u, syms, _NOPINS;
                                                                release_rate=r)
            end
        end
    end
    @testset "an earlier evaluation leaves nothing behind" begin
        d = _with_wrong_signs(_corpus(:G6PD)); syms = cha_coords(:G6PD)
        lo, hi = cha_coord_bounds(:G6PD)
        o = _cha_objective(:G6PD, d, syms, _NOPINS)
        a, b = _points(lo, hi, 2)
        la = o(a); lb = o(b)
        @test o(a) == la
        @test o(b) == lb
    end
end

@testset "public loss functions equal the frozen loss bit for bit" begin
    function dict_points(enzyme, variant; scale=:relative, n=8)
        syms = cha_coords(enzyme, variant; scale=scale)
        lo, hi = cha_coord_bounds(enzyme, variant; scale=scale)
        [Dict(syms .=> 10 .^ u) for u in _points(lo, hi, n)]
    end
    @testset "centered: $enzyme $variant" for (enzyme, variant) in _RELATIVE_CASES
        d = _corpus(enzyme)
        for c in dict_points(enzyme, variant)
            @test ChaFit.cha_centered_logratio_loss(enzyme, nothing, d, c; variant=variant) ==
                  FrozenChaLoss.cha_centered_logratio_loss(enzyme, nothing, d, c; variant=variant)
        end
    end
    @testset "centered with explicit kf, Et, release_rate, release_eq, kr" begin
        d = _corpus(:G6PD)
        for c in dict_points(:G6PD, :_deploy; n=3), kw in (
                (; kf=2.0), (; Et=3.0), (; release_rate=5.0), (; release_eq=1e-5),
                (; kr=0.25), (; keq=13.7), (; kf=2.0, Et=3.0, release_rate=5.0, release_eq=1e-5, kr=0.25))
            @test ChaFit.cha_centered_logratio_loss(:G6PD, nothing, d, c; kw...) ==
                  FrozenChaLoss.cha_centered_logratio_loss(:G6PD, nothing, d, c; kw...)
        end
    end
    @testset "absolute" begin
        d = _abs_data()
        for c in dict_points(:G6PD, :_deploy)
            @test ChaFit.cha_absolute_logratio_loss(:G6PD, nothing, d, c; kf=178.0) ==
                  FrozenChaLoss.cha_absolute_logratio_loss(:G6PD, nothing, d, c; kf=178.0)
        end
    end
    @testset "row log-ratios and groups" begin
        d = _with_wrong_signs(_corpus(:PGD))
        for c in dict_points(:PGD, :cha_base; n=3)
            a = fill(NaN, nrows(d)); b = fill(NaN, nrows(d))
            pa, ga = ChaFit._cha_row_logratios!(a, :PGD, nothing, d, c)
            pb, gb = FrozenChaLoss._cha_row_logratios!(b, :PGD, nothing, d, c)
            @test pa == pb
            @test ga == gb
            @test isequal(a, b)
        end
    end
    @testset "macro tuple: names and values" begin
        for (enzyme, variant) in _RELATIVE_CASES, c in dict_points(enzyme, variant; n=2)
            @test ChaFit.cha_macro_tuple(enzyme, c; keq=1.7, variant=variant) ===
                  FrozenChaLoss.cha_macro_tuple(enzyme, c; keq=1.7, variant=variant)
        end
    end
    @testset "non-uniform keq within a figure still throws" begin
        d = _corpus(:G6PD)
        k = copy(d.keq); k[1] *= 2
        bad = Dataset(d.concs, d.rate, d.group, k, d.Et)
        c = first(dict_points(:G6PD, :_deploy; n=1))
        @test_throws ErrorException ChaFit.cha_centered_logratio_loss(:G6PD, nothing, bad, c)
        @test_throws ErrorException _cha_objective(:G6PD, bad, cha_coords(:G6PD), _NOPINS)
    end
end

@testset "seeded fit follows the frozen CMA-ES path" begin
    # A small iteration budget and a huge time budget: maxtime never binds, so the stopping
    # point depends only on the values the objective returns.
    for (enzyme, variant, pins) in ((:G6PD, :SS_NADPH_release_rate_eq, _NOPINS),
                                    (:PGD, :cha_base, resolve_cha_pins(:PGD, :cha_base, :mode2)))
        d = _corpus(enzyme)
        syms = cha_coords(enzyme, variant)
        lo, hi = cha_coord_bounds(enzyme, variant)
        frozen = u -> FrozenChaLoss._cha_loss_with_pins(enzyme, nothing, d, u, syms, pins; variant=variant)
        fit = ChaFit.cha_fit_candidate(enzyme, nothing, d; n_restarts=2, maxiter=40, maxtime=1e4,
                                       seed=7, pins=pins, variant=variant)
        losses = Float64[]; endpoints = Vector{Float64}[]
        for r in 1:2
            s = hash((7, r)) % UInt64
            x0 = lo .+ (hi .- lo) .* rand(MersenneTwister(s), length(syms))
            o = minimize(frozen, x0, 0.1; lower=lo, upper=hi, maxiter=40, maxtime=1e4, seed=s,
                         verbosity=0)
            push!(losses, fbest(o)); push!(endpoints, collect(xbest(o)))
        end
        best = argmin(losses)
        @test fit.loss == losses[best]
        @test [x.loss for x in fit.restarts] == losses
        @test [x.dist for x in fit.restarts] == [norm(e .- endpoints[best]) for e in endpoints]
        final = [haskey(pins, s) ? pins[s] : endpoints[best][i] for (i, s) in enumerate(syms)]
        @test fit.coords == Dict(syms .=> 10 .^ final)
    end
end

# ##########################################################################################
#                                 Gate: it must stay fast                                    #
# ##########################################################################################
# Allocation and inference catch most ways of losing the speed. Run-time dispatch that does not
# allocate gets past both, so the last check times the objective against the frozen code in this
# process, which cancels out machine speed. Measured ratio at the rewrite: about 11x.

_allocs(o, u) = (o(u); @allocated o(u))

@testset "objective: no allocation, concrete result" begin
    @testset "$enzyme $variant" for (enzyme, variant) in _RELATIVE_CASES
        d = _corpus(enzyme)
        syms = cha_coords(enzyme, variant)
        lo, hi = cha_coord_bounds(enzyme, variant)
        u = first(_points(lo, hi, 1))
        o = _cha_objective(enzyme, d, syms, _NOPINS; variant=variant)
        @test _allocs(o, u) == 0
        @test (@inferred o(u)) isa Float64
    end
    @testset "with pins" begin
        d = _corpus(:G6PD); syms = cha_coords(:G6PD)
        lo, hi = cha_coord_bounds(:G6PD)
        o = _cha_objective(:G6PD, d, syms, resolve_cha_pins(:G6PD, :SS_NADPH_release_rate_eq, :mode2))
        @test _allocs(o, first(_points(lo, hi, 1))) == 0
    end
    @testset "absolute scale" begin
        d = _abs_data()
        syms = cha_coords(:G6PD, :_deploy; scale=:absolute)
        lo, hi = cha_coord_bounds(:G6PD, :_deploy; scale=:absolute)
        o = _cha_objective(:G6PD, d, syms, _NOPINS; scale=:absolute)
        u = first(_points(lo, hi, 1))
        @test _allocs(o, u) == 0
        @test (@inferred o(u)) isa Float64
    end
end

@testset "objective: at least 4x faster than the frozen loss" begin
    d = _corpus(:G6PD); syms = cha_coords(:G6PD)
    lo, hi = cha_coord_bounds(:G6PD)
    u = first(_points(lo, hi, 1))
    o = _cha_objective(:G6PD, d, syms, _NOPINS)
    frozen = x -> FrozenChaLoss._cha_loss_with_pins(:G6PD, nothing, d, x, syms, _NOPINS)
    best(f) = (f(u); minimum(@elapsed(for _ in 1:200; f(u); end) for _ in 1:5))
    # Coverage counters cost time on every executed line, which swamps a loss this cheap: the
    # ratio measured about 1.9x under Pkg.test(coverage=true), as CI runs it. The time says nothing
    # there, so only a run without coverage checks it; the allocation gate above runs in both.
    if Base.JLOptions().code_coverage == 0
        @test best(frozen) >= 4 * best(o)
    else
        @test_skip best(frozen) >= 4 * best(o)
    end
end

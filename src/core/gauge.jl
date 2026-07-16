################################################################################
#                  MANDATORY SCALE GAUGE + ANALYTIC kcat                       #
################################################################################
#
# The first steady-state forward rate constant (k<i>f) of each mechanism is
# pinned to 1.0 (the "k1f=1" scale gauge). This removes the single overall
# rate-scale degeneracy that every mechanism carries, so the free parameter
# vector has one fewer dimension than `fitted_params`. The gauge is ALWAYS on
# (both Mode 1 and Mode 2) and is implemented here, IN-PIPELINE, on EnzymeRates
# `main` — we do NOT use the gauged-fitting machinery (`loss_gauged!`,
# `params_with_gauge`, `*_kcat`, `can_gauge_kcat`, `VmaxRatioBox`).
#
# `build_params` is the SINGLE place that assembles the `rate_equation` params
# NamedTuple from a free log10 vector. It is AD-friendly: the returned values
# track `eltype(logθ)` (e.g. ForwardDiff.Dual) so sensitivities can be taken
# through it (Task 3.1).

"The SS catalysis step's forward rate constant — the kf=1 scale gauge. Upstream EnzymeRates
 names the chemistry forward rate `k_<from>_to_<to>` (SS releases are `koff_*`/`kon_*`,
 rapid-equilibrium constants `K_*`), so it is the unique fitted param of the `k_*_to_*` shape."
function gauge_param(mech)
    hits = [s for s in EnzymeRates.fitted_params(mech)
            if (t = String(s); startswith(t, "k_") && occursin("_to_", t))]
    isempty(hits) && error("no SS catalysis forward rate (k_*_to_*) to gauge in $(mech)")
    length(hits) > 1 && error("ambiguous gauge: multiple k_*_to_* forward rates $(hits) in $(mech)")
    first(hits)
end

"Fittable parameters minus the gauge parameter (the free vector's symbols)."
free_params(mech) = filter(!=(gauge_param(mech)), collect(EnzymeRates.fitted_params(mech)))

"Build the rate_equation params NamedTuple: gauge=1, free from 10^logθ, plus Keq/E_total."
function build_params(mech, logθ::AbstractVector; keq::Real, etotal::Real=1.0)
    free = free_params(mech)
    g = gauge_param(mech)
    # Value container element type tracks eltype(logθ) so AD (ForwardDiff.Dual)
    # propagates through build_params -> rate_equation (Task 3.1).
    T = promote_type(eltype(logθ), typeof(float(keq)), typeof(float(etotal)))
    vals = Dict{Symbol,T}(g => one(T))
    for (s, lv) in zip(free, logθ)
        vals[s] = T(10.0)^lv
    end
    fps = EnzymeRates.fitted_params(mech)
    names = (fps..., :Keq, :E_total)
    NamedTuple{names}((map(s -> vals[s], fps)..., T(keq), T(etotal)))
end

"Apparent kcat = rate per E_total at saturating substrate, zero product."
function analytic_kcat(mech, logθ; keq::Real, sat::Real=1.0)
    mets = EnzymeRates.metabolites(mech)
    # Saturating substrates = the mechanism's declared substrates (enzyme-agnostic).
    # `EnzymeRates.substrates` returns the substrate tuple (e.g. (:NADP,:G6P) for
    # G6PD, (:NADP,:PGA) for PGD); analytic_kcat saturates these and zeroes products.
    subs = EnzymeRates.substrates(mech)
    # Fail loud if the saturating substrates are absent from the mechanism: a silent
    # `s in subs` miss would zero every concentration and yield kcat ≈ 0 with no error.
    @assert all(in(mets), subs) "analytic_kcat substrates $(subs) not all present in mechanism metabolites $(mets)"
    concs = NamedTuple{Tuple(mets)}(map(s -> s in subs ? sat : 0.0, mets))
    p = build_params(mech, logθ; keq=keq)
    EnzymeRates.rate_equation(mech, concs, p) / p.E_total
end

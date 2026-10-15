################################################################################
#                 ENZYMERATES NAME LAYER (FRE vocabulary <-> upstream)           #
################################################################################
#
# FRE's vocabulary is EnzymeRates 0.5.0's composition-semantic parameter naming:
#   K_<ligand>_<form>          RE binding (dissociation) constant, e.g. K_G6P_ENADP
#   koff_<ligand>_<form>       SS release (off) rate,                e.g. koff_NADPH_E
#   kon_<ligand>_<form>        SS binding (on) rate,                 e.g. kon_NADPH_E
#   k_<from>_to_<to>           SS chemistry rate,                    e.g. k_EG6PNADP_to_ENADPHPGLn
# The Cha derivations, the deploy maps, micro_parameters.jl keys and every downstream
# consumer (PentosePhosphatePathway.jl param names) are written in it. EnzymeRates >= 0.6
# names every constant structurally after its reaction (`K_<complex>_to_<form>_<ligand>`,
# `k_<a>_to_<b>`), so this layer translates at the boundary: every FRE call that hands
# EnzymeRates a parameter NamedTuple, or reads parameter symbols back from it, goes through
# the `er_*` functions below. Parameter order is ASCII-sorted FRE names (0.5.0's order), so
# free vectors and fixtures keep their positions. test/test_enzymerates_names.jl replays a
# 0.5.0-generated oracle (names, order, step view, rate values) against this layer.

"Translate one upstream (structural) parameter name to FRE's (0.5.0) name."
function legacy_param_name(s::Symbol)
    s in (:Keq, :E_total) && return s
    mt = match(r"^(K|k)_(.+)_to_(.+)$", String(s))
    mt === nothing && error("legacy_param_name: unrecognized parameter name $s")
    kind, a, b = mt.captures
    as, bs = split(a, '_'), split(b, '_')
    if kind == "K"
        # RE constants are named in the release direction: complex -> form + ligand.
        (length(as) == 1 && length(bs) == 2) ||
            error("legacy_param_name: $s is not a binding constant (no 0.5.0 form)")
        return Symbol("K_", bs[2], "_", bs[1])
    end
    length(as) == 1 && length(bs) == 1 && return s                   # chemistry
    length(as) == 1 && length(bs) == 2 && return Symbol("koff_", bs[2], "_", bs[1])
    length(as) == 2 && length(bs) == 1 && return Symbol("kon_", as[2], "_", as[1])
    error("legacy_param_name: unrecognized rate-constant shape $s")
end

# Per-mechanism-type translation, cached: upstream fitted names, their FRE names, and the
# FRE fitted tuple in FRE (sorted) order.
const _ER_NAME_CACHE = Dict{Any,Any}()
const _ER_NAME_LOCK = ReentrantLock()

function _er_names(m)
    lock(_ER_NAME_LOCK) do
        get!(_ER_NAME_CACHE, typeof(m)) do
            up = Tuple(EnzymeRates.fitted_params(m))
            fre = map(legacy_param_name, up)
            allunique(fre) || error("legacy name collision in $(typeof(m)): $(fre)")
            dep, _ = EnzymeRates._dependent_param_exprs(typeof(m))
            rename = Dict{Symbol,Symbol}(s => legacy_param_name(s)
                                         for s in (up..., keys(dep)...))
            (up = up, fre = fre, sorted = Tuple(sort(collect(fre); by = String)),
             rename = rename)
        end
    end
end

"Fitted parameter names in FRE vocabulary, ASCII-sorted (EnzymeRates 0.5.0's order)."
er_fitted_params(m) = _er_names(m).sorted

"All rate_equation parameter names in FRE vocabulary: fitted, then Keq and E_total."
function er_parameters(m)
    up = Tuple(EnzymeRates.parameters(m))
    up[end-1:end] == (:Keq, :E_total) ||
        error("er_parameters: expected (..., :Keq, :E_total), got $(up)")
    (er_fitted_params(m)..., :Keq, :E_total)
end

"Re-key an FRE-vocabulary params NamedTuple into upstream names and order (AD-safe)."
function _to_upstream(m, p::NamedTuple)
    n = _er_names(m)
    NamedTuple{(n.up..., :Keq, :E_total)}((map(s -> getfield(p, s), n.fre)..., p.Keq, p.E_total))
end

"EnzymeRates.rate_equation with an FRE-vocabulary params NamedTuple."
er_rate_equation(m, concs, p::NamedTuple) = EnzymeRates.rate_equation(m, concs, _to_upstream(m, p))

_rename_syms(e::Symbol, d) = get(d, e, e)
_rename_syms(e::Expr, d) = Expr(e.head, map(a -> _rename_syms(a, d), e.args)...)
_rename_syms(e, d) = e

"Haldane/Wegscheider closures `(dependent => expr, independent)` in FRE vocabulary."
function er_dependent_param_exprs(m)
    d = _er_names(m).rename
    dep, _ = EnzymeRates._dependent_param_exprs(typeof(m))
    (Dict{Symbol,Any}(d[s] => _rename_syms(ex, d) for (s, ex) in dep), er_fitted_params(m))
end

"EnzymeRates' raw symbolic (num, den, ...) POLYs with rate-constant symbols in FRE vocabulary."
function er_raw_rate_polys(m)
    d = _er_names(m).rename
    map(EnzymeRates._raw_symbolic_rate_polys(typeof(m))) do poly
        Dict([get(d, s, s) => e for (s, e) in mono] => c for (mono, c) in poly)
    end
end

"Elementary steps as 0.5.0's `reactions`: ((form, consumed...), (form, released...), is_eq, group)."
function er_reactions(m)
    mm = EnzymeRates.Mechanism(m)
    out = Any[]
    for (g, group) in enumerate(EnzymeRates.steps(mm)), s in group
        push!(out, ((EnzymeRates.name(EnzymeRates.from_species(s)),
                     map(EnzymeRates.name, EnzymeRates.consumed(s))...),
                    (EnzymeRates.name(EnzymeRates.to_species(s)),
                     map(EnzymeRates.name, EnzymeRates.released(s))...),
                    EnzymeRates.is_equilibrium(s), g))
    end
    Tuple(out)
end

"Substrate names, as 0.5.0's `substrates(mechanism)` (upstream now takes the reaction)."
er_substrates(m) = Tuple(map(EnzymeRates.name,
    EnzymeRates.substrates(EnzymeRates.reaction(EnzymeRates.Mechanism(m)))))

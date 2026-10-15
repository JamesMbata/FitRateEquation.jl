# Macroscopic kinetic constants read from the Haldane-reduced symbolic rate law.
#
# EnzymeRates._raw_symbolic_rate_polys returns the RAW (num, den) POLYs whose
# rate-constant symbols include DEPENDENT ones — Wegscheider-tied binding K's
# (e.g. K3) and the Haldane-derived catalysis reverse k (k5r) — which are NOT in
# build_params/fitted_params. We rebuild the full symbol->value map with the SAME
# closures the generated rate_equation uses (EnzymeRates._dependent_param_exprs),
# so num/den evaluate to exactly rate_equation. Gauged denominator coefficients
# (each denominator term divided by the free-enzyme constant term) are the apparent
# association constants; macro constants are read from them via the alias map (5.2).

# Recursive numeric evaluation of a dependent-param / arithmetic Expr over a
# Symbol=>value map (handles Number, Symbol, and :call Exprs with + - * / ^).
function _eval_expr(e, vals)
    e isa Number && return float(e)
    e isa Symbol && return float(vals[e])
    if e isa Expr && e.head === :call
        op = e.args[1]
        a = [_eval_expr(x, vals) for x in e.args[2:end]]
        op === :+ && return sum(a)
        op === :* && return prod(a)
        op === :- && return length(a) == 1 ? -a[1] : a[1] - sum(a[2:end])
        op === :/ && return a[1] / a[2]
        op === :^ && return a[1] ^ a[2]
        error("unhandled operator $op in dependent-param expr")
    end
    error("unhandled expr node $e ($(typeof(e)))")
end

# Complete rate-constant value map at logθ: independent params (+gauge, Keq,
# E_total) from build_params, then the dependent params from the Haldane/
# Wegscheider closures. Element type Float64 (closed-form, not AD).
function _micro_values(mech, logθ; keq::Real)
    p = build_params(mech, logθ; keq=keq)
    vals = Dict{Symbol,Float64}()
    for s in keys(p); vals[s] = float(getfield(p, s)); end
    dep_exprs, _ = EnzymeRates._dependent_param_exprs(typeof(mech))
    for (sym, expr) in dep_exprs
        vals[sym] = _eval_expr(expr, vals)
    end
    vals
end

# Evaluate a raw POLY (Dict{monomial=>rational}) at concentrations `concs`
# (NamedTuple) and a COMPLETE rate-constant value map `vals` (from _micro_values).
function _eval_poly(poly, concs::NamedTuple, vals::AbstractDict)
    metset = keys(concs)
    total = 0.0
    for (mono, coeff) in poly
        term = float(coeff)
        for (sym, e) in mono
            base = sym in metset ? getfield(concs, sym) :
                   haskey(vals, sym) ? vals[sym] :
                   error("unknown POLY symbol $sym")
            term *= float(base)^e
        end
        total += term
    end
    total
end

# Split a monomial into (concentration-part, rate-constant-part), each sorted.
function _split_mono(mono, metset)
    conc = Pair{Symbol,Int}[]; kpart = Pair{Symbol,Int}[]
    for (s, e) in mono
        push!(s in metset ? conc : kpart, s => e)
    end
    (sort!(conc; by=first), sort!(kpart; by=first))
end

"Gauged denominator coefficients keyed by concentration-monomial (Vector{Pair}).
 Each value is the numeric coefficient at θ divided by the constant (gauge) term."
function gauged_denominator_coeffs(mech, logθ; keq::Real)
    _, den = EnzymeRates._raw_symbolic_rate_polys(typeof(mech))
    metset = Set(EnzymeRates.metabolites(mech))
    vals = _micro_values(mech, logθ; keq=keq)
    bins = Dict{Vector{Pair{Symbol,Int}}, Float64}()
    for (mono, coeff) in den
        conc, kpart = _split_mono(mono, metset)
        val = float(coeff) * prod((vals[s]^e for (s, e) in kpart); init=1.0)
        bins[conc] = get(bins, conc, 0.0) + val
    end
    g0 = get(bins, Pair{Symbol,Int}[], 0.0)
    @assert g0 != 0 "denominator has no constant (free-enzyme) term to gauge"
    Dict(k => v / g0 for (k, v) in bins)
end

_monolabel(conc) = join(["$(s)^$(e)" for (s,e) in conc], "_")

# Resolve the free micro-K of a rapid-equilibrium DEAD-END binding step where `lig`
# binds the single free `form`. The dead-end product form (e.g. E_dN) is what
# distinguishes this from the productive release step — under RE canonicalization the
# released NADPH also lands on lhs_mets, so a bare "lig on lhs" test is NOT enough.
# Returns `nothing` if the mechanism has no such dead-end step (e.g. V1/V2 of a PGD
# whose enzyme-wide `ki_micro_direct` only applies to the V3 dead-end variant).
function _micro_direct_sym(mech, lig::Symbol, form::Symbol)
    free = Set(free_params(mech))
    de = _deadend_forms(mech)
    for st in _mechanism_steps(mech)
        st.is_eq || continue
        (lig in st.lhs_mets) || continue            # ligand participates
        (form in st.lhs_forms) || continue          # binds this free form
        any(f -> f in de, st.rhs_forms) || continue # product is a DEAD-END form (not release)
        s = Symbol("K$(st.param_index)")
        s in free && return s
    end
    nothing
end

"Macroscopic constants at θ for `enzyme` (default: inferred from the mechanism).
 Named constants Ki_X / Kd_X (1/coeff of [X]^1) and Km_X (coeff[cosub]/coeff[pair]);
 everything else an anonymous denominator-monomial lump."
function macro_constants(mech, logθ; keq::Real, enzyme::Symbol=_enzyme_of(mech))
    w = _wiring(enzyme)
    g = gauged_denominator_coeffs(mech, logθ; keq=keq)
    pair = w.substrate_pair
    out = NamedTuple[]
    for (conc, coeff) in g
        isempty(conc) && continue                  # the gauge term itself
        key = conc
        if haskey(w.ki_map, conc)                  # dead-end inhibitor: Ki = 1/coeff
            push!(out, (name=w.ki_map[conc], role=:named, value=1.0/coeff,
                        recipe = gg -> 1.0/gg[key], micro=nothing))
        elseif haskey(w.kd_map, conc)              # substrate: emit BOTH binary Kd and apparent Km
            push!(out, (name=w.kd_map[conc], role=:named, value=1.0/coeff,
                        recipe = gg -> 1.0/gg[key], micro=nothing))
            if haskey(w.km_ratio, conc) && haskey(g, pair)
                kmname, cosub = w.km_ratio[conc]
                push!(out, (name=kmname, role=:named, value=g[cosub]/g[pair],
                            recipe = gg -> gg[cosub]/gg[pair], micro=nothing))
            end
        else                                       # anonymous denominator lump
            push!(out, (name=Symbol("lump_", _monolabel(conc)), role=:anonymous,
                        value=coeff, recipe = gg -> gg[key], micro=nothing))
        end
    end
    # Ki-ratio named constants: noncompetitive dead-end Ki read as the cross-term ratio
    # g[carrier]/g[carrier·ligand] (e.g. G6PD Ki_NADPH on E·G6P), NOT the bare [ligand]
    # term (which is the reverse release Km).
    for (nm, (carrier, cross)) in w.ki_ratio
        (haskey(g, carrier) && haskey(g, cross)) || continue
        push!(out, (name=nm, role=:named, value=g[carrier]/g[cross],
                    recipe = gg -> gg[carrier]/gg[cross], micro=nothing))
    end
    # Micro-direct named constants (e.g. V3's decoupled forward Ki_NADPH_fwd): read
    # 1/K of the dead-end binding step directly, bypassing the (lumped) gauged coeff.
    if !isempty(w.ki_micro_direct)
        vals = _micro_values(mech, logθ; keq=keq)
        for (nm, (lig, form)) in w.ki_micro_direct
            sym = _micro_direct_sym(mech, lig, form)   # nothing if no dead-end (V1/V2)
            (sym === nothing || sym ∉ keys(vals)) && continue
            push!(out, (name=nm, role=:named, value=1.0/vals[sym],
                        recipe = (gg -> NaN), micro=sym))
        end
    end
    sort!(out; by=x->String(x.name))
    out
end

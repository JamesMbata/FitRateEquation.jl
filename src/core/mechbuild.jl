# Mechanism builder + step introspection, vendored from fitting/mechanism_id/
# {panel/engine.jl, anchor.jl} during the EnzymeRates upstream migration. Only the
# builder helpers consensus_macro uses are kept (the ranking-panel machinery —
# _subsets/_patterns/PanelSpec/_deadend_combos/build_panel — is dropped).
#
# The upstream @enzyme_mechanism DSL rejects opaque form names (E_NB) and requires
# decomposed call notation E(NADP,G6P). `_mech` derives each form's bound-metabolite
# composition from the step graph (`_form_compositions`) and emits call notation; the
# terse opaque form names in enzymes/<enzyme>.jl are unchanged.

# Per-row penalty for a predicted/observed rate-sign mismatch in the centered-log-ratio
# loss. Vendored from the consensus_macro core/loss.jl const (the coeff-space
# centered_logratio_loss it accompanied is superseded by ChaFit.cha_centered_logratio_loss
# and was not carried, but this const is still consumed by the Cha loss and test_cha_fit.jl).
const _SIGN_PENALTY = 10.0

"""
    _form_compositions(metset, steps) -> Dict{Symbol,Vector{Symbol}}

Bound-metabolite composition of every enzyme form, derived from the `(lhs, rhs, op)`
step list. Conservation across a BINDING/RELEASE step (a step with ≥1 free metabolite):
`bound(to) = (bound(from) ∪ free_lhs) \\ free_rhs`. CHEMISTRY/iso steps (no free
metabolite on either side) CONVERT bound substrates→products and are EXCLUDED from
propagation — their endpoints get composition from the binding (forward) and release
(backward) steps instead, so an SS super-node resolves without a substrate→product
pairing. Iterates to a fixpoint from `bound(:E) = []`.
"""
function _form_compositions(metset::Set{Symbol}, steps)
    isform(x) = !(x in metset)
    comp = Dict{Symbol,Vector{Symbol}}(:E => Symbol[])
    changed = true
    while changed
        changed = false
        for (lhs, rhs, _op) in steps
            lf = filter(isform, lhs); lm = filter(in(metset), lhs)
            rf = filter(isform, rhs); rm = filter(in(metset), rhs)
            (length(lf) == 1 && length(rf) == 1) || continue   # single enzyme form each side
            (isempty(lm) && isempty(rm)) && continue           # chemistry/iso step: skip
            L, R = lf[1], rf[1]
            if haskey(comp, L) && !haskey(comp, R)
                comp[R] = sort(unique(setdiff(union(comp[L], lm), rm))); changed = true
            elseif haskey(comp, R) && !haskey(comp, L)
                comp[L] = sort(unique(setdiff(union(comp[R], rm), lm))); changed = true
            end
        end
    end
    comp
end

# Build one `lhs op rhs` step Expr; `termexpr` maps each term Symbol to its emitted
# form (default identity = opaque names; `_mech` passes a call-notation mapper).
function _step(lhs::Vector{Symbol}, rhs::Vector{Symbol}, op::Symbol, termexpr=identity)
    mkside(syms) = length(syms) == 1 ? termexpr(syms[1]) : Expr(:call, :+, map(termexpr, syms)...)
    Expr(:call, op, mkside(lhs), mkside(rhs))
end

"""
    _mech(subs, prods, steps; regs=Symbol[]) -> EnzymeMechanism

Assemble the `@enzyme_mechanism` body AST from substrate/product symbol lists
and a vector of `(lhs::Vector{Symbol}, rhs::Vector{Symbol}, op::Symbol)` steps,
then evaluate it through the surface macro.

`regs` declares extra metabolites (e.g. the dead-end ligand `ATP`) so they are
recognized as species — not enzyme forms — in any dead-end binding step. NADPH
is already a declared product, so it needs no regulator entry to dead-end-bind.
A `regulators:` line is emitted only when `regs` is non-empty.
"""
function _mech(subs::Vector{Symbol}, prods::Vector{Symbol}, steps::Vector;
               regs::Vector{Symbol}=Symbol[])
    metset = Set{Symbol}(vcat(subs, prods, regs))
    comp = _form_compositions(metset, steps)
    # form symbol -> call-notation Expr (E / E(m1,m2,...)); metabolite -> bare symbol.
    termexpr(x) = x in metset ? x :
        (haskey(comp, x) ?
            (isempty(comp[x]) ? :E : Expr(:call, :E, comp[x]...)) :
            error("_mech: enzyme form $x has no derivable composition " *
                  "(not reachable from :E via any binding/release step)"))
    stepblock = Expr(:block)
    for (lhs, rhs, op) in steps
        push!(stepblock.args, _step(lhs, rhs, op, termexpr))
    end
    sub_line = Expr(:tuple, Expr(:call, :(:), :substrates, subs[1]), subs[2:end]...)
    prod_line = Expr(:tuple, Expr(:call, :(:), :products, prods[1]), prods[2:end]...)
    steps_line = Expr(:call, :(:), :steps, stepblock)
    body = Expr(:block, sub_line, prod_line)
    if !isempty(regs)
        push!(body.args, Expr(:tuple, Expr(:call, :(:), :regulators, regs[1]), regs[2:end]...))
    end
    push!(body.args, steps_line)
    eval(:(@enzyme_mechanism $body))
end

"""
    _deadend_step(form::Symbol, ligand::Symbol) -> (lhs, rhs, op)

A rapid-equilibrium dead-end binding step `form + ligand ⇌ form_<ligand>`, with
a fresh uniquely-tagged product form (`E_ATP`, `E_G_ATP`, `E_dN`, `E_G_dN`, …)
so it cannot alias an existing enzyme form. NADPH uses the short `dN`
("dead-end NADPH") tag.
"""
function _deadend_step(form::Symbol, ligand::Symbol)
    tag = ligand === :NADPH ? "dN" : String(ligand)
    newform = Symbol(String(form), "_", tag)
    ([form, ligand], [newform], :(⇌))
end

# Structural view of a built mechanism's elementary steps, in source order, with
# metabolites/forms split out per side and each step's kinetic-group param index.
function _mechanism_steps(mech)
    rxns = EnzymeRates.reactions(mech)
    met_set = Set(EnzymeRates.metabolites(mech))
    # Representative step index for each kinetic group (the index its rate
    # constants are named after; see _group_param_symbols in EnzymeRates).
    group_rep = Dict{Int,Int}()
    for (i, step) in enumerate(rxns)
        g = step[4]
        haskey(group_rep, g) || (group_rep[g] = i)
    end
    steps = NamedTuple[]
    for (i, (lhs, rhs, is_eq, group)) in enumerate(rxns)
        lhs_v = collect(Symbol, lhs)
        rhs_v = collect(Symbol, rhs)
        push!(steps, (
            lhs       = lhs_v,
            rhs       = rhs_v,
            lhs_mets  = Symbol[s for s in lhs_v if s in met_set],
            rhs_mets  = Symbol[s for s in rhs_v if s in met_set],
            lhs_forms = Symbol[s for s in lhs_v if s ∉ met_set],
            rhs_forms = Symbol[s for s in rhs_v if s ∉ met_set],
            is_eq     = is_eq,
            group     = group,
            index     = i,
            param_index = group_rep[group],
        ))
    end
    steps
end

# (`step_for_param`, which mapped `k<i>f`/`K<i>` names to steps via a numeric index, was
# dropped in the upstream migration: upstream param names are composition-semantic, not
# step-indexed, so bounds classification now reads the name directly — see core/bounds.jl.)

"Enzyme forms that appear in exactly one elementary step (graph degree 1) — the
 dead-end forms (e.g. E_dN, E_ATP). All cyclic forms have degree >= 2."
function _deadend_forms(mech)
    deg = Dict{Symbol,Int}()
    for st in _mechanism_steps(mech)
        for f in (st.lhs_forms..., st.rhs_forms...)
            deg[f] = get(deg, f, 0) + 1
        end
    end
    Set(f for (f, d) in deg if d == 1)
end

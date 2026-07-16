# Generic consensus-mechanism builder helper shared across enzymes. The
# enzyme-specific topology (substrate/product lists, step wiring, dead-end forms)
# lives in enzymes/<enzyme>.jl and calls this helper + _mech directly.

# Rapid-equilibrium dead-end binding steps: for each (forms, ligand) spec, one
# `_deadend_step(form, ligand)` per open enzyme form the ligand abortively binds.
function _deadends(specs)
    extra = Tuple[]
    for (forms, lig) in specs, f in forms
        push!(extra, _deadend_step(f, lig))
    end
    extra
end

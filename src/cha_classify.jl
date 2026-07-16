# =========================================================================================
#         Identifiability + classification in MACRO-COORDINATE space (the Cha fit)
# =========================================================================================
#
# The Cha twin of coeff_identifiability.jl. The Cha fit varies the named forward shape
# constants DIRECTLY through cha_rate_* (cha_coords: substrate Kd's, alpha, dead-end/regulator
# Ki's, the bare-[NADPH] reverse-release Km). As in coeff space, the data_identified vs
# unconstrained decision is made HERE, over cha_coords (in log10), via the Hessian of
# cha_centered_logratio_loss. The cha_coords are FIXED named constants (no fiber / null-space:
# the gauge kf/Et, the swept release_rate fiber, and the Haldane kr are all held outside the
# coords), so every flat eigendirection is a genuinely unconstrained coordinate.
#
# This mirrors coeff_identifiability.jl::{coeff_identifiable_functions, classify_coords}
# EXACTLY -- same eigen/rank cut (tol=2e-2), same calibrated pseudo-inverse covariance and CI
# math, same class SYMBOLS (:data_identified / :literature_pinned / :unconstrained) -- only
# the parameter space (cha_coords) and the loss (cha_centered_logratio_loss) differ.
#
# The Hessian uses the SAME central-difference stencil as coeff_identifiability.jl (copied
# from identifiability.jl::_fd_hessian). Although cha_rate_* is pure arithmetic, the loss
# `cha_centered_logratio_loss` (cha_fit.jl) writes intermediate residuals into a hard
# `Vector{Float64}` buffer, so ForwardDiff Duals cannot flow through it -- exactly the same
# reason coeff space uses finite differences over its Float64-typed coeff bundle. The FD
# stencil at h=1e-3 reproduces the curvature without needing AD through that buffer.
#
# A coord's macro VALUE depends ONLY on its own log-parameter (value = 10^x_j), so its
# value-gradient is a scaled unit vector (ln(10)*value * e_j); `frac` is then just the
# projection of that unit direction onto the stiff eigen-subspace -- identical to coeff space.
# =========================================================================================

module ChaClassify

export cha_identifiable_functions, classify_cha

# Parent-relative imports so ChaClassify composes BOTH as a Main-level module (tests) AND as a
# FitRateEquation submodule (pipeline). The test header includes cha_fit.jl before this file.
using ..ChaFit
using ..ChaLaws
using ..FitRateEquation
using LinearAlgebra
using Statistics: median

# Indices of cha_coords that are NOT pinned (pinned coords are held fixed during the fit and
# excluded from the Hessian -- mirror coeff_identifiability.jl::_coord_unpinned_idx).
_cha_unpinned_idx(coord_syms, pins) = [j for (j, s) in enumerate(coord_syms) if !haskey(pins, s)]

# Central-difference Hessian of f at x over the given coordinate indices. Bit-identical to
# identifiability.jl::_fd_hessian (the stencil coeff_identifiability.jl reuses) -- copied here
# because that helper is a private FitRateEquation-internal not exported to this Main-level module.
function _cha_fd_hessian(f, x::Vector{Float64}, idx::Vector{Int}; h::Float64=1e-3)
    n = length(idx); H = zeros(n, n)
    e(k) = (v = zeros(length(x)); v[idx[k]] = h; v)
    for a in 1:n, b in a:n
        fpp = f(x .+ e(a) .+ e(b)); fpm = f(x .+ e(a) .- e(b))
        fmp = f(x .- e(a) .+ e(b)); fmm = f(x .- e(a) .- e(b))
        H[a,b] = H[b,a] = (fpp - fpm - fmp + fmm) / (4h^2)
    end
    H
end

"""
    cha_identifiable_functions(enzyme, mech, d, coords_dict; keq, pins=Dict(),
                               release_rate=<default>, release_eq=<default>, kr=nothing,
                               tol=2e-2) -> (eigvals, eigvecs, rank, idx)

Macro-coordinate identifiability at a Cha fit optimum: central-difference Hessian of
`cha_centered_logratio_loss` over the UNPINNED free cha_coords (in log10), eigendecomposed.
Returns the same shape as `coeff_identifiable_functions` (eigvals desc, eigvecs, rank, idx).
`tol` is the SAME practical-identifiability cutoff (stiff iff λ > tol·λmax).
"""
function cha_identifiable_functions(enzyme::Symbol, mech, d::Dataset, coords_dict;
                                    keq::Union{Nothing,Real}=nothing,
                                    pins::AbstractDict=Dict{Symbol,Float64}(),
                                    release_rate::Union{Nothing,Real}=nothing,
                                    release_eq::Union{Nothing,Real}=nothing,
                                    kr::Union{Nothing,Real}=nothing, tol::Real=2e-2,
                                    variant::Symbol=:_deploy)
    coord_syms = cha_coords(enzyme, variant)
    x = [log10(coords_dict[s]) for s in coord_syms]
    idx = _cha_unpinned_idx(coord_syms, pins)
    f = xv -> begin
        cd = Dict(coord_syms[k] => 10.0^xv[k] for k in eachindex(coord_syms))
        _cha_loss(enzyme, mech, d, cd; keq=keq,
                  release_rate=release_rate, release_eq=release_eq, kr=kr, variant=variant)
    end
    H = _cha_fd_hessian(f, collect(float.(x)), idx)
    F = eigen(Symmetric(H))
    order = sortperm(F.values; rev=true)
    vals = F.values[order]; vecs = F.vectors[:, order]
    λmax = maximum(abs, vals)
    rank = λmax == 0 ? 0 : count(>(tol*λmax), abs.(vals))
    (eigvals=vals, eigvecs=vecs, rank=rank, idx=idx)
end

# Thread the optional release_rate/release_eq/kr through to cha_centered_logratio_loss only
# when provided; otherwise let the loss apply its per-enzyme defaults (which depend on the
# coords dict, so they must NOT be precomputed against a stale dict).
function _cha_loss(enzyme, mech, d, cd; keq, release_rate, release_eq, kr, variant::Symbol=:_deploy)
    kwargs = Dict{Symbol,Any}(:keq => keq, :kr => kr, :variant => variant)
    release_rate !== nothing && (kwargs[:release_rate] = release_rate)
    release_eq   !== nothing && (kwargs[:release_eq]   = release_eq)
    cha_centered_logratio_loss(enzyme, mech, d, cd; kwargs...)
end

"""
    classify_cha(enzyme, mech, d, coords_dict, pins, idf; keq, sigma2=1.0,
                 variant=nothing, mode=nothing, stiff_frac=0.8, ci_rel_tol=1.0,
                 release_rate=<default>, release_eq=<default>, kr=nothing) -> Vector{NamedTuple}

Classify each cha_coord at a Cha fit optimum. Fields `(name, value, class, ci)`. Mirrors
`classify_coords`'s CI/stiff-subspace logic and class SYMBOLS, parameter space = cha_coords:
  - `:literature_pinned` -- `s` is in `pins` (held at its literature value during the fit);
  - `:data_identified` iff the coord's value-gradient lies mostly in the stiff eigen-subspace
    (`frac >= stiff_frac`) AND its relative CI is within `ci_rel_tol`; CI = the calibrated
    `sqrt(2·σ̂²·grad' H⁺ grad)`. Otherwise `:unconstrained`. CI reported only for
    `:data_identified`, NaN otherwise.

`variant`/`mode` are accepted (so the Task-15 call site can pass them). `classify_cha` itself
does NOT re-evaluate the Cha loss (it consumes the precomputed `idf`), so `variant` is inert in
THIS function — the variant-dependent curvature (HK1 alpha per H1/H3) must be applied upstream by
passing the same `variant` to `cha_identifiable_functions`. `mode` likewise carries no
Cha-coord machinery (override/conflated/product-side-literature are coeff-space-only).
"""
function classify_cha(enzyme::Symbol, mech, d::Dataset, coords_dict, pins, idf;
                      keq::Real=median(d.keq), sigma2::Real=1.0, variant=:_deploy, mode=nothing,
                      stiff_frac::Real=0.8, ci_rel_tol::Real=1.0,
                      release_rate::Union{Nothing,Real}=nothing,
                      release_eq::Union{Nothing,Real}=nothing, kr::Union{Nothing,Real}=nothing)
    coord_syms = cha_coords(enzyme, variant)
    idx   = idf.idx
    stiff = idf.eigvecs[:, 1:idf.rank]
    # Calibrated pseudo-inverse covariance on the UNPINNED coords (same form as classify_coords).
    Σ = idf.rank == 0 ? zeros(length(idx), length(idx)) :
        idf.eigvecs[:, 1:idf.rank] * Diagonal(1 ./ idf.eigvals[1:idf.rank]) * idf.eigvecs[:, 1:idf.rank]'
    # Position of each unpinned coord within the idx-restricted (Hessian/eigvec) subspace.
    pos_in_idx = Dict(j => p for (p, j) in enumerate(idx))

    out = NamedTuple[]
    for (j, s) in enumerate(coord_syms)
        value = coords_dict[s]
        if haskey(pins, s)
            push!(out, (name=s, value=value, class=:literature_pinned, ci=NaN)); continue
        end
        # Data-driven coord: value = 10^x_j, so the value-gradient w.r.t. the log10-coords is
        # ln(10)·value in the s-direction and zero elsewhere. Restricted to idx, gsub is zero
        # except at this coord's position within idx.
        gsub = zeros(length(idx))
        p = get(pos_in_idx, j, 0)
        if p != 0
            gsub[p] = log(10) * value
        end
        frac = norm(gsub) == 0 ? 0.0 : norm(stiff' * gsub) / norm(gsub)
        ci   = sqrt(2 * sigma2 * max(gsub' * Σ * gsub, 0.0))
        rel  = (isfinite(ci) && value != 0) ? ci / abs(value) : Inf
        cls  = (frac >= stiff_frac && rel <= ci_rel_tol) ? :data_identified : :unconstrained
        push!(out, (name=s, value=value, class=cls,
                    ci=(cls === :data_identified ? ci : NaN)))
    end
    out
end

end # module ChaClassify

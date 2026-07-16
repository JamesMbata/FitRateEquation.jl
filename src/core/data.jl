# Corpus data loading, vendored from fitting/mechanism_id/{types,run}.jl during the
# EnzymeRates upstream migration (self-containment). `Dataset`/`nrows`/`load_dataset`/
# `_to_float` are the only data-layer names consensus_macro uses.

# Dataset: one figure-grouped kinetic corpus, units already in Molar.
# `concs` is parameterized on the row NamedTuple type `T` so that `load_dataset` can build
# a CONCRETELY-typed vector (all rows of a corpus share one metabolite key set) -- indexing
# `d.concs[i]` then returns a concrete NamedTuple instead of forcing dynamic dispatch on
# every row of every loss evaluation. `Dataset` (unparameterized) still matches any `T` as a
# type annotation (`d::Dataset`), so existing call sites are unaffected.
struct Dataset{T<:NamedTuple}
    concs::Vector{T}           # per row, keyed by metabolite symbol -> concentration (M)
    rate::Vector{Float64}      # measured rate (arbitrary per-figure units)
    group::Vector{String}      # (Article, Fig) group key per row
    keq::Vector{Float64}       # per-row apparent Keq
end
nrows(d::Dataset) = length(d.rate)

# Coerce a CSV cell to Float64; missing/blank/non-numeric -> `default`.
function _to_float(x, default::Float64=0.0)
    x === missing && return default
    if x isa Real
        return Float64(x)
    end
    if x isa AbstractString
        s = strip(x)
        isempty(s) && return default
        v = tryparse(Float64, s)
        return v === nothing ? default : v
    end
    return default
end

"Load a Dataset from a config: read CSV, drop zero/blank/non-finite rates, convert
 µM->M, build per-row concs/group/keq. Group key is Article|Fig."
function load_dataset(cfg)
    df = CSV.read(cfg.data_csv, DataFrame)
    metsyms = collect(keys(cfg.metabolites))
    # One concrete row type per config's metabolite key set (G6PD/PGD/HK1 etc. each get
    # their own `T` -- this is NOT a single global type).
    T = NamedTuple{Tuple(metsyms), NTuple{length(metsyms),Float64}}
    concs = T[]; rate = Float64[]; grp = String[]; keq = Float64[]
    for row in eachrow(df)
        r = _to_float(row[cfg.rate_col], NaN)
        (isfinite(r) && r != 0.0) || continue   # drop zero / blank / non-finite rates
        vals = map(metsyms) do s
            col, unit = cfg.metabolites[s]
            x = _to_float(row[col], 0.0)         # missing concentrations -> 0.0
            unit === :uM ? x / 1e6 : x
        end
        push!(concs, T(Tuple(vals)))
        push!(rate, r)
        push!(grp, string(row[cfg.article_col], "|", row[cfg.fig_col]))
        push!(keq, _to_float(row[cfg.keq_col], NaN))
    end
    Dataset(concs, rate, grp, keq)
end

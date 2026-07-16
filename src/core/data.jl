# Corpus data loading, vendored from fitting/mechanism_id/{types,run}.jl during the
# EnzymeRates upstream migration (self-containment). `Dataset`/`nrows`/`load_dataset`/
# `_to_float` are the only data-layer names consensus_macro uses.

# Dataset: one figure-grouped kinetic corpus, units already in Molar.
struct Dataset
    concs::Vector{NamedTuple}  # per row, keyed by metabolite symbol -> concentration (M)
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
    concs = NamedTuple[]; rate = Float64[]; grp = String[]; keq = Float64[]
    metsyms = collect(keys(cfg.metabolites))
    for row in eachrow(df)
        r = _to_float(row[cfg.rate_col], NaN)
        (isfinite(r) && r != 0.0) || continue   # drop zero / blank / non-finite rates
        vals = map(metsyms) do s
            col, unit = cfg.metabolites[s]
            x = _to_float(row[col], 0.0)         # missing concentrations -> 0.0
            unit === :uM ? x / 1e6 : x
        end
        push!(concs, NamedTuple{Tuple(metsyms)}(Tuple(vals)))
        push!(rate, r)
        push!(grp, string(row[cfg.article_col], "|", row[cfg.fig_col]))
        push!(keq, _to_float(row[cfg.keq_col], NaN))
    end
    Dataset(concs, rate, grp, keq)
end

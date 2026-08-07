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

"Read the corpus CSV into the canonical fit DataFrame: one Molar column per metabolite
 symbol, plus Rate / source (Article|Fig) / Apparent_Keq, with zero, blank and non-finite
 rate rows dropped. `X_axis_label` is carried through when the corpus has it (G6PD/PGD do;
 HK1 does not) — the per-figure plot renderer needs it, the fit does not.

 This is the SINGLE corpus loader: `load_dataset` and the plotter's snapshot both come from
 it, so they cannot drift apart."
function read_corpus(cfg)
    raw = CSV.read(cfg.data_csv, DataFrame)
    df  = DataFrame()
    for s in collect(keys(cfg.metabolites))
        col, unit = cfg.metabolites[s]
        vals = _to_float.(raw[!, col], 0.0)      # missing concentrations -> 0.0
        df[!, s] = unit === :uM ? vals ./ 1e6 : vals
    end
    df.Rate         = _to_float.(raw[!, cfg.rate_col], NaN)
    df.source       = string.(raw[!, cfg.article_col], "|", raw[!, cfg.fig_col])
    df.Apparent_Keq = _to_float.(raw[!, cfg.keq_col], NaN)
    hasproperty(raw, :X_axis_label) && (df.X_axis_label = string.(raw[!, "X_axis_label"]))
    filter!(r -> isfinite(r.Rate) && r.Rate != 0.0, df)   # same drop as load_dataset
    return df
end

"Build a Dataset from a canonical corpus DataFrame (the output of `read_corpus`, optionally
 row-filtered). Rebuilds the same concretely-typed row NamedTuple vector `load_dataset`
 built, so `d.concs` stays concrete and loss evaluation does not fall back to dynamic
 dispatch."
function dataset_from_corpus(df::DataFrame, cfg)
    metsyms = collect(keys(cfg.metabolites))
    T = NamedTuple{Tuple(metsyms), NTuple{length(metsyms),Float64}}
    concs = T[T(Tuple(Float64(row[s]) for s in metsyms)) for row in eachrow(df)]
    Dataset(concs, Vector{Float64}(df.Rate), Vector{String}(df.source),
            Vector{Float64}(df.Apparent_Keq))
end

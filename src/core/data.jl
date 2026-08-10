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
 µM->M, build per-row concs/group/keq. Group key is Article|Fig.

 Thin composition of `read_corpus` (CSV -> canonical DataFrame) and `dataset_from_corpus`.
 Signature and values are unchanged; the split exists so `fit_consensus_equation` can snapshot the exact
 rows it fit (`fit_corpus.csv`) without a second, drift-prone loader."
load_dataset(cfg) = dataset_from_corpus(read_corpus(cfg), cfg)

# Columns `read_corpus` writes itself, after the metabolite loop. A config whose metabolite
# key collides with one of these would have its data silently overwritten, so the collision
# is rejected up front — where the offending key can still be named.
const _RESERVED_CORPUS_COLS = (:Rate, :source, :Apparent_Keq, :X_axis_label)

# THE metabolite ordering. `cfg.metabolites` is an unordered Dict, and three places
# independently derived an order from it: the corpus column order, the `d.concs` NamedTuple
# TYPE PARAMETER, and the renderer's adapter sweep order. They agreed only because they all
# called keys() on the same Dict in one process, and none was stable across Julia versions.
# Sorting once here makes all three stable and provably identical to each other.
metabolite_syms(cfg) = sort(collect(keys(cfg.metabolites)))

"Read the corpus CSV into the canonical fit DataFrame: one Molar column per metabolite
 symbol, plus Rate / source (Article|Fig) / Apparent_Keq, with zero, blank and non-finite
 rate rows dropped. `X_axis_label` is carried through when the corpus has it (G6PD/PGD do;
 HK1 does not) — the per-figure plot renderer needs it, the fit does not.

 This is the SINGLE corpus loader: `load_dataset` and the plotter's snapshot both come from
 it, so they cannot drift apart."
function read_corpus(cfg)
    raw = CSV.read(cfg.data_csv, DataFrame)
    for s in metabolite_syms(cfg)
        s in _RESERVED_CORPUS_COLS && error(
            "read_corpus: config metabolite key :$s collides with a reserved corpus " *
            "column ($(join(_RESERVED_CORPUS_COLS, ", "))). Rename the metabolite symbol; " *
            "read_corpus writes these four itself and would overwrite it.")
    end
    df  = DataFrame()
    for s in metabolite_syms(cfg)
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
function dataset_from_corpus(df::AbstractDataFrame, cfg)
    metsyms = metabolite_syms(cfg)
    T = NamedTuple{Tuple(metsyms), NTuple{length(metsyms),Float64}}
    concs = T[T(Tuple(Float64(row[s]) for s in metsyms)) for row in eachrow(df)]
    Dataset(concs, Vector{Float64}(df.Rate), Vector{String}(df.source),
            Vector{Float64}(df.Apparent_Keq))
end

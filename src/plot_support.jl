# ##########################################################################################
#      Non-Makie helpers for the consensus-macro Cha fit-vs-data plotter (Task 8 split)      #
# ##########################################################################################
#
# Ported from PPP_Experiments/fitting/consensus_macro/plot_consensus_fit.jl. These helpers
# are pure string/DataFrames/CSV/EnzymeRates logic -- none touch CairoMakie -- so they live
# in the main module and are reachable (as FitRateEquation.detect_enzyme, etc.) without
# loading the CairoMakie extension. The rendering loop itself (the only CairoMakie-dependent
# code, `plot_fit_on_data` + `FitRateEquation.plot_consensus_fit`) lives in
# ext/FitRateEquationMakieExt.jl and calls these helpers by qualified name.
# ##########################################################################################

# ------------------------------------------------------------------------------------------
#                                  ENZYME / CONFIG LOOKUP
# ------------------------------------------------------------------------------------------

function config_for(enzyme::Symbol)
    enzyme === :G6PD ? g6pd_config() :
    enzyme === :PGD  ? pgd_config()  :
    enzyme === :HK1  ? hk1_config()  :
    error("config_for: unknown enzyme $enzyme (expected :G6PD, :PGD, or :HK1)")
end

"Detect the enzyme from a results path: the component directly under `fitting/`."
function detect_enzyme(results_dir::AbstractString)
    m = match(r"(?:^|/)fitting/([^/]+)/", normpath(abspath(results_dir)) * "/")
    m === nothing && error("detect_enzyme: no `fitting/<ENZYME>/` segment in $results_dir")
    enz = Symbol(m.captures[1])
    enz in (:G6PD, :PGD, :HK1) ||
        error("detect_enzyme: unrecognized enzyme `$enz` under fitting/ in $results_dir")
    return enz
end

# ------------------------------------------------------------------------------------------
#                          READ FITTED cha_coords FROM macro_constants
# ------------------------------------------------------------------------------------------
# macro_constants.csv holds one row per (variant, mode, name): the classed cha_coords PLUS
# derived readoffs (Km_G6P / Km_PGA / kcatKm_*). Keep only the actual cha_coords for this
# (enzyme, variant); the derived readoffs are NOT law inputs.

function read_coords(mc::DataFrame, enzyme::Symbol, variant::Symbol, mode::Symbol)
    want = ChaFit.cha_coords(enzyme, variant)
    sub  = mc[(mc.variant .== String(variant)) .& (mc.mode .== String(mode)), :]
    coords = Dict{Symbol,Float64}()
    for s in want
        rows = sub[sub.name .== String(s), :]
        nrow(rows) == 1 ||
            error("read_coords: expected exactly 1 row for coord $s in $variant/$mode, got $(nrow(rows))")
        coords[s] = Float64(rows.value[1])
    end
    return coords
end

# ------------------------------------------------------------------------------------------
#                        CHA ADAPTER OVER THE RENDERER'S 2-CALL INTERFACE
# ------------------------------------------------------------------------------------------
# plot_fit_on_data(mech, ...) touches the mechanism only through metabolites(mech) and
# rate_equation(mech, concs, params). ChaAdapter satisfies both from the deployed Cha law
# at CHA_DEPLOY_RELEASE_RATE. It honors the renderer's PER-FIGURE keq: the renderer is
# called with Keq = :Apparent_Keq, so it passes each figure's apparent keq as params.Keq;
# rate_equation rebuilds the macro tuple's Haldane kr for that keq. The Cha FORWARD rate is
# keq-independent (kr enters only when products are present), so forward panels are
# unaffected; reverse panels use the figure's own apparent keq. `default_keq` is a fallback
# only for callers that pass no Keq (the renderer always does via the :Apparent_Keq path).

struct ChaAdapter
    enzyme::Symbol
    coords::Dict{Symbol,Float64}
    variant::Symbol
    metab_syms::Vector{Symbol}
    default_keq::Float64
end

EnzymeRates.metabolites(a::ChaAdapter) = Tuple(a.metab_syms)

# The deploy-fiber macro tuple for one figure's apparent keq.
_cha_adapter_tuple(a::ChaAdapter, keq::Real) =
    ChaFit.cha_macro_tuple(a.enzyme, a.coords;
        keq = keq,
        release_rate = ChaFit.CHA_DEPLOY_RELEASE_RATE,
        variant = a.variant)

function EnzymeRates.rate_equation(a::ChaAdapter, concs, params)
    ratefn = a.enzyme === :G6PD ? ChaLaws.cha_rate_G6PD :
             a.enzyme === :PGD  ? ChaLaws.cha_rate_PGD  :
             a.enzyme === :HK1  ? ChaLawsHK1.cha_rate_HK1 :
             error("ChaAdapter rate_equation: unknown enzyme $(a.enzyme)")
    keq = hasproperty(params, :Keq) ? Float64(params.Keq) : a.default_keq
    m   = _cha_adapter_tuple(a, keq)
    return ratefn(m; ChaFit._cha_row_kwargs(a.enzyme, concs)...)
end

"Build a ChaAdapter for one (variant, mode). Curves are drawn at the DEPLOY release rate;
keq is supplied per figure by the renderer (:Apparent_Keq path). `default_keq` is only a
fallback for a caller that passes no Keq."
function build_cha_adapter(enzyme::Symbol, coords::AbstractDict, variant::Symbol,
                           default_keq::Real)
    metab_syms = collect(keys(config_for(enzyme).metabolites))
    return ChaAdapter(enzyme, Dict{Symbol,Float64}(coords), variant, metab_syms,
                      Float64(default_keq))
end

# ------------------------------------------------------------------------------------------
#                       BUILD THE PLOTTING DATAFRAME FROM THE CORPUS CSV
# ------------------------------------------------------------------------------------------
# Mirror load_dataset: read the corpus CSV, convert µM->M, drop zero/blank/
# non-finite rates. Add the renderer's required columns (Rate, source, X_axis_label,
# Apparent_Keq). Column names/units come from the same config the fit used. `cfg.data_csv`
# is already an absolute path (resolved via pkgdir(FitRateEquation) by the config
# constructors), so -- unlike the source script -- no repo-root prefix is needed here.

_to_float(x) = x === missing ? NaN :
    (x isa AbstractString ? something(tryparse(Float64, x), NaN) : Float64(x))

function build_plot_df(cfg)
    raw = CSV.read(cfg.data_csv, DataFrame)
    hasproperty(raw, :X_axis_label) || error(
        "build_plot_df: $(cfg.name) corpus ($(cfg.data_csv)) has no `X_axis_label` column, " *
        "which the per-figure panel renderer requires. $(cfg.name) is not supported yet.")
    df = DataFrame()
    for (sym, (col, unit)) in cfg.metabolites
        vals = _to_float.(raw[!, col])
        vals[isnan.(vals)] .= 0.0                       # missing concentration -> absent (0)
        df[!, sym] = unit === :uM ? vals ./ 1e6 : vals
    end
    df.Rate         = _to_float.(raw[!, cfg.rate_col])
    df.source       = string.(raw[!, cfg.article_col], "|", raw[!, cfg.fig_col])
    df.X_axis_label = string.(raw[!, "X_axis_label"])
    df.Apparent_Keq = _to_float.(raw[!, cfg.keq_col])
    filter!(row -> isfinite(row.Rate) && row.Rate != 0.0, df)   # same drop as load_dataset
    return df
end

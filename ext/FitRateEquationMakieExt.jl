# ##########################################################################################
#         CairoMakie rendering for the consensus-macro Cha fit-vs-data plotter               #
# ##########################################################################################
#
# Ported from PPP_Experiments/fitting/consensus_macro/plot_consensus_fit.jl (the render loop)
# and fitting/G6PD/rate_eq/plot_fit_on_data.jl (the per-figure panel renderer). This is a
# package extension: it loads ONLY when CairoMakie is loaded alongside FitRateEquation, and
# it implements the `FitRateEquation.plot_consensus_fit` stub (defined unconditionally in the
# main module). All non-Makie helpers (detect_enzyme, config_for, read_coords,
# build_cha_adapter, build_plot_df) live in the main module (src/plot_support.jl); this file
# calls them by qualified name and contains ONLY CairoMakie-dependent code.
#
# Predictions come from the DEPLOYED Cha law: the fitted macro constants are read from
# <run_dir>/macro_constants.csv, assembled with FitRateEquation.ChaFit.cha_macro_tuple at
# CHA_DEPLOY_RELEASE_RATE, and evaluated through FitRateEquation.ChaLaws.cha_rate_*, via the
# ChaAdapter defined in plot_support.jl.
#
# HK1 is scoped out: its corpus has no X_axis_label column (needed to pick each panel's
# swept metabolite), so detecting :HK1 and calling build_plot_df raises a clear error.
# ##########################################################################################

module FitRateEquationMakieExt

using FitRateEquation, EnzymeRates, CSV, DataFrames, Printf, Statistics, CairoMakie

# ------------------------------------------------------------------------------------------
#                    PER-FIGURE PANEL RENDERER (ported from plot_fit_on_data.jl)
# ------------------------------------------------------------------------------------------
# Adapted from DataDrivenEnzymeRateEqs.jl for use with EnzymeRates.jl mechanisms.
# One panel per source figure, scatter for data points, lines for fitted curves.
#
# Arguments:
#   mechanism     -- EnzymeMechanism-like instance (here: a FitRateEquation.ChaAdapter)
#   fitted_params -- NamedTuple of fitted parameters (real space)
#   data          -- DataFrame with columns: Rate, source, X_axis_label,
#                    and one column per metabolite in metabolites(mechanism)
#   Keq           -- equilibrium constant passed through to rate_equation.
#                    A scalar Real applies one value to every point; a Symbol
#                    names a per-row column in `data` (e.g. :Apparent_Keq) so a
#                    pH-corrected apparent K'eq can vary per figure.
#
# Keyword arguments:
#   enzyme_name     -- label for y-axis (default "")
#   num_col         -- number of columns in the figure grid (default 5)
#   scaler          -- global size scaling factor (default 4.0)
#   absolute_rates  -- if true, plot raw predicted rates without per-figure
#                      Vmax correction (default false)

function plot_fit_on_data(
    mechanism,
    fitted_params::NamedTuple,
    data::DataFrame,
    Keq::Union{Real,Symbol};
    enzyme_name::String = "",
    num_col::Int = 5,
    scaler = 4.0,
    absolute_rates::Bool = false,
)
    # Scalar Keq: one params NamedTuple for every point. Symbol Keq: build params
    # per point/series from the named per-row column. `_params_for` centralizes
    # both paths so the scalar behaviour is byte-for-byte unchanged.
    keq_is_col = Keq isa Symbol
    _params_for(keq_val) = merge(fitted_params, (Keq = keq_val, E_total = 1.0))
    params_full = keq_is_col ? nothing : _params_for(Keq)
    MetNames = EnzymeRates.metabolites(mechanism)

    fontsize  = scaler * 5
    markersize = scaler * 3
    linewidth  = scaler * 1
    set_theme!(Theme(
        fontsize = fontsize,
        Axis = (
            titlefont = :regular,
            xticksize = scaler * 1,
            yticksize = scaler * 1,
            xlabelfont = :bold,
            ylabelfont = :bold,
            yticklabelpad = scaler * 1,
            ylabelpadding = scaler * 3,
        ),
        Legend = (titlefont = :regular,),
    ))

    sources  = unique(data.source)
    num_rows = ceil(Int, length(sources) / num_col)
    size_pt  = 72 .* (scaler .* (7, num_rows))
    fig = Figure(size = size_pt)

    for (i, source) in enumerate(sources)
        row_idx = (i - 1) ÷ num_col + 1
        col_idx = (i - 1) % num_col + 1
        gl = fig[row_idx, col_idx] = GridLayout()

        fig_data = data[data.source .== source, :]
        x_axis_metabolite = Symbol(fig_data.X_axis_label[1])
        other_metabolites = [m for m in MetNames if m != x_axis_metabolite]

        # Per-figure Vmax correction: geometric mean of predicted/measured ratio
        if absolute_rates
            fig_vmax = 1.0
        else
            log_ratios = Float64[]
            for row in eachrow(fig_data)
                concs = NamedTuple{MetNames}(Tuple(row[col] for col in MetNames))
                params_pt = keq_is_col ? _params_for(row[Keq]) : params_full
                pred = EnzymeRates.rate_equation(mechanism, concs, params_pt)
                if row.Rate != 0.0 && sign(pred) == sign(row.Rate) && pred != 0.0
                    push!(log_ratios, log(pred / row.Rate))
                end
            end
            fig_vmax = isempty(log_ratios) ? 1.0 : exp(mean(log_ratios))
        end

        ax = Axis(gl[1, 1],
            title = replace(source, "_" => " "),
            xticks = LinearTicks(3),
            yticks = LinearTicks(3),
            limits = begin
                maximum(fig_data.Rate) > 0.0 ?
                    (nothing, (-0.05 * maximum(fig_data.Rate), nothing)) :
                    (nothing, (nothing, -0.05 * minimum(fig_data.Rate)))
            end,
            ytickformat = ys -> ["$(round(x/maximum(abs.(ys)), sigdigits=1))" for x in ys],
            xtickformat = xs -> ["$(round(x/1e-3, sigdigits=2))" for x in xs],
            yticklabelrotation = pi / 2,
            xlabelpadding = scaler * 0,
            ylabelpadding = scaler * 0,
            xticklabelpad = 0,
            yticklabelpad = 0,
            xlabel = _crop_to_vowel_end(string(x_axis_metabolite), 5) * ", mM",
            ylabel = col_idx == 1 ? "$(enzyme_name) Rate" : "",
        )

        # Identify constant vs changing metabolites across series in this panel
        const_metabs = [m for m in other_metabolites
                        if length(unique(fig_data[!, m])) == 1 &&
                           unique(fig_data[!, m])[1] != 0.0]
        changing_metabs = [m for m in other_metabolites
                          if length(unique(fig_data[!, m])) > 1]

        # One series per unique combination of non-x-axis metabolite concentrations
        unique_combos = unique(fig_data[!, other_metabolites])
        for combo_row in eachrow(unique_combos)
            # Filter data matching this combination
            mask = trues(nrow(fig_data))
            for m in other_metabolites
                mask .&= fig_data[!, m] .== combo_row[m]
            end
            series_data = fig_data[mask, :]

            # Per-row apparent K'eq is constant within a figure/series, so the
            # continuous prediction line uses this series' value (scalar path
            # reuses the single params NamedTuple).
            series_params = keq_is_col ? _params_for(series_data[1, Keq]) : params_full

            # Build prediction function for line plot
            function predict(x)
                concs = NamedTuple{MetNames}(Tuple(
                    m == x_axis_metabolite ? x : combo_row[m] for m in MetNames
                ))
                EnzymeRates.rate_equation(mechanism, concs, series_params) / fig_vmax
            end

            # Legend label from changing metabolite concentrations
            label = join([_unit_conc_format(combo_row[m]) for m in changing_metabs], ", ")
            isempty(label) && (label = " ")

            scatter!(ax, series_data[!, x_axis_metabolite], series_data.Rate,
                markersize = markersize, label = label)
            lines!(ax, 0 .. maximum(fig_data[!, x_axis_metabolite]),
                predict, linewidth = linewidth, label = label)
        end

        # Legend with constant metabolites in title, changing metabolites as labels
        legend_title = ""
        for m in const_metabs
            legend_title *= _crop_to_vowel_end(string(m), 5) * "=" *
                            _unit_conc_format(unique(fig_data[!, m])[1]) * "\n"
        end
        for (j, m) in enumerate(changing_metabs)
            legend_title *= _crop_to_vowel_end(string(m), 5)
            j < length(changing_metabs) && (legend_title *= ", ")
        end
        isempty(legend_title) && (legend_title = "No var metabs")

        Legend(gl[1, 2], ax, legend_title;
            merge = true, unique = true,
            labelsize = fontsize, titlesize = fontsize,
            patchsize = scaler .* (6.0f0, 5.0f0),
            patchlabelgap = scaler * 2,
            padding = scaler .* (0.5f0, 0.0f0, 0.0f0, 0.0f0),
            framevisible = false, rowgap = 0, titlegap = 0,
            titlevalign = :top, titlehalign = :left,
            valign = :top, halign = :left,
        )
        colgap!(gl, 1)
    end

    colgap!(fig.layout, scaler * 5)
    rowgap!(fig.layout, scaler * 5)
    return fig
end

function _crop_to_vowel_end(str, max_length::Int = 5)
    vowels = Set("aeiouyAEIOUY")
    length(str) > max_length && (str = str[1:max_length])
    while !isempty(str) && last(str) in vowels
        str = str[1:end-1]
    end
    return str
end

function _unit_conc_format(conc::Number)
    conc == 0.0 && return "----"
    conc >= 0.1   && return @sprintf("%.2f", round(conc, sigdigits=2)) * "M"
    conc >= 10e-3 && return @sprintf("%.0f", round(conc, sigdigits=2)/1e-3) * "mM"
    conc >= 1e-3  && return @sprintf("%.1f", round(conc, sigdigits=2)/1e-3) * "mM"
    conc >= 1e-6  && return @sprintf("%.0f", round(conc, sigdigits=2)/1e-6) * "uM"
    conc >= 1e-9  && return @sprintf("%.0f", round(conc, sigdigits=2)/1e-9) * "nM"
    conc >= 0.1e-9 && return @sprintf("%.0f", round(conc, sigdigits=2)/1e-12) * "pM"
    return "Number Out of Scale"
end

# ------------------------------------------------------------------------------------------
#                                          MAIN
# ------------------------------------------------------------------------------------------
# Implements the stub declared in the main module (`function plot_consensus_fit end`). Called
# as FitRateEquation.plot_consensus_fit(run_dir) -- active once CairoMakie is loaded, since
# Julia resolves extension methods onto the parent's generic function.

function FitRateEquation.plot_consensus_fit(run_dir::AbstractString)
    isdir(run_dir) || error("Results dir not found: $run_dir")
    mc_file = joinpath(run_dir, "macro_constants.csv")
    isfile(mc_file) || error("Missing macro_constants.csv in $run_dir")

    enzyme = FitRateEquation.detect_enzyme(run_dir)
    cfg    = FitRateEquation.config_for(enzyme)
    println("Results dir: $run_dir")
    println("Enzyme: $enzyme")

    df = FitRateEquation.build_plot_df(cfg)   # errors here for HK1 (no X_axis_label)
    println("Data: $(nrow(df)) rows, $(length(unique(df.source))) source figures")

    mc = CSV.read(mc_file, DataFrame)
    plots_dir = joinpath(run_dir, "plots")
    mkpath(plots_dir)

    cells = unique(collect(zip(mc.variant, mc.mode)))
    for (variant_s, mode_s) in cells
        variant = Symbol(variant_s); mode = Symbol(mode_s)
        out_png = joinpath(plots_dir, "$(variant_s)_$(mode_s)_fit.png")
        @printf("\n%s / %s\n", variant_s, mode_s)
        try
            coords  = FitRateEquation.read_coords(mc, enzyme, variant, mode)
            adapter = FitRateEquation.build_cha_adapter(enzyme, coords, variant, cfg.deploy_keq)
            fig = plot_fit_on_data(adapter, (;), df, :Apparent_Keq;
                                   enzyme_name = cfg.name, absolute_rates = false)
            save(out_png, fig, px_per_unit = 4)
            println("  saved: $out_png")
        catch e
            @warn "Plot failed for $variant_s/$mode_s" exception = (e, catch_backtrace())
        end
    end
    println("\nDone. Plots written to: $plots_dir")
    return plots_dir
end

end # module

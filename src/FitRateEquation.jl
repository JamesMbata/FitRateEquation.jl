module FitRateEquation

using EnzymeRates
using LinearAlgebra, Statistics, Random
using Distributed
using ForwardDiff
using CSV, DataFrames
using Dates
using Pkg, ClusterManagers

# Vendored core (data loading, mechanism builder, gauge). bounds/loss/fit/structural
# are NOT carried — the Cha path uses cha_coord_bounds / cha_centered_logratio_loss /
# cha_classify instead.
include("core/data.jl")
include("core/er_names.jl")
include("core/mechbuild.jl")
include("core/gauge.jl")

export Dataset, load_dataset, read_corpus, dataset_from_corpus, nrows
export gauge_param, free_params, build_params, analytic_kcat
export legacy_param_name, er_fitted_params, er_parameters, er_rate_equation
export er_dependent_param_exprs, er_raw_rate_polys, er_reactions, er_substrates
export _mechanism_steps, _deadend_forms, _deadend_step, _mech, _SIGN_PENALTY

include("enzyme_wiring.jl")
include("mechanisms.jl")
include("enzymes/g6pd.jl")
include("enzymes/pgd.jl")
include("cv.jl")
include("promotable.jl")
include("macro_collect.jl")
include("cha_laws.jl")
include("cha_invert.jl")
include("cha_fit.jl")
include("cha_classify.jl")
include("cha_deploy.jl")
include("cha_koffq_report.jl")
include("worker_setup.jl")
include("run.jl")
include("configs/G6PD.jl")
include("configs/PGD.jl")
include("plot_support.jl")
include("cli.jl")

export EnzymeWiring, register_enzyme!
export consensus_variants
export macro_constants
export mode_agreement
export fit_consensus_equation, write_outputs
export setup_workers
export run_g6pd, run_pgd
export cli_main

# plot_consensus_fit stub — the real method lives in the CairoMakie package extension.
function plot_consensus_fit end
export plot_consensus_fit

end # module

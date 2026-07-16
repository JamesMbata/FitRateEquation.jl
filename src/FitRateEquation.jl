module FitRateEquation

using EnzymeRates
using LinearAlgebra, Statistics, Random
using Distributed
using ForwardDiff
using CSV, DataFrames

# Vendored core (data loading, mechanism builder, gauge). bounds/loss/fit/structural
# are NOT carried — the Cha path uses cha_coord_bounds / cha_centered_logratio_loss /
# cha_classify instead.
include("core/data.jl")
include("core/mechbuild.jl")
include("core/gauge.jl")

export Dataset, load_dataset, nrows
export gauge_param, free_params, build_params, analytic_kcat
export _mechanism_steps, _deadend_forms, _deadend_step, _mech, _SIGN_PENALTY

include("enzyme_wiring.jl")
include("mechanisms.jl")
include("enzymes/g6pd.jl")
include("enzymes/pgd.jl")
# HK1 stays guarded: its mechanisms use bespoke low-level EnzymeMechanism construction +
# the allosteric DSL reworked upstream; the include is guarded so G6PD/PGD load, and HK1
# auto-re-enables once ported.
const HK1_AVAILABLE = try
    include("enzymes/hk1.jl"); true
catch err
    @warn "FitRateEquation: HK1 wiring disabled on this EnzymeRates (deferred port)" exception=(err, catch_backtrace())
    false
end
include("cv.jl")
include("macro_collect.jl")
include("cha_laws.jl")
include("cha_invert.jl")
include("cha_fit.jl")
include("cha_classify.jl")
include("cha_deploy.jl")
include("cha_koffq_report.jl")
include("run.jl")
include("configs/G6PD.jl")
include("configs/PGD.jl")
include("configs/HK1.jl")

export EnzymeWiring, register_enzyme!
export consensus_variants
export macro_constants
export mode_agreement
export run_variant, run_all, write_outputs
export g6pd_config, pgd_config, hk1_config

# plot_consensus_fit stub — the real method lives in the CairoMakie package extension.
function plot_consensus_fit end
export plot_consensus_fit

end # module

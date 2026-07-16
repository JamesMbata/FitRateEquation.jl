# ##########################################################################################
#                  Render check for the CairoMakie plotting extension                        #
# ##########################################################################################
# NOT included in runtests.jl / the default suite: `using CairoMakie` cold-precompiles for
# minutes, and this runs a real (smoke) G6PD fit. Run manually / in a dedicated CI job:
#   julia --project test/test_plot_render.jl
#
# Confirms FitRateEquation.plot_consensus_fit is ACTIVE once CairoMakie is loaded (the
# extension mechanism resolves the method onto the main-module stub) and writes at least one
# PNG into <run_dir>/plots/.

using Test, FitRateEquation, CairoMakie

@testset "plot render (CairoMakie)" begin
    out = mktempdir()
    run_g6pd(smoke=true, nprocs=1, outdir=out)
    plot_consensus_fit(out)
    @test isdir(joinpath(out, "plots"))
    @test !isempty(readdir(joinpath(out, "plots")))
end

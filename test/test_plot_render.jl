# ##########################################################################################
#                  Render check for the CairoMakie plotting extension                       #
# ##########################################################################################
# NOT included in runtests.jl / the default suite: `using CairoMakie` cold-precompiles for
# minutes, and this runs a real (smoke) G6PD fit. CairoMakie is a [weakdeps] entry, so it is
# not loadable from the main project env either — run this through the dedicated sub-env:
#   julia --project=test/render -e 'using Pkg; Pkg.instantiate()'   # first run on a checkout
#   julia --project=test/render test/test_plot_render.jl
# (that env commits no Manifest, so a fresh clone needs the instantiate step first)
#
# Confirms FitRateEquation.plot_consensus_fit is ACTIVE once CairoMakie is loaded (the
# extension mechanism resolves the method onto the main-module stub) and writes at least one
# PNG into <run_dir>/plots/.
#
# The fit deliberately uses a NON-DEFAULT corpus. The bug this file exists to guard against
# was the plotter re-deriving its corpus from a config instead of reading the run's own
# snapshot: against the default corpus that is invisible, because the re-derived corpus and
# the real one are identical. Fitting a subset makes the two distinguishable, so a
# regression shows up as wrong row counts here rather than as a silently wrong plot.

using Test, FitRateEquation, CairoMakie, CSV, DataFrames

@testset "plot render (CairoMakie, custom corpus)" begin
    bundled = CSV.read(g6pd_config().data_csv, DataFrame)
    # Three of the seven articles — three, not one, so leave-one-article-out CV still has
    # non-empty train and test sets per fold.
    keep = Set(["Wang2002", "Ozer2002", "Beutler1986"])
    sub  = bundled[in.(string.(bundled.Article), Ref(keep)), :]
    @test nrow(sub) < nrow(bundled)

    tmpcsv = joinpath(mktempdir(), "subset_corpus.csv")
    CSV.write(tmpcsv, sub)

    out = mktempdir()
    # `run_g6pd` takes no `data_csv` (only `run_g6pd_noatp` does), so the corpus override
    # goes in through the config, exactly as run_g6pd's own body would then run it: the
    # deploy variants, anchor_reverse and seed all keep their defaults, and the explicit
    # budget is the smoke budget (`_budget(true)`). No workers are added, so this is the
    # serial path — the `nprocs=1` of a `run_g6pd(smoke=true, nprocs=1)` call.
    fit_consensus_equation(g6pd_config(; data_csv=tmpcsv); outdir=out,
            n_restarts=2, maxiter=150, maxtime=120.0)

    # The plotter's own reader sees the run's corpus, not the bundled default.
    fc = FitRateEquation.read_fit_corpus(out)
    @test nrow(fc) > 0                       # `all`/`issubset` over empty is vacuously true
    @test nrow(fc) < nrow(bundled)
    @test Set(first.(split.(string.(fc.source), "|"))) ⊆ keep

    # …and so does plot_consensus_fit ITSELF. The three assertions above only exercise
    # read_fit_corpus; a plotter that went back to re-deriving the corpus from the config
    # would still satisfy them, and would still write a full plots/ dir — just one drawn
    # over the wrong 565 rows. Its announced row/figure counts are the only externally
    # visible evidence of which corpus it actually drew, so assert on those.
    logfile = joinpath(mktempdir(), "render_stdout.log")   # a file, not an IOBuffer:
    open(logfile, "w") do io                               # redirect_stdout needs a real fd
        redirect_stdout(io) do
            plot_consensus_fit(out)
        end
    end
    rendered = read(logfile, String)
    print(rendered)                                   # keep the per-cell log visible
    @test occursin("Data: $(nrow(fc)) rows, $(length(unique(fc.source))) source figures",
                   rendered)

    @test isdir(joinpath(out, "plots"))
    @test !isempty(readdir(joinpath(out, "plots")))
end

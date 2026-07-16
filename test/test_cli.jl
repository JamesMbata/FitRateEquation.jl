# In-process CLI: parsing + dispatch. No subprocess, no preflight, no real fits — dispatch
# is exercised via a RECORDING stub table (kwargs captured into a Ref), never the real
# run_g6pd/run_pgd/run_hk1/run_g6pd_noatp runners.
using Test, FitRateEquation
const CLI = FitRateEquation

@testset "cli parsing" begin
    @test CLI.parse_cli(String[])[1] == "help"
    @test CLI.parse_cli(["g6pd", "--smoke", "--nprocs", "4"]) ==
        ("g6pd", (smoke=true, nprocs=4, outdir=nothing, rundir=nothing, data=nothing))
    @test CLI.parse_cli(["pgd", "--outdir", "/tmp/x"])[2].outdir == "/tmp/x"
    @test CLI.parse_cli(["plot", "some/dir"]) ==
        ("plot", (smoke=false, nprocs=nothing, outdir=nothing, rundir="some/dir", data=nothing))
    @test_throws ErrorException CLI.parse_cli(["bogus"])
    @test_throws ErrorException CLI.parse_cli(["g6pd", "--nprocs", "0"])
    @test_throws ErrorException CLI.parse_cli(["plot"])

    # g6pd-noatp + --data (the Task-7 correction over the brief: a genuinely distinct
    # ATP-free fit with its own data-CSV override, not a `run_g6pd` alias).
    @test CLI.parse_cli(["g6pd-noatp"]) ==
        ("g6pd-noatp", (smoke=false, nprocs=nothing, outdir=nothing, rundir=nothing, data=nothing))
    @test CLI.parse_cli(["g6pd-noatp", "--data", "/tmp/noatp.csv", "--smoke"]) ==
        ("g6pd-noatp", (smoke=true, nprocs=nothing, outdir=nothing, rundir=nothing, data="/tmp/noatp.csv"))

    # --data is gated to g6pd-noatp: every other subcommand must reject it clearly.
    @test_throws ErrorException CLI.parse_cli(["g6pd", "--data", "/tmp/x.csv"])
    @test_throws ErrorException CLI.parse_cli(["pgd", "--data", "/tmp/x.csv"])
    @test_throws ErrorException CLI.parse_cli(["hk1", "--data", "/tmp/x.csv"])
    @test_throws ErrorException CLI.parse_cli(["plot", "some/dir", "--data", "/tmp/x.csv"])
    # --data with no value still errors (missing-value check), even though sub is gated OK.
    @test_throws ErrorException CLI.parse_cli(["g6pd-noatp", "--data"])
end

@testset "cli dispatch (recording)" begin
    rec = Ref{Any}(nothing)
    disp = Dict(
        :g6pd               => (; kw...) -> (rec[] = (:g6pd, kw); nothing),
        :pgd                => (; kw...) -> (rec[] = (:pgd, kw); nothing),
        :hk1                => (; kw...) -> (rec[] = (:hk1, kw); nothing),
        Symbol("g6pd-noatp") => (; kw...) -> (rec[] = (Symbol("g6pd-noatp"), kw); nothing),
        :plot               => (; kw...) -> (rec[] = (:plot, kw); nothing),
    )

    @test CLI.cli_main(["help"]; dispatch=disp) == 0
    @test rec[] === nothing   # help never touches dispatch

    @test CLI.cli_main(["g6pd", "--smoke"]; dispatch=disp) == 0
    @test rec[][1] == :g6pd
    @test rec[][2][:smoke] == true
    @test rec[][2][:nprocs] === nothing
    @test rec[][2][:outdir] === nothing
    @test !haskey(rec[][2], :data_csv)   # plain g6pd never gets a data_csv kwarg

    @test CLI.cli_main(["pgd", "--outdir", "/tmp/pgd_out"]; dispatch=disp) == 0
    @test rec[][1] == :pgd
    @test rec[][2][:outdir] == "/tmp/pgd_out"

    @test CLI.cli_main(["hk1", "--nprocs", "2"]; dispatch=disp) == 0
    @test rec[][1] == :hk1
    @test rec[][2][:nprocs] == 2

    @test CLI.cli_main(["plot", "some/run/dir"]; dispatch=disp) == 0
    @test rec[][1] == :plot
    @test rec[][2][:rundir] == "some/run/dir"
    @test length(rec[][2]) == 1

    # g6pd-noatp routes to its own dispatch entry with data_csv threaded through (nothing
    # when --data is omitted, the CSV path when given).
    @test CLI.cli_main(["g6pd-noatp", "--smoke"]; dispatch=disp) == 0
    @test rec[][1] == Symbol("g6pd-noatp")
    @test rec[][2][:smoke] == true
    @test rec[][2][:data_csv] === nothing

    @test CLI.cli_main(["g6pd-noatp", "--data", "/tmp/noatp.csv", "--outdir", "/tmp/o"]; dispatch=disp) == 0
    @test rec[][1] == Symbol("g6pd-noatp")
    @test rec[][2][:data_csv] == "/tmp/noatp.csv"
    @test rec[][2][:outdir] == "/tmp/o"
end

@testset "cli default dispatch wires run_g6pd_noatp" begin
    # No real fit runs here — just confirm the DEFAULT dispatch table (used by bin/fitrateequation)
    # maps g6pd-noatp to the actual exported run_g6pd_noatp function, not a placeholder alias.
    d = CLI._default_dispatch()
    @test d[Symbol("g6pd-noatp")] === FitRateEquation.run_g6pd_noatp
    @test d[:g6pd] === FitRateEquation.run_g6pd
    @test d[:pgd] === FitRateEquation.run_pgd
    @test d[:hk1] === FitRateEquation.run_hk1
end

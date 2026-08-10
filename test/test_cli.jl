# In-process CLI: argument parsing only. `cli_main` dispatches straight to
# `fit_consensus_equation`, so there is no injection seam to stub — the parser tests below
# cover the whole CLI surface without spawning a fit.
using Test
using FitRateEquation: parse_cli

@testset "cli parsing" begin
    @test parse_cli(String[])[1] == "help"
    @test parse_cli(["g6pd", "--smoke", "--nprocs", "4"]) ==
        ("g6pd", (smoke=true, nprocs=4, outdir=nothing, rundir=nothing, data=nothing, variant=nothing))
    @test parse_cli(["pgd", "--outdir", "/tmp/x"])[2].outdir == "/tmp/x"
    @test parse_cli(["plot", "some/dir"]) ==
        ("plot", (smoke=false, nprocs=nothing, outdir=nothing, rundir="some/dir", data=nothing, variant=nothing))
    @test_throws ErrorException parse_cli(["bogus"])
    @test_throws ErrorException parse_cli(["g6pd", "--nprocs", "0"])
    @test_throws ErrorException parse_cli(["plot"])
end

@testset "cli --variant" begin
    sub, o = parse_cli(["g6pd", "--variant", "no_atp", "--smoke"])
    @test sub == "g6pd"
    @test o.variant == "no_atp"
    @test o.smoke

    # g6pd-noatp is gone
    @test_throws ErrorException parse_cli(["g6pd-noatp"])

    # --data is valid on any enzyme now
    sub, o = parse_cli(["pgd", "--data", "/tmp/x.csv"])
    @test o.data == "/tmp/x.csv"

    # --variant is valid on every enzyme subcommand, rejected elsewhere
    @test parse_cli(["hk1", "--variant", "full_re"])[2].variant == "full_re"
    @test_throws ErrorException parse_cli(["plot", "some/dir", "--variant", "no_atp"])

    # --data is valid on every enzyme subcommand, rejected on plot
    @test parse_cli(["hk1", "--data", "/tmp/h.csv"])[2].data == "/tmp/h.csv"
    @test_throws ErrorException parse_cli(["plot", "some/dir", "--data", "/tmp/x.csv"])

    # missing-value checks for the new flags
    @test_throws ErrorException parse_cli(["g6pd", "--variant"])
    @test_throws ErrorException parse_cli(["g6pd", "--data"])
end

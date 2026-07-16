# In-process CLI dispatcher. No subprocess / --project / preflight — the package is
# already loaded, so subcommands call the exported runners directly.
const _CLI_SUBS = ("g6pd", "pgd", "g6pd-noatp", "hk1", "plot", "help")

const CLI_USAGE = """
FitRateEquation — consensus rate-equation fitter (G6PD / PGD / HK1)

Usage: fitrateequation <subcommand> [flags]
  g6pd | pgd | hk1                Fit an enzyme (writes artifacts to --outdir)
  g6pd-noatp                      Fit the ATP-free G6PD variant (writes artifacts to --outdir)
  plot <run_dir>                  Render the fitted law over the corpus (needs CairoMakie)
  help                            Show this message
Flags: --smoke  --nprocs N  --outdir DIR  --data CSV (g6pd-noatp only, overrides the corpus)
"""

# `data` is only ever set for `g6pd-noatp` (see the `--data` gate below); every other
# subcommand carries `data=nothing` so all subs share one option NamedTuple shape.
_EMPTY_OPTS = (smoke=false, nprocs=nothing, outdir=nothing, rundir=nothing, data=nothing)

function parse_cli(argv::AbstractVector{<:AbstractString})
    isempty(argv) && return ("help", _EMPTY_OPTS)
    sub = String(argv[1])
    (sub in ("-h", "--help", "help")) && return ("help", _EMPTY_OPTS)
    sub in _CLI_SUBS || error("unknown subcommand: $sub\n\n$CLI_USAGE")
    smoke = false; nprocs = nothing; outdir = nothing; rundir = nothing; data = nothing
    i = 2
    while i <= length(argv)
        tok = String(argv[i])
        if tok == "--smoke"
            smoke = true; i += 1
        elseif tok == "--nprocs"
            i < length(argv) || error("--nprocs requires a value")
            n = tryparse(Int, argv[i+1])
            (n === nothing || n < 1) && error("--nprocs must be a positive integer")
            nprocs = n; i += 2
        elseif tok == "--outdir"
            i < length(argv) || error("--outdir requires a value")
            outdir = String(argv[i+1]); i += 2
        elseif tok == "--data"
            # `--data` overrides the fitted corpus and only makes sense for the ATP-free
            # variant (run_g6pd_noatp's data_csv kwarg) — reject it clearly elsewhere
            # rather than silently ignoring it on g6pd/pgd/hk1/plot.
            sub == "g6pd-noatp" || error("--data is only valid with g6pd-noatp\n\n$CLI_USAGE")
            i < length(argv) || error("--data requires a value")
            data = String(argv[i+1]); i += 2
        elseif startswith(tok, "-")
            error("unknown flag: $tok\n\n$CLI_USAGE")
        elseif sub == "plot" && rundir === nothing
            rundir = tok; i += 1
        else
            error("unexpected argument: $tok")
        end
    end
    sub == "plot" && rundir === nothing && error("plot requires a <run_dir>\n\n$CLI_USAGE")
    return (sub, (smoke=smoke, nprocs=nprocs, outdir=outdir, rundir=rundir, data=data))
end

_default_dispatch() = Dict(
    :g6pd => run_g6pd, :pgd => run_pgd, Symbol("g6pd-noatp") => run_g6pd_noatp,
    :hk1 => run_hk1, :plot => (; rundir, kw...) -> plot_consensus_fit(rundir),
)

function cli_main(argv::AbstractVector{<:AbstractString}; dispatch=_default_dispatch())
    sub, o = parse_cli(argv)
    sub == "help" && (print(CLI_USAGE); return 0)
    fn = dispatch[Symbol(sub)]
    if sub == "plot"
        fn(; rundir=o.rundir)
    elseif sub == "g6pd-noatp"
        fn(; smoke=o.smoke, nprocs=o.nprocs, outdir=o.outdir, data_csv=o.data)
    else
        fn(; smoke=o.smoke, nprocs=o.nprocs, outdir=o.outdir)
    end
    return 0
end

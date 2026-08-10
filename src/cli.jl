# In-process CLI dispatcher. No subprocess / --project / preflight — the package is
# already loaded, so subcommands call fit_consensus_equation directly.
const _CLI_SUBS = ("g6pd", "pgd", "hk1", "plot", "help")

const CLI_USAGE = """
FitRateEquation — consensus rate-equation fitter (G6PD / PGD / HK1)

Usage: fitrateequation <subcommand> [flags]
  g6pd | pgd | hk1                Fit an enzyme (writes artifacts to --outdir)
  plot <run_dir>                  Render the fitted law over the corpus (needs CairoMakie)
  help                            Show this message
Flags: --smoke  --nprocs N  --outdir DIR  --data CSV  --variant NAME
  --variant NAME  Fit an alternative rate law (e.g. no_atp, full_re, no_g6p_atp_deadend)
  --data CSV      Fit your own corpus (canonical columns required)
"""

_EMPTY_OPTS = (smoke=false, nprocs=nothing, outdir=nothing, rundir=nothing, data=nothing, variant=nothing)

function parse_cli(argv::AbstractVector{<:AbstractString})
    isempty(argv) && return ("help", _EMPTY_OPTS)
    sub = String(argv[1])
    (sub in ("-h", "--help", "help")) && return ("help", _EMPTY_OPTS)
    sub in _CLI_SUBS || error("unknown subcommand: $sub\n\n$CLI_USAGE")
    smoke = false; nprocs = nothing; outdir = nothing; rundir = nothing; data = nothing; variant = nothing
    i = 2
    while i <= length(argv)
        tok = String(argv[i])
        if tok == "--smoke"
            smoke = true; i += 1
        elseif tok == "--nprocs"
            i < length(argv) || error("--nprocs requires a value")
            n = tryparse(Int, argv[i+1]); (n === nothing || n < 1) && error("--nprocs must be a positive integer")
            nprocs = n; i += 2
        elseif tok == "--outdir"
            i < length(argv) || error("--outdir requires a value")
            outdir = String(argv[i+1]); i += 2
        elseif tok == "--data"
            sub in ("g6pd", "pgd", "hk1") || error("--data is only valid with an enzyme subcommand\n\n$CLI_USAGE")
            i < length(argv) || error("--data requires a value")
            data = String(argv[i+1]); i += 2
        elseif tok == "--variant"
            sub in ("g6pd", "pgd", "hk1") || error("--variant is only valid with an enzyme subcommand\n\n$CLI_USAGE")
            i < length(argv) || error("--variant requires a value")
            variant = String(argv[i+1]); i += 2
        elseif startswith(tok, "-")
            error("unknown flag: $tok\n\n$CLI_USAGE")
        elseif sub == "plot" && rundir === nothing
            rundir = tok; i += 1
        else
            error("unexpected argument: $tok")
        end
    end
    sub == "plot" && rundir === nothing && error("plot requires a <run_dir>\n\n$CLI_USAGE")
    return (sub, (smoke=smoke, nprocs=nprocs, outdir=outdir, rundir=rundir, data=data, variant=variant))
end

function cli_main(argv::AbstractVector{<:AbstractString})
    sub, o = parse_cli(argv)
    sub == "help" && (print(CLI_USAGE); return 0)
    if sub == "plot"
        plot_consensus_fit(o.rundir)
    else
        enz = Symbol(sub)                       # :g6pd / :pgd / :hk1
        variants = o.variant === nothing ? nothing : [Symbol(o.variant)]
        fit_consensus_equation(enz; smoke=o.smoke, nprocs=o.nprocs, outdir=o.outdir,
                               data_csv=o.data, variants=variants)
    end
    return 0
end

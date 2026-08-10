# FitRateEquation.jl

## 1. What this does

This package takes published enzyme-kinetics measurements (how fast an enzyme
converts its substrates into products, under different concentrations) and fits
them to a single, literature-consistent "consensus" rate equation for that
enzyme. It reports back the key numbers that describe the enzyme's behavior —
things like how tightly it binds each substrate, and how strongly related
molecules inhibit it — along with plots and a written report so you can see how
well the equation matches the data. It currently supports two enzymes fully
(G6PD and PGD, both central to red-blood-cell metabolism), with a third (HK1)
planned but not yet available.

You do not need to know any biochemistry or write any fitting code yourself —
the package ships with the enzyme data already built in, and one function call
runs the whole pipeline.

## 2. Install Julia

FitRateEquation.jl is written in the [Julia](https://julialang.org/)
programming language, so you need Julia installed first (version 1.11 or
later).

Download and install it from **<https://julialang.org/downloads/>**, following
the instructions for your operating system.

Once installed, open a terminal and start the Julia console (the "REPL",
short for Read-Eval-Print Loop) by typing:

```sh
julia
```

You'll see a `julia>` prompt. From there, typing the character `]` switches
you into **package mode**, where the prompt changes to `pkg>` — this is
Julia's built-in tool for installing and managing packages (like `pip` for
Python or `npm` for JavaScript). Press Backspace at an empty `pkg>` prompt to
return to the normal `julia>` prompt.

You can also run the `Pkg.add(...)` commands below directly at the `julia>`
prompt (without switching to `]` mode) — that's what the commands in this
README use, since it's copy-pasteable either way.

## 3. Install the package

Run the following two commands at the `julia>` prompt:

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
Pkg.add(url="https://github.com/JamesMbata/FitRateEquation.jl")
```

**Why two commands?** FitRateEquation.jl depends on another package,
EnzymeRates.jl, that does the low-level rate-equation algebra. EnzymeRates.jl
isn't published in Julia's central package registry (the way most packages
are), so Julia can't find it automatically — you have to tell it exactly
where to get it. A "git URL" here is just the web address of the code
repository; `Pkg.add(url="...")` tells Julia "download and install the
package that lives at this address" instead of looking it up by name in the
registry. Install EnzymeRates.jl first, then FitRateEquation.jl, in that
order.

## 4. Quickstart

Once both packages are installed, run a fast "smoke test" fit for G6PD:

```julia
using FitRateEquation
run_g6pd(smoke=true)
```

`smoke=true` uses a tiny fitting budget so this finishes quickly (a couple of
minutes) — it's meant to confirm everything is wired up correctly, not to
produce publication-quality numbers. Drop `smoke=true` (or set it to `false`)
for the full, slower fit.

When it finishes, you'll find a new folder under `./results/` (something like
`results/G6PD_2026-07-16_smoke/`) containing **seven files** — the fitted
constants, goodness-of-fit statistics, and a human-readable report. See
[Understanding the outputs](#6-understanding-the-outputs) below for what each
one means.

## 5. Run on your own data

If you have your own enzyme-kinetics measurements, you can fit them the same
way, as long as your CSV file uses the right columns. For G6PD, the columns
FitRateEquation.jl expects are:

| Column | What it is | Units |
|---|---|---|
| `[NADP] (uM)` | NADP⁺ concentration in that measurement | µM |
| `[G6P] (uM)` | Glucose-6-phosphate concentration | µM |
| `[NADPH] (uM)` | NADPH concentration | µM |
| `[PGLn] (uM)` | 6-phosphogluconolactone concentration | µM |
| `[ATP] (uM)` | ATP concentration (a regulator of G6PD) | µM |
| `Rate_V` | the measured reaction rate for that row | your assay's rate units |
| `Article` | which paper/experiment the row came from | text label |
| `Fig` | which figure within that paper | text label |
| `Apparent_Keq` | the apparent equilibrium constant for that row's conditions | dimensionless |

`Article` and `Fig` are used to group rows from the same source figure
together (so the fit knows which points were measured under the same
calibration) — they don't need to mean anything beyond "this row came from
this figure of this paper."

Your CSV must use the **exact** column names in the table above — this is a
fixed, canonical schema. A file with differently named columns is rejected with
an error that prints the header it expected, so you can rename your columns to
match. Column *remapping* (pointing the fitter at differently named columns) is
not supported. Extra annotation columns you don't need are allowed and ignored.

Point FitRateEquation.jl at your CSV like this:

```julia
using FitRateEquation
fit_consensus_equation(:g6pd; data_csv="/path/to/my_corpus.csv")
```

Add `outdir="my_results"` to choose where the output folder is written.

### Which rate law you're fitting

Each enzyme has one **deployed consensus rate law** — the equation wired into the
downstream `PentosePhosphatePathway.jl` simulation — plus a few alternative laws
you can fit for comparison. They all fit the same corpus and write the same seven
output files; they differ only in the structure of the rate equation (which
product-release steps are treated as steady-state, and which substrate/product
"dead-end" inhibitions are included).

There is one entry point for every fit: `fit_consensus_equation(:g6pd | :pgd | :hk1)`.
By default it fits that enzyme's deployed consensus law; to fit an alternative law
instead, name the variant with the `variants` keyword. `run_g6pd()` / `run_pgd()` /
`run_hk1()` are thin aliases for the deployed-law call, kept for discoverability.

| Enzyme | Rate law | How to run |
|---|---|---|
| G6PD | **Consensus law (deployed)** — random-order Bi-Bi, steady-state NADPH release | `run_g6pd()` or `fit_consensus_equation(:g6pd)` |
| G6PD | ATP-free — ATP dropped as a metabolite (ATP-bearing rows filtered out), for data measured without ATP | `fit_consensus_equation(:g6pd; variants=[:no_atp])` |
| G6PD | Ablation — drops the E·G6P·ATP dead-end | `fit_consensus_equation(:g6pd; variants=[:no_g6p_atp_deadend])` |
| G6PD | Ablation — drops the E·G6P·NADPH dead-end | `fit_consensus_equation(:g6pd; variants=[:no_g6p_nadph_deadend])` |
| G6PD | Ablation — drops both G6P dead-ends | `fit_consensus_equation(:g6pd; variants=[:no_g6p_both_deadends])` |
| G6PD | Pure rapid-equilibrium form (conflation diagnostic; not deployable) | `fit_consensus_equation(:g6pd; variants=[:RE_rate_eq])` |
| PGD | **Consensus law (deployed)** — Topham Bi-Ter, steady-state Ru5P release | `run_pgd()` or `fit_consensus_equation(:pgd)` |
| PGD | Fully rapid-equilibrium (evaluation only) | `fit_consensus_equation(:pgd; variants=[:full_re])` |

`variants` takes a Julia list `[...]`; each name is a Symbol (a leading colon, e.g.
`:no_atp`). Selecting a variant automatically applies its data filter (e.g. the
ATP-free variant drops the ATP-bearing rows) and names its output folder after it.

Every call takes the same `smoke`, `nprocs`, and `outdir` keywords as in the
[Quickstart](#4-quickstart), and any of these fits can be launched from the
[command line](#8-command-line) with `--variant NAME`. The `:full_re` PGD variant
has its own write-up in [§10](#10-an-alternative-pgd-rate-law-full_re); for the full
technical detail on every variant — what each dead-end means and why the ablations
exist — see [`AGENTS.md`](AGENTS.md).

## 6. Understanding the outputs

Every run writes seven files to its output folder:

- **`macro_constants.csv`** — the headline result: the fitted binding/inhibition
  constants (Km, Kd, Ki) for the enzyme, with each one labeled as
  data-derived, taken from the literature, or unconstrained by the data.
- **`fit_corpus.csv`** — the exact measurements this run was fitted to, saved
  alongside the result. If you pointed the fitter at your own data file, or had
  it skip some rows, this is what it actually used — and it's what the plots are
  drawn against.
- **`goodness_of_fit.csv`** — how well the fitted equation matches the data,
  including a cross-validation score computed by holding out one paper's data
  at a time.
- **`identifiable_functions.csv`** — which combinations of constants the data
  can actually pin down versus which are underdetermined.
- **`micro_parameters.jl`** — the fitted result rewritten as a ready-to-use
  Julia parameter block, for anyone wiring this rate law into a larger
  simulation.
- **`report.md`** — a plain-language written summary of the fit, readable
  without opening any of the other files.
- **`provenance.toml`** — a record of exactly how the run was produced
  (package version, random seed, fitting budget), so the result can be
  reproduced later.

## 7. Plots

To see the fitted curve plotted against the actual data points, install
CairoMakie (a plotting package) and call `plot_consensus_fit`:

```julia
using Pkg
Pkg.add("CairoMakie")

using CairoMakie, FitRateEquation
plot_consensus_fit("results/G6PD_2026-07-16_smoke")
```

This writes one PNG image per source figure into a `plots/` subfolder inside
your results folder. Plotting only works after `CairoMakie` has been loaded
(`using CairoMakie`) — FitRateEquation.jl deliberately doesn't require it for
fitting, since it's a large dependency you only need if you want pictures.
Plotting currently works for G6PD and PGD results only.

The plots are drawn against the `fit_corpus.csv` saved inside the results folder,
so they always show the data the fit actually used. Results folders produced by
versions before 0.2.0 don't contain that file and can't be plotted — re-run the
fit if you have an old one.

## 8. Command line

If you'd rather run a fit from your terminal than from inside Julia,
FitRateEquation.jl ships a small command-line interface. It does the same thing
as the `run_*` functions above — it just takes its instructions as command-line
arguments instead.

The simplest way to reach it (works from any project where the package is
installed) is:

```sh
julia --project -e 'using FitRateEquation; FitRateEquation.cli_main(ARGS)' -- g6pd --smoke
```

Everything after the `--` is passed to the tool. The subcommands are:

| Subcommand | What it does |
|---|---|
| `g6pd` | Fit G6PD |
| `pgd` | Fit PGD |
| `hk1` | Fit HK1 *(not yet available — see Troubleshooting)* |
| `plot <run_dir>` | Render the fitted law over the data for a finished run |
| `help` | Print usage |

And the flags:

| Flag | What it does |
|---|---|
| `--smoke` | Use the fast, low-budget fit (same as `smoke=true`) |
| `--nprocs N` | Use `N` worker processes for the fit |
| `--outdir DIR` | Write the run's seven output files into `DIR` |
| `--data CSV` | Fit your own CSV (canonical columns required) instead of the built-in corpus; valid with any enzyme subcommand |
| `--variant NAME` | Fit an alternative rate law instead of the deployed one (e.g. `no_atp`, `full_re`, `no_g6p_atp_deadend`) |

To fit the ATP-free G6PD variant, pass `g6pd --variant no_atp`; the older
`g6pd-noatp` subcommand no longer exists. For example, at smoke budget:

```sh
julia --project -e 'using FitRateEquation; FitRateEquation.cli_main(ARGS)' -- g6pd --variant no_atp --smoke
```

For example, to fit PGD at full budget and put the results in a folder called
`pgd_run`:

```sh
julia --project -e 'using FitRateEquation; FitRateEquation.cli_main(ARGS)' -- pgd --outdir pgd_run
```

The package also includes a ready-made launcher script at
`bin/fitrateequation` that wraps the line above, so if you have the package
checked out you can run it directly:

```sh
julia --project bin/fitrateequation g6pd --smoke
```

## 9. Troubleshooting

**`Pkg.add(url=...)` for FitRateEquation.jl fails, or complains it can't find
EnzymeRates:** make sure you ran the EnzymeRates.jl `Pkg.add(url=...)`
command *first*, as a separate step (see [Install the
package](#3-install-the-package)). EnzymeRates.jl is not in the general
registry, so Julia cannot resolve it automatically as a dependency.

**Julia complains about your Julia version, or something doesn't precompile:**
FitRateEquation.jl requires **Julia 1.11 or later**. Run `julia --version` in
your terminal to check, and download a newer release from
<https://julialang.org/downloads/> if needed.

**`run_hk1()` throws an error:** this is expected — HK1 support is not yet
available in this release (the underlying mechanism hasn't been ported over
yet). `run_g6pd` and `run_pgd` (and their general form,
`fit_consensus_equation(:g6pd | :pgd)`, including every variant) are fully
available.

## 10. An alternative PGD rate law (`:full_re`)

Alongside the default PGD fit, the package includes an alternative,
experimental rate equation for PGD called `:full_re` (a "fully
rapid-equilibrium" form). You run it by naming the variant:

```julia
using FitRateEquation
fit_consensus_equation(:pgd; variants=[:full_re])   # add smoke=true for a fast plumbing check
```

It writes the same seven output files as `run_pgd()`, into a folder labeled
`PGD_fullre_...`. On the built-in data it generalizes to unseen papers
noticeably better than the default PGD law and better captures how the product
NADPH slows the enzyme down. It is an **evaluation variant only** — it is not
wired into the downstream simulation model. The full write-up of what it does
and how it compares is in
[`docs/pgd_fullre_evaluation.md`](docs/pgd_fullre_evaluation.md).

## 11. Going deeper

This README covers everyday use. For the full model details — the exact rate
equation being fit, what's held fixed versus what's fit from data, the
per-enzyme reaction mechanisms, and the reasoning behind them — see
[`AGENTS.md`](AGENTS.md).

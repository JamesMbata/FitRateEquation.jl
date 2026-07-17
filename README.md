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
`results/G6PD_2026-07-16_smoke/`) containing **six files** — the fitted
constants, goodness-of-fit statistics, and a human-readable report. See
[Understanding the outputs](#6-understanding-the-outputs) below for what each
one means.

## 5. Run on your own data

If you have your own enzyme-kinetics measurements, you can fit them the same
way, as long as your CSV file has the right columns. For G6PD, the columns
FitRateEquation.jl expects (by default; you can point it at differently named
columns too — see below) are:

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

Point FitRateEquation.jl at your CSV like this:

```julia
using FitRateEquation
cfg = g6pd_config(data_csv="/path/to/my_corpus.csv")
run_all(cfg, outdir="my_results")
```

If your columns are named differently, build the config yourself and pass a
custom `metabolites` mapping (symbol → `(csv_column_name, :uM)`) plus
`rate_col`, `article_col`, `fig_col`, and `keq_col` — see
`src/configs/G6PD.jl` in the package source for the exact fields, or the
deeper reference in [`AGENTS.md`](AGENTS.md).

## 6. Understanding the outputs

Every run writes six files to its output folder:

- **`macro_constants.csv`** — the headline result: the fitted binding/inhibition
  constants (Km, Kd, Ki) for the enzyme, with each one labeled as
  data-derived, taken from the literature, or unconstrained by the data.
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

## 8. Troubleshooting

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
yet). `run_g6pd`, `run_pgd`, and `run_g6pd_noatp` are fully available.

## 9. Going deeper

This README covers everyday use. For the full model details — the exact rate
equation being fit, what's held fixed versus what's fit from data, the
per-enzyme reaction mechanisms, and the reasoning behind them — see
[`AGENTS.md`](AGENTS.md).

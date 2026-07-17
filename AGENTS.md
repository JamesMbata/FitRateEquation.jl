# FitRateEquation — machine / LLM reference

> **Audience:** agents and engineers who need the full model, gauge, and result
> detail. If you just want to *run* the pipeline, start with `README.md` (the human
> user guide). This file is the deep reference: the Cha-form law, the gauge, the
> per-enzyme mechanisms, the de-conflation findings, mode semantics, outputs, and
> layout.

Consensus macro-constant extractor for the G6PD and PGD rate equations, packaged
as a standalone, installable Julia package.

FitRateEquation does **not** select a mechanism (unlike the now-removed
`fitting/mechanism_id/` and `fitting/apparent_id/` ranking pipelines in the sibling
`PPP_Experiments` repo, which *enumerated and ranked* hundreds of candidate
skeletons). It **imposes** the literature-consensus mechanism for each enzyme and
extracts the identifiable forward **macro constants** (Km / Kd / Ki) from it. This
sidesteps — rather than resolves — the family ties that a ranking pipeline keeps
producing on this forward-saturation-dominated corpus (the binding *order* is not
identifiable from the data; see `docs/G6PD_session_context.md` and
`docs/PGD_session_context.md` in the sibling `PPP_Experiments` repo).

> **EnzymeRates dependency:** this package runs on **upstream EnzymeRates**
> (`DenisTitovLab/EnzymeRates.jl`, pinned in `Project.toml`'s `[sources]` for
> dev/CI; a downstream consumer installs it separately by URL — see `README.md`
> §3). It is self-contained — it vendors the data loader / mechanism builder /
> gauge helpers under `src/core/` (formerly borrowed from `fitting/mechanism_id/`).
> Mechanism topologies are written with terse opaque form names in
> `src/enzymes/<enzyme>.jl`; the builder auto-derives the decomposed call notation
> `E(NADP,G6P)` the upstream DSL requires. Rate-constant references use upstream's
> composition-semantic names (`K_NADP_E`, `koff_NADPH_E`, `k_EG6PNADP_to_ENADPHPGLn`,
> …). **HK1 is guarded**: its bespoke low-level `EnzymeMechanism` construction +
> allosteric DSL need a separate port from the upstream migration that G6PD/PGD
> already went through. `src/enzymes/hk1.jl` is included inside a `try/catch` in
> `src/FitRateEquation.jl`; if it fails to load on the installed EnzymeRates build,
> `HK1_AVAILABLE` is set to `false`, G6PD/PGD load normally, HK1 tests auto-skip,
> and `run_hk1()` raises a clear error pointing back here rather than silently
> doing nothing. There is currently no bundled EnzymeRates release on which HK1
> wiring succeeds — treat `run_hk1` as **not yet available** until this note is
> updated.

## What it fits — the Cha-form macro-direct law

The package fits each enzyme's consensus law as a **Cha-form partial-equilibrium
(macro/micro) hybrid**, optimizing the **named forward macro constants directly**
through a closed-form rate law (`cha_rate_G6PD` / `cha_rate_PGD` in
`src/cha_laws.jl`). This is the Phase-A/B reparameterization; it **retired** the
older coefficient-space path (`coeff_fit`/`coeff_params`/`coeff_eval`/
`macro_bridge`/`coeff_identifiability`) that predated this package.

The Cha law is the **promoted** law: one product-release step besides catalysis is
forced to steady state, the rest stay rapid-equilibrium. The pure rapid-equilibrium
variant is simply its `release_rate → ∞` limit, so the package runs **only the
deploy (promoted) variant per enzyme** — it is not a multi-variant panel:

- **G6PD → `SS_NADPH_release_rate_eq`** — catalysis (step 5) *and* NADPH-release
  (step 7) both SS, with PGLn-release (step 6) rapid-equilibrium between them. `E_C`
  and `E_H` form one SS **super-node** `S = {E_C ⇌ E_H}`, internally RE-distributed
  by `r = PGLn/Kd_6PGLn`. The Cha C-factor is `C = 1 + kf/koffQ`.
- **PGD → `:cha_base`** — Topham asymmetric-random Bi-Ter with the one
  gauge-mandated SS release: **Ru5P-release** (Topham 1986 rate-limiting). CO₂ and
  NADPH release stay rapid-equilibrium (PGD's NADPH release is fast, >800 s⁻¹, so no
  silent fiber is promoted). Super-node `S = {E_C ⇌ E_1}`, RE-linked by
  `s = CO2/Kd_CO2`.

The numeric Cha laws are exactness-anchored against `EnzymeRates.rate_equation` for
the real micro mechanism at **120/120 @ rtol 1e-10** — this exactness anchor is
carried forward as the package's `test/test_cha_laws.jl` / `test_cha_pgd_laws.jl`
suite (the original derivation notebooks, `derivations/cha_derive_{g6pd,pgd}.py`,
live in the sibling `PPP_Experiments` repo and are not part of this package).

### What is fit vs. fixed (the gauge)

Under the per-(Article, Fig) **mean-centered log-ratio loss** the overall Vmax is
gauged out, so the kinetic-scale parameters are **fixed, not fit**:

- `kf = 1.0`, `Et = 1.0` (kcat = 1 / unit-enzyme gauge);
- the promoted SS-release rate (`release_rate`: `koffQ` for G6PD, the Ru5P-release
  `koff` for PGD) is a **swept fiber** held at a healthy default. The fit *loss* is
  insensitive to it because the only fiber-**invariant** forward observable is the
  specificity constant `kcat/Km = kf/(alpha·Kd)`; the apparent `kcat` and `Km`
  *individually* slide along the fiber (`kcat = kf·r/(kf+r)`, `Km = alpha·Kd/C`,
  `C = 1+kf/r`), holding their ratio fixed. So `release_rate` is unidentifiable from the
  forward corpus, **but any Km/kcat readoff must use the same `release_rate` the law is
  deployed at** (`CHA_DEPLOY_RELEASE_RATE = 1e3`) or it reports a constant the deployed
  model does not exhibit;
- the reverse catalysis `kr` is **Haldane-determined** from `keq`.

So the **free fit coordinates** (`cha_coords`, in `src/cha_fit.jl`) are exactly the
data-identifiable forward shape constants:

| Enzyme | `cha_coords` (fit directly) |
|---|---|
| G6PD | `Kd_NADP`, `Kd_G6P`, `Kd_6PGLn`, `alpha`, `Ki_NADPH`, `Ki_ATP`, `Ki_ATP_EG`, `Km_NADPH_rev` |
| PGD  | `Kd_NADP`, `Kd_PGA`, `alpha`, `Kd_CO2`, `Ki_NADPH`, `Ki_ATP`, `Ki_ATP_EN`, `Km_NADPH_rev` |

The **apparent Michaelis constants** are NOT coords — they are **derived** readoffs
`Km = alpha·Kd / C`, `C = 1+kf/release_rate` (`cha_apparent_km`): `Km_G6P` (G6PD);
`Km_PGA`, `Km_NADP` (PGD). They are read at `CHA_DEPLOY_RELEASE_RATE` (so they match the
deployed law) and surfaced with class `:derived`, each paired with its fiber-invariant
specificity `kcatKm_* = kf/(alpha·Kd)` (`cha_specificity`). The Mode-2/3 `Km_PGA` *anchor*
also reads at `CHA_DEPLOY_RELEASE_RATE` (`_cha_anchor_penalty`), so the penalty pulls the
**deployed** apparent `Km_PGA` onto the literature target — a PGD consensus *deploy* at
`release_rate = 1e3` lands `Km_PGA` on the 38–80 µM band, not `2×` it.

Semantics: **`Km`** = ratio-type Michaelis constant, **`Kd`** = binary free-enzyme
dissociation constant, **`Ki`** = inhibitor / dead-end constant, **`alpha`** =
dimensionless ternary-interaction factor.

### G6PD — random-order sequential Bi-Bi (Wang 2002 consensus)

`NADP + G6P ⇌ NADPH + PGLn` (6-phosphogluconolactone). King-Altman skeleton
(`src/enzymes/g6pd.jl`):

```
E + NADP   ⇌ E_N        ┐ random substrate binding
E + G6P    ⇌ E_G        │
E_N + G6P  ⇌ E_NB       │
E_G + NADP ⇌ E_NB       ┘
E_NB → E_C   (SS — catalysis, gauge anchor)
E_C  ⇌ E_H + PGLn       (RE — PGLn released first)
E_H  ⇌ E + NADPH        (SS — promoted NADPH release; super-node {E_C ⇌ E_H})
```

- **Dead-ends (rapid-equilibrium):** ATP competitive vs NADP (binds open-
  dinucleotide forms `[E, E_G]`); NADPH dead-end on `[E_G]` only — `E·G6P·NADPH`,
  noncompetitive vs G6P (Wang 2002 Fig 4 / kinetics-literature §15.1). The `E_N` and
  the free-E `E·NADPH` dead-ends were dropped: dropping the free-E `E·NADPH`
  (the reverse of productive NADPH release) leaves the forward `Ki_NADPH` as a clean
  single cross-term constant. See `docs/G6PD_session_context.md` (sibling
  `PPP_Experiments` repo). **Regulator:** ATP.
- **Forward `Ki_NADPH`** is read as the **cross-term** `E·G6P·NADPH` dead-end
  (`g[G6P]/g[G6P·NADPH]`), distinct from the bare-[NADPH] reverse-release Km
  (`Km_NADPH_rev`). Anchoring `Km_NADPH_rev` (3.9 µM, all modes) de-conflates the two.

### PGD — random-binding / ordered-release Bi-Ter (Topham 1986)

`NADP + PGA (6PG) ⇌ CO2 + Ru5P + NADPH`. King-Altman skeleton (`src/enzymes/pgd.jl`):

```
E + NADP   ⇌ E_N        ┐ random substrate binding
E + PGA    ⇌ E_G        │
E_N + PGA  ⇌ E_NB       │
E_G + NADP ⇌ E_NB       ┘
E_NB → E_C   (SS — catalysis, gauge anchor)
E_C ⇌ E_1 + CO2         (RE — CO2 out first; super-node {E_C ⇌ E_1})
E_1 → E_2 + Ru5P        (SS — promoted Ru5P release, Topham rate-limiting)
E_2 ⇌ E + NADPH         (RE — NADPH released last)
```

- **Dead-ends:** ATP competitive vs PGA, plus a distinct ATP dead-end on `E·NADP`
  (`Ki_ATP`, `Ki_ATP_EN`); the E·PGA NADPH dead-end carries the decoupled
  forward `Ki_NADPH` cross term (`g[PGA]/g[PGA·NADPH]`), distinct from the
  bare-[NADPH] reverse Km (`Km_NADPH_rev`). **Regulator:** ATP.
- **`KdRu`** (Ru5P-release equilibrium `koff/kon`) is a Dalziel **nuisance fiber** —
  unidentifiable on forward + product-inhibition data — held at `CHA_KDRU_DEFAULT`,
  never pinned. It is distinct from `Km_NADPH_rev` (the NADPH-release equilibrium);
  both enter the Haldane formula independently.

### The recurring finding (why the de-conflation matters)

Both enzymes share one structural problem: a **single `Ki_NADPH` symbol** would
otherwise serve both the forward product-inhibition (E·NADPH) term and the reverse
productive-release Km. The Cha law separates them into the **cross-term dead-end**
(`Ki_NADPH`) and the **bare-[NADPH] reverse Km** (`Km_NADPH_rev`). Outcomes at full
budget:

- **G6PD:** the de-conflation **works** — Mode-1 data-identifies the cross-term
  `Ki_NADPH ≈ 37 µM (±7 µM)`, on the forward scale with a finite CI (not the reverse
  ~0.2–2.4 µM). But it reads **above** the 9–24 µM literature band; Mode 2 pins it to
  15 µM.
- **PGD:** the forward cross-term de-conflation was **refuted** on this corpus
  (Mode 1 reports `Ki_NADPH` diagnostic/unconstrained), so it is **pin-only** 17 µM
  (Cottreau) in Modes 2/3. `Km_PGA` is **not data-identified** — `alpha` rails to its
  bound (the corpus lacks sub-Km [6PG] coverage), so the apparent Km is anchored, not
  recovered. The high/unstable PGD leave-one-article-out CV traces to a
  **forward-CO₂ product-inhibition coverage gap** (the leave-Weisz1985-out fold);
  trimmed CV ≈ G6PD's.

> **Topology is FROZEN.** Do not add SS release steps to chase the forward
> `Ki_NADPH` — proven inert (the shared symbol + corpus geometry is the blocker, not
> the SS-step count). The only de-conflation routes are a decoupled forward symbol
> (the cross-term) or bench data. See `test/test_rec4_topology_freeze.jl`.

## Modes

Each enzyme is fit in **per-enzyme modes** (`modes_for`: G6PD = `mode1, mode2`;
PGD = `mode1, mode2, mode3`), written to separate, never-cross-ranked leaderboards.
Pins are resolved structurally by `ChaFit.resolve_cha_pins`; the PGD `Km_PGA` anchor
goes through the apparent-Km penalty (`cha_anchors` in `src/run.jl`), never a hard
pin (it is a derived constant, not a coord).

- **Mode 1** — free fit. The data-only, family-agnostic headline. (G6PD still
  anchors `Km_NADPH_rev` in *all* modes to protect the cross-term identifiability.)
- **Mode 2** — literature pins applied (`Ki_ATP`, `Ki_NADPH`; PGD soft-anchors
  `Km_PGA` to the 38–80 µM band midpoint, 59 µM, weight 1.0 — provisional). This is
  the ODE drop-in.
- **Mode 3** (PGD only) — the explicit `Km_PGA` **hard override** (38 µM) realized as
  a high-weight (100) apparent-Km anchor.

Pinned constants are tagged `:literature_pinned` **structurally** (a coordinate that
has a pin), never from profile curvature. `report.md` surfaces a **mode-agreement
check** on the forward shape constants: if a "pinned" constant was not actually flat
(it coupled into the forward fit) the point estimates diverge by `> tol_dex` and a
distortion warning fires — never hidden. Pinned/derived rows never enter the check,
so an intended data-vs-physiology contrast is not mislabeled a flatness warning.

**Convention (mode 1 is deploy for new laws):** for any law adopting it, **mode 1 is
the deploy source of truth; mode 2/3 are diagnostics/validation only**, not deploy
candidates. This is a convention for *future* re-derivations — the currently
deployed full-G6PD law (mode2, literature-pinned `Ki_ATP`/`Ki_NADPH`, "the ODE
drop-in" above) and PGD (mode2/3) are unchanged and keep their existing deploy
status until they are themselves re-derived. The first law to adopt the convention
is the ATP-free G6PD variant (`:no_atp`, `src/enzymes/g6pd.jl`): with the ATP
dead-ends dropped entirely, its mode1 free fit is fully data-identified, so mode1 is
its deploy candidate and mode2 (`Ki_NADPH` literature pin only) is diagnostics-only.

## Running

This is an installed Julia package, not an in-repo script pipeline: there is no
required working directory and no `--project=...` flag to remember; once `using
FitRateEquation` has loaded, every runner works from any current directory. See
`README.md` for the novice-level walkthrough; this section is the complete
reference.

**Exported runners** (each writes the six artifacts below to `outdir` and returns
the `run_all` results):

```julia
using FitRateEquation

run_g6pd(; outdir=nothing, smoke=false, nprocs=nothing)
run_pgd(;  outdir=nothing, smoke=false, nprocs=nothing)
run_hk1(;  outdir=nothing, smoke=false, nprocs=nothing)         # errors: HK1 guarded, see above
run_g6pd_noatp(; outdir=nothing, smoke=false, nprocs=nothing, data_csv=nothing)
```

- `smoke=true` → `n_restarts=2, maxiter=150` (vs full `n_restarts=48, maxiter=1000`)
  — fast plumbing check, not a fit to trust.
- `nprocs` → local worker-count override; defaults to the `FRE_NPROCS` env var if
  set, else `max(1, min(3, Sys.CPU_THREADS - 1))` (the guarded 3-worker default in
  `src/worker_setup.jl`, capped low to avoid the RAM oversubscription that
  OOM-crashes full-budget fits — each worker holds its own copy of the model/data,
  so memory, not CPU, is the binding constraint). Inside a SLURM allocation,
  `setup_workers` detects it (`SLURM_JOB_CPUS_PER_NODE`) and uses `SlurmManager`
  instead, sized to the allocation — `nprocs` is ignored in that case.
- Default `outdir`: `./results/<ENZYME>_<YYYY-MM-DD>[_smoke]/` (relative to the
  caller's `pwd()`), e.g. `results/G6PD_2026-07-16_smoke/`. `run_g6pd_noatp`'s
  default is labeled `G6PD_noatp_...` so it never collides with plain `run_g6pd`.
- **Own data:** build a config pointing at your CSV (same schema as the bundled
  corpus), then run it directly:
  ```julia
  cfg = g6pd_config(; data_csv="/path/to/my_corpus.csv")   # also pgd_config, hk1_config
  run_all(cfg; outdir="my_results", n_restarts=48, maxiter=1000)
  ```
- **Fit variants per enzyme:** G6PD `:SS_NADPH_release_rate_eq` (deploy) + `:no_atp`
  (ATP-free, via `run_g6pd_noatp`); PGD `:cha_base`; HK1 `:H1`, `:H4` (not runnable
  while guarded). `run_all(cfg; variants=[…], row_filter=…)` exposes a custom
  variant set / row filter directly for advanced use.
- Configs (data CSV, `deploy_keq`, metabolite columns/units) live in
  `src/configs/G6PD.jl`, `src/configs/PGD.jl`, `src/configs/HK1.jl`. The bundled
  corpora resolve via `pkgdir(FitRateEquation)` so they load correctly regardless
  of installation location; `data_csv` overrides them. The fit itself uses each
  figure's own apparent Keq, read per-row from the corpus `keq_col` column; the
  config's `deploy_keq` is the single readout/deploy value — G6PD `13.655` (apparent
  Keq at cellular pH 7.2 / 37 °C), PGD `0.17` (pH-flat Bi-Ter CO₂ aq, 37 °C).

**CLI** — a thin shim (`bin/fitrateequation`, `using FitRateEquation;
cli_main(ARGS)`) dispatching the same runners in-process, no subprocess:

```sh
bin/fitrateequation <g6pd|pgd|g6pd-noatp|hk1> [--smoke] [--nprocs N] [--outdir DIR]
bin/fitrateequation g6pd-noatp [--data CSV]
bin/fitrateequation plot <run_dir>
bin/fitrateequation help
```

### Plotting the fit over the data

`plot_consensus_fit(run_dir)` (loaded from the `CairoMakie` package extension,
`ext/FitRateEquationMakieExt.jl` — `using CairoMakie` first) renders the fitted law
against the digitized corpus — one panel per source figure (`Article|Fig`), data as
scatter with the fitted prediction as a line. It is a **standalone post-hoc step**:
it does not run as part of a fit, and it is not auto-generated.

```julia
using CairoMakie, FitRateEquation
plot_consensus_fit("results/G6PD_2026-07-16_smoke")
```

- **`run_dir`** is a results dir (absolute or relative). The **enzyme is
  auto-detected** from `macro_constants.csv`'s content (not the path, so a
  relocated/renamed run dir still resolves); the fitted constants are read back
  from that file. Output is one `<run_dir>/plots/<variant>_<mode>_fit.png` per row.
- Each prediction is the **deployed** Cha law: the macro tuple is assembled at
  `CHA_DEPLOY_RELEASE_RATE` and evaluated through `cha_rate_*`, so the curve matches
  the law written to `micro_parameters.jl`. Each figure's curve uses **that
  figure's own apparent Keq** (the reverse arm's Haldane `kr`), read per-row from
  the corpus. A per-figure Vmax gauge normalizes each panel, so the plots show
  **shape agreement**, not absolute scale.
- **G6PD and PGD only.** HK1 is not yet supported — its corpus lacks the
  `X_axis_label` column the per-figure renderer needs, so pointing the plotter at
  an HK1 run dir raises a clear error rather than rendering.
- The non-Makie logic (enzyme detection, config lookup, coordinate readback, plot
  dataframe assembly) lives in `src/plot_support.jl` and is reachable without
  loading CairoMakie; only the render loop itself lives in the extension. The first
  invocation after `using CairoMakie` precompiles it (~a few minutes; cached after).
  Because it only reads `macro_constants.csv`, it never re-runs or perturbs a fit.

## Outputs (per run dir)

| File | Contents |
|---|---|
| `macro_constants.csv` | Resolved `cha_coords` per mode (class `data_identified` / `literature_pinned` / `unconstrained`) **plus** the derived apparent Michaelis constants (`Km_G6P` / `Km_PGA` / `Km_NADP`, class `derived`). |
| `goodness_of_fit.csv` | In-sample loss + leave-one-article-out CV (`mean ± se`) per mode. |
| `identifiable_functions.csv` | Per-mode eigen-spectrum of the identifiable macro-coordinate subspace. |
| `micro_parameters.jl` | ODE-ready deploy block — the **closed-form** `ChaDeploy.cha_deploy_micro` inverse of the fitted macro coords (no LM / root-finding; reproduces the law to machine precision). |
| `report.md` | Per-mode macro tables, cross-mode agreement, the G6PD `koffQ` hybrid block, the PGD `Km_PGA` gap warning, and the de-conflation caveat. |
| `provenance.toml` | Package version / EnzymeRates SHA, budget, row count, seed, smoke flag — reproducibility record. `_git_sha` reads `"unknown"` in an installed depot (no `.git` to inspect) — this is expected and already `try/catch`-guarded; the package version is the durable provenance handle for an installed copy. |

`koffQ` is reported as a **hybrid** (`src/cha_koffq_report.jl`): deployed at a
healthy swept default (flux-neutral), and *additionally* a reverse-weighted
diagnostic value + wide CI + deploy↔data gap — report-only, never fed to deploy.

A `plots/` subdirectory is **not** written by the fit run itself — it is produced
on demand by `plot_consensus_fit` (see [Plotting the fit over the
data](#plotting-the-fit-over-the-data) above).

## Determinism: the byte-identity smoke fixtures

`test/fixtures/{g6pd,pgd}_smoke_macro_constants.csv` are committed reference
`macro_constants.csv` outputs from a fixed-seed smoke run
(`test/test_byte_identity.jl`). Every test run re-runs `run_g6pd(smoke=true,
nprocs=1)` / `run_pgd(smoke=true, nprocs=1)` into a temp dir and diffs the result
against the fixture: the `variant`/`mode`/`name`/`value`/`class` columns must match
**exactly** (this is the strict gate on the fit output itself — the seeded CMA-ES
budget is `maxiter`-bound, not wall-clock, so results are reproducible regardless of
worker count or machine); the `ci` column (a finite-difference-Hessian confidence
interval from `cha_classify`) is compared with `rtol=1e-6` instead, since FD-Hessian
evaluation picks up last-bit floating-point noise from accumulated execution state
when run mid-suite — a fresh-process run reproduces the fixture exactly on every
column including `ci`. If you change the fit path (loss, bounds, CMA-ES config,
gauge) and these fixtures fail, that is the regression gate doing its job: either
the change is a bug, or the fixtures need a deliberate, reviewed regeneration — do
not silently update them to make the test pass.

## Layout

```
src/FitRateEquation.jl              module entry point (includes below, in load order)
src/run.jl                          Cha-path fit/CV + reduction + write_outputs; exported runners
src/worker_setup.jl                 shared local/Savio addprocs setup (guarded 3-worker default)
src/cli.jl                          in-process CLI dispatcher (cli_main; see bin/fitrateequation)

src/cha_laws.jl       closed-form Cha rate laws (cha_rate_G6PD / cha_rate_PGD)
src/cha_fit.jl        cha_coords, gauge/fiber/Haldane assembly, centered-logratio loss,
                       CMA-ES fit, apparent-Km readoff, resolve_cha_pins
src/cha_classify.jl   macro-coordinate identifiability (FD Hessian) + class assignment
src/cha_invert.jl     closed-form macro↔micro readoff + koffQ fiber map (exactness anchor)
src/cha_deploy.jl     closed-form macro→micro deploy inverse (micro_parameters.jl)
src/cha_koffq_report.jl   G6PD koffQ hybrid report (swept deploy + reverse diagnostic)
src/promotable.jl     declarative registry of promotable slow steps (G6PD nadph_release)

src/enzymes/g6pd.jl, src/enzymes/pgd.jl, src/enzymes/hk1.jl
                       per-enzyme topology + lit values + alias/ki-ratio maps (hk1.jl
                       guarded — see the EnzymeRates dependency note above)
src/configs/G6PD.jl, src/configs/PGD.jl, src/configs/HK1.jl
                       data CSV path (pkgdir-resolved), deploy_keq, metabolite columns
src/enzyme_wiring.jl   EnzymeWiring registry struct + accessors + modes_for
src/mechanisms.jl      generic dead-end / King-Altman builder helper

src/core/data.jl, src/core/mechbuild.jl, src/core/gauge.jl
                       vendored data loader / mechanism builder / gauge (formerly borrowed
                       from fitting/mechanism_id/ in the sibling repo)
src/macro_collect.jl   macroscopic-constant readback from the Haldane-reduced symbolic rate
                       law; retained as the exactness-anchor test spine (test_pgd_macro_collect.jl)
src/cv.jl              LIVE: leave-one-article-out fold helpers (_article, _article_folds,
                       _subset) used by run.jl's Cha-space LOOCV (_cha_loocv)

src/plot_support.jl    non-Makie plot helpers (enzyme detection, config lookup, coordinate
                       readback, plot dataframe assembly) — reachable without CairoMakie
ext/FitRateEquationMakieExt.jl
                       CairoMakie package extension: the actual render loop
                       (plot_consensus_fit), loads only when CairoMakie is loaded

data/                  bundled corpora: G6PD_all_EnzymeData.csv, PGD_EnzymeData_with_CO2.csv,
                       Choe_HK1_kinetic_data.csv
bin/fitrateequation    CLI entry point (using FitRateEquation; cli_main(ARGS))
test/                  regression suite (runtests.jl); test/fixtures/ holds the byte-identity
                       smoke references (see above)
```

`identifiability.jl` and `pins.jl` — the coefficient-space identifiability
classifier and Mode-2 coord-pin resolver from the pre-package pipeline — were the
genuinely dead pair (superseded by `cha_classify.jl` and `ChaFit.resolve_cha_pins`
respectively) and are **not carried into this package**. `cv.jl`, by contrast, is
correctly labeled above as **live**: its fold helpers are the ones `run.jl` calls
for every CV computation, Cha-path included — do not treat it as legacy.

## See also

- `docs/G6PD_session_context.md`, `docs/PGD_session_context.md` (sibling
  `PPP_Experiments` repo) — full run history and the de-conflation diagnoses.
- `~/projects/PentosePhosphatePathway.jl/docs/{g6pd,pgd}_kinetics_literature.md` —
  literature anchors.
- `fitting/mechanism_id/`, `fitting/apparent_id/` (sibling `PPP_Experiments` repo)
  — the ranking pipelines this package complements rather than replaces.

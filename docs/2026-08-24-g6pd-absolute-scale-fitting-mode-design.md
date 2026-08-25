# G6PD Absolute-Scale Fitting Mode — Design Spec

**Date:** 2026-08-24
**Branch:** `absolute-scale-mode`
**Status:** Design — awaiting review before implementation planning
**Scope:** `FitRateEquation.jl` (G6PD Cha fitting path)

---

## 1. Summary

Add an orthogonal `scale` axis to the G6PD Cha fitting path:

- `scale=:relative` (default) — the current per-(Article,Fig) mean-centered
  log-ratio loss. Bit-for-bit unchanged.
- `scale=:absolute` — an *uncentered* log-ratio loss that scores absolute rate
  magnitudes. The enzyme turnover `kcat` becomes a free, dimensionally-honest
  fit coordinate, and each reaction's measured enzyme concentration `Et` (a new
  required corpus column) enters the prediction as a linear prefactor.

The `scale` axis is orthogonal to the existing pin-set axis (`mode1`/`mode2`)
and mechanism axis (`variant`) and composes freely with both.

## 2. Motivation

The centered loss lets every measurement group float its own overall scale
(per-group mean subtraction). That freedom is exactly what allows structurally
different mechanisms to collapse onto the same fit — under the centered loss the
only fiber-invariant observable the forward corpus constrains is the specificity
constant `kcat/Km = kf/(α·Kd)` (see `cha_fit.jl` header at
`_default_release_rate`). Uncentering removes the per-group scale freedom, so
**between-condition relative magnitudes become discriminating signal**.

**Primary objective (user):** improve the ability to discriminate between the
existing G6PD dead-end variants (dead-end forms; also tightens the forward shape
constants). Recovering an absolute Vmax/kcat is a *byproduct*, not the
deliverable.

**Known caveat carried into the design:** the prior consensus_macro RE-vs-SS work
found that *forward-only* data goes machine-zero-blind to the release-mechanism
(RE vs SS) axis. Dead-end **forms** are a different axis, so absolute forward
data plausibly adds signal there — but this mode is built to *test* that
empirically, not to assume it. The RE-vs-SS release axis is explicitly out of
scope.

## 3. Goals / Non-Goals

### Goals
- `scale=:relative | :absolute` keyword threaded end-to-end (API + CLI).
- Uncentered absolute loss, sharing the per-row arithmetic core with the
  centered loss (no duplication, no drift).
- `kcat` as a first-class, dimensionally-honest fit coordinate (s⁻¹), single
  source of truth via a scale-parameterized `cha_coords`.
- Fiber-free (`C = 1`) absolute mode — no release-rate back-solve.
- Per-row enzyme concentration `Et` read from a required `[G6PD] (nM)` column.
- Explicit `pins::Dict{Symbol,Float64}` value-override for the reverse/degenerate
  constants, driving the escalation ladder by hand, documented in the README.
- Dual model-selection output: leave-one-group-out CV **and** in-sample loss +
  identifiability.
- `scale=:relative` reproduces the current pipeline bit-for-bit (regression lock).

### Non-Goals
- No new mechanisms or binding-order variants (the six existing G6PD variants
  only).
- No corpus article-filtering: absolute mode runs on a supplied single-scale,
  forward-only CSV (the user's forthcoming dataset), not a slice of the mixed
  corpus.
- No attempt to identify koffQ / the release-rate fiber from forward data (known
  unidentifiable).
- PGD and HK1 remain relative-only. The machinery is left generalizable, but
  `scale=:absolute` is guarded to `:G6PD` for now (explicit "not yet wired"
  error).
- Absolute Vmax is emitted (a byproduct) but is not the acceptance criterion.

## 4. Background: the scale/fiber relationship

The Cha law carries a promoted SS-release fiber (koffQ for G6PD). Along that
fiber, with catalysis forward rate `kf` and release rate `r`:

```
kcat = kf·r/(kf+r)          Km = α·Kd·r/(kf+r) = α·Kd/C,   C = 1 + kf/r
kcat/Km = kf/(α·Kd)         (fiber-invariant)
```

Two consequences drive the design:

1. **`kcat` is capped at `r`.** A finite `r` with the user's `kcat ≈ 178 s⁻¹`
   (from SA 181 U · ~59 kDa / 60 000) and the deploy `r = 1e3` would give
   `kf = kcat·r/(r−kcat) ≈ 216`, hence `C ≈ 1.22` — a ~22% fiber factor silently
   distorting the apparent `Km`. Unacceptable.
2. **The fiber is forward-unidentifiable** and the RE-vs-SS axis is out of scope,
   so we *fix* it — at `C = 1` exactly.

**Decision: absolute mode is fiber-free (`C = 1`).** This mirrors how
`cha_apparent_km` already treats fully-RE PGD and HK1 (`cha_fit.jl:431-432`).
Then `kcat ≡ kf` directly, `Km = α·Kd` exactly (matching the deployed law's
catalysis-limited `r=1e3` readoff to within ~0.1%), and there is no back-solve
and no release-rate units reconciliation. The Haldane `kr` for G6PD
(`cha_haldane_kr`, `cha_fit.jl:133-137`) consumes `kf` and `Km_NADPH_rev` and is
**independent of the release-rate magnitude**, so thermodynamic consistency is
preserved unchanged.

## 5. API surface

New keyword `scale::Symbol = :relative` on:
- `fit_consensus_equation(enzyme; …, scale=:relative)` — `run.jl`.
- The `run_g6pd(...)` alias forwards it.
- CLI: `--scale relative|absolute` in `cli.jl`.

Threading: `fit_consensus_equation` → `_fit_consensus`/`run_variants`
(`run.jl`) → `cha_fit_candidate` → `_cha_loss_with_pins` → the loss and
`cha_coords`/`cha_coord_bounds`.

Guard: `scale=:absolute` with `enzyme != :G6PD` raises an explicit
"absolute scale not yet wired for <enzyme>" error (no silent fallback).

Example:

```julia
fit_consensus_equation(:g6pd;
    variant = :no_atp,
    mode    = :mode2,
    scale   = :absolute,                     # new orthogonal axis
    data_csv = "my_single_scale_forward.csv",
    pins    = Dict(:Kd_6PGLn => log10(2.1e-4),
                   :Km_NADPH_rev => log10(3.9e-6)))
```

## 6. Scale-parameterized coordinates (single source of truth)

`cha_coords` and `cha_coord_bounds` gain `scale::Symbol = :relative`. For
`:G6PD, scale=:absolute` they append **`:kcat`** to the binding-constant coords:

```julia
cha_coords(:G6PD, variant; scale=:absolute) ==
    [<existing binding coords for variant>…, :kcat]
```

Because the whole fit pipeline keys off `cha_coords` (bounds, the CMA-ES vector,
the pin-overwrite, `_assert_pin_is_coord`), `:kcat` then flows everywhere
automatically. This is the honest single-source-of-truth choice: the coord set
*is* the parameter-vector definition, and `:kcat` genuinely is a free parameter
of the absolute-mode fit.

**Bound:** `:kcat` is the one rate-dimensioned coord (as `:alpha` and
`:split_ratio` are already non-`M` special cases in `cha_coord_bounds`). Range
`[10, 1000] s⁻¹` (log10 `[1, 3]`). The user's SA-derived `~178 s⁻¹` is a
**report-time validation check** against the literature `150–250 s⁻¹` band, not
a hard constraint — a fit that lands outside the band must be visible, not
clamped.

**Downstream consumers** (`classify_cha`, `cha_invert`/`cha_deploy`) get a small
`:kcat`-aware branch: in absolute mode `kcat` is the fitted, identifiable scale
and feeds `analytic_kcat`/the micro map directly (no gauge-1 assumption).

## 7. The loss (shared core + two aggregators)

Refactor `cha_centered_logratio_loss` into a shared core plus two thin
aggregators:

- **`_cha_row_logratios!(…)`** — the shared per-row core: the group loop,
  per-figure `keq` resolution, `cha_macro_tuple` assembly, `cha_rate_*`
  evaluation, and the sign/finite `_SIGN_PENALTY` bookkeeping. Fills the
  `logratio` vector, accumulates the penalty, and returns the per-group index
  sets. **This is the arithmetic that must never drift between modes.**
- **`cha_centered_logratio_loss`** — thin wrapper: aggregates per-group variance
  (mean-centered), preserving the current byte-identical fold order (the
  existing float-associativity-sensitive structure is retained exactly).
- **`cha_absolute_logratio_loss`** — thin wrapper: aggregates `sum(logratio²)`
  uncentered.

`_cha_loss_with_pins` dispatches on `scale`. Per-figure `keq` is still resolved
in the core (it varies across figures); in the uncentered path `keq` no longer
cancels, but for forward data far from equilibrium the reverse term is
negligible.

### 7.1 Per-row `Et`

`Et` is a **linear prefactor** on the Cha rate (`v = Et · kcat · f(concs,
params)`; enzyme conservation gives `v ∝ Et_total`). So the core evaluates the
unit-`Et` rate once and scales by the row's `Et`:

```
absolute mode:  logratio_i = log(Et_i · v_unit_i) − log(v_obs_i)
```

No per-row tuple rebuild — `cha_macro_tuple` stays at `Et = 1`, and `Et_i`
multiplies the prediction. This is what realizes "single global scale": differing
enzyme amounts across reactions are absorbed by the recorded `Et`, while the
single shared `:kcat` coord is jointly constrained by all rows. `Et` is **data,
not a fit parameter** — the parameterization is unchanged. In relative mode `Et`
is unused (centering discards it).

## 8. Data: the `[G6PD] (nM)` requirement

When `scale=:absolute`, the loader **asserts** the corpus contains a
`[G6PD] (nM)` column and reads it into a per-row `Et` vector; a missing column is
a hard error that names the column. Relative mode ignores the column (optional
there).

Config: `g6pd_config()` gains
`enzyme_conc_col = "[G6PD] (nM)"`, `enzyme_conc_unit = :nM`. This is an enzyme
concentration, not a metabolite, so it lives **outside** the `metabolites` map.

**Units:** the loader converts `[G6PD]` nM → M (÷1e9) and Rate → M·s⁻¹, so
`kcat` comes out in s⁻¹ and the `[10, 1000] s⁻¹` bound holds. The exact rate
unit conversion pins to the new dataset's rate column and is a data-prep
assumption to confirm at implementation time.

**Data assumptions for the absolute corpus** (the user's forthcoming dataset):
single rate scale (one purification workflow), forward-only (no PGLn/reverse
rows), same EnzymeData schema as the mixed corpus **plus** the `[G6PD] (nM)`
column.

## 9. Pins & the escalation ladder

Reverse and near-degenerate constants are pinned via an explicit
`pins::Dict{Symbol,Float64}` (coord ⇒ log10 value) keyword on the fit entry
points, merged over the mode-derived pins and validated by the existing
`_assert_pin_is_coord` guard (errors on a non-coord or silent no-op). Explicit
values let the user paste **data-determined** constants from a prior relative fit
on the mixed corpus rather than being forced to literature values.

Division of labor between the two modes:

- **Relative + mixed corpus** → data-determine `Kd_6PGLn` and `Km_NADPH_rev`
  (informed by the literature reverse rows in the mixed corpus).
- **Absolute + new forward-only CSV** → discriminate variants, with those
  reverse constants pinned in.

Manual escalation ladder (run at **`mode1`** — forward `Ki`s free — so the
explicit `pins` are the sole source of clamping and each rung is honest; using
`mode2` here would pre-pin `Ki_ATP`/`Ki_NADPH` and defeat rungs 1–4):

```
rung 1:  scale=:absolute, anchor_reverse=false                              # nothing pinned
rung 2:  scale=:absolute, anchor_reverse=false, pins=Dict(:Kd_6PGLn=>…)     # + PGLn (K_PGLn)
rung 3:  … pins=Dict(:Kd_6PGLn=>…, :Km_NADPH_rev=>…)                        # + reverse NADPH
rung 4:  … additionally :Ki_NADPH                                           # + forward dead-end Ki
```

Note on naming: the 6-phosphogluconolactone release step is rapid-equilibrium, so
its dissociation constant `Kd_6PGLn` *is* the 6PGLn product-inhibition constant;
there is no separate `Ki_6PGLn` coord. `Km_NADPH` maps to `Km_NADPH_rev`.

**README:** a dedicated section documents the `pins` override, the ladder, and
where the pinned values come from (user request).

## 10. Model selection & outputs

Absolute runs report **both** (user choice):

- **In-sample:** absolute loss + free-param count + Hessian identifiability
  (`cha_identifiable_functions` / `classify_cha`).
- **Leave-one-group-out CV:** fold by `Fig` (single article ⇒ the existing
  leave-one-article-out is degenerate). `cv.jl` gains a group-column fold helper
  paralleling `_article_folds`.

The seven output artifacts are retained. `report.md` and `provenance.toml`
record `scale=:absolute`, the fitted `kcat` (with the 150–250 s⁻¹ validation
verdict), and the CV fold basis. `micro_parameters.jl` carries a genuine absolute
scale (byproduct, labeled as such). The CairoMakie extension plots absolute
predicted-vs-measured rates (no per-figure recentering) when `scale=:absolute`.

## 11. Files touched

- `src/cha_fit.jl` — loss refactor (core + two aggregators); `cha_coords` /
  `cha_coord_bounds` scale-parameterized; `:kcat` handling; per-row `Et` in the
  core; `pins` merge; `_cha_loss_with_pins`/`cha_fit_candidate` scale threading.
- `src/run.jl` — `scale` through `fit_consensus_equation` / `_fit_consensus` /
  `run_variants`; guard; `pins` keyword; CV group-fold selection; report /
  provenance fields.
- `src/configs/G6PD.jl` — `enzyme_conc_col` / `enzyme_conc_unit`.
- `src/core/data.jl` — absolute-mode `[G6PD] (nM)` validation + `Et` column.
- `src/cv.jl` — leave-one-group-out (by `Fig`) fold helper.
- `src/cha_invert.jl` / `src/cha_deploy.jl` — `:kcat`-aware readoff/deploy.
- `src/cli.jl` — `--scale`.
- `ext/FitRateEquationMakieExt.jl` — absolute predicted-vs-measured plotting.
- `README.md` — `pins` override + absolute-mode section.
- `test/` — see §12.

## 12. Testing

- **Relative invariance:** `scale=:relative` reproduces the current fit
  bit-for-bit (regression lock over the loss refactor).
- **Uncentered arithmetic:** hand-computed small-case check of
  `cha_absolute_logratio_loss`.
- **`Et` prefactor:** doubling a row's `Et` shifts its predicted-rate log-ratio
  by exactly `log 2`.
- **`kcat ≡ kf` (fiber-free):** with `C = 1`, the fitted `:kcat` round-trips
  through `cha_macro_tuple`/`analytic_kcat`.
- **Scale recovery (synthetic):** generate single-scale forward data (with a
  `[G6PD] (nM)` column) from a known `kcat` and a known generating variant;
  assert absolute mode recovers `kcat` within tolerance and that the generating
  variant wins on CV — the end-to-end discrimination smoke test.
- **`[G6PD] (nM)` requirement:** absolute mode on a corpus lacking the column
  errors with a message naming the column.
- **Pins guard:** a bogus `pins` key errors; a valid one clamps.

## 13. Open items (confirm at implementation)

- Exact rate-column unit conversion for the new dataset (→ M·s⁻¹).
- Final `:kcat` bound width (default `[10, 1000] s⁻¹`).
- CV fold granularity: by `Fig` (default) vs `Experiment_date` — both are single-
  scale; `Fig` chosen for finer folds.
- Worktree Julia env: EnzymeRates resolution needs the symlink + copied Manifest
  fix before the baseline/tests run (per prior session notes).

## 14. Ethos note

Per the user's global CLAUDE.md, development cost is explicitly **not** weighted;
quality, simplicity-as-clarity, robustness, and long-term maintainability are.
This design chose the shared-core + two-aggregator loss and `kcat` as an honest
named coordinate (over a boolean-switched megafunction with `kf` smuggled into
the binding-constant dict) on those grounds, accepting more code for a clearer,
drift-free, more generalizable structure.

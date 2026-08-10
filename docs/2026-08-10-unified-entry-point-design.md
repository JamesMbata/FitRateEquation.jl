# Design: collapse the fitting API to a single entry point

- **Date:** 2026-08-10
- **Status:** Approved (design); implementation not started
- **Target repo:** `FitRateEquation.jl`
- **Version bump:** `0.4.0 → 0.5.0` (breaking)

## 1. Problem

The `0.4.0` rename of `run_all → fit_consensus_equation` left the public fitting
surface split across two layers that read as peers:

- `fit_consensus_equation(cfg; …)` — the general engine (takes a config object, requires
  `outdir`, does *not* set up workers, has its own raw budget knobs).
- `run_g6pd` / `run_pgd` / `run_hk1` / `run_g6pd_noatp` / `run_pgd_fullre` — thin-ish
  wrappers that each add real convenience logic on top (budget preset from `smoke`,
  `setup_workers`, default outdir naming, and — for the last two — a variant + row-filter
  preset).

The symptom that surfaced this: README §5 documents alternative rate laws as
`fit_consensus_equation(g6pd_config(); variants=[:variant_name], …)`, where `:variant_name`
is an unexplained placeholder and the general engine is exposed directly to novices
alongside the friendly `run_*` wrappers. Two documented ways in, an easy-to-misuse general
function, and an over-promised "bring your own columns" config path.

## 2. Goals / non-goals

**Goals**

- One documented entry point for fitting: `fit_consensus_equation`.
- Keep a small, discoverable novice front door (`run_g6pd`/`run_pgd`/`run_hk1`).
- Fold every piece of convenience logic (budget, workers, outdir, variant presets) into the
  one function so calling it directly is never the "expert-only, easy-to-misuse" path.
- Close the `:no_atp` row-filter footgun.
- Make the corpus column contract explicit and enforced.

**Non-goals**

- No change to the fitting math, Cha law, gauge, modes, or output artifacts.
- No renaming of `fit_consensus_equation` (the EnzymeRates `fit_rate_equation` collision that
  motivated the `0.4.0` name still stands — see `reference_enzymerates_export_collisions`).
- No deprecation shims: removed names are removed outright (this is a pre-1.0 breaking bump).
- No new custom-column / schema-remapping capability — the opposite: it is being removed.

## 3. Confirmed decisions

1. **Public API is symbol-only.** `fit_consensus_equation(:g6pd; …)`. Config builders
   (`g6pd_config`/`pgd_config`/`hk1_config`) stay as internal helpers but leave the exported
   surface and all user docs.
2. **Canonical columns enforced.** Corpus columns are validated on load against the enzyme's
   canonical schema; a non-conforming CSV errors with the exact expected header. **Extra
   annotation columns are allowed (ignored), not rejected** — "no custom columns" means "no
   remapping / bring-your-own-schema," not "the file may contain nothing else."
3. **`:no_atp` / `:full_re` become variant arguments.** `run_g6pd_noatp` and `run_pgd_fullre`
   are **deleted outright** (no deprecation). `run_g6pd`/`run_pgd`/`run_hk1` (one per enzyme,
   deploy law) survive as thin discoverability aliases.
4. **`smoke=false` defaults to the full runner budget** `(n_restarts=48, maxiter=1000,
   maxtime=300)`, superseding the old implicit `fit_consensus_equation` defaults
   `(8, 1e6, 20)`.

## 4. Public surface: before → after

| | Before (`0.4.0`) | After (`0.5.0`) |
|---|---|---|
| Fit fns | `fit_consensus_equation(cfg)` + `run_g6pd` + `run_pgd` + `run_hk1` + `run_g6pd_noatp` + `run_pgd_fullre` | `fit_consensus_equation(:enzyme)` + `run_g6pd` / `run_pgd` / `run_hk1` |
| Config builders | `g6pd_config` / `pgd_config` / `hk1_config` **exported** | same fns, **un-exported** |
| Deleted | — | `run_g6pd_noatp`, `run_pgd_fullre` |

Exported fitting surface: 6 → 4, only one of which is the documented path.

## 5. The unified signature

```julia
fit_consensus_equation(
    enzyme::Symbol;                 # :g6pd | :pgd | :hk1  (case-normalized)
    variants = <deploy default for enzyme>,   # e.g. [:no_atp], [:full_re], [:no_g6p_both_deadends]
    data_csv = nothing,             # your own file; must use the canonical columns
    smoke    = false,               # false → full budget, true → fast plumbing budget
    outdir   = nothing,             # nothing → results/<LABEL>_<date>[_smoke]/
    nprocs   = nothing,             # forwarded to setup_workers (idempotent)
    anchor_reverse = <variant-aware default>,   # G6PD-only; semantics unchanged
    # power-user escape hatch — overrides the smoke→budget mapping when non-nothing:
    n_restarts = nothing, maxiter = nothing, maxtime = nothing,
    seed = 1,
)
```

- The current `fit_consensus_equation(cfg; outdir::REQUIRED, …)` becomes an **internal**
  method (working name `_fit_consensus(cfg; …)`), carrying the existing load → filter →
  build-dataset → pmap → reduce → `write_outputs` pipeline unchanged.
- The public symbol method: normalizes the enzyme symbol, builds the config internally
  (`g6pd_config(; data_csv)` etc.), resolves the variant profiles (§6), maps `smoke`+explicit
  overrides to a budget (§7), calls `setup_workers(nprocs)`, computes the default outdir
  (§7), then forwards to `_fit_consensus`.
- HK1 stays guarded: `fit_consensus_equation(:hk1)` raises the same clear
  `HK1_AVAILABLE == false` error the current `run_hk1` does.

## 6. Variant registry (folding `:no_atp` / `:full_re` in correctly)

`:no_atp` cannot naively become "`variants=[:no_atp]`": it also requires `drop_atp_rows` and a
distinct outdir label. Each variant therefore gets a **profile** in the enzyme registry
(`src/enzyme_wiring.jl`, backed by the per-enzyme files), the single source of truth:

```
variant_profile(:g6pd, :SS_NADPH_release_rate_eq) → (mechanism=…, row_filter=identity,      label="")
variant_profile(:g6pd, :no_atp)                    → (mechanism=…, row_filter=drop_atp_rows, label="noatp")
variant_profile(:g6pd, :no_g6p_atp_deadend)        → (mechanism=…, row_filter=identity,      label="")
variant_profile(:pgd,  :cha_base)                  → (mechanism=…, row_filter=identity,      label="")
variant_profile(:pgd,  :full_re)                   → (mechanism=…, row_filter=identity,      label="fullre")
```

The profile carries the **mechanism**, the **row_filter**, and the **outdir label** —
centralizing what is currently scattered across `run_g6pd_noatp` (variant + filter + label),
`_mech_for`, and the per-runner outdir strings.

The `anchor_reverse` default is **out of scope for the profile**: it stays in the existing
variant-aware logic (`_default_anchor_reverse` / `_G6PD_ANCHOR_OPTIONAL_VARIANTS`), unchanged.
Folding it in would add risk to working G6PD-specific behavior for no user-facing benefit.

**Correctness win:** selecting `:no_atp` now *always* pulls its `drop_atp_rows` filter. Today,
`fit_consensus_equation(cfg; variants=[:no_atp])` without the filter silently fits an ATP-blind
law to ATP-bearing rows — a latent footgun this closes.

**Row-filter conflict rule:** one run builds one `Dataset`, so it can apply only one filter.

- All selected variants declaring `identity` (or the same non-identity filter) → fine
  (unchanged multi-variant behavior; variants land in separate leaderboards in one dir).
- Two selected variants declaring *different* non-identity filters → **error**. In practice
  `:no_atp` is always fit alone, so this never bites normal use; the error prevents a silently
  wrong corpus.

## 7. Budget / workers / outdir

- **Budget:** `smoke=false` → `(48, 1000, 300)`; `smoke=true` → `(2, 150, 120)`. Any of
  `n_restarts` / `maxiter` / `maxtime` passed non-`nothing` overrides the corresponding preset
  field (power-user escape). This supersedes the old implicit `(8, 1e6, 20)` defaults, which
  were reachable only by internal/test callers.
- **Workers:** the unified fn calls `setup_workers(nprocs)`. `setup_workers` is already
  idempotent (no-ops when workers exist, is serial for `n ≤ 1`), so tests that manage workers
  themselves are unaffected. Under SLURM it uses `SlurmManager` and ignores `nprocs`, as today.
- **Outdir label** (`LABEL` in `results/<LABEL>_<date>[_smoke]/`):
  - deploy variant → `<ENZYME>` (e.g. `G6PD`, `PGD`).
  - single non-deploy variant selected → `<ENZYME>_<variant-label>` (reproduces the old
    `G6PD_noatp_…`, `PGD_fullre_…`).
  - multiple variants → `<ENZYME>` (they share one dir; separate leaderboards, as today).

## 8. Schema enforcement

On load (in / immediately after `read_corpus`, `src/core/data.jl`), validate corpus columns
against the enzyme's canonical schema derived from the internal config:

- **Required:** the enzyme's metabolite columns (e.g. G6PD `[NADP] (uM)`, `[G6P] (uM)`,
  `[NADPH] (uM)`, `[PGLn] (uM)`, `[ATP] (uM)`) plus `Rate_V`, `Article`, `Fig`, `Apparent_Keq`.
- **Missing / misnamed required column →** error naming the offending column(s) and printing
  the exact expected header for that enzyme.
- **Extra columns →** allowed and ignored (the bundled corpora carry annotation columns such
  as `Experiment_date`).

## 9. `run_*` wrappers + CLI

- `run_g6pd(; outdir, smoke, nprocs, anchor_reverse) = fit_consensus_equation(:g6pd; outdir,
  smoke, nprocs, anchor_reverse)`; same one-liner for `run_pgd`, `run_hk1` (HK1 guard kept).
  Retained purely for tab-completion discoverability (`run_<TAB>`).
- **CLI (`src/cli.jl`):**
  - subcommands become `g6pd | pgd | hk1 | plot | help` (drop `g6pd-noatp`).
  - add `--variant NAME` (maps to `variants=[Symbol(NAME)]`).
  - `--data CSV` becomes valid for any enzyme subcommand (no longer gated to `g6pd-noatp`).
  - `g6pd-noatp` → `g6pd --variant no_atp`.

## 10. Blast radius / touch list

| Area | Change |
|---|---|
| `src/run.jl` | public symbol method + internal `_fit_consensus`; budget/outdir/worker folding; delete `run_g6pd_noatp`, `run_pgd_fullre` |
| `src/enzyme_wiring.jl` + `src/enzymes/*.jl` | `variant_profile` registry (mechanism + row_filter + label); `anchor_reverse` default logic unchanged |
| `src/core/data.jl` | column-schema validation + error message |
| `src/configs/*.jl` | keep; drop from exports |
| `src/FitRateEquation.jl` | export list: drop `g6pd_config`/`pgd_config`/`hk1_config`, `run_g6pd_noatp`, `run_pgd_fullre` |
| `src/cli.jl` | `--variant` flag; generic `--data`; drop `g6pd-noatp` sub |
| `test/` | ~16 files reference `*_config` (47 hits) → symbol API or path-qualified; update `test_cli.jl`, `test_runners.jl`, `test_cha_noatp.jl`, `test_cha_g6pd_deadend_variants.jl`, `test_plot_render.jl` for the deleted runners/budget default |
| `README.md` §5, `AGENTS.md` | rewrite to the single documented path; delete the custom-columns over-promise; document the schema contract |
| `Project.toml` | version `0.4.0 → 0.5.0` |

## 11. Testing strategy

- **Byte-identity smoke fixtures** (`test/test_byte_identity.jl`): must still pass — the
  variant×mode×name structure is unchanged. Verify `run_g6pd(smoke=true, nprocs=1)` /
  `run_pgd(smoke=true, nprocs=1)` still resolve through the new wrapper to the same structure.
- **Parallel equivalence** (`test/test_parallel_equivalence.jl`): confirm the auto
  `setup_workers` call does not perturb the serial-vs-pmap bit-identity guarantee (order-
  preserving per-cell seed makes results worker-count-invariant; the test controls `nprocs`).
- **New tests:**
  - schema validation: a CSV missing a required column errors with the expected header; a CSV
    with an extra column fits fine.
  - `:no_atp` auto-filter: `fit_consensus_equation(:g6pd; variants=[:no_atp])` drops ATP rows
    (fit_corpus row count matches the old `run_g6pd_noatp`).
  - row-filter conflict: two different non-identity filters in one run errors.
  - outdir labeling: deploy vs `:no_atp` vs `:full_re` default dirs.
  - CLI: `g6pd --variant no_atp` equals the old `g6pd-noatp`; unknown `--variant` errors.
- **Full local + CI green** before tagging `0.5.0` (registered package — see the
  `verify-tests-local-and-ci` skill).

## 12. Migration notes (for the CHANGELOG / breaking-change body)

- `run_g6pd_noatp()` → `fit_consensus_equation(:g6pd; variants=[:no_atp])` (or CLI
  `g6pd --variant no_atp`).
- `run_pgd_fullre()` → `fit_consensus_equation(:pgd; variants=[:full_re])`.
- `fit_consensus_equation(g6pd_config(); …)` → `fit_consensus_equation(:g6pd; …)`;
  `g6pd_config`/`pgd_config`/`hk1_config` are no longer exported.
- Custom column names / metabolite remapping are no longer supported; conform your CSV to the
  canonical header (printed by the new validation error).

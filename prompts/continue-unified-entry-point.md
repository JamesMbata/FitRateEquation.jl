# Handoff: execute the FitRateEquation single-entry-point refactor

You are continuing an approved, planned refactor of the **`FitRateEquation.jl`** package
(`/home/james/projects/FitRateEquation.jl`). Brainstorming, spec, and the implementation plan
are all done and committed. Your job is to **execute the plan, task by task, using
subagent-driven development.**

## First actions (do these before anything else)

1. Invoke the **`superpowers:subagent-driven-development`** skill and follow it: one fresh
   subagent per task, a combined spec+quality review between tasks (one review pass, not two —
   per the maintainer's standing preference), then move to the next task.
2. Read these two files in full — they are the source of truth:
   - Spec: `docs/2026-08-10-unified-entry-point-design.md`
   - Plan: `docs/2026-08-10-unified-entry-point-plan.md`  ← the task-by-task checklist you execute
3. Confirm you are on branch **`feat/unified-entry-point`** (`git -C /home/james/projects/FitRateEquation.jl branch --show-current`). Do NOT create a new branch. All work lands here; a PR to `main` comes at the end (Task 9).

> Note: `docs/superpowers/` is **gitignored** in this repo, so the spec and plan live in the
> tracked `docs/` directory (flat), not under `docs/superpowers/`. Read the paths above.

## What the refactor does (one paragraph)

Collapse the fitting API to one documented entry point, `fit_consensus_equation(:enzyme; …)`,
which folds in enzyme selection, the `smoke`→budget mapping, worker setup, auto-outdir, and
per-variant presets (row filter + outdir label). `run_g6pd`/`run_pgd`/`run_hk1` survive as thin
discoverability aliases; `run_g6pd_noatp`/`run_pgd_fullre` are deleted (their `:no_atp`/`:full_re`
laws become `variants=[…]` arguments). The config builders (`g6pd_config`/`pgd_config`/`hk1_config`)
are demoted to un-exported internals, and corpus columns are validated against a canonical schema
on load. Breaking bump `0.4.0 → 0.5.0`.

## The nine tasks (all in the plan, with full code)

1. `variant_profile` registry (row_filter + outdir label) — additive
2. Corpus column-schema validation in `read_corpus` — additive
3. New public `fit_consensus_equation(::Symbol; …)`, coexisting with the cfg-method
4. Thin `run_*` aliases; delete `run_g6pd_noatp`/`run_pgd_fullre` + `_run_enzyme`
5. Rename the cfg-method → internal `_fit_consensus`; repoint callers
6. Un-export `g6pd_config`/`pgd_config`/`hk1_config` (16 test files, 47 refs)
7. CLI `--variant` rework; drop `g6pd-noatp` subcommand
8. Docs: README §5 + §8, AGENTS.md
9. Version bump `0.5.0` + full local & CI verification

Execute them **in order** — the sequence is deliberately chosen to keep the test suite green
between tasks (the new symbol method lands alongside the old cfg-method in Task 3; the breaking
rename is Task 5, after the new path is proven).

## Hard constraints (carry into every task)

- **Do NOT rename `fit_consensus_equation`.** The name avoids a collision with EnzymeRates'
  exported `fit_rate_equation` under `using FitRateEquation, EnzymeRates`. Task 9 step 3 smoke-tests
  this co-import — keep it.
- **No deprecation shims.** Deleted names are gone outright.
- **Public API is symbol-only** (`:g6pd`/`:pgd`/`:hk1`, case-normalized). Never surface a config
  object in user-facing signatures or docs.
- **Determinism gate stays green:** `test/test_byte_identity.jl` (variant×mode×name structure) and
  `test/test_parallel_equivalence.jl` (serial ≡ pmap bit-identity). Verify both in Tasks 3 and 5.
- **No change to the fitting math**, Cha law, gauge, modes, or the seven output artifacts.
- **Commits:** Conventional Commits with a scope; end each body with:
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`

## Execution gotchas specific to this repo (important)

- **Julia runs are slow and must not run concurrently.** Concurrent `julia` precompilation
  deadlocks on the precompile lock. Run each task's test step **synchronously and one at a time**
  (a subagent running a foreground `julia --project test/<file>.jl` is fine; do not launch two Julia
  processes at once). The per-task smoke-budget tests are ~1–3 min each.
- **The full suite (Task 9, `Pkg.test()`) is long.** Run it from the **main session** with
  `run_in_background` and wait for the completion notification — do not run it inside a background
  subagent (background subagents get SIGTERM'd when the parent yields).
- **`smoke=true, nprocs=1`** for every in-task fit — keeps runs fast and deterministic.
- **Plan open item (Task 3):** confirm the `fit_corpus.csv` reader name in `src/plot_support.jl`
  (the plan assumes `read_fit_corpus`); if it differs, use the actual name. The assertion only needs
  the `ATP` column of `fit_corpus.csv`.
- **`variant_profile` placement:** the plan puts it in `src/run.jl` (next to `drop_atp_rows`), not
  `enzyme_wiring.jl` as the spec sketched — because `drop_atp_rows` is defined in `run.jl` and
  `enzyme_wiring.jl` loads earlier. Follow the plan.

## Definition of done

All nine tasks committed on `feat/unified-entry-point`; `julia --project -e 'using Pkg; Pkg.test()'`
green locally; the co-import smoke passes; then push and open a PR to `main` and confirm the full CI
run is green (use the `verify-tests-local-and-ci` skill — this is a registered package version bump,
so local and CI must agree). Report the PR link.

Start now: invoke `superpowers:subagent-driven-development`, read the spec + plan, then dispatch the
Task 1 subagent.

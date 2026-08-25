# Handoff: Execute the G6PD absolute-scale fitting mode (inline execution)

## Your task

Implement the **G6PD absolute-scale fitting mode** for `FitRateEquation.jl` by
executing an existing, approved implementation plan **task-by-task, inline in
this session, using the `superpowers:executing-plans` skill** (batch execution
with review checkpoints — no subagents). The design and plan are already written,
reviewed, and committed — you are in the *execution* phase, not design.

**Invoke `superpowers:executing-plans` first**, then work the plan's 12 tasks in
order in this session.

## Where everything is

- **Repo:** `FitRateEquation.jl` (sibling of `PentosePhosphatePathway.jl` and
  `PPP_Experiments`).
- **Worktree (work here, do NOT create a new one):**
  `/home/james/projects/FitRateEquation.jl/.worktrees/absolute-scale-mode`
- **Branch:** `absolute-scale-mode` (already checked out in that worktree).
- **Spec:** `docs/2026-08-24-g6pd-absolute-scale-fitting-mode-design.md`
- **Plan (your task list):** `docs/2026-08-24-g6pd-absolute-scale-fitting-mode-plan.md`
- Both are committed (spec `2bd4953`, plan `111e766`). Read the plan in full
  before starting; read the spec for the "why" behind any task.

## What is being built (one-paragraph context)

`FitRateEquation.jl` fits G6PD kinetic data to a Cha-form rate law. Today it uses
a per-(Article,Fig) **mean-centered** log-ratio loss that discards each group's
overall rate scale — which is what lets structurally different candidate
mechanisms fit equally well. This feature adds an orthogonal
`scale=:relative|:absolute` axis. `scale=:absolute` uses an **uncentered** loss
so absolute rate magnitudes discriminate the dead-end variants, and makes the
turnover `kcat` a free coordinate. The user is building a new single-scale,
forward-only dataset (one improved purification; specific activity 181 U ⇒
`kcat ≈ 178 s⁻¹`) with a per-reaction `[G6PD] (nM)` enzyme-concentration column;
this mode is being built in preparation for that data.

## Key design decisions already locked (do not re-litigate)

- **Orthogonal axis:** `scale::Symbol=:relative` (default) / `:absolute`, threaded
  `fit_consensus_equation → _fit_consensus → cha_fit_candidate → loss`. Composes
  with the existing `mode1`/`mode2` (pins) and `variant` (mechanism) axes.
- **Loss = shared core + two aggregators.** Extract `_cha_row_logratios!`
  (per-row arithmetic) and keep a thin `cha_centered_logratio_loss` (variance) +
  new `cha_absolute_logratio_loss` (uncentered sum-of-squares). Relative output
  MUST stay **byte-identical** (locked by `test/test_byte_identity.jl`).
- **`kcat` is a first-class coord** via a scale-parameterized
  `cha_coords(enzyme, variant; scale)` (single source of truth — bounds, pins,
  overwrite all key off it). Bound `[10,1000] s⁻¹`. Expected 150–250 s⁻¹ is a
  report verdict, never a clamp.
- **Fiber-free `C=1`.** In absolute mode `kcat ≡ kf`, realized by fitting `kf=kcat`
  at `CHA_ABS_RELEASE_RATE = 1e8`. This is exact and numerically safe *because
  forward-only data has `PGLn = 0`*, so in `cha_rate_G6PD` the fiber term
  `kf·gAB/koffQ → 0` and `konQ` only ever appears as `konQ/koffQ = 1/Km_NADPH_rev`.
- **Per-row `Et`** (Molar) from a new required `[G6PD] (nM)` column, applied as a
  linear prefactor on the prediction (`v = Et · kcat · f`). New 5th `Dataset`
  field with a 4-arg back-compat constructor (existing tests construct
  `Dataset(concs, rate, group, keq)` positionally — must keep working).
- **Explicit `pins::Dict{Symbol,Float64}` override** (guarded by
  `_assert_pin_is_coord`) so the user pins reverse constants (`Kd_6PGLn`,
  `Km_NADPH_rev`, then `Ki_NADPH`) to data-determined values from a prior
  relative fit. README-documented. The escalation ladder runs at `mode1`.
- **Model selection:** both leave-one-group-out CV (fold by `Fig` via new
  `_group_folds`) AND in-sample loss + identifiability.
- **Guards:** `scale=:absolute` errors for any enzyme ≠ `:G6PD`; and errors if the
  corpus lacks a finite `[G6PD] (nM)` column.

## Execution constraints for THIS repo (read before starting)

1. **All Julia runs stay in this (main) session — that's the point of inline
   execution.** Do NOT offload Julia to subagents (background subagents get
   SIGTERM'd on yield). Run tasks **sequentially**; never launch two `julia`
   processes at once (they deadlock on the precompile lock).
2. **Per-task vs full-suite runs.** During a task, run only the single new/changed
   test file — `julia --project test/test_xxx.jl` — for a fast loop. For the final
   full-suite gate (`julia --project -e 'using Pkg; Pkg.test()'`, can take many
   minutes) use the harness `run_in_background` and wait for the completion
   notification (it completes + notifies reliably from the main session — verified
   for ~24-min runs). Do not hand long runs off for manual checking.
3. **Worktree Julia env.** Manifest.toml is gitignored, so the worktree needs
   `Pkg.instantiate()` first (Task 0). EnzymeRates resolves from a `[sources]`
   git URL; if that fails offline, copy the main checkout's
   `/home/james/projects/FitRateEquation.jl/Manifest.toml` into the worktree and
   re-instantiate. **Baseline must be green before Task 1** (it was deferred at
   plan-writing time).
4. **Bit-identity lock.** After every task that touches `src/cha_fit.jl` or
   `src/run.jl`, run `test/test_byte_identity.jl` — relative-mode output must not
   drift. If it drifts, the refactor changed float fold order; fix before moving on.
5. **Highest-uncertainty task: Task 6.** Threading `scale`/`pins` through the
   `_build_tasks` / `_run_fit_task` / `_reduce_cells` machinery in `run.jl` was
   only partially mapped when the plan was written. **Grep all `resolve_cha_pins`
   and `cha_fit_candidate` call sites in `run.jl`** and thread through each; check
   that diff carefully.
6. **Two plan spots say "match the actual shape in the file"** (not placeholders):
   `classify_cha`'s exact signature (Task 8) and the CLI arg-parser name (Task 10).
   Read the current code and match it.

## Inline-execution rhythm (from executing-plans)

Work in batches with checkpoints. For each task: write the failing test → run it
(foreground) and confirm it fails → implement → run it and confirm it passes →
run `test/test_byte_identity.jl` if `cha_fit.jl`/`run.jl` changed → commit.
Pause at natural checkpoints (e.g. after Tasks 1–3 loss/Et core, after Task 6
wiring, after Task 8 deploy) to let the user review before continuing. Keep the
plan's checkboxes updated as you go.

## Repo conventions (from CLAUDE.md)

- Run scripts: `julia --project <path>`. Full tests: `julia --project -e 'using Pkg; Pkg.test()'`.
- Units: Molar model; µM/nM converted on import.
- 6-phosphogluconolactone is `PGLn`/`6PGLn`; 6-phosphogluconate is `PGA` (never `6-PG`).
- Commits: Conventional Commits `type(scope): summary`; end body with
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. Ask before new branches.
- Do NOT weight development cost in design calls — prefer clarity, robustness,
  maintainability (this drove the shared-core + honest-`kcat` choices).
- Design docs / plans live in `docs/` with a `YYYY-MM-DD-` prefix (NOT
  `docs/superpowers/`, which is gitignored here).

## How to start

1. `cd /home/james/projects/FitRateEquation.jl/.worktrees/absolute-scale-mode`
2. Read the plan and spec.
3. Invoke `superpowers:executing-plans`.
4. Do **Task 0** (instantiate + green baseline) before any implementation.
5. Execute Tasks 1–12 in order with the rhythm above, pausing at checkpoints.
6. When all 12 pass, run the full suite (background) from this session and report.

## Definition of done

All 12 tasks committed on `absolute-scale-mode`; `julia --project -e 'using Pkg;
Pkg.test()'` green; `test/test_byte_identity.jl` confirms relative mode is
byte-identical to `main`; the end-to-end synthetic test (Task 12) recovers the
planted `kcat` within 10%. Do NOT merge — leave the branch for the user to review
(finishing-a-development-branch is a separate, later step).

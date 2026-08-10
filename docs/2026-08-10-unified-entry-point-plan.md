# Unified Entry Point Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse the FitRateEquation fitting API to one documented entry point — `fit_consensus_equation(:enzyme; …)` — with `run_g6pd`/`run_pgd`/`run_hk1` kept as thin discoverability aliases.

**Architecture:** A new public method `fit_consensus_equation(enzyme::Symbol; …)` absorbs the budget/worker/outdir/variant-preset logic the `run_*` wrappers hold today, then forwards to the existing pipeline renamed to the internal `_fit_consensus(cfg; …)`. Variant presets (`:no_atp` filter, outdir labels) move into a `variant_profile` lookup; corpus columns are validated on load; the config builders are demoted to un-exported internals.

**Tech Stack:** Julia 1.11+, `Distributed`/`ClusterManagers` (worker setup), `CSV`/`DataFrames` (corpus), `Test` (suite). Depends on upstream `DenisTitovLab/EnzymeRates.jl` (pinned in `[sources]`).

## Global Constraints

- **Package version bump:** `0.4.0 → 0.5.0` (breaking) — set in `Project.toml`.
- **Do not rename `fit_consensus_equation`** — the EnzymeRates `fit_rate_equation` export collision that motivated the name still stands. Verify exported-name behavior under `using FitRateEquation, EnzymeRates`.
- **No deprecation shims:** removed names (`run_g6pd_noatp`, `run_pgd_fullre`, exported `*_config`) are removed outright.
- **Public API is symbol-only:** `:g6pd`/`:pgd`/`:hk1` (case-normalized to internal `:G6PD`/`:PGD`/`:HK1`). Never expose a config object in user-facing signatures/docs.
- **Determinism gate stays green:** `test/test_byte_identity.jl` (variant×mode×name structure) and `test/test_parallel_equivalence.jl` (serial≡pmap bit-identity) must pass unchanged.
- **No change to the fitting math**, Cha law, gauge, modes, or the seven output artifacts.
- **Budget presets (verbatim):** `smoke=true` → `(n_restarts=2, maxiter=150, maxtime=120.0)`; `smoke=false` → `(n_restarts=48, maxiter=1000, maxtime=300.0)`.
- **Enzyme display names** (for outdir/config `name`): `"G6PD"`, `"PGD"`, `"HK1"`.
- **Commit style:** Conventional Commits with scope; end body with the repo's `Co-Authored-By` trailer.

---

## File Structure

- `src/run.jl` — add `variant_profile`, the enzyme/config/label/filter helpers, and the public `fit_consensus_equation(::Symbol; …)`; rename the current `fit_consensus_equation(cfg; …)` to `_fit_consensus`; rewrite `run_g6pd`/`run_pgd`/`run_hk1`; delete `run_g6pd_noatp`, `run_pgd_fullre`, `_run_enzyme`.
- `src/core/data.jl` — add `_validate_corpus_columns` and call it in `read_corpus`.
- `src/FitRateEquation.jl` — trim exports (drop `g6pd_config`/`pgd_config`/`hk1_config`, `run_g6pd_noatp`, `run_pgd_fullre`).
- `src/cli.jl` — subcommands `g6pd|pgd|hk1|plot|help`; add `--variant`; generalize `--data`; drop `g6pd-noatp`.
- `test/` — new tests for the symbol method, schema validation, `:no_atp` auto-filter, outdir labels, filter conflict, CLI `--variant`; update tests that call the old cfg-method or the deleted runners.
- `README.md` §5 + `AGENTS.md` — rewrite to the one path; document the schema contract; delete the custom-columns over-promise.
- `Project.toml` — version bump.

---

### Task 1: `variant_profile` registry (row_filter + outdir label)

Additive; no existing behavior changes. Lives in `src/run.jl` next to `drop_atp_rows` (not `enzyme_wiring.jl` as the spec sketched, because `drop_atp_rows` is defined in `run.jl` and `enzyme_wiring.jl` loads earlier — keeping them together respects load order and cohesion).

**Files:**
- Modify: `src/run.jl` (after `drop_atp_rows`, ~line 170)
- Test: `test/test_variant_profile.jl` (new)

**Interfaces:**
- Consumes: `drop_atp_rows` (existing, `src/run.jl:167`).
- Produces: `variant_profile(enzyme::Symbol, variant::Symbol) -> NamedTuple{(:row_filter,:label)}` where `row_filter` is a `DataFrame->DataFrame` function (default `identity`) and `label::String` (default `""`). `enzyme` is the internal upper-case symbol (`:G6PD`/`:PGD`/`:HK1`).

- [ ] **Step 1: Write the failing test**

Create `test/test_variant_profile.jl`:

```julia
using Test, FitRateEquation
using FitRateEquation: variant_profile, drop_atp_rows

@testset "variant_profile" begin
    # :no_atp carries the ATP row filter and the "noatp" outdir label
    p = variant_profile(:G6PD, :no_atp)
    @test p.row_filter === drop_atp_rows
    @test p.label == "noatp"

    # :full_re labels the outdir but does not filter rows
    p = variant_profile(:PGD, :full_re)
    @test p.row_filter === identity
    @test p.label == "fullre"

    # deploy / any other variant: no filter, no label
    @test variant_profile(:G6PD, :SS_NADPH_release_rate_eq) == (row_filter=identity, label="")
    @test variant_profile(:PGD, :cha_base) == (row_filter=identity, label="")
    @test variant_profile(:G6PD, :no_g6p_atp_deadend) == (row_filter=identity, label="")
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_variant_profile.jl`
Expected: FAIL — `UndefVarError: variant_profile not defined`.

- [ ] **Step 3: Write minimal implementation**

In `src/run.jl`, immediately after `drop_atp_rows` (after line 170):

```julia
# Per-(enzyme, variant) run profile: the row filter to apply to the corpus and the outdir
# label suffix. Centralizes what run_g6pd_noatp/run_pgd_fullre used to hardcode. Selecting a
# variant now ALWAYS pulls its filter, closing the footgun where a bare variants=[:no_atp]
# call fit an ATP-blind law to ATP-bearing rows. Default: no filtering, no label.
function variant_profile(enzyme::Symbol, variant::Symbol)
    enzyme === :G6PD && variant === :no_atp  && return (row_filter=drop_atp_rows, label="noatp")
    enzyme === :PGD  && variant === :full_re && return (row_filter=identity,      label="fullre")
    return (row_filter=identity, label="")
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_variant_profile.jl`
Expected: PASS (5 tests).

- [ ] **Step 5: Register the test file**

In `test/runtests.jl`, add `include("test_variant_profile.jl")` alongside the other `include(...)` lines.

- [ ] **Step 6: Commit**

```bash
git add src/run.jl test/test_variant_profile.jl test/runtests.jl
git commit -m "feat(run): variant_profile registry for row-filter + outdir label"
```

---

### Task 2: Corpus column-schema validation

Additive guard in the single loader. Errors early with the exact expected header when a required column is missing/misnamed; extra annotation columns are allowed.

**Files:**
- Modify: `src/core/data.jl` (`read_corpus`, ~line 61)
- Test: `test/test_schema_validation.jl` (new)

**Interfaces:**
- Consumes: `cfg.metabolites` (Dict sym → `(csv_col, unit)`), `cfg.rate_col`, `cfg.article_col`, `cfg.fig_col`, `cfg.keq_col` (all existing config fields).
- Produces: `_validate_corpus_columns(raw::DataFrame, cfg)` — throws `ErrorException` naming missing columns + expected header; returns `nothing` on success. Called first inside `read_corpus`.

- [ ] **Step 1: Write the failing test**

Create `test/test_schema_validation.jl`:

```julia
using Test, FitRateEquation, CSV, DataFrames
using FitRateEquation: g6pd_config, read_corpus

# Minimal 1-row G6PD corpus with the canonical columns.
function _write_canonical_csv(path; drop=nothing, extra=false)
    cols = Dict(
        "[NADP] (uM)" => 100.0, "[G6P] (uM)" => 200.0, "[NADPH] (uM)" => 0.0,
        "[PGLn] (uM)" => 0.0, "[ATP] (uM)" => 0.0,
        "Rate_V" => 1.0, "Article" => "T", "Fig" => "1", "Apparent_Keq" => 13.0,
    )
    drop === nothing || delete!(cols, drop)
    extra && (cols["Notes"] = "annotation")
    CSV.write(path, DataFrame(cols))
    path
end

@testset "schema validation" begin
    dir = mktempdir()

    # conforming corpus loads
    ok = _write_canonical_csv(joinpath(dir, "ok.csv"))
    @test read_corpus(g6pd_config(; data_csv=ok)) isa DataFrame

    # extra annotation column is allowed (ignored)
    ex = _write_canonical_csv(joinpath(dir, "extra.csv"); extra=true)
    @test read_corpus(g6pd_config(; data_csv=ex)) isa DataFrame

    # missing required column errors, and the message names it + the expected header
    bad = _write_canonical_csv(joinpath(dir, "bad.csv"); drop="[G6P] (uM)")
    err = try read_corpus(g6pd_config(; data_csv=bad)); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("[G6P] (uM)", err.msg)
    @test occursin("Rate_V", err.msg)   # expected header is printed
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_schema_validation.jl`
Expected: FAIL — the missing-column case currently throws an `ArgumentError` from `raw[!, col]`, not an `ErrorException` naming the header (assertion mismatch).

- [ ] **Step 3: Write minimal implementation**

In `src/core/data.jl`, add before `read_corpus` (after line 52):

```julia
# Validate that a raw corpus DataFrame carries every column the config requires (metabolite
# source columns + rate/article/fig/keq). Extra annotation columns are allowed and ignored.
# "no custom columns" here means no remapping/bring-your-own-schema — the required columns
# must be present under their canonical names. Errors early with the full expected header.
function _validate_corpus_columns(raw::AbstractDataFrame, cfg)
    required = String[cfg.metabolites[s][1] for s in metabolite_syms(cfg)]
    append!(required, String[cfg.rate_col, cfg.article_col, cfg.fig_col, cfg.keq_col])
    have    = Set(names(raw))
    missing_cols = [c for c in required if !(c in have)]
    isempty(missing_cols) && return nothing
    error("read_corpus: corpus $(cfg.data_csv) is missing required column(s): " *
          join(missing_cols, ", ") * ".\nExpected columns (extra columns are allowed): " *
          join(required, ", ") * ".")
end
```

Then make it the first statement inside `read_corpus`, right after the CSV read (`src/core/data.jl:62`):

```julia
function read_corpus(cfg)
    raw = CSV.read(cfg.data_csv, DataFrame)
    _validate_corpus_columns(raw, cfg)
    for s in metabolite_syms(cfg)
        # ... unchanged ...
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_schema_validation.jl`
Expected: PASS.

- [ ] **Step 5: Verify the bundled corpora still load**

Run: `julia --project -e 'using FitRateEquation; using FitRateEquation: read_corpus, g6pd_config, pgd_config; read_corpus(g6pd_config()); read_corpus(pgd_config()); println("ok")'`
Expected: prints `ok` (bundled corpora carry all required columns).

- [ ] **Step 6: Register + commit**

Add `include("test_schema_validation.jl")` to `test/runtests.jl`, then:

```bash
git add src/core/data.jl test/test_schema_validation.jl test/runtests.jl
git commit -m "feat(data): validate corpus columns against canonical schema on load"
```

---

### Task 3: New public `fit_consensus_equation(::Symbol; …)` (coexists with the cfg method)

Add the symbol method that folds in enzyme selection, variant profiles, budget, workers, and auto-outdir, then forwards to the existing cfg-method (still named `fit_consensus_equation(cfg; …)` at this stage — renamed in Task 5). Keeping both methods live means existing tests stay green until Task 5.

**Files:**
- Modify: `src/run.jl` (helpers near the runners ~line 635; new method after the cfg-method ~line 293)
- Test: `test/test_symbol_entry.jl` (new)

**Interfaces:**
- Consumes: `run_variants` (`src/run.jl:37`), `_default_anchor_reverse` (`:59`), `variant_profile` (Task 1), `_budget` (`:635`), `_default_outdir` (`:640`), `setup_workers` (`src/worker_setup.jl:39`), `g6pd_config`/`pgd_config`/`hk1_config` (configs), `HK1_AVAILABLE` (`src/FitRateEquation.jl:29`), and the cfg-method `fit_consensus_equation(cfg; …)`.
- Produces:
  - `fit_consensus_equation(enzyme::Symbol; variants=nothing, data_csv=nothing, smoke::Bool=false, outdir=nothing, nprocs=nothing, anchor_reverse=nothing, n_restarts=nothing, maxiter=nothing, maxtime=nothing, seed::Int=1)` → results vector.
  - Helpers: `_canonical_enzyme(::Symbol)->Symbol`, `_enzyme_config(enz::Symbol, data_csv)->cfg`, `_combined_row_filter(enz::Symbol, variants)->Function`, `_labeled(name::AbstractString, label::AbstractString)->String`.

- [ ] **Step 1: Write the failing test**

Create `test/test_symbol_entry.jl`:

```julia
using Test, FitRateEquation
using FitRateEquation: _canonical_enzyme, _combined_row_filter, _labeled, drop_atp_rows

@testset "symbol entry helpers" begin
    @test _canonical_enzyme(:g6pd) === :G6PD
    @test _canonical_enzyme(:G6PD) === :G6PD
    @test _canonical_enzyme(:pgd)  === :PGD
    @test_throws ErrorException _canonical_enzyme(:nope)

    # single :no_atp variant -> its filter; deploy -> identity
    @test _combined_row_filter(:G6PD, [:no_atp]) === drop_atp_rows
    @test _combined_row_filter(:G6PD, [:SS_NADPH_release_rate_eq]) === identity

    @test _labeled("G6PD", "")      == "G6PD"
    @test _labeled("G6PD", "noatp") == "G6PD_noatp"
end

@testset "symbol entry smoke fit (G6PD)" begin
    dir = mktempdir()
    res = fit_consensus_equation(:g6pd; smoke=true, nprocs=1, outdir=dir)
    @test !isempty(res)
    @test isfile(joinpath(dir, "macro_constants.csv"))
    @test isfile(joinpath(dir, "fit_corpus.csv"))
end

@testset "symbol entry :no_atp auto-filters ATP rows" begin
    dir = mktempdir()
    fit_consensus_equation(:g6pd; variants=[:no_atp], smoke=true, nprocs=1, outdir=dir)
    corpus = FitRateEquation.read_fit_corpus(joinpath(dir, "fit_corpus.csv"))
    @test all(<=(0.0), corpus.ATP)   # ATP-bearing rows dropped
end

@testset "conflicting row filters error" begin
    # two DIFFERENT non-identity filters in one run is rejected. Simulate by a variant list
    # whose profiles disagree; only :no_atp has a non-identity filter today, so pair it with a
    # hand-rolled second filter via _combined_row_filter's contract.
    @test_throws ErrorException FitRateEquation._combined_row_filter_check(
        [drop_atp_rows, df -> df])
end
```

> Note: `read_fit_corpus` exists in `src/plot_support.jl` (reads a run's `fit_corpus.csv`). If its exact name differs, read that file and use the correct reader; the `fit_corpus.csv` has an `ATP` column in Molar.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_symbol_entry.jl`
Expected: FAIL — helpers/method undefined.

- [ ] **Step 3: Write minimal implementation**

In `src/run.jl`, add the helpers near the runner section (after `_default_outdir`, ~line 642):

```julia
# Map a user enzyme symbol (:g6pd or :G6PD, any case) to the internal upper-case symbol.
function _canonical_enzyme(e::Symbol)
    u = Symbol(uppercase(String(e)))
    u in (:G6PD, :PGD, :HK1) ||
        error("fit_consensus_equation: unknown enzyme :$e (expected :g6pd, :pgd, or :hk1)")
    u
end

# Build the (internal) config for an enzyme, optionally overriding the bundled corpus path.
function _enzyme_config(enz::Symbol, data_csv)
    if enz === :G6PD
        data_csv === nothing ? g6pd_config() : g6pd_config(; data_csv=data_csv)
    elseif enz === :PGD
        data_csv === nothing ? pgd_config() : pgd_config(; data_csv=data_csv)
    else
        data_csv === nothing ? hk1_config() : hk1_config(; data_csv=data_csv)
    end
end

# One run builds one Dataset, so it can apply only one row filter. Collect the distinct
# non-identity filters across the selected variants' profiles; >1 is a hard error.
function _combined_row_filter_check(filters)
    nonid = unique(f -> objectid(f), [f for f in filters if f !== identity])
    length(nonid) > 1 && error("fit_consensus_equation: selected variants declare conflicting " *
        "row filters; a single run can apply only one. Fit them in separate runs.")
    isempty(nonid) ? identity : nonid[1]
end

_combined_row_filter(enz::Symbol, variants) =
    _combined_row_filter_check([variant_profile(enz, v).row_filter for v in variants])

_labeled(name::AbstractString, label::AbstractString) =
    isempty(label) ? String(name) : string(name, "_", label)
```

Then add the public symbol method after the cfg-method (`fit_consensus_equation(cfg; …)` ends at line 293):

```julia
"""
    fit_consensus_equation(enzyme::Symbol; variants, data_csv, smoke, outdir, nprocs,
                           anchor_reverse, n_restarts, maxiter, maxtime, seed)

The single entry point. `enzyme` is `:g6pd`, `:pgd`, or `:hk1` (case-insensitive). Selects the
enzyme's deploy variant by default; pass `variants=[:no_atp]`, `[:full_re]`, `[:no_g6p_atp_deadend]`,
… to fit an alternative law (its row filter + outdir label are applied automatically). `smoke=true`
uses the fast plumbing budget. `data_csv` fits your own corpus (canonical columns required).
Writes the seven artifacts to `outdir` (default `results/<LABEL>_<date>[_smoke]/`) and returns the
results. `n_restarts`/`maxiter`/`maxtime` override the smoke→budget mapping when given.
"""
function fit_consensus_equation(enzyme::Symbol; variants=nothing, data_csv=nothing,
        smoke::Bool=false, outdir=nothing, nprocs=nothing, anchor_reverse=nothing,
        n_restarts=nothing, maxiter=nothing, maxtime=nothing, seed::Int=1)
    enz = _canonical_enzyme(enzyme)
    enz === :HK1 && !HK1_AVAILABLE &&
        error("HK1 is not available on this EnzymeRates build (deferred port). See AGENTS.md.")
    vars = variants === nothing ? run_variants(enz) : Vector{Symbol}(variants)
    ar   = anchor_reverse === nothing ? _default_anchor_reverse(enz, vars) : anchor_reverse
    cfg  = _enzyme_config(enz, data_csv)
    rf   = _combined_row_filter(enz, vars)
    label = length(vars) == 1 ? variant_profile(enz, vars[1]).label : ""
    b  = _budget(smoke)
    nr = n_restarts === nothing ? b.n_restarts : n_restarts
    mi = maxiter    === nothing ? b.maxiter    : maxiter
    mt = maxtime    === nothing ? b.maxtime    : maxtime
    od = outdir === nothing ? _default_outdir(_labeled(String(cfg.name), label), smoke) : outdir
    setup_workers(nprocs)
    @info "FitRateEquation run starting" enzyme=enz nworkers=nworkers() smoke outdir=od anchor_reverse=ar variants=vars
    fit_consensus_equation(cfg; outdir=od, n_restarts=nr, maxiter=mi, maxtime=mt, seed=seed,
                           variants=vars, row_filter=rf, anchor_reverse=ar)
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_symbol_entry.jl`
Expected: PASS (helpers + smoke fit + auto-filter + conflict error).

- [ ] **Step 5: Register + commit**

Add `include("test_symbol_entry.jl")` to `test/runtests.jl`, then:

```bash
git add src/run.jl test/test_symbol_entry.jl test/runtests.jl
git commit -m "feat(run): symbol-keyed fit_consensus_equation entry point"
```

---

### Task 4: Rewrite `run_*` wrappers; delete `run_g6pd_noatp`/`run_pgd_fullre`

Point the surviving wrappers at the symbol method and remove the two variant-runners and the now-unused `_run_enzyme`.

**Files:**
- Modify: `src/run.jl` (`run_g6pd`/`run_pgd`/`run_hk1` ~lines 660-728; delete `run_g6pd_noatp`, `run_pgd_fullre`, `_run_enzyme` ~line 644)
- Modify: `src/FitRateEquation.jl:59` (export line — drop the two names)
- Modify: `test/test_runners.jl`, `test/test_cha_noatp.jl`, `test/test_cha_g6pd_deadend_variants.jl`, `test/test_plot_render.jl` (drop deleted-runner references / switch to symbol API)

**Interfaces:**
- Consumes: `fit_consensus_equation(::Symbol; …)` (Task 3).
- Produces: `run_g6pd(; outdir, smoke, nprocs, anchor_reverse=true)`, `run_pgd(; outdir, smoke, nprocs)`, `run_hk1(; outdir, smoke, nprocs)` — each returns the results vector, unchanged semantics.

- [ ] **Step 1: Update the runner-behavior test first**

In `test/test_runners.jl`, remove any `run_g6pd_noatp`/`run_pgd_fullre` cases and assert the new equivalences. Add:

```julia
@testset "run_* wrappers forward to the symbol entry" begin
    dir = mktempdir()
    r1 = run_g6pd(; smoke=true, nprocs=1, outdir=joinpath(dir, "a"))
    r2 = fit_consensus_equation(:g6pd; smoke=true, nprocs=1, outdir=joinpath(dir, "b"))
    @test length(r1) == length(r2)   # same cell set

    # noatp is now a variant, not a runner
    @test !isdefined(FitRateEquation, :run_g6pd_noatp)
    @test !isdefined(FitRateEquation, :run_pgd_fullre)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project test/test_runners.jl`
Expected: FAIL — `run_g6pd_noatp`/`run_pgd_fullre` still defined.

- [ ] **Step 3: Rewrite the wrappers and delete the two runners**

Replace `run_g6pd`/`run_pgd`/`run_hk1` bodies (keep docstrings, trimmed) with:

```julia
run_g6pd(; outdir=nothing, smoke=false, nprocs=nothing, anchor_reverse=true) =
    fit_consensus_equation(:g6pd; outdir, smoke, nprocs, anchor_reverse)

run_pgd(; outdir=nothing, smoke=false, nprocs=nothing) =
    fit_consensus_equation(:pgd; outdir, smoke, nprocs)

run_hk1(; outdir=nothing, smoke=false, nprocs=nothing) =
    fit_consensus_equation(:hk1; outdir, smoke, nprocs)
```

Delete the entire `function run_g6pd_noatp(...) … end`, `run_pgd_fullre(...) = …`, and the now-unused `function _run_enzyme(...) … end` block. In `src/FitRateEquation.jl:59`, change the export line to:

```julia
export run_g6pd, run_pgd, run_hk1
```

- [ ] **Step 4: Fix remaining references to the deleted runners**

- `test/test_cha_noatp.jl` — replace `run_g6pd_noatp(...)` calls with `fit_consensus_equation(:g6pd; variants=[:no_atp], …)`.
- `test/test_cha_g6pd_deadend_variants.jl` — this already uses the cfg-method for ablation variants; leave the cfg-method calls for Task 5, but replace any `run_g6pd_noatp` reference with the symbol API.
- `test/test_plot_render.jl` — replace `run_pgd_fullre`/`run_g6pd_noatp` with `fit_consensus_equation(:pgd; variants=[:full_re], …)` / `fit_consensus_equation(:g6pd; variants=[:no_atp], …)`.

- [ ] **Step 5: Run the affected tests**

Run: `julia --project test/test_runners.jl && julia --project test/test_cha_noatp.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/run.jl src/FitRateEquation.jl test/test_runners.jl test/test_cha_noatp.jl test/test_cha_g6pd_deadend_variants.jl test/test_plot_render.jl
git commit -m "refactor(run)!: run_* thin aliases; drop run_g6pd_noatp/run_pgd_fullre"
```

---

### Task 5: Rename the cfg-method to internal `_fit_consensus`; repoint callers

Make the config-object path internal. Mechanical: same signature, new name.

**Files:**
- Modify: `src/run.jl` (the `function fit_consensus_equation(cfg; …)` def at line 263 → `_fit_consensus`; its call inside the symbol method from Task 3)
- Modify: cfg-method call sites in tests: `test/test_parallel_equivalence.jl`, `test/test_outputs.jl`, `test/test_pgd_outputs.jl`, `test/test_cha_g6pd_deadend_variants.jl`, `test/test_plot_consensus_fit.jl`, `test/test_plot_render.jl`

**Interfaces:**
- Consumes: nothing new.
- Produces: `_fit_consensus(cfg; outdir, n_restarts=8, maxiter=1_000_000, maxtime=20.0, seed=1, variants=run_variants(Symbol(cfg.name)), row_filter=identity, anchor_reverse=_default_anchor_reverse(...))` — the current cfg-method body verbatim, renamed. Not exported.

- [ ] **Step 1: Rename the definition**

In `src/run.jl:263`, change `function fit_consensus_equation(cfg; …)` to `function _fit_consensus(cfg; …)` (leave the body untouched). Update the docstring above it to describe the internal role.

- [ ] **Step 2: Repoint the symbol method**

In the Task 3 symbol method, change the final forwarding call from `fit_consensus_equation(cfg; …)` to `_fit_consensus(cfg; …)`.

- [ ] **Step 3: Repoint test call sites**

For each cfg-method caller, choose the target:
- **Low-level determinism tests** (`test_parallel_equivalence.jl`) — they need the config path with explicit worker control. Replace `fit_consensus_equation(cfg; kw…)` with `FitRateEquation._fit_consensus(cfg; kw…)` (path-qualified; identical signature).
- **Output/plot tests** (`test_outputs.jl`, `test_pgd_outputs.jl`, `test_plot_consensus_fit.jl`, `test_plot_render.jl`, `test_cha_g6pd_deadend_variants.jl`) — prefer the public symbol API where they just want a run: e.g. `fit_consensus_equation(:g6pd; variants=[:no_g6p_atp_deadend], smoke=true, nprocs=1, outdir=dir, anchor_reverse=true)`. Where a test asserts an exact low-budget number tied to the old cfg defaults, keep it on `FitRateEquation._fit_consensus(cfg; …)` to preserve the `(8, 1e6, 20)` defaults.

- [ ] **Step 4: Run the suite for the touched files**

Run: `julia --project test/test_parallel_equivalence.jl && julia --project test/test_outputs.jl && julia --project test/test_pgd_outputs.jl`
Expected: PASS. Also run `julia --project test/test_byte_identity.jl` — the variant×mode×name structure must be unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/run.jl test/test_parallel_equivalence.jl test/test_outputs.jl test/test_pgd_outputs.jl test/test_cha_g6pd_deadend_variants.jl test/test_plot_consensus_fit.jl test/test_plot_render.jl
git commit -m "refactor(run)!: rename cfg-method to internal _fit_consensus"
```

---

### Task 6: Un-export the config builders

Demote `g6pd_config`/`pgd_config`/`hk1_config` to internal. They remain callable as `FitRateEquation.g6pd_config(...)`.

**Files:**
- Modify: `src/FitRateEquation.jl:57` (export line)
- Modify: every test using bare `g6pd_config`/`pgd_config`/`hk1_config` (16 files; 47 refs)

**Interfaces:**
- Consumes: nothing new.
- Produces: no exported `*_config`; internal definitions unchanged.

- [ ] **Step 1: Remove from exports**

In `src/FitRateEquation.jl:57`, delete the line `export g6pd_config, pgd_config, hk1_config`.

- [ ] **Step 2: Repoint test references**

In each test file that used a bare `*_config`, add it to the file's `using` line:

```julia
using FitRateEquation: g6pd_config, pgd_config, hk1_config   # (only those actually used)
```

Affected files (from `grep -rln 'g6pd_config\|pgd_config\|hk1_config' test/`): `test_cha_classify.jl`, `test_plot_render.jl`, `test_cha_koffq.jl`, `test_data_path.jl`, `test_parallel_equivalence.jl`, `test_run_fit.jl`, `test_pgd_outputs.jl`, `test_data.jl`, `test_cha_pgd_fullre_deploy_classify.jl`, `test_cv.jl`, `test_cha_fit.jl`, `test_plot_consensus_fit.jl`, `test_cha_pgd_fullre_fit.jl`, `test_pgd_data.jl`, `test_outputs.jl`, `test_config_deploy_keq.jl`. (Some may already import them explicitly — check before adding to avoid a duplicate-import warning.)

- [ ] **Step 3: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (or only failures owned by later tasks — CLI/docs). Fix any missed `*_config` import.

- [ ] **Step 4: Commit**

```bash
git add src/FitRateEquation.jl test/
git commit -m "refactor(api)!: un-export g6pd_config/pgd_config/hk1_config"
```

---

### Task 7: CLI — `--variant`, generic `--data`, drop `g6pd-noatp`

**Files:**
- Modify: `src/cli.jl`
- Modify: `test/test_cli.jl`

**Interfaces:**
- Consumes: `fit_consensus_equation(::Symbol; …)`, `plot_consensus_fit`.
- Produces: subcommands `g6pd|pgd|hk1|plot|help`; flags `--smoke`, `--nprocs N`, `--outdir DIR`, `--data CSV`, `--variant NAME`. `cli_main(argv)` dispatches `g6pd/pgd/hk1` to `fit_consensus_equation(:g6pd|:pgd|:hk1; smoke, nprocs, outdir, data_csv, variants)`.

- [ ] **Step 1: Update the CLI test first**

In `test/test_cli.jl`, replace `g6pd-noatp` expectations. Add parser assertions:

```julia
using FitRateEquation: parse_cli
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
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project test/test_cli.jl`
Expected: FAIL — `o.variant` undefined; `g6pd-noatp` still accepted.

- [ ] **Step 3: Rewrite `src/cli.jl`**

```julia
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
```

> The old `cli_main(argv; dispatch=...)` injection is dropped — `test/test_cli.jl` should test `parse_cli` directly (the dispatch table existed only to stub runners; parser tests cover the surface without spawning fits). If a dispatch seam is still wanted, keep a `dispatch` kwarg mapping `Symbol(sub) => fit_consensus_equation`-closures, but the simpler form above is preferred.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project test/test_cli.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cli.jl test/test_cli.jl
git commit -m "refactor(cli)!: --variant flag, generic --data, drop g6pd-noatp subcommand"
```

---

### Task 8: Documentation — README §5 + AGENTS.md

**Files:**
- Modify: `README.md` (§5 "Run on your own data" + "Which rate law you're fitting"; §8 command line)
- Modify: `AGENTS.md` (running section, config note, CLI block)

**Interfaces:** none (prose).

- [ ] **Step 1: Rewrite README §5**

- Replace the custom-columns paragraph with the canonical-schema contract: "Your CSV must use the exact column names in the table above; a non-conforming file errors with the expected header. Column *remapping* is not supported."
- Replace the "own data" example with:
  ```julia
  using FitRateEquation
  fit_consensus_equation(:g6pd; data_csv="/path/to/my_corpus.csv")
  ```
- Rewrite the "Which rate law you're fitting" table so every row is a full, self-contained call (no `…` placeholder, no `g6pd_config()`):

  | Enzyme | Rate law | How to run |
  |---|---|---|
  | G6PD | Consensus (deployed) | `run_g6pd()` or `fit_consensus_equation(:g6pd)` |
  | G6PD | ATP-free | `fit_consensus_equation(:g6pd; variants=[:no_atp])` |
  | G6PD | Drop E·G6P·ATP dead-end | `fit_consensus_equation(:g6pd; variants=[:no_g6p_atp_deadend])` |
  | G6PD | Drop E·G6P·NADPH dead-end | `fit_consensus_equation(:g6pd; variants=[:no_g6p_nadph_deadend])` |
  | G6PD | Drop both G6P dead-ends | `fit_consensus_equation(:g6pd; variants=[:no_g6p_both_deadends])` |
  | G6PD | Pure RE (diagnostic) | `fit_consensus_equation(:g6pd; variants=[:RE_rate_eq])` |
  | PGD | Consensus (deployed) | `run_pgd()` or `fit_consensus_equation(:pgd)` |
  | PGD | Fully RE (evaluation) | `fit_consensus_equation(:pgd; variants=[:full_re])` |

- Add one line explaining the syntax: "`variants` takes a Julia list `[...]`; each name is a Symbol (leading colon). Selecting a variant automatically applies its data filter and names its output folder."

- [ ] **Step 2: Update README §8 (CLI)**

Drop the `g6pd-noatp` row; document `--variant NAME` and that `--data` works for any enzyme. Update the example to `… -- g6pd --variant no_atp --smoke`.

- [ ] **Step 3: Update AGENTS.md**

- "Running" section: replace the runner list — `run_g6pd`/`run_pgd`/`run_hk1` remain; note `fit_consensus_equation(:enzyme; variants=[…])` is the general form; remove `run_g6pd_noatp`/`run_pgd_fullre`.
- Config note: state builders are internal (`FitRateEquation.g6pd_config`), not exported; canonical columns enforced.
- CLI block: mirror the new subcommands/flags.

- [ ] **Step 4: Verify docs test (if any)**

Run: `julia --project test/test_docs.jl`
Expected: PASS (this test checks README/AGENTS invariants; update its expected strings if it pins the old names).

- [ ] **Step 5: Commit**

```bash
git add README.md AGENTS.md test/test_docs.jl
git commit -m "docs: single entry point; canonical-schema contract; CLI --variant"
```

---

### Task 9: Version bump + full local & CI verification

**Files:**
- Modify: `Project.toml` (version)

**Interfaces:** none.

- [ ] **Step 1: Bump the version**

In `Project.toml`, set `version = "0.5.0"`.

- [ ] **Step 2: Full local suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all green. Investigate any failure before proceeding — do not mark complete on a red suite.

- [ ] **Step 3: Co-import smoke (name-collision guard)**

Run: `julia --project -e 'using FitRateEquation, EnzymeRates; r = fit_consensus_equation(:g6pd; smoke=true, nprocs=1, outdir=mktempdir()); println(length(r), " cells; ok")'`
Expected: prints a nonzero cell count + `ok` (confirms `fit_consensus_equation` resolves unambiguously under co-import).

- [ ] **Step 4: Push branch and open PR; confirm CI green**

Follow the `verify-tests-local-and-ci` skill: push `feat/unified-entry-point`, open the PR, and confirm the full CI run (not just one job) is green before requesting merge. Because this is a registered package version bump, both local and CI must agree.

- [ ] **Step 5: Commit**

```bash
git add Project.toml
git commit -m "chore(release): bump version to 0.5.0 (breaking: single entry point)"
```

---

## Self-Review

**Spec coverage:** symbol-only API (Task 3, 6) ✓; un-export configs (Task 6) ✓; schema enforcement + extra columns allowed (Task 2) ✓; `:no_atp`/`:full_re` → variants with auto-filter/label (Tasks 1, 3, 4) ✓; delete `run_g6pd_noatp`/`run_pgd_fullre` (Task 4) ✓; keep `run_g6pd`/`run_pgd`/`run_hk1` (Task 4) ✓; budget/workers/outdir folded in (Task 3) ✓; row-filter conflict rule (Tasks 1, 3) ✓; CLI `--variant`, drop `g6pd-noatp` (Task 7) ✓; docs rewrite (Task 8) ✓; `0.5.0` bump + CI (Task 9) ✓; `anchor_reverse` default logic unchanged (Task 3 reuses `_default_anchor_reverse`) ✓; determinism gate preserved (Tasks 3, 5 verify `test_byte_identity`/`test_parallel_equivalence`) ✓.

**Placeholder scan:** no TBD/TODO; every code step carries real code; test transformations name exact files. The one `…` in the README table is explicitly removed by Task 8.

**Type consistency:** `variant_profile` returns `(row_filter, label)` used identically in Tasks 1/3; `_fit_consensus` keeps the exact signature of the renamed cfg-method (Task 5); the symbol method's kwargs match what `run_*` pass (Task 4) and what `cli_main` passes (Task 7).

**Open verification note:** Task 3's `read_fit_corpus` reference must be confirmed against `src/plot_support.jl`; if the reader has a different name, use the actual one (the assertion only needs the `ATP` column of `fit_corpus.csv`).

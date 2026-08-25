# G6PD Absolute-Scale Fitting Mode — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Add an orthogonal `scale=:relative|:absolute` axis to the G6PD Cha fitting path, so an uncentered loss over single-scale forward-only data can discriminate the dead-end variants and recover `kcat`.

**Architecture:** `scale` threads from `fit_consensus_equation` down to the loss. A shared per-row loss core feeds two aggregators (centered variance / uncentered sum-of-squares). `kcat` becomes a first-class coord via a scale-parameterized `cha_coords`; absolute mode is fiber-free (`C=1`) by fitting `kf=kcat` at a large release rate. Per-row enzyme concentration `Et` (new required `[G6PD] (nM)` column) enters as a linear prefactor. Reverse constants are pinned via an explicit `pins` override.

**Tech Stack:** Julia ≥ 1.11, EnzymeRates.jl (git dep), CMAEvolutionStrategy, DataFrames, CSV, Test.

**Spec:** `docs/2026-08-24-g6pd-absolute-scale-fitting-mode-design.md`

## Global Constraints

- Units: model is Molar. `[G6PD] (nM)` → M (÷1e9). Absolute-mode rate column is M·s⁻¹; `kcat` reported in s⁻¹.
- Fiber-free absolute mode: `C=1` realized by `kf=kcat`, `release_rate = CHA_ABS_RELEASE_RATE = 1e8`. Forward-only data only (`P=PGLn=0`), which is what makes the large-`release_rate` limit numerically safe.
- `scale=:relative` must reproduce the current pipeline **bit-for-bit** (regression lock in `test_byte_identity.jl`).
- `scale=:absolute` is guarded to `:G6PD`; any other enzyme errors explicitly.
- `kcat` coord bound: `[10, 1000]` s⁻¹ (log10 `[1, 3]`). Expected 150–250 s⁻¹ is a report-time validation verdict, never a clamp.
- Naming: `Kd_6PGLn` is the 6PGL (K_PGLn) constant; `Km_NADPH_rev` is the reverse NADPH constant; there is no `Ki_6PGLn` coord.
- Pin values in the escalation ladder are user-supplied (from a prior relative fit); the ladder is run at `mode1`.
- Section headers in scripts: ASCII `#` banners, 90 cols, centered title (per user CLAUDE.md).
- Commits: Conventional Commits `type(scope): summary`; end body with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

---

## Task 0: Worktree environment + green baseline

**Files:** none (environment only).

- [x] **Step 1: Instantiate the worktree project**

Run:
```bash
cd /home/james/projects/FitRateEquation.jl/.worktrees/absolute-scale-mode
julia --project -e 'using Pkg; Pkg.instantiate()'
```
Expected: resolves EnzymeRates from the `[sources]` git URL and precompiles. If it fails to resolve EnzymeRates (offline / depth-sensitive source), copy the main checkout's Manifest and retry:
```bash
cp /home/james/projects/FitRateEquation.jl/Manifest.toml .   # gitignored; local only
julia --project -e 'using Pkg; Pkg.instantiate()'
```

- [x] **Step 2: Run the full baseline suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (clean baseline). If any test fails on a fresh `main` baseline, STOP and report — do not build on a red baseline.

- [x] **Step 3: No commit** (environment only; Manifest is gitignored).

---

## Task 1: `Et` column — Dataset field, config, loader, validation

Adds the per-row enzyme concentration carried by a new `[G6PD] (nM)` corpus column. Relative mode ignores it; absolute mode requires it (validated in Task 6).

**Files:**
- Modify: `src/core/data.jl` (Dataset struct + 4-arg back-compat constructor; `read_corpus`; `dataset_from_corpus`)
- Modify: `src/configs/G6PD.jl` (`enzyme_conc_col`, `enzyme_conc_unit`)
- Test: `test/test_data.jl` (append), and a small fixture CSV under `test/`

**Interfaces:**
- Produces: `Dataset` gains field `Et::Vector{Float64}` (5th, M units; `NaN` where absent). Back-compat `Dataset(concs, rate, group, keq)` fills `Et` with `NaN`. `g6pd_config()` gains `enzyme_conc_col::String`, `enzyme_conc_unit::Symbol`.

- [x] **Step 1: Write the failing test**

Create a tiny fixture `test/fixtures/g6pd_abs_mini.csv` (forward-only, 4 rows, includes `[G6PD] (nM)`):
```csv
[NADP] (uM),[G6P] (uM),[NADPH] (uM),[PGLn] (uM),[ATP] (uM),Rate_V,Specific_Activity,Article,Fig,X_axis_label,Experiment_date,pH,Temperature,Apparent_Keq,[G6PD] (nM)
50,200,0,0,0,1.0e-6,1,Mbata2026,1a,NADP,26-09-01,7.6,25,13.655,5.0
100,200,0,0,0,1.5e-6,1,Mbata2026,1a,NADP,26-09-01,7.6,25,13.655,5.0
200,200,0,0,0,1.8e-6,1,Mbata2026,1a,NADP,26-09-01,7.6,25,13.655,5.0
400,200,0,0,0,2.0e-6,1,Mbata2026,1a,NADP,26-09-01,7.6,25,13.655,5.0
```

Append to `test/test_data.jl`:
```julia
@testset "Et column: [G6PD] (nM) -> M, present in Dataset" begin
    cfg = FitRateEquation.g6pd_config(
        data_csv = joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    @test cfg.enzyme_conc_col == "[G6PD] (nM)"
    @test cfg.enzyme_conc_unit == :nM
    d = FitRateEquation.load_dataset(cfg)
    @test length(d.Et) == FitRateEquation.nrows(d)
    @test all(d.Et .== 5.0e-9)            # 5 nM -> 5e-9 M
end

@testset "Et back-compat: 4-arg Dataset fills Et with NaN" begin
    d0 = FitRateEquation.load_dataset(FitRateEquation.g6pd_config())
    d1 = Dataset(d0.concs, d0.rate, d0.group, d0.keq)   # 4-arg positional (existing idiom)
    @test length(d1.Et) == FitRateEquation.nrows(d1)
    @test all(isnan, d1.Et)
end
```

- [x] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_data.jl`
Expected: FAIL — `enzyme_conc_col` not in config / `Dataset` has no field `Et`.

- [x] **Step 3: Implement**

In `src/core/data.jl`, extend the struct and add a back-compat constructor:
```julia
struct Dataset{T<:NamedTuple}
    concs::Vector{T}
    rate::Vector{Float64}
    group::Vector{String}
    keq::Vector{Float64}
    Et::Vector{Float64}        # per-row enzyme concentration (M); NaN where absent
end
# Back-compat: existing 4-arg positional callers (tests, direct constructions) get NaN Et.
Dataset(concs, rate, group, keq) = Dataset(concs, rate, group, keq, fill(NaN, length(rate)))
```

In `read_corpus`, after the `Apparent_Keq` line and before `filter!`, add the opportunistic Et column:
```julia
    if hasproperty(cfg, :enzyme_conc_col) && cfg.enzyme_conc_col in names(raw)
        ev = _to_float.(raw[!, cfg.enzyme_conc_col], NaN)
        df.Et = cfg.enzyme_conc_unit === :nM ? ev ./ 1e9 :
                cfg.enzyme_conc_unit === :uM ? ev ./ 1e6 : ev
    else
        df.Et = fill(NaN, nrow(raw))
    end
```
(Place the `df.Et` assignment before `filter!` so the drop applies to Et too.)

In `dataset_from_corpus`, pass Et through:
```julia
    Dataset(concs, Vector{Float64}(df.Rate), Vector{String}(df.source),
            Vector{Float64}(df.Apparent_Keq), Vector{Float64}(df.Et))
```

In `src/configs/G6PD.jl`, add the two fields to the returned NamedTuple:
```julia
        enzyme_conc_col = "[G6PD] (nM)",
        enzyme_conc_unit = :nM,
```

- [x] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_data.jl`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/core/data.jl src/configs/G6PD.jl test/test_data.jl test/fixtures/g6pd_abs_mini.csv
git commit -m "feat(g6pd): carry per-row Et from [G6PD] (nM) column

Adds a fifth Dataset field Et (Molar), populated by the loader from an
optional [G6PD] (nM) column, with a 4-arg back-compat constructor so
existing positional constructions still work. G6PD config gains
enzyme_conc_col/enzyme_conc_unit.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Refactor the loss into a shared core (relative unchanged)

Extracts the per-row arithmetic so centered and (Task 3) absolute aggregators share it with zero drift. Relative output must be byte-identical.

**Files:**
- Modify: `src/cha_fit.jl` (`cha_centered_logratio_loss` → shared `_cha_row_logratios!` + thin centered wrapper)
- Test: `test/test_byte_identity.jl` (already locks relative output) + `test/test_cha_fit.jl`

**Interfaces:**
- Produces: `ChaFit._cha_row_logratios!(logratio, enzyme, mech, d, coords; keq, kf, Et, release_rate, release_eq, kr, variant) -> (penalty::Float64, groups::Vector{Vector{Int}})` fills `logratio` (per-row `log(pred)-log(obs)`, `NaN` on sign/finite penalty) and returns accumulated penalty + per-group row-index sets. `cha_centered_logratio_loss` keeps its exact signature and output.

- [x] **Step 1: Write the failing test**

Append to `test/test_cha_fit.jl`:
```julia
@testset "shared loss core: centered wrapper == direct core aggregation" begin
    using Statistics: median
    d = load_dataset(g6pd_config()); keq = median(d.keq)
    m = FitRateEquation.v2_mechanism()
    coords = Dict(s => getfield(cha_macro_readoffs_G6PD(m, -3 .+ 2 .* rand(length(free_params(m))); keq=keq), s)
                  for s in cha_coords(:G6PD))
    L = ChaFit.cha_centered_logratio_loss(:G6PD, m, d, coords; keq=keq)
    lr = fill(NaN, FitRateEquation.nrows(d))
    pen, groups = ChaFit._cha_row_logratios!(lr, :G6PD, m, d, coords; keq=keq)
    total = pen
    for idx in groups
        vals = filter(isfinite, lr[idx])
        isempty(vals) && continue
        μ = sum(vals)/length(vals)
        total += sum(x -> (x-μ)^2, vals)
    end
    @test L ≈ total / FitRateEquation.nrows(d)
end
```

- [x] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cha_fit.jl`
Expected: FAIL — `_cha_row_logratios!` undefined.

- [x] **Step 3: Implement**

In `src/cha_fit.jl`, extract the group/keq/rate loop from `cha_centered_logratio_loss` into `_cha_row_logratios!` (preserving the group iteration order and the `_SIGN_PENALTY`/finite handling exactly). The core does NOT do mean-centering — it only fills `logratio` and returns `(penalty, groups)`. Then rewrite `cha_centered_logratio_loss` to call the core and aggregate per-group variance in the SAME fold order as today (accumulate `penalty` first, then add group variances in `groups` order — preserving float associativity):
```julia
function _cha_row_logratios!(logratio, enzyme::Symbol, mech, d::Dataset, coords::AbstractDict;
        keq::Union{Nothing,Real}=nothing, kf::Real=1.0, Et::Real=1.0,
        release_rate::Real=_default_release_rate(enzyme),
        release_eq::Real=_default_release_eq(enzyme, coords),
        kr::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy)
    cha_rate_enz = enzyme === :G6PD ? ChaLaws.cha_rate_G6PD :
                   enzyme === :PGD  ? (variant === :full_re ? ChaLaws.cha_rate_PGD_fullRE :
                                                              ChaLaws.cha_rate_PGD) :
                   enzyme === :HK1  ? ChaLawsHK1.cha_rate_HK1 :
                   error("_cha_row_logratios!: unknown enzyme $enzyme")
    penalty = 0.0
    groups = [findall(==(g), d.group) for g in unique(d.group)]
    for idx in groups
        keq_g = keq === nothing ? (only(unique(d.keq[idx]))) : keq
        m = cha_macro_tuple(enzyme, coords; keq=keq_g, kf=kf, Et=Et,
                            release_rate=release_rate, release_eq=release_eq, kr=kr, variant=variant)
        for i in idx
            v = cha_rate_enz(m; _cha_row_kwargs(enzyme, d.concs[i])...)
            o = d.rate[i]
            if !isfinite(v) || v == 0 || sign(v) != sign(o)
                penalty += _SIGN_PENALTY
                logratio[i] = NaN
            else
                logratio[i] = log(abs(v)) - log(abs(o))
            end
        end
    end
    (penalty, groups)
end

function cha_centered_logratio_loss(enzyme::Symbol, mech, d::Dataset, coords::AbstractDict;
        keq::Union{Nothing,Real}=nothing, kf::Real=1.0, Et::Real=1.0,
        release_rate::Real=_default_release_rate(enzyme),
        release_eq::Real=_default_release_eq(enzyme, coords),
        kr::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy)
    n = nrows(d)
    logratio = fill(NaN, n)
    penalty, groups = _cha_row_logratios!(logratio, enzyme, mech, d, coords; keq=keq, kf=kf,
        Et=Et, release_rate=release_rate, release_eq=release_eq, kr=kr, variant=variant)
    total = penalty
    for idx in groups
        vals = filter(isfinite, logratio[idx])
        isempty(vals) && continue
        μ = sum(vals) / length(vals)
        total += sum(x -> (x - μ)^2, vals)
    end
    total / n
end
```
Note: the original stashed group variances to preserve `penalty`-first ordering; the rewrite adds `penalty` first then the group variances in `groups` order — same order. If `test_byte_identity.jl` shows any last-bit drift, restore the original two-phase `group_variances` stash inside the wrapper.

- [x] **Step 4: Run tests to verify they pass (incl. bit-identity)**

Run: `julia --project test/test_cha_fit.jl && julia --project test/test_byte_identity.jl`
Expected: PASS both. `test_byte_identity.jl` confirms relative output unchanged.

- [x] **Step 5: Commit**

```bash
git add src/cha_fit.jl test/test_cha_fit.jl
git commit -m "refactor(cha): extract shared per-row loss core

Splits cha_centered_logratio_loss into a shared _cha_row_logratios! core
plus a thin centered aggregator, so the absolute-mode aggregator can reuse
the exact per-row arithmetic. Relative output is byte-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Absolute (uncentered) aggregator + per-row `Et`

Adds the uncentered aggregator and the linear per-row `Et` prefactor, and wires the `scale`-based dispatch into `_cha_loss_with_pins`.

**Files:**
- Modify: `src/cha_fit.jl` (`cha_absolute_logratio_loss`; per-row `Et` in the core; `_cha_loss_with_pins` gains `scale`)
- Test: `test/test_cha_absolute.jl` (new; add to `runtests.jl`)

**Interfaces:**
- Consumes: `_cha_row_logratios!` (Task 2).
- Produces: `ChaFit.cha_absolute_logratio_loss(enzyme, mech, d, coords; keq, kf, release_rate, release_eq, kr, variant) -> Float64` (uncentered `sum(logratio²)/n`, per-row `Et` from `d.Et`). `_cha_loss_with_pins(...; scale::Symbol=:relative)`.

- [ ] **Step 1: Write the failing test**

Create `test/test_cha_absolute.jl`:
```julia
using FitRateEquation
using FitRateEquation: g6pd_config
using FitRateEquation.ChaFit
using FitRateEquation.ChaLaws
using Test

@testset "absolute loss: Et is a linear prefactor (log-shift)" begin
    cfg = g6pd_config(data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    m = FitRateEquation.v2_mechanism()
    # planted coords (any physically valid point); kcat as :kf via scale-absolute coords (Task 4)
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    n = FitRateEquation.nrows(d)
    lr = fill(NaN, n)
    ChaFit._cha_row_logratios!(lr, :G6PD, m, d, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE, Et=1.0)   # Et=1 baseline
    # Doubling every row's Et must shift each finite log-ratio by exactly log(2).
    d2 = Dataset(d.concs, d.rate, d.group, d.keq, 2 .* d.Et)
    lr2 = fill(NaN, n)
    ChaFit._cha_row_logratios!(lr2, :G6PD, m, d2, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
    for i in 1:n
        (isfinite(lr[i]) && isfinite(lr2[i])) || continue
        @test lr2[i] - lr[i] ≈ log(2) atol=1e-10
    end
end

@testset "absolute loss: uncentered sum of squares" begin
    cfg = g6pd_config(data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    m = FitRateEquation.v2_mechanism()
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    L = ChaFit.cha_absolute_logratio_loss(:G6PD, m, d, coords; keq=13.655, kf=178.0)
    lr = fill(NaN, FitRateEquation.nrows(d))
    pen, _ = ChaFit._cha_row_logratios!(lr, :G6PD, m, d, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
    expect = (pen + sum(x -> x^2, filter(isfinite, lr))) / FitRateEquation.nrows(d)
    @test L ≈ expect
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cha_absolute.jl`
Expected: FAIL — `CHA_ABS_RELEASE_RATE` / `cha_absolute_logratio_loss` undefined; `Et` kwarg not applied per-row.

- [ ] **Step 3: Implement**

In `src/cha_fit.jl`:
1. Add the constant near `CHA_DEPLOY_RELEASE_RATE`:
```julia
# Absolute-mode release rate: large enough that the SS-release fiber factor C=1+kf/koffQ ≈ 1
# across the kcat bound (kf ≤ 1000), so kcat ≡ kf. Forward-only data (P=0) makes this exact
# and numerically safe (konQ appears only as konQ/koffQ = 1/Km_NADPH_rev).
const CHA_ABS_RELEASE_RATE = 1.0e8
```
2. Change the per-row prediction in `_cha_row_logratios!` to apply `d.Et` when finite. Replace the `v = cha_rate_enz(...)` line body with:
```julia
            vunit = cha_rate_enz(m; _cha_row_kwargs(enzyme, d.concs[i])...)
            eti = isnan(d.Et[i]) ? 1.0 : d.Et[i]
            v = eti * vunit
```
   (In relative mode `d.Et` is `NaN` ⇒ `eti=1.0` ⇒ prediction unchanged; the tuple keeps `Et=1.0`, so this is a strict no-op for the centered path — `test_byte_identity.jl` still passes.)
3. Add the aggregator:
```julia
function cha_absolute_logratio_loss(enzyme::Symbol, mech, d::Dataset, coords::AbstractDict;
        keq::Union{Nothing,Real}=nothing, kf::Real=1.0,
        release_rate::Real=CHA_ABS_RELEASE_RATE,
        release_eq::Real=_default_release_eq(enzyme, coords),
        kr::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy)
    n = nrows(d)
    logratio = fill(NaN, n)
    penalty, _ = _cha_row_logratios!(logratio, enzyme, mech, d, coords; keq=keq, kf=kf,
        Et=1.0, release_rate=release_rate, release_eq=release_eq, kr=kr, variant=variant)
    total = penalty
    for x in logratio
        isfinite(x) && (total += x * x)
    end
    total / n
end
```

- [ ] **Step 4: Run tests to verify pass (incl. bit-identity no-op)**

Run: `julia --project test/test_cha_absolute.jl && julia --project test/test_byte_identity.jl`
Expected: PASS both.

- [ ] **Step 5: Add to runtests and commit**

Add `include("test_cha_absolute.jl")` to `test/runtests.jl` (after `test_cha_fit.jl`).
```bash
git add src/cha_fit.jl test/test_cha_absolute.jl test/runtests.jl
git commit -m "feat(cha): uncentered absolute loss + per-row Et prefactor

Adds cha_absolute_logratio_loss (uncentered sum-of-squares) and applies
d.Et as a linear prefactor in the shared core. Relative mode is a strict
no-op (Et is NaN -> 1.0).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Scale-parameterized `cha_coords` + `:kcat` + loss dispatch

Makes `:kcat` a first-class coord in absolute mode and realizes fiber-free `C=1` by mapping `:kcat → kf` with `release_rate=CHA_ABS_RELEASE_RATE`.

**Files:**
- Modify: `src/cha_fit.jl` (`cha_coords`, `cha_coord_bounds`, `_cha_loss_with_pins`, `cha_fit_candidate` gain `scale`)
- Test: `test/test_cha_absolute.jl` (append)

**Interfaces:**
- Produces: `cha_coords(enzyme, variant=:_deploy; scale::Symbol=:relative)` appends `:kcat` for `(:G6PD, :absolute)`. `cha_coord_bounds(...; scale)` gives `:kcat` bound `[1,3]` (log10). `_cha_loss_with_pins(...; scale)` and `cha_fit_candidate(...; scale)`. In absolute mode the loss pops `:kcat` from the coord dict and passes `kf=coords[:kcat]`, `release_rate=CHA_ABS_RELEASE_RATE`.

- [ ] **Step 1: Write the failing test**

Append to `test/test_cha_absolute.jl`:
```julia
@testset "cha_coords appends :kcat only for (:G6PD, :absolute)" begin
    @test !(:kcat in cha_coords(:G6PD))
    cs = cha_coords(:G6PD, :_deploy; scale=:absolute)
    @test :kcat in cs
    @test cs[end] == :kcat
    lo, hi = ChaFit.cha_coord_bounds(:G6PD, :_deploy; scale=:absolute)
    @test length(lo) == length(cs)
    @test (lo[end], hi[end]) == (1.0, 3.0)   # kcat in [10, 1000] s^-1
end

@testset "fiber-free C=1: kcat == fitted kf, Km == alpha*Kd" begin
    # At release_rate = CHA_ABS_RELEASE_RATE the apparent Km loses its fiber factor.
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.3,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    km = ChaFit.cha_apparent_km(:G6PD, coords, :Km_G6P;
                                kf=178.0, release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
    @test km ≈ coords[:alpha]*coords[:Kd_G6P] rtol=1e-5   # C -> 1
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cha_absolute.jl`
Expected: FAIL — `cha_coords` has no `scale` kwarg.

- [ ] **Step 3: Implement**

In `cha_coords`, add `scale::Symbol=:relative` and, at the end of the `:G6PD` branch, append `:kcat` when `scale === :absolute`:
```julia
function cha_coords(enzyme::Symbol, variant::Symbol=:_deploy; scale::Symbol=:relative)
    base = <existing body returning the binding-constant coords>
    if enzyme === :G6PD && scale === :absolute
        return vcat(base, :kcat)
    end
    return base
end
```
(Refactor the existing return-chain into `base` first, then the append.)

In `cha_coord_bounds`, add `scale` and thread it to `cha_coords`; add the `:kcat` case:
```julia
function cha_coord_bounds(enzyme::Symbol, variant::Symbol=:_deploy; scale::Symbol=:relative)
    coords = cha_coords(enzyme, variant; scale=scale)
    lo = Float64[]; hi = Float64[]
    for s in coords
        if s === :alpha
            push!(lo, -2.0); push!(hi, 2.0)
        elseif s === :split_ratio
            push!(lo, log10(2.0)); push!(hi, 3.0)
        elseif s === :kcat
            push!(lo, 1.0); push!(hi, 3.0)      # 10 .. 1000 s^-1
        else
            push!(lo, -9.0); push!(hi, 0.0)
        end
    end
    lo, hi
end
```

In `_cha_loss_with_pins`, add `scale::Symbol=:relative`; when absolute, pop `:kcat` and dispatch:
```julia
function _cha_loss_with_pins(enzyme, mech, d, u, coords_syms, pins, anchors;
        keq::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy, scale::Symbol=:relative)
    if !isempty(pins)
        u = collect(u)
        for (idx, k) in enumerate(coords_syms)
            haskey(pins, k) && (u[idx] = pins[k])
        end
    end
    coords_dict = Dict(coords_syms .=> 10 .^ u)
    if scale === :absolute
        kcat = pop!(coords_dict, :kcat)
        L = cha_absolute_logratio_loss(enzyme, mech, d, coords_dict;
                keq=keq, kf=kcat, release_rate=CHA_ABS_RELEASE_RATE, variant=variant)
    else
        L = cha_centered_logratio_loss(enzyme, mech, d, coords_dict; keq=keq, variant=variant)
    end
    L += _cha_anchor_penalty(enzyme, coords_dict, anchors; variant=variant)
    L
end
```

In `cha_fit_candidate`, add `scale::Symbol=:relative`, thread it to `cha_coords`/`cha_coord_bounds`/`_cha_loss_with_pins`:
```julia
function cha_fit_candidate(enzyme::Symbol, mech, d::Dataset; n_restarts::Int=8,
        maxiter::Int=1_000_000, maxtime::Real=20.0, seed::Int=1,
        keq::Union{Nothing,Real}=nothing, pins::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
        anchors=nothing, variant::Symbol=:_deploy, scale::Symbol=:relative)
    coords_syms = cha_coords(enzyme, variant; scale=scale)
    lo, hi = cha_coord_bounds(enzyme, variant; scale=scale)
    objective = u -> _cha_loss_with_pins(enzyme, mech, d, u, coords_syms, pins, anchors;
                                         keq=keq, variant=variant, scale=scale)
    <unchanged restart loop / pin-overwrite / return>
end
```

- [ ] **Step 4: Run tests to verify pass (incl. bit-identity)**

Run: `julia --project test/test_cha_absolute.jl && julia --project test/test_byte_identity.jl && julia --project test/test_cha_fit.jl`
Expected: PASS all (relative default path unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/cha_fit.jl test/test_cha_absolute.jl
git commit -m "feat(cha): scale-parameterized cha_coords with :kcat coord

Absolute mode appends :kcat (bound 10..1000 s^-1) to the G6PD coords and
maps it to kf at CHA_ABS_RELEASE_RATE (fiber-free C=1). _cha_loss_with_pins
and cha_fit_candidate gain a scale kwarg; relative default is unchanged.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Explicit `pins` value-override merge

Lets the caller pin any coord to an explicit log10 value on top of the mode-derived pins, validated by the existing guard. (`cha_fit_candidate` already accepts `pins`; this task adds the guard-checked *merge* at the resolve layer and the entry-point plumbing spot.)

**Files:**
- Modify: `src/cha_fit.jl` (`resolve_cha_pins` gains an `extra::Dict` merge with guard)
- Test: `test/test_cha_fit.jl` (append)

**Interfaces:**
- Produces: `resolve_cha_pins(enzyme, variant, mode; anchor_reverse=true, extra::Dict{Symbol,Float64}=Dict(), scale::Symbol=:relative)` merges `extra` over the mode pins after asserting each key is a coord for `(enzyme,variant,scale)`.

- [ ] **Step 1: Write the failing test**

Append to `test/test_cha_fit.jl`:
```julia
@testset "resolve_cha_pins merges explicit extra pins (guarded)" begin
    p = ChaFit.resolve_cha_pins(:G6PD, :_deploy, :mode1;
            extra=Dict(:Kd_6PGLn=>log10(2.1e-4)), scale=:absolute)
    @test p[:Kd_6PGLn] == log10(2.1e-4)
    # bogus coord errors
    @test_throws ErrorException ChaFit.resolve_cha_pins(:G6PD, :_deploy, :mode1;
            extra=Dict(:NotACoord=>0.0), scale=:absolute)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cha_fit.jl`
Expected: FAIL — `resolve_cha_pins` has no `extra`/`scale` kwargs.

- [ ] **Step 3: Implement**

Extend `resolve_cha_pins` signature and, before returning `pins`, merge `extra` with the guard (thread `scale` into the coord checks, since `:kcat` is only a coord in absolute mode):
```julia
function resolve_cha_pins(enzyme::Symbol, variant::Symbol, mode::Symbol;
        anchor_reverse::Bool=true, extra::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
        scale::Symbol=:relative)
    <existing body, but pass scale to _assert_pin_is_coord via cha_coords>
    for (k, v) in extra
        _assert_pin_is_coord(enzyme, k, variant; scale=scale)
        pins[k] = v
    end
    pins
end
```
Update `_assert_pin_is_coord` to accept `scale`:
```julia
function _assert_pin_is_coord(enzyme::Symbol, name::Symbol, variant::Symbol=:_deploy;
                              scale::Symbol=:relative)
    name in cha_coords(enzyme, variant; scale=scale) && return nothing
    error(<existing message>)
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_cha_fit.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cha_fit.jl test/test_cha_fit.jl
git commit -m "feat(cha): explicit guarded pins override in resolve_cha_pins

Adds extra::Dict merge (with _assert_pin_is_coord guard) so callers pin
coords to data-determined values from a prior relative fit. scale is
threaded so :kcat is recognized as a coord only in absolute mode.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Thread `scale` + `pins` through `run.jl`; guard non-G6PD; assert `[G6PD] (nM)`

Wires the new axis end-to-end and enforces the absolute-mode data requirement.

**Files:**
- Modify: `src/run.jl` (`fit_consensus_equation`, `_fit_consensus`, task build/reduce plumbing)
- Test: `test/test_run_fit.jl` (append)

**Interfaces:**
- Consumes: `cha_fit_candidate(...; scale)`, `resolve_cha_pins(...; extra, scale)`, `d.Et`.
- Produces: `fit_consensus_equation(enzyme; …, scale::Symbol=:relative, pins::Dict{Symbol,Float64}=Dict())`. `_fit_consensus(cfg; …, scale, pins)`. Absolute + non-G6PD errors; absolute with a corpus lacking finite `Et` errors naming `[G6PD] (nM)`.

- [ ] **Step 1: Write the failing test**

Append to `test/test_run_fit.jl`:
```julia
@testset "absolute mode: guards" begin
    # non-G6PD is rejected
    @test_throws ErrorException fit_consensus_equation(:pgd; scale=:absolute, smoke=true,
        outdir=mktempdir())
    # G6PD absolute on a corpus WITHOUT [G6PD] (nM) errors naming the column
    err = try
        fit_consensus_equation(:g6pd; scale=:absolute, smoke=true, outdir=mktempdir())
        nothing
    catch e; e end
    @test err isa ErrorException
    @test occursin("[G6PD] (nM)", sprint(showerror, err))
end

@testset "absolute mode: smoke fit on mini forward corpus runs" begin
    out = mktempdir()
    res = fit_consensus_equation(:g6pd; scale=:absolute, smoke=true, outdir=out,
        data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"),
        pins=Dict(:Kd_6PGLn=>log10(2.1e-4), :Km_NADPH_rev=>log10(3.9e-6)))
    @test !isempty(res)
    @test isfile(joinpath(out, "macro_constants.csv"))
end
```
(If the mini fixture is too small for the CV fold count, this test asserts only that the in-sample fit + outputs run; CV specifics are Task 7.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_run_fit.jl`
Expected: FAIL — `fit_consensus_equation` has no `scale`/`pins`; no guard.

- [ ] **Step 3: Implement**

In `fit_consensus_equation`, add `scale::Symbol=:relative` and `pins::Dict{Symbol,Float64}=Dict{Symbol,Float64}()`; after `enz = _canonical_enzyme(...)` add the guard:
```julia
    scale === :absolute && enz !== :G6PD &&
        error("absolute scale is not yet wired for $enz (G6PD only).")
```
Pass `scale` and `pins` into `_fit_consensus`.

In `_fit_consensus`, add `scale` and `pins` kwargs. After building `d`, enforce the Et requirement in absolute mode:
```julia
    if scale === :absolute
        (hasproperty(cfg, :enzyme_conc_col) && any(isfinite, d.Et)) ||
            error("absolute scale requires a \"[G6PD] (nM)\" column with finite values; " *
                  "none found in $(cfg.data_csv).")
    end
```
Thread `scale` and `pins` into the task machinery: `_build_tasks(...; scale)` and each `resolve_cha_pins(...; extra=pins, scale=scale)` call site (there are three: the tasks builder, the reduce path, and any inline single-fit path — grep `resolve_cha_pins` in `run.jl`), and pass `scale` into every `cha_fit_candidate` call inside `_run_fit_task`/`_reduce_cells`. Add `scale` to the `meta` NamedTuple and to the `write_outputs(...; scale=scale)` call (consumed in Task 9).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_run_fit.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/run.jl test/test_run_fit.jl
git commit -m "feat(run): thread scale + pins end-to-end; guard absolute mode

fit_consensus_equation gains scale and pins. Absolute mode is guarded to
G6PD and requires a finite [G6PD] (nM) column. scale/pins flow through the
task build/reduce path into cha_fit_candidate and resolve_cha_pins.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Leave-one-group-out CV (by `Fig`) for absolute mode

Adds a single-article CV fold helper and selects it when `scale=:absolute`.

**Files:**
- Modify: `src/cv.jl` (group-column fold helper paralleling `_article_folds`)
- Modify: `src/run.jl` (use group folds when absolute)
- Test: `test/test_cv.jl` (append)

**Interfaces:**
- Produces: `_group_folds(d::Dataset) -> Vector{NamedTuple{(:train,:test)}}` — one held-out fold per unique `d.group`.

- [ ] **Step 1: Write the failing test**

Append to `test/test_cv.jl`:
```julia
@testset "_group_folds: one fold per unique group, train/test partition" begin
    cfg = FitRateEquation.g6pd_config(
        data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    folds = FitRateEquation._group_folds(d)
    @test length(folds) == length(unique(d.group))
    for f in folds
        @test sort(vcat(f.train, f.test)) == collect(1:FitRateEquation.nrows(d))
        @test isempty(intersect(f.train, f.test))
    end
end
```
(The mini fixture is a single group `Mbata2026|1a`, so `folds` has length 1 with an empty complement — assert length and partition still hold; a multi-group fixture can be added if richer coverage is wanted.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cv.jl`
Expected: FAIL — `_group_folds` undefined.

- [ ] **Step 3: Implement**

In `src/cv.jl`, mirror `_article_folds` but fold on `d.group`:
```julia
function _group_folds(d::Dataset)
    groups = unique(d.group)
    all_idx = collect(1:length(d.group))
    [ (train = [i for i in all_idx if d.group[i] != g],
       test  = [i for i in all_idx if d.group[i] == g]) for g in groups ]
end
```
In `run.jl`, where LOO-article folds are chosen (`_cha_loocv` / `_article_folds` call site), select `_group_folds(d)` when `scale === :absolute`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_cv.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cv.jl src/run.jl test/test_cv.jl
git commit -m "feat(cv): leave-one-group-out folds for absolute mode

Single-article absolute-mode data makes leave-one-article-out degenerate;
_group_folds folds by Fig group instead. run.jl selects it under scale=:absolute.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: `:kcat`-aware classify / invert / deploy

Makes the readoff/classification/deploy path handle the extra `:kcat` coord (identifiable scale) instead of assuming the `kf=1` gauge.

**Files:**
- Modify: `src/cha_classify.jl` (accept/label `:kcat`), `src/cha_invert.jl` / `src/cha_deploy.jl` (use fitted `kcat` as `kf`)
- Test: `test/test_cha_classify.jl` (append)

**Interfaces:**
- Consumes: absolute-mode fit `coords` containing `:kcat`.
- Produces: `classify_cha` returns a row for `:kcat` (labeled data-identified in absolute mode); `cha_deploy_micro` uses `kf = coords[:kcat]` when present (else the gauge `kf=1`).

- [ ] **Step 1: Write the failing test**

Append to `test/test_cha_classify.jl`:
```julia
@testset "classify handles :kcat coord in absolute mode" begin
    using Statistics: median
    d = FitRateEquation.load_dataset(FitRateEquation.g6pd_config())
    m = FitRateEquation.v2_mechanism()
    coords = Dict{Symbol,Float64}(s => 1e-4 for s in cha_coords(:G6PD))
    coords[:alpha] = 1.0; coords[:kcat] = 178.0
    classed = ChaClassify.classify_cha(:G6PD, m, d, coords, Dict{Symbol,Float64}(), nothing;
                                       scale=:absolute)
    @test any(c -> c.name == :kcat, classed)
end
```
(Match the actual `classify_cha` positional/kwarg shape used elsewhere in this test file; add a `scale` kwarg defaulting to `:relative`.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cha_classify.jl`
Expected: FAIL — `classify_cha` lacks `scale`/`:kcat` handling.

- [ ] **Step 3: Implement**

- `classify_cha`: add `scale::Symbol=:relative`; when absolute, include `:kcat` in the classified coords (report its fitted value; identifiability via the existing Hessian machinery over the absolute loss). Thread `scale` to any internal `cha_coords`/`cha_identifiable_functions` calls.
- `cha_deploy_micro` / `cha_invert`: read `kf = get(coords, :kcat, 1.0)` so the deployed micro map carries the fitted absolute scale (byproduct) instead of the `kf=1` gauge; everything else unchanged.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_cha_classify.jl && julia --project test/test_cha_deploy.jl`
Expected: PASS both (relative path unchanged: `get(coords, :kcat, 1.0)` → 1.0).

- [ ] **Step 5: Commit**

```bash
git add src/cha_classify.jl src/cha_invert.jl src/cha_deploy.jl test/test_cha_classify.jl
git commit -m "feat(cha): kcat-aware classify/invert/deploy

Absolute-mode fits carry a :kcat coord; classify reports it and deploy uses
it as kf. Relative mode is unchanged (kf defaults to the gauge 1.0).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 9: Outputs — report/provenance scale fields + kcat verdict; absolute plotting

Surfaces the new axis in the artifacts and plots absolute rates.

**Files:**
- Modify: `src/run.jl` (`write_outputs`, report/provenance writers)
- Modify: `ext/FitRateEquationMakieExt.jl` (absolute predicted-vs-measured)
- Test: `test/test_outputs.jl` (append)

**Interfaces:**
- Consumes: `meta.scale`, absolute-mode `coords[:kcat]`.
- Produces: `report.md` records `scale`, fitted `kcat`, and a `150–250 s⁻¹` in-band verdict; `provenance.toml` records `scale`.

- [ ] **Step 1: Write the failing test**

Append to `test/test_outputs.jl`:
```julia
@testset "absolute outputs record scale + kcat verdict" begin
    out = mktempdir()
    fit_consensus_equation(:g6pd; scale=:absolute, smoke=true, outdir=out,
        data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"),
        pins=Dict(:Kd_6PGLn=>log10(2.1e-4), :Km_NADPH_rev=>log10(3.9e-6)))
    prov = read(joinpath(out, "provenance.toml"), String)
    @test occursin("scale", prov) && occursin("absolute", prov)
    report = read(joinpath(out, "report.md"), String)
    @test occursin("kcat", lowercase(report))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_outputs.jl`
Expected: FAIL — provenance/report lack `scale`/`kcat`.

- [ ] **Step 3: Implement**

- Thread `scale` into `write_outputs` and add a `scale = "..."` line to the provenance TOML writer.
- In the report writer, when `meta.scale === :absolute`, add a line with the fitted `kcat` and a verdict: `in-band (150–250 s⁻¹)` vs `OUT OF BAND`.
- In `ext/FitRateEquationMakieExt.jl`, when the run is absolute, plot predicted-vs-measured without per-figure recentering (multiply prediction by row `Et`). Follow the existing per-figure plotting entry (`plot_consensus_fit`), branching on the run's recorded scale.

- [ ] **Step 4: Run tests to verify pass**

Run: `julia --project test/test_outputs.jl && julia --project test/test_plot_render.jl`
Expected: PASS both.

- [ ] **Step 5: Commit**

```bash
git add src/run.jl ext/FitRateEquationMakieExt.jl test/test_outputs.jl
git commit -m "feat(run): record scale + kcat verdict; absolute plotting

report.md/provenance.toml record scale and the fitted kcat with a
150-250 s^-1 in-band verdict; the Makie extension plots absolute
predicted-vs-measured rates when scale=:absolute.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 10: CLI `--scale`

**Files:**
- Modify: `src/cli.jl`
- Test: `test/test_cli.jl` (append)

**Interfaces:**
- Produces: `--scale relative|absolute` parsed into the `fit_consensus_equation` call.

- [ ] **Step 1: Write the failing test**

Append to `test/test_cli.jl` (mirror the file's existing arg-parse test idiom):
```julia
@testset "CLI parses --scale" begin
    opts = FitRateEquation.CLI._parse_args(["g6pd", "--scale", "absolute"])  # match actual parser name
    @test opts.scale == :absolute
end
```
(Use the real parser function/return shape in `cli.jl`; if the CLI parses directly into a call, assert via a dry-run flag instead.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_cli.jl`
Expected: FAIL — `--scale` unknown.

- [ ] **Step 3: Implement**

Add `--scale` handling to `cli_main`/the arg parser, defaulting to `:relative`, forwarded to `fit_consensus_equation(...; scale=...)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project test/test_cli.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cli.jl test/test_cli.jl
git commit -m "feat(cli): add --scale relative|absolute

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 11: README — absolute mode + `pins` override

**Files:**
- Modify: `README.md`
- Test: `test/test_docs.jl` (if it link-checks/section-checks README; else manual)

**Interfaces:** none (docs).

- [ ] **Step 1: Write/adjust the doc check**

If `test/test_docs.jl` asserts README sections, add an assertion that a `## Absolute-scale fitting` section and a `pins` mention exist:
```julia
@testset "README documents absolute mode + pins" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    @test occursin("Absolute-scale", readme)
    @test occursin("pins", readme)
    @test occursin("[G6PD] (nM)", readme)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project test/test_docs.jl`
Expected: FAIL — section absent.

- [ ] **Step 3: Write the README section**

Add a `## Absolute-scale fitting (G6PD)` section covering: what absolute mode does and when to use it; the `[G6PD] (nM)` column requirement and units; the `scale=:absolute` keyword; the explicit `pins` override with the escalation ladder (rungs 1–4 at `mode1`, values from a prior relative fit); and the `kcat` band verdict. Include the four-rung code example from the spec.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project test/test_docs.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add README.md test/test_docs.jl
git commit -m "docs(readme): document absolute-scale mode and pins override

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 12: End-to-end synthetic discrimination smoke test

Proves the pipeline recovers a planted `kcat` and that the generating variant wins on the dual metrics — the acceptance test for the whole feature.

**Files:**
- Test: `test/test_absolute_discrimination.jl` (new; add to `runtests.jl`)

**Interfaces:** consumes the full absolute pipeline.

- [ ] **Step 1: Write the failing test**

Create `test/test_absolute_discrimination.jl`:
```julia
using FitRateEquation
using FitRateEquation: g6pd_config
using FitRateEquation.ChaFit
using FitRateEquation.ChaLaws
using Statistics: median
using Test

@testset "absolute mode recovers planted kcat" begin
    # Build synthetic single-scale forward data from a known kcat + deployed variant.
    d0 = load_dataset(g6pd_config())
    m  = FitRateEquation.v2_mechanism()
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    kcat_true = 178.0; Et = fill(5e-9, FitRateEquation.nrows(d0))
    rates = Float64[]
    for (i, cc) in enumerate(d0.concs)
        mtup = ChaFit.cha_macro_tuple(:G6PD, coords; keq=13.655, kf=kcat_true, Et=Et[i],
                    release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
        push!(rates, ChaLaws.cha_rate_G6PD(mtup; ChaFit._cha_row_kwargs(:G6PD, cc)...))
    end
    dsyn = Dataset(d0.concs, rates, d0.group, d0.keq, Et)
    fit = cha_fit_candidate(:G6PD, m, dsyn; n_restarts=8, maxiter=400, maxtime=90.0,
                            seed=1, keq=13.655, scale=:absolute,
                            pins=Dict(:Kd_6PGLn=>log10(2e-4), :Km_NADPH_rev=>log10(3.9e-6)))
    @test isapprox(fit.coords[:kcat], kcat_true; rtol=0.1)   # within 10%
end
```
(If forward-only synthetic rows must exclude reverse, all `d0.concs` already have `PGLn=0` for the deployed G6PD corpus figures; if any figure carries `PGLn>0`, filter those rows out before planting.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project test/test_absolute_discrimination.jl`
Expected: FAIL initially only if earlier tasks incomplete; otherwise it should pass once the pipeline is wired. (Written last, it is the integration gate.)

- [ ] **Step 3: (No new implementation)** — this test exercises Tasks 1–8. If it fails, debug the wired path, not the test.

- [ ] **Step 4: Run to verify it passes; then the full suite**

Run: `julia --project test/test_absolute_discrimination.jl && julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS, full suite green.

- [ ] **Step 5: Add to runtests and commit**

Add `include("test_absolute_discrimination.jl")` to `test/runtests.jl`.
```bash
git add test/test_absolute_discrimination.jl test/runtests.jl
git commit -m "test(cha): end-to-end absolute-mode kcat recovery

Plants a known kcat into synthetic forward data and asserts absolute mode
recovers it within 10% — the integration gate for the feature.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:** §5 API → Task 6/10; §6 scale-parameterized coords/`:kcat` → Task 4; §7 loss core+aggregators+`Et` → Tasks 2/3; §4 fiber-free `C=1` → Task 4 (`CHA_ABS_RELEASE_RATE`); §8 `[G6PD] (nM)` requirement → Task 1 (loader) + Task 6 (assertion); §9 pins+ladder → Task 5 (+ README Task 11); §10 dual CV+in-sample → Task 7 (CV) + Task 8/9 (identifiability/report); §11 files → all tasks; §12 testing → Tasks 2 (bit-identity), 3 (arithmetic/Et), 4 (fiber), 12 (synthetic), 6 (requirement/guard), 5 (pins guard). All spec sections map to a task.

**Placeholder scan:** No TBD/TODO. Two soft spots flagged for the executor to match live code, not invent: the exact `classify_cha` signature (Task 8) and the CLI parser name (Task 10) — both instruct "match the actual shape in the file", with concrete surrounding assertions. These are integration-with-existing-code notes, not missing content.

**Type consistency:** `Dataset` is 5-field everywhere after Task 1 (4-arg constructor preserved for existing callers). `cha_coords`/`cha_coord_bounds`/`_cha_loss_with_pins`/`cha_fit_candidate`/`resolve_cha_pins`/`_assert_pin_is_coord` all gain a consistent `scale::Symbol=:relative` kwarg. `CHA_ABS_RELEASE_RATE` is defined once (Task 3) and used in Tasks 4/12. `_cha_row_logratios!` signature is stable across Tasks 2/3.

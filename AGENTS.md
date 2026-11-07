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

> **EnzymeRates dependency:** this package runs on **upstream EnzymeRates 0.9.x**
> (`DenisTitovLab/EnzymeRates.jl`, registered in General; `[compat]` pins `0.9.1`).
> EnzymeRates >= 0.6 names constants structurally (`K_EG6P_to_E_G6P`,
> `k_ENADPH_to_E_NADPH`); FRE keeps the 0.5.0 composition-semantic names as its own
> vocabulary and translates at the boundary in `src/core/er_names.jl`. **Never call
> `EnzymeRates.rate_equation` / `fitted_params` / `parameters` / `substrates` /
> `_dependent_param_exprs` / `_raw_symbolic_rate_polys` directly — use the `er_*`
> wrappers**, or names will not match. `test/test_enzymerates_names.jl` replays a
> 0.5.0-generated oracle (`test/fixtures/er050_oracle.toml`) through the layer.
> It is self-contained — it vendors the data loader / mechanism builder /
> gauge helpers under `src/core/` (formerly borrowed from `fitting/mechanism_id/`).
> Mechanism topologies are written with terse opaque form names in
> `src/enzymes/<enzyme>.jl`; the builder auto-derives the decomposed call notation
> `E(NADP,G6P)` the upstream DSL requires. Rate-constant references use FRE's (EnzymeRates 0.5.0)
> composition-semantic names (`K_NADP_E`, `koff_NADPH_E`, `k_EG6PNADP_to_ENADPHPGLn`,
> …).

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
  gauge-mandated SS release: **Ru5P-release**. CO₂ and
  NADPH release stay rapid-equilibrium (PGD's NADPH release is fast, >800 s⁻¹, so no
  silent fiber is promoted). Super-node `S = {E_C ⇌ E_1}`, RE-linked by
  `s = CO2/Kd_CO2`. Topham et al. 1986 (*Eur J Biochem* 156) supply the topology (random
  binding, ordered release), not the RE/SS assignment: they place a slow, possibly
  rate-limiting step after hydride transfer and before Ru5P release (CO₂ dissociation in
  their Scheme 1, or decarboxylation, or isomerization of E·NADPH·Ru5P·CO₂), with
  half-of-the-sites reactivity as the alternative reading (p. 566), and they reject
  rapid-equilibrium binary complexes in TEA (p. 565). The >800 s⁻¹ is a lower bound for
  binary E·NADPH only. *(Corrected 2026-11-06: this said "Ru5P-release (Topham 1986
  rate-limiting)".)*

An **additive** fully-rapid-equilibrium PGD variant, **`:full_re`**, is also
registered (run it with `fit_consensus_equation(:pgd; variants=[:full_re])`). It keeps only catalysis at steady state
and makes *every* product release a rapid-equilibrium dissociation (ordered
CO₂ → Ru5P → NADPH), so there is **no promoted SS-release fiber** (`C = 1`,
apparent Km = α·Kd) and NADPH becomes a **competitive free-E ligand**. The law is
`ChaLaws.cha_rate_PGD_fullRE` with readoff `ChaInvert.cha_macro_readoffs_PGD_fullRE`;
every code path is guarded on `variant === :full_re`, so `:cha_base` (and G6PD)
are byte-identical when it is not selected. It is an **evaluation-only** variant —
not deployed into `PentosePhosphatePathway.jl`. See
[`docs/pgd_fullre_evaluation.md`](docs/pgd_fullre_evaluation.md) for the fit / CV /
Holten findings and the α = 1 recommendation.

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
  `koff` for PGD) is **fixed at `CHA_DEPLOY_RELEASE_RATE = 1e3`** (kf = 1 gauge units).
  The fit, the `micro_parameters.jl` deploy block, the Km readoffs and the plots all use
  this one value (`_default_release_rate` resolves to it), so the emitted constants are the
  fitted law. The value is a **convention** (release fast relative to catalysis,
  `C = 1+kf/r ≈ 1`), not a measurement. `release_rate` is free **only if `alpha` moves with
  it**: in the forward law (P = 0) the shape depends on `alpha/C`, so `alpha/C`, and hence
  `Km = alpha·Kd/C`, is what the data pin. Evaluating the same `alpha` at another rate is a
  different law. Under the centered loss the overall scale is gauged out, so the
  `kcat/Km = kf/(alpha·Kd)` value is meaningful only on the kf = 1 gauge; a refit at another
  rate moves it along with `alpha`. On the G6PD corpus a full refit fits equally well at any
  `koffQ` from 0.1 to 1e6 (see the report block below). **PGD's `koff` is not free the same
  way:** with `KdRu` held at 1, a full refit of `:cha_base` Mode 1 gave loss 0.118 at `koff` = 1,
  0.053 at 100, 0.025 at 1e3, 0.021 at 1e4 and 0.039 at 1e6, and `Km_PGA` moved from railed
  values to about 50-60 µM. On that corpus `koff = 1e3` and `KdRu = 1` shape the fit, so they are
  choices to revisit before deploying the PGD consensus law. **Run dirs with `package_version` < 0.10.0** were fit at
  `release_rate = 1` (`C = 2`) but emitted and plotted the same `alpha` at 1e3, so their
  `micro_parameters.jl`, reported Km's and plots describe a law with twice the fitted Km's and
  a worse loss than `goodness_of_fit.csv` reports. Absolute mode keeps its own
  `CHA_ABS_RELEASE_RATE = 1e8` and `:full_re` has no promoted release (`C = 1`); neither changed;
- the reverse catalysis `kr` is **Haldane-determined** from `keq`.

So the **free fit coordinates** (`cha_coords`, in `src/cha_fit.jl`) are exactly the
data-identifiable forward shape constants:

| Enzyme | `cha_coords` (fit directly) |
|---|---|
| G6PD | `Kd_NADP`, `Kd_G6P`, `Kd_6PGLn`, `alpha`, `Ki_NADPH`, `Ki_ATP`, `Ki_ATP_EG`, `Km_NADPH_rev` |
| PGD  | `Kd_NADP`, `Kd_PGA`, `alpha`, `Kd_CO2`, `Ki_NADPH`, `Ki_ATP`, `Ki_ATP_EN`, `Km_NADPH_rev` |
| PGD `:no_atp` | as `:cha_base` without `Ki_ATP`, `Ki_ATP_EN` |
| PGD `:full_re` | `Kd_NADP`, `Kd_PGA`, `alpha`, `Kd_NADPH`, `Kd_Ru5P`, `Kd_CO2` |

The **apparent Michaelis constants** are NOT coords — they are **derived** readoffs
`Km = alpha·Kd / C`, `C = 1+kf/release_rate` (`cha_apparent_km`): `Km_G6P` (G6PD);
`Km_PGA`, `Km_NADP` (PGD). They are read at `CHA_DEPLOY_RELEASE_RATE`, the rate the law is
fit and deployed at, and surfaced with class `:derived`, each paired with its specificity
`kcatKm_* = kf/(alpha·Kd)` on the kf = 1 gauge (`cha_specificity`). No mode pins or pulls an
apparent Km: PGD's Mode-2 soft `Km_PGA` anchor and Mode-3 `Km_PGA` override (the
`_cha_anchor_penalty` / `cha_anchors` machinery) were retired 2026-10-21.

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
  (`Km_NADPH_rev`). The opt-in anchor (`anchor_reverse=true`: `Km_NADPH_rev` = 3.9 µM,
  all modes) de-conflates the two when the data cannot; by default both are fit.
- **Lactone binding is fit, not pinned.** The Cha coord `Kd_6PGLn` is free in every
  variant and mode. The literature table's `Ki_6PGLn` (1027 µM, the true δ-lactone
  reverse Km, Beutler & Kuhl 1986 Fig 2a; was the 200 µM γ-lactone-analogue mean before
  0.9.0) is reference only: `resolve_cha_pins` never reads it, and pinning it would count
  the bundled corpus's Beutler Fig 2a rows twice.

### PGD — random-binding / ordered-release Bi-Ter (Topham 1986)

`NADP + PGA ⇌ CO2 + Ru5P + NADPH`. King-Altman skeleton (`src/enzymes/pgd.jl`):

```
E + NADP   ⇌ E_N        ┐ random substrate binding
E + PGA    ⇌ E_G        │
E_N + PGA  ⇌ E_NB       │
E_G + NADP ⇌ E_NB       ┘
E_NB → E_C   (SS — catalysis, gauge anchor)
E_C ⇌ E_1 + CO2         (RE — CO2 out first; super-node {E_C ⇌ E_1})
E_1 → E_2 + Ru5P        (SS — promoted Ru5P release; not Topham's slow step, see above)
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
  15 µM. *(This 37 µM figure predates 0.9.0, when `Km_NADPH_rev` was still anchored by
  default, and 0.10.0, which fits at the deploy `koffQ`; `report.md` states the current
  fitted value and its band position, so read it there.)*
- **PGD (0.10.0 fit at the deploy `koff` = 1e3, `KdRu` = 1):** the earlier conclusion,
  that the forward cross-term de-conflation was **refuted** (Mode 1 `Ki_NADPH`
  diagnostic/unconstrained, `alpha` railed, `Km_PGA` not data-identified), came from fits at
  `koff` = 1 and does **not** hold on this branch. Mode 1 now data-identifies `Ki_NADPH`
  (about 22 µM), `Kd_PGA` is data-identified, and `alpha` is about 12, off its limit; the loss
  fell 4.7-fold (0.118 to 0.025). **These new PGD values are not established:** the loss depends
  on the fixed `koff` and `KdRu` (see the `koff` series above), so the identification may
  reflect that choice. Mode 2 still pins `Ki_NADPH` to 17 µM (Cottreau). No mode anchors
  `Km_PGA` (2026-10-21). The high/unstable PGD leave-one-article-out CV traces to a
  **forward-CO₂ product-inhibition coverage gap** (the leave-Weisz1985-out fold);
  trimmed CV ≈ G6PD's. *This diagnosis predates the 2026-10-21 corpus correction,
  which raised Weisz 1985 Fig 6A–6D `[CO2]` 10× (and 6C/6D `[6PGA]` 10×), and the
  2026-10-22 removal of Ceyhan 2005 (a whole article, so leave-one-article-out CV now
  has 6 folds, not 7); it has not been re-run on the current corpus.*

> **Topology is FROZEN.** Do not add SS release steps to chase the forward
> `Ki_NADPH` — proven inert (the shared symbol + corpus geometry is the blocker, not
> the SS-step count). The only de-conflation routes are a decoupled forward symbol
> (the cross-term) or bench data. See `test/test_rec4_topology_freeze.jl`.

## Modes

**G6PD's rate law is under active development.** The mode/anchor/ablation machinery
below (`anchor_reverse`, mode1 vs mode2, the dead-end ablation variants) is
scaffolding for a not-yet-settled law, not a permanent API surface — expect these
gates to simplify as the mechanism converges (e.g. dropping the `anchor_reverse`
distinction if the anchor turns out unnecessary on the full corpus, or promoting an
ablation variant to be the deployed law outright).

Each (enzyme, variant) is fit in **modes** (`modes_for`: `mode1, mode2`), written to
separate, never-cross-ranked leaderboards. Pins are resolved structurally by
`ChaFit.resolve_cha_pins`. `modes_for` drops Mode 2 for a variant it would leave
identical to Mode 1 (no extra pin) — PGD `:full_re` has no `Ki_NADPH` coord, so it runs
`mode1` only.

- **Mode 1** — free fit. The data-only, family-agnostic headline. (G6PD pins
  `Km_NADPH_rev` in *all* modes only with the opt-in `anchor_reverse=true`.)
- **Mode 2** — literature pins applied. G6PD: `Ki_ATP` + `Ki_NADPH` (15 µM). PGD:
  **`Ki_NADPH` (17 µM) only**, so Mode 2 isolates the effect of that one literature
  constant.

**PGD mode changes (2026-10-21).** Mode 3 (the 38 µM `Km_PGA` override) and Mode 2's
soft 59 µM `Km_PGA` anchor were removed, along with the `anchors` penalty machinery that
served only them; PGD Mode 2 also stopped pinning `Ki_ATP` (the PGD corpus has no ATP, so
the pin was fit-inert). Run dirs written before this change (`package_version` < 0.8.0 in
`provenance.toml`) still carry `mode3` rows and plot normally; their Mode-2 values were fit
under the old anchor.

**`anchor_reverse` (G6PD switch, default `false` since 0.9.0).** `fit_consensus_equation`/
`run_g6pd`/`resolve_cha_pins` accept `anchor_reverse`; `true` pins `Km_NADPH_rev` to the
3.9 µM literature value (Gordon 1995, measured against the γ-lactone analogue) in every
mode. The default was `true` (variant-aware, `false` only for the three ablation variants)
until 2026-11-03, when it became `false` for every variant: the anchor is now opt-in and
`Km_NADPH_rev` is fit from the data.

The risk the anchor guarded against is real and data-dependent: at a single G6P the
`E·G6P·NADPH` dead-end (`Ki_NADPH`) and the bare-[NADPH] release term (`Km_NADPH_rev`)
collapse into one apparent Ki, and on the full 565-row forward-only corpus the deployed law
left `Ki_NADPH` railed without the anchor (an identifiability failure, not a parsimony
trade-off). Separating the two from data needs NADPH series at ≥ 2 G6P levels (or
NADPH-varied reverse data). So the **`NOT DEPLOYABLE` banner is data-driven**
(`_nadph_unidentified`): an anchor-off G6PD fit is tagged in `micro_parameters.jl`/
`report.md` (per variant × mode) exactly when its forward-NADPH coord classified
`:unconstrained` — `Ki_NADPH` where the variant has the dead-end, else `Km_NADPH_rev` (the
dead-end-dropped ablations, where the release term carries all forward NADPH inhibition).
Mode 2 literature-pins `Ki_NADPH`, so it passes by construction. The anchor state +
variants are always recorded in `provenance.toml`. No-op for PGD.

Pinned constants are tagged `:literature_pinned` **structurally** (a coordinate that
has a pin), never from profile curvature. `report.md` surfaces a **mode-agreement
check** on the forward shape constants: if a "pinned" constant was not actually flat
(it coupled into the forward fit) the point estimates diverge by `> tol_dex` and a
distortion warning fires — never hidden. Pinned/derived rows never enter the check,
so an intended data-vs-physiology contrast is not mislabeled a flatness warning.

**Standing rule (parsimony-first): mode 1 is deploy for future re-derivations.**
Mode 1 is the deploy source of truth; mode 2 is a literature-sensitivity diagnostic
only, never an independent deploy path. The gate is mode 1 passing coupled-flux/
simulation-accuracy validation downstream in `PPP_Experiments` — not literature
agreement; a simpler law disagreeing with a historical assay value may just mean the
assay isn't what the live parameter should be. Generalizes the prior `:no_atp`-only
convention to every future re-derivation — **not retroactive**: the deployed
full-G6PD law (mode2-pinned) is grandfathered until re-derived. PGD has nothing to
grandfather from this package: PentosePhosphatePathway.jl's `rate_PGD` is still the C2 law
with the `PGD_K1`…`PGD_k7r` constants of 2026-04-07 (its `4829776`), which predate this
package and the mode2/3 anchors, and no consensus PGD fit has been deployed. *(Corrected
2026-11-06: this said the deployed PGD was "fit under the retired mode2/3 `Km_PGA` anchors".
Those anchors applied only to the undeployed consensus PGD fits.)*
`:no_atp` (`src/enzymes/g6pd.jl`) is the first law under the rule: its mode1 fit is
fully data-identified, so mode1 is its deploy candidate.

**Dead-end-dropped ablation variants (2026-07-20; anchor default updated 2026-07-21, then
superseded by the global `anchor_reverse=false` default on 2026-11-03):**
three further exact `Ki -> Inf` limits of the deployed law, each dropping one or both
of the two abortive-ternary G6P-cross-term dead ends (`src/enzymes/g6pd.jl`):
`:no_g6p_nadph_deadend` (drops `E·G6P·NADPH`, `Ki_NADPH -> Inf`),
`:no_g6p_atp_deadend` (drops `E·G6P·ATP`, `Ki_ATP_EG -> Inf`),
`:no_g6p_both_deadends` (drops both) — all three defaulted to `anchor_reverse=false`
before every variant did.
`:no_g6p_atp_deadend` originally needed `anchor_reverse=true` to de-conflate
`Ki_NADPH` (2026-07-20 report below), but on the smaller, more-physiological 134-row
`mydata` corpus it jointly data-identifies BOTH `Km_NADPH_rev` and `Ki_NADPH` without
the anchor (`results/no_g6p_atp_deadend_anchor_false_mydata`, 2026-07-21) — not yet
confirmed on the full 565-row historical corpus (`SS_anchor_reverse_false`, same day,
still shows `Ki_NADPH` unconstrained there), and possibly reflecting `mydata`'s more
physiological assay conditions rather than a like-for-like comparison. Unresolved;
the exemption stayed scoped to these three variants until the 2026-11-03 global default
change, which replaced the variant list with the per-fit identifiability check above.

Model-selection verdict on the original (anchor-appropriate) fits
(`results/G6PD_deadend_variants_report.md`): `:no_g6p_atp_deadend` is a clean
parsimony win (forward constants unchanged, CV improves an order of magnitude in
mode2) — `Ki_ATP_EG` was an unconstrained nuisance dimension with zero forward-fit
benefit, consistent with its ~30 mM literature Ki (Ninfali 1996, superphysiological).
The NADPH-dead-end drops are NOT recommended: `alpha` loses `data_identified` status
and LOAO CV degrades or destabilizes — the bare-`[NADPH]` reverse-release channel
does not cleanly substitute for the dedicated forward cross-term. None of the three
has been deployed; only `:no_g6p_atp_deadend` is a plausible future deploy candidate,
pending the coupled-oxidative-flux check.

**Reading `unconstrained`:** it flags a term worth testing by ablation, not proof it's
droppable — the finding above is the counterexample: `Ki_ATP_EG` (`unconstrained`)
was safe to drop, `Ki_NADPH` (`data_identified`) was not. Classification alone
doesn't predict drop-safety; only an actual ablation does.

## Running

This is an installed Julia package, not an in-repo script pipeline: there is no
required working directory and no `--project=...` flag to remember; once `using
FitRateEquation` has loaded, every runner works from any current directory. See
`README.md` for the novice-level walkthrough; this section is the complete
reference.

**One documented entry point** — `fit_consensus_equation(enzyme::Symbol; …)` (writes
the eight artifacts below to `outdir` and returns the results):

```julia
using FitRateEquation

fit_consensus_equation(:g6pd | :pgd;
    variants=nothing, data_csv=nothing, smoke=false, outdir=nothing, nprocs=nothing,
    anchor_reverse=false,                      # G6PD-only; true pins Km_NADPH_rev (below)
    n_restarts=nothing, maxiter=nothing, maxtime=nothing, seed=nothing)  # power-user budget
```

The enzyme selector is symbol-only and case-insensitive. `variants` names alternative
laws (`[:no_atp]`, `[:full_re]`, `[:no_g6p_atp_deadend]`, …); omitting it fits the
deployed consensus law. Selecting a variant applies its bundled `row_filter` and names
the output folder after it.

**Thin discoverability aliases** (deployed-law calls; each forwards to
`fit_consensus_equation`):

```julia
run_g6pd(; outdir=nothing, smoke=false, nprocs=nothing, anchor_reverse=false)
run_pgd(;  outdir=nothing, smoke=false, nprocs=nothing)
```

There is no `run_g6pd_noatp` / `run_pgd_fullre`: use
`fit_consensus_equation(:g6pd; variants=[:no_atp])` and
`fit_consensus_equation(:pgd; variants=[:full_re])`.

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
  caller's `pwd()`), e.g. `results/G6PD_2026-07-16_smoke/`. A variant run is labeled
  after the variant (e.g. `G6PD_noatp_...`) so it never collides with the deployed-law
  run.
- **Own data:** point the entry point at your CSV via `data_csv`. The file must use
  the enzyme's **canonical column schema** (a non-conforming header errors, printing
  the header it expected); column remapping is not supported.
  ```julia
  fit_consensus_equation(:g6pd; data_csv="/path/to/my_corpus.csv", outdir="my_results")
  ```
- **Fit variants per enzyme:** G6PD `:SS_NADPH_release_rate_eq` (deploy) + `:no_atp`
  (ATP-free) + the dead-end-dropped ablations `:no_g6p_nadph_deadend`,
  `:no_g6p_atp_deadend`, `:no_g6p_both_deadends`; PGD `:cha_base` (deploy) + `:full_re`
  (fully-RE evaluation variant) + `:no_atp` (ATP-free `:cha_base`, evaluation variant;
  keeps every row). All are
  selected the same way — `fit_consensus_equation(:enzyme; variants=[…])` — with no
  dedicated runners; each variant carries its own `row_filter`, so e.g. G6PD `:no_atp` drops
  the ATP-bearing rows while the G6P dead-end ablations keep ATP a metabolite. The
  internal `fit_consensus_equation(cfg; variants=[…], row_filter=…)` method (a config
  object plus explicit filter) exposes a custom variant set / row filter for advanced
  in-package use. A `row_filter` is a `DataFrame -> DataFrame` function applied to the
  canonical corpus (`read_corpus`'s output) *before* the `Dataset` is built, so the
  rows it keeps are exactly the rows snapshotted to `fit_corpus.csv`
  (`drop_atp_rows` is the bundled example).
- **Config builders are internal, not exported.** `FitRateEquation.g6pd_config` /
  `pgd_config` (data CSV, `deploy_keq`, metabolite columns/units) live in
  `src/configs/G6PD.jl`, `src/configs/PGD.jl` and are the objects
  the symbol entry point resolves internally — user code does not build them; it passes
  `data_csv=` to `fit_consensus_equation(:enzyme; …)`. The bundled corpora resolve via
  `pkgdir(FitRateEquation)` so they load correctly regardless of installation location;
  `data_csv` overrides them, and its header must match the canonical schema. The fit
  itself uses each figure's own apparent Keq, read per-row from the corpus `keq_col`
  column; the config's `deploy_keq` is the single readout/deploy value — G6PD `13.655`
  (apparent Keq at cellular pH 7.2 / 37 °C), PGD `0.17` (pH-flat Bi-Ter CO₂ aq, 37 °C).

**CLI** — a thin shim (`bin/fitrateequation`, `using FitRateEquation;
cli_main(ARGS)`) dispatching `fit_consensus_equation` in-process, no subprocess:

```sh
bin/fitrateequation <g6pd|pgd> [--smoke] [--nprocs N] [--outdir DIR] [--data CSV] [--variant NAME]
bin/fitrateequation plot <run_dir>
bin/fitrateequation help
```

Subcommands are `g6pd | pgd | plot | help`. `--data` and `--variant` are valid
with any enzyme subcommand (`--variant NAME` fits an alternative law, e.g. `no_atp`,
`full_re`, `no_g6p_atp_deadend`); the old `g6pd-noatp` subcommand is gone — use
`g6pd --variant no_atp`.

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
  from that file, and the data points come from that run's own `fit_corpus.csv`.
  Output is one `<run_dir>/plots/<variant>_<mode>_fit.png` per row.
- **Run dirs written before 0.2.0 are not plottable** and raise a clear error. They
  predate `fit_corpus.csv`, and a custom `data_csv` or `row_filter` cannot be
  recovered from outside the run — earlier versions silently fell back to the bundled
  default corpus, plotting the fitted law against data it was never fit on. Re-run
  the fit to get a plottable dir.
- Each prediction is the **deployed** Cha law: the macro tuple is assembled at
  `CHA_DEPLOY_RELEASE_RATE` and evaluated through `cha_rate_*`, so the curve matches
  the law written to `micro_parameters.jl`. Each figure's curve uses **that
  figure's own apparent Keq** (the reverse arm's Haldane `kr`), read per-row from
  the corpus. A per-figure Vmax gauge normalizes each panel, so the plots show
  **shape agreement**, not absolute scale.
- **G6PD and PGD only.**
- The non-Makie logic (enzyme detection, config lookup, coordinate readback,
  `read_fit_corpus`) lives in `src/plot_support.jl` and is reachable without
  loading CairoMakie; only the render loop itself lives in the extension. The first
  invocation after `using CairoMakie` precompiles it (~a few minutes; cached after).
  Because it only reads the run dir (`macro_constants.csv` + `fit_corpus.csv`), it
  never re-runs or perturbs a fit.
- **Running the render check.** `test/test_plot_render.jl` is not part of the default
  suite (Makie cold-precompiles for minutes and the file runs a real smoke fit), and
  CairoMakie is a `[weakdeps]` entry, so it is not loadable from the main project env
  either. Run it through the committed sub-environment — that env is a `Project.toml`
  only (no Manifest), so a fresh checkout needs the instantiate step once:
  ```sh
  julia --project=test/render -e 'using Pkg; Pkg.instantiate()'   # first run only
  julia --project=test/render test/test_plot_render.jl
  ```
  It fits a deliberately non-default corpus so a regression to config-derived corpora is
  visible as a row-count failure.

## Outputs (per run dir)

| File | Contents |
|---|---|
| `macro_constants.csv` | Resolved `cha_coords` per mode (class `data_identified` / `literature_pinned` / `unconstrained`) **plus** the derived apparent Michaelis constants (`Km_G6P` / `Km_PGA` / `Km_NADP`, class `derived`). |
| `fit_corpus.csv` | The exact rows this run fit — Molar units, post-`row_filter`, with `Rate` / `source` (`Article\|Fig`) / `Apparent_Keq` / `X_axis_label`. Makes the run dir self-describing about its corpus; the plotter reads this, never a config default. |
| `goodness_of_fit.csv` | In-sample loss + CV under **both** fold schemes per variant × mode: `cv_article_{mean,se,n_folds}` (leave-one-article-out, the **headline**: figures from one article share a lab, enzyme prep and assay, so holding out a figure leaks) then `cv_figure_{mean,se,n_folds}` (leave-one-figure-out, optimistic; for comparison with EnzymeRates.jl, whose LOOCV holds out one `Article_Fig` group, includes 1-row groups and uses SE = std/√n_groups, all reproduced here). Since 0.9.0; earlier run dirs have a single article-CV pair `cv_mean`/`cv_se`. |
| `cv_folds.csv` | Per-fold losses behind both schemes: `variant, mode, scheme (article\|figure), held_out (article or Article\|Fig key), n_test_rows, loss`; article folds first, each scheme in corpus order. For paired per-fold comparisons (map `\|`→`_` to join EnzymeRates' `cv_fold_<group>` columns). Since 0.9.0. |
| `identifiable_functions.csv` | Per-mode eigen-spectrum of the identifiable macro-coordinate subspace. |
| `micro_parameters.jl` | ODE-ready deploy block — the **closed-form** `ChaDeploy.cha_deploy_micro` inverse of the fitted macro coords (no LM / root-finding; reproduces the law to machine precision). |
| `report.md` | Per-mode macro tables, cross-mode agreement, the G6PD `koffQ` hybrid block, and the de-conflation caveat. |
| `provenance.toml` | Package version / EnzymeRates SHA, budget, source corpus path (`data_csv`), row count, seed, smoke flag — reproducibility record. `_git_sha` reads `"unknown"` in an installed depot (no `.git` to inspect) — this is expected and already `try/catch`-guarded; the package version is the durable provenance handle for an installed copy. |

The `koffQ` block (`src/cha_koffq_report.jl`) states the deploy convention and tabulates
a profile: every free constant refit at each `koffQ`, with the PGLn-bearing rows up-weighted.
Report-only, never fed to deploy. On the current corpus the profile is flat, because `alpha`
moves with `C` and `Kd_6PGLn` absorbs the reverse-row effect (those rows sit at a single
[NADPH]), so the corpus does not identify `koffQ`. Since 0.10.0 this replaces the old "data-identified
koffQ" row, a scan with the constants held fixed that just returned the rate they were fit at.

A `plots/` subdirectory is **not** written by the fit run itself — it is produced
on demand by `plot_consensus_fit` (see [Plotting the fit over the
data](#plotting-the-fit-over-the-data) above).

Since **0.3.0**, the exported `write_outputs` **requires `corpus=`** — the exact
(post-`row_filter`) DataFrame that was fit. This is a **breaking change**: a caller
that omits it now errors instead of writing a run dir with no `fit_corpus.csv`,
which the plotter cannot render. `fit_consensus_equation` always passes it, so only direct
`write_outputs` callers are affected. 0.3.0 also changed the `fit_corpus.csv`
artifact itself: its metabolite columns are now emitted **sorted alphabetically**
(G6PD: `ATP, G6P, NADP, NADPH, PGLn`, previously the config `Dict`'s hash order
`NADP, ATP, G6P, PGLn, NADPH`), so read that file by column name — any consumer
reading it positionally breaks silently.

## Determinism: the byte-identity smoke fixtures

`test/fixtures/{g6pd,pgd}_smoke_macro_constants.csv` are committed reference
`macro_constants.csv` outputs from a fixed-seed smoke run
(`test/test_byte_identity.jl`). Every test run re-runs `run_g6pd(smoke=true,
nprocs=1)` / `run_pgd(smoke=true, nprocs=1)` into a temp dir and compares
**variant/mode/name** against the fixture: which coordinates each variant x mode
produces, and in what order. That property is machine-independent, so the gate is
always on, on every machine and in CI.

The **fitted values are not compared**, and this is deliberate. At smoke budget the
seeded CMA-ES trajectory is under-converged, so different BLAS/RNG reduction
ordering puts a different machine on a materially different point — measured spread
is ~30-58% on `data_identified` coordinates and ~90-100% on `unconstrained` ones.
`class` is not exempt: it is decided by `cha_classify` from the Hessian at that same
drifted optimum. No single tolerance both survives that spread and still catches a
real regression.

An exact value/class comparison did live here as an opt-in tier. It was retired
because it gated on Julia `VERSION` while the property it actually required was
reference-*machine* identity: it could not pass in CI (rotating `ubuntu-latest`
images) and did not pass on the maintainer's machine either, so it never gated
anything. Behavioural regressions in the fit path are covered by
`test_cha_deploy`, `test_cha_classify`, `test_mode_agreement` and the
goodness-of-fit assertions. Exact reproduction of fit output *is* still gated,
just not across machines: `test/test_parallel_equivalence.jl` asserts the
`fit_consensus_equation` pmap path is **bit-identical** to the serial baseline — two runs in the
same process, which is a determinism guarantee the package can actually make.

If a structural comparison fails, that is the gate doing its job: the coordinate set
or its ordering changed. Either the change is a bug, or the fixtures need a
deliberate, reviewed regeneration — do not silently update them to make the test
pass.

## Fit speed: the loss must not allocate or dispatch

CMA-ES calls the objective thousands of times per restart, so its cost sets the cost of a
run. `_cha_objective` (`src/cha_fit.jl`) works out each figure's rows and keq once per fit
(`_ChaLossData`) and picks the rate law once (`_cha_law`). Everything enzyme-specific then
dispatches on that law's type, so no row chooses anything at run time. Do not look up
`enzyme::Symbol` or build a `Dict` inside the row loop: that cost about 10x on 0.10.0.

`test/test_loss_speed.jl` compares the objective, the public loss functions and a seeded fit
with `==` against `test/frozen_cha_loss.jl`, a verbatim copy of the 0.10.0 loss. It also
requires zero allocation, an inferred `Float64`, and at least 4x the frozen loss's speed. Do not
edit the frozen copy to make a test pass; a deliberate change to the arithmetic means updating
it in the same reviewed change.

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
src/cha_koffq_report.jl   G6PD koffQ report (deploy convention + all-constants-refit profile)
src/promotable.jl     declarative registry of promotable slow steps (G6PD nadph_release)

src/enzymes/g6pd.jl, src/enzymes/pgd.jl
                       per-enzyme topology + lit values + alias/ki-ratio maps
src/configs/G6PD.jl, src/configs/PGD.jl
                       data CSV path (pkgdir-resolved), deploy_keq, metabolite columns
src/enzyme_wiring.jl   EnzymeWiring registry struct + accessors + modes_for
src/mechanisms.jl      generic dead-end / King-Altman builder helper

src/core/data.jl, src/core/mechbuild.jl, src/core/gauge.jl
                       vendored data loader / mechanism builder / gauge (formerly borrowed
                       from fitting/mechanism_id/ in the sibling repo)
src/macro_collect.jl   macroscopic-constant readback from the Haldane-reduced symbolic rate
                       law; retained as the exactness-anchor test spine (test_pgd_macro_collect.jl)
src/cv.jl              LIVE: CV fold schemes (CV_SCHEMES = article, figure; _cv_folds,
                       _article_folds, _group_folds, _subset, _cv_summary) used by run.jl's
                       Cha-space CV (_cha_loocv and the pmap task path)

src/plot_support.jl    non-Makie plot helpers (enzyme detection, config lookup, coordinate
                       readback, plot dataframe assembly) — reachable without CairoMakie
ext/FitRateEquationMakieExt.jl
                       CairoMakie package extension: the actual render loop
                       (plot_consensus_fit), loads only when CairoMakie is loaded

data/                  bundled corpora: G6PD_all_EnzymeData.csv (byte copy of
                       PPP_Experiments fitting/G6PD/rate_eq/, 477 rows / 20 figures /
                       6 articles as of 2026-11-06 (436 / 18 before Mbata2026 Figs 1B,
                       1C were added); δ-lactone-only reverse data since 2026-11-03 — Beutler
                       1986 Fig 2a in, Beutler Fig 3 and Gordon 1994 Fig 1A/1B γ-lactone
                       analogue rows out),
                       PGD_EnzymeData_with_CO2.csv (PGD: byte copy of
                       PPP_Experiments fitting/PGD/PGD_EnzymeData_with_CO2.csv, 375 rows /
                       20 figures / 6 articles; the 400-row corpus minus Ceyhan 2005 Fig 6,
                       which was assayed in 100 mM phosphate + 10 mM Mg2+ and has a 5-10x
                       higher apparent Km_PGA than every other figure. Bundled 2026-10-22 as
                       PGD_EnzymeData_with_CO2_no_Ceyhan.csv; renamed 2026-10-31 to match
                       PPP_Experiments 9eab16c, contents unchanged)
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

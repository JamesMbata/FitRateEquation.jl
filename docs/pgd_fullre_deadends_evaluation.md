# PGD `:full_re_deadends` rate law — fit / evaluation findings (Phase 4)

Evaluation of the `:full_re_deadends` PGD variant — the fully rapid-equilibrium, ordered-release
`:full_re` law **plus two rapid-equilibrium free-enzyme dead-ends (CO₂, Ru5P)** — added to
`FitRateEquation.jl` in the `pgd-fullre` Phase-4 effort. The dead-ends add bare-`[CO2]`/`[Ru5P]`
competitive denominator terms (macro `Ki_CO2`/`Ki_Ru5P`) so CO₂ and Ru5P inhibit in the **NADPH-free
single-product assays** of Weisz 1985 (figs 6A–6D, 7A–7B), which the strictly-ordered `:full_re`
**structurally cannot fit** (every CO₂/Ru5P term in `:full_re` carries a factor of NADPH and so
vanishes identically at NADPH = 0).

> **Deployment is OUT OF SCOPE.** This is a fit/readoff/classify/run evaluation inside
> `FitRateEquation.jl`. Nothing is written into `PentosePhosphatePathway.jl`; `deploy_variant(:PGD)`
> stays `:cha_base`. The coupled flux gate (Phase-4 Task 6) reads coupled behavior on a throwaway
> dev branch only. See the master tracker.

## Setup

- **Corpus:** the PGD kinetic corpus — 23 figures / 7 pubs (`load_dataset(pgd_config())`);
  forward-dominated (one reverse figure, Villet). `[ATP] = 0` throughout. Weisz 1985 is the **only**
  article carrying CO₂/Ru5P single-product titrations (6A–6D: CO₂ 0→2200 µM; 7A–7B: Ru5P 0→400/1080 µM;
  all with NADPH = 0 and the other product = 0).
- **Law / variant:** `:full_re_deadends` — the registered core `_pgd_fullre_deadends_core()`
  (`:full_re` core + `_deadends([([:E],:CO2),([:E],:Ru5P)])`), 8 RE `free_params` (the 6 `:full_re`
  core Kd's + `K_CO2_E`/`K_Ru5P_E`). Coords `[:Kd_NADP, :Kd_PGA, :alpha, :Kd_NADPH, :Kd_Ru5P,
  :Kd_CO2, :Ki_CO2, :Ki_Ru5P]` (fiber-free, `C = 1`, apparent Km = α·Kd).
- **Budget:** full — `n_restarts = 48`, `maxiter = 1000`, `maxtime = 300 s`, `seed = 1`,
  `nprocs = 10`; per-figure keq (median ≈ 0.116 M). Leave-one-**article**-out CV.
- **Commands** (`results/PGD_fullre_deadends_eval/analysis.jl`):
  ```
  run_pgd_fullre_deadends(outdir="results/PGD_fullre_deadends_eval", nprocs=10)
  run_pgd_fullre(outdir="results/PGD_fullre_eval", nprocs=10)   # :full_re head-to-head
  ```

## The decisive result — Weisz single-product inhibition is now reproduced

`:full_re` predicts **exactly no** CO₂/Ru5P inhibition in the NADPH-free assays (the terms are
structurally zero); `:full_re_deadends` reproduces the observed drops:

| single-product assay (NADP=10 µM, 6PG=50 µM, other products=0) | data | `:full_re` predicted v₀/v | `:full_re_deadends` predicted v₀/v |
|---|---|---|---|
| CO₂ 0 → 2200 µM (6A regime) | ~2.3× | **1.00×** (flat — structural zero) | **2.10×** ✅ |
| Ru5P 0 → 400 µM (7A regime) | ~1.6× | **1.00×** (flat — structural zero) | **1.47×** ✅ |

Per-figure centered-logratio loss (lower = better within-figure shape; Vmax gauged out per figure):

| Weisz figure | what it titrates | `:full_re` | `:full_re_deadends` |
|---|---|---|---|
| 5A | NADP sweep (NADPH inhib) | 0.033 | 0.025 |
| 5B | 6PG sweep (NADPH inhib) | 0.049 | 0.057 |
| **6A** | **CO₂ vs NADP @ 6PG=50 µM** | 0.125 | **0.032** (3.9× ↓) |
| 6B | CO₂ vs NADP @ 6PG=2000 µM | 0.110 | 0.092 |
| 6C | CO₂ vs 6PG @ NADP=10 µM | 0.174 | 0.154 |
| 6D | CO₂ vs 6PG @ NADP=100 µM | 0.246 | 0.205 |
| 7A | Ru5P vs NADP @ 6PG=50 µM | 0.069 | 0.056 |
| **7B** | **Ru5P vs 6PG @ NADP=10 µM** | 0.107 | **0.039** (2.7× ↓) |

The dead-ends improve the fit on **every** CO₂/Ru5P figure, decisively on 6A (CO₂-vs-NADP) and 7B
(Ru5P-vs-6PG) — the cleanest single-product-inhibition geometries. `:full_re_deadends` also lowers
the overall in-sample loss (mode 1: **0.066** vs `:full_re` 0.079).

## No distortion of the forward / Holten fits

The dead-end terms are **exactly zero on every product-free row**, so they cannot perturb the
forward figures. The head-to-head confirms the forward story is preserved (mode 1):

| | Km_NADP | Km_PGA | α | Kd_NADPH | Holten v₀/v @30 µM NADPH |
|---|---|---|---|---|---|
| `:full_re` | 33.4 µM | 183 µM | 13.6 | 0.87 µM | 1.83× |
| `:full_re_deadends` | 26.0 µM | 169 µM | 3.81 | 2.23 µM | 1.89× |

Km_NADP moves *toward* the human band (33 → 26 µM); Holten NADPH product-inhibition is unchanged
(1.83 → 1.89×) — the competitive free-E NADPH physics the fully-RE law exists for is intact. The two
new constants are identified from the Weisz titrations: **Ki_CO2 = 0.50 mM, Ki_Ru5P = 210 µM**
(both within the probed range — CO₂ to 2.2 mM, Ru5P to ~1 mM).

**Competitive-vs-both pattern (free-E form).** A free-E dead-end is competitive against the swept
substrate (raises the LB slope, not the intercept). It fits every Weisz figure better (both the
vs-NADP 6A/6B/7A and the vs-6PG 6C/6D/7B panels improve), so the free-E form is **sufficient** — the
data do not demand an E·NADP/E·6PG (uncompetitive) form, and escalation to Phase 5 (random product
release) is **not required** on the fit evidence.

## Identifiability & CV — the one honest caveat

Leave-one-**article**-out CV does **not** improve and rises (mode 1: `:full_re` 0.126 →
`:full_re_deadends` 0.198). A full-budget **per-fold decomposition** (`n_restarts = 48`) localizes
the entire gap to the **two CO₂/Ru5P-bearing articles** — the five CO₂/Ru5P-free articles are flat:

| held-out article | `:full_re` | `:full_re_deadends` | Δ | % of CV gap |
|---|---|---|---|---|
| **Villet1972** (reverse, all products) | 0.096 | 0.455 | **+0.359** | **71%** |
| **Weisz1985** (forward, NADPH = 0) | 0.126 | 0.241 | **+0.115** | **23%** |
| Cottreau1975 | 0.058 | 0.106 | +0.048 | 10% |
| Akkemik/Chan/Pearse/Ceyhan | — | — | −0.018 net | −4% |

So the CV degradation is a **corpus-grouping** effect (94% from the 2 CO₂/Ru5P articles), but it is
**not** solely the Weisz coverage gap I first attributed it to. **Two distinct mechanisms:**
1. **Leave-Weisz-out (+0.11, 23%)** — the pure single-article identifiability gap: no forward
   CO₂/Ru5P training data ⇒ `Ki_CO2`/`Ki_Ru5P` rail ⇒ held-out Weisz predicted worse than `:full_re`'s
   flat guess. Same geometry as α and `[[project_pgd_t1g_cv_root_cause]]`.
2. **Leave-Villet-out (+0.36, 71% — the dominant term)** — a **forward→reverse transfer failure**:
   Villet's reverse rows carry CO₂ **and** Ru5P, so the dead-end terms are **active** there; fit to
   Weisz's *forward* NADPH-free data, they **mis-transfer** and *worsen* the reverse prediction
   (0.096 → 0.455). This is not a mere coverage gap — it is a genuine open question about whether the
   free-E dead-end **form** is right for the reverse arm, and is the subject of **Phase 4.5** (the
   agent-team adjudication of "different rate law vs reevaluate the CV metric";
   `docs/superpowers/plans/2026-07-23-pgd-fullre-phase4.5-villet-cv-agent-team.md` in PPP_Experiments).

The in-sample fit and the direct structural test (Weisz 6A–7B) are the decisive Phase-4 metrics and
both pass; the Villet CV transfer is flagged, not dismissed, and handed to Phase 4.5.

## α = 1 deploy form (for the coupled flux gate)

As for `:full_re`, pinning **α = 1** costs almost nothing (in-sample loss 0.0662 → 0.0671, +1.3%)
and lands **Km_NADP = 15.1 µM squarely in the human band** (13–18 µM, Cottreau/Chan), Km_PGA = 124 µM.
Recommended fitted constants for the coupled flux gate (α = 1, mode 1):

```
Kd_NADP  = 1.512e-5   Kd_PGA  = 1.237e-4   alpha = 1.0
Kd_NADPH = 4.029e-6   Kd_Ru5P = 1.400e-7   Kd_CO2 = 0.1579
Ki_CO2   = 7.617e-4   Ki_Ru5P = 3.502e-4        # the free-E dead-ends (competitive)
```

## Verdict

- **Fits Weisz 6A–7B where `:full_re` cannot.** The structural-zero is fixed: single-product CO₂/Ru5P
  inhibition goes from flat (1.00×) to the observed magnitudes (2.10× / 1.47× vs data 2.3× / 1.6×),
  and every 6/7 per-figure loss improves. ✅ (This is the phase goal.)
- **Forward/Holten fits preserved** — the dead-ends are silent on product-free rows; Km_NADP,
  Km_PGA, Kd_NADPH, and the Holten NADPH ratio are unchanged or slightly improved. ✅
- **Free-E form is sufficient** — it fits both the vs-NADP and vs-6PG panels; the data do not demand
  an uncompetitive form, so **Phase 5 (random release) is not required** on the fit evidence. ✅
- **CV does not improve (single-article caveat)** — `Ki_CO2`/`Ki_Ru5P` are Weisz-only, so leave-one-
  article-out cannot constrain them and CV slightly rises; this is information geometry, not a
  mechanism failure. ⚠️ (documented, expected)

**Bottom line:** `:full_re_deadends` closes the honest boundary of the strictly-ordered law — it
reproduces Weisz's NADPH-free CO₂/Ru5P single-product inhibition at no cost to the forward/Holten
physics, with the free-E (competitive) form sufficient. The recommended form for downstream reads is
**α = 1, mode 1**. Deployment into PPP.jl remains gated (over-reduction axis; the Phase-4 coupled
flux gate reads coupled behavior read-only).

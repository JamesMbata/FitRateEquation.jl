# PGD fully-RE (`:full_re`) rate law — fit / evaluation findings

Evaluation of the fully rapid-equilibrium, ordered-product-release PGD law (`:full_re`) added
to `FitRateEquation.jl` in the `pgd-fullre` effort (Phases 1–3). The law is: **random RE
substrate binding (NADP, 6PG), a single SS catalytic gauge, ordered RE product release
CO₂ → Ru5P → NADPH**. Because NADPH releases last it re-binds free E first in reverse ⇒ it is a
**competitive free-E ligand** (`Q/Kd_NADPH` in the denominator) — the term the deployed
`cha_base` law dilutes through its SS-release super-node.

> **Deployment is OUT OF SCOPE.** This is a fit/readoff/classify/run evaluation inside
> `FitRateEquation.jl`. Nothing is written into `PentosePhosphatePathway.jl`; `deploy_variant(:PGD)`
> stays `:cha_base`. See the Deployment note.

## Setup

- **Corpus:** the PGD kinetic corpus — **436 rows / 23 articles** (`load_dataset(pgd_config())`);
  forward-dominated (one reverse figure). `[ATP] = 0` across all rows, so the effectors-off core is
  well-posed (V1's ATP dead-ends would be silent anyway).
- **Law / variant:** `:full_re` — the registered effectors-off core (`_pgd_fullre_core()`), 6 RE
  binding `free_params`, no `K_ATPinh_*`. Core coords
  `[:Kd_NADP, :Kd_PGA, :alpha, :Kd_NADPH, :Kd_Ru5P, :Kd_CO2]` (fiber-free, `C = 1`,
  apparent Km = α·Kd).
- **Keq:** per-figure during the fit (scalar median ≈ 0.116 M); deploy anchor `keq = 0.17 M`
  (Villet & Dalziel, 38 °C) in `provenance.toml`.
- **Budget:** full — `n_restarts = 48`, `maxiter = 1000`, `maxtime = 300 s`, `seed = 1`,
  `nprocs = 10`. Three PGD modes (mode 1 unanchored; modes 2/3 apply the Km_PGA anchor
  59 µM / 38 µM). Leave-one-**article**-out CV.
- **Commands:**
  ```
  julia --project=. -e 'using FitRateEquation; run_pgd_fullre(outdir="results/PGD_fullre_eval", nprocs=10)'
  julia --project=. -e 'using FitRateEquation; run_pgd(outdir="results/PGD_cha_base_eval", nprocs=10)'   # cha_base head-to-head
  julia --project=. results/PGD_fullre_eval/analysis.jl     # α=1 vs α-free, Holten
  julia --project=. results/PGD_fullre_eval/prodinhib.jl    # product-inhibition fingerprint
  ```

## Fit & CV (leave-one-article-out) — `:full_re` vs `cha_base`

| mode | law | in-sample loss | CV mean ± se |
|---|---|---|---|
| mode 1 | **full_re** | 0.079 | **0.117 ± 0.043** |
| mode 1 | cha_base | 0.107 | 1.262 ± 1.077 |
| mode 2 | **full_re** | 0.085 | **0.153 ± 0.053** |
| mode 2 | cha_base | 0.187 | 1.702 ± 1.315 |
| mode 3 | **full_re** | 0.089 | **0.120 ± 0.036** |
| mode 3 | cha_base | 0.190 | 1.739 ± 1.335 |

**`:full_re` is decisively better on generalization** — CV ≈ **8–11× lower** than `cha_base` in
every mode, with a far tighter standard error (0.04–0.05 vs 1.1–1.3, i.e. no single article
dominates the leave-one-out for `:full_re`, whereas `cha_base`'s CV is one-fold-dominated).
In-sample loss is also lower. This is with **fewer** free parameters (6 vs `cha_base`'s 8).

## Identifiability & α

`:full_re` mode-1 identifiability spectrum (6 core coords): 5 stiff eigenvalues + 1 soft
(`data_identified` × 5, one `unconstrained`), consistent across all three modes. The single soft
direction is **the α–Kd ridge**, and it is real:

| quantity | mode 1 | mode 2 | mode 3 |
|---|---|---|---|
| α | 13.6 (identified) | 1.41 (identified) | 0.58 (unconstrained) |
| Kd_Ru5P | **unconstrained** | 9.9 µM | 6.6 µM |

α does not agree across modes (0.58 → 13.6): it trades off against Kd_NADP/Kd_PGA and, in
mode 1, leaves Kd_Ru5P unconstrained.

**α = 1 vs α-free (mode 1), direct fits:**

| | loss | apparent Km_NADP | apparent Km_PGA | spectrum |
|---|---|---|---|---|
| α-free | 0.0793 | 33.4 µM | 183 µM | Kd_Ru5P **unconstrained** (rank 5/6) |
| **α = 1** | 0.0826 (+4%) | **14.1 µM** | 107 µM | **all 5 free coords identified** (rank 5/5) |

Pinning **α = 1** costs only ~4% in loss but (i) **rescues Kd_Ru5P** (unconstrained →
data_identified), (ii) pulls apparent **Km_NADP into the human band** (33 → 14 µM; lit 13–18 µM,
Cottreau/Chan), and (iii) removes the one soft direction. The corpus does not support a distinct
α (no crossed 2-D NADP×6PG grids at the needed resolution) — the extra DOF is spent on a
non-generalizing α. **Recommendation: deploy α = 1** (a single-site RE law); keep α-free only as a
diagnostic. (Km_PGA stays high in mode 1 — 107–183 µM — because the corpus lacks sub-Km [6PG]
coverage; the 38–80 µM band is imposed by the mode-2/3 anchor, which lands Km_PGA at 60/38 µM.)

## Holten NADPH product-inhibition — the law's raison d'être

The whole point of the fully-RE law is that NADPH re-binds free E first in reverse, giving a
**competitive `Q/Kd_NADPH` free-E term**. The corpus **data-identifies** this constant:

| | competitive NADPH constant | how obtained |
|---|---|---|
| **full_re** | **Kd_NADPH = 0.87 µM** (mode 1; 2.1 µM mode 2) | **data-identified** from the forward corpus |
| cha_base | Km_NADPH_rev = 22.3 µM (mode 2) | diluted through the SS-release super-node |
| cha_base | Ki_NADPH = 17 µM | **literature-PINNED** (Cottreau); rails to the 1 M bound when unanchored (mode 1) |

The ratio **Km_NADPH_rev / Kd_NADPH ≈ 25.5×** directly reproduces the **~28× SS-release
dilution** that the `hos_cha` session identified. The decisive difference is **identifiability**:
`:full_re` constrains a **sub-µM competitive Kd_NADPH straight from the forward data**, whereas
`cha_base` **cannot identify forward NADPH inhibition at all** — it rails Ki_NADPH to the bound and
only inhibits once Ki_NADPH is pinned to the literature 17 µM (which then "collapses the coupled
ODE at high NADPH demand", per `hos_cha`).

Intrinsic NADPH inhibition produced by the fitted `:full_re` law (v₀/v at [NADP] = 30 µM,
[6PG] = 200 µM): **1.83× at 30 µM NADPH, 3.78× at 100 µM** — physiological product inhibition
carried by a data-constrained constant, no pin required.

## Product-inhibition fingerprint (structural, from the closed form)

Denominator `D = 1 + A/Kd_NADP + B/Kd_PGA + AB/(α·Kd_NADP·Kd_PGA) + Q/Kd_NADPH +
Q·R/(Kd_NADPH·Kd_Ru5P) + Q·R·C/(Kd_NADPH·Kd_Ru5P·Kd_CO2)` (A = NADP, B = 6PG, Q = NADPH,
R = Ru5P, C = CO₂). Numeric spot-checks on `cha_rate_PGD_fullRE` (mode-1 coords):

| product | prediction (ordered RE release) | numeric check |
|---|---|---|
| **NADPH** | competitive vs **both** substrates (free-E term) | apparent Km_NADP 28.6 → 37.4 µM at 20 µM NADPH (**1.31×**, competitive) |
| **Ru5P** | inhibits **only** when NADPH present (Q·R term) | v₀/v = **1.000 alone**, **1.37×** with 5 µM NADPH |
| **CO₂** | inhibits **only** with NADPH **and** Ru5P present (Q·R·C term) | v₀/v = **1.000 alone**, **1.25×** with NADPH + Ru5P |

This is the exact ordered-release signature: NADPH competitive against both substrates; Ru5P and
CO₂ are **conditional** inhibitors gated by the later-released product(s). It departs from the
classic graphical assignments that treat each product's inhibition independently. The config-gated
dead-ends (`:Ki_ATP`/`:Ki_NADPH`) are OFF by default and were not needed: `[ATP] = 0` throughout
the corpus, and the NADPH E·6PG dead-end is the LOOCV-refuted term (PGD Problem 5).

### Known limitation — Weisz 1985 6A–7B CO₂/Ru5P single-product inhibition (structural)

The same ordered-release structure that makes NADPH the sole free-E competitor also makes the law
**unable to reproduce the Weisz 1985 CO₂ (6A–6D) and Ru5P (7A–7B) product-inhibition figures**, and
this is **structural, not a parameter artifact**. Those assays titrate CO₂ (0→2200 µM) or Ru5P
(0→400/1080 µM) with **NADPH = 0** (and the other product = 0) — a single-product inhibition design.
In the law, CO₂ appears *only* in `Q·R·C/(…)` and Ru5P *only* in `Q·R/(…)` and `Q·R·C/(…)`; every
CO₂/Ru5P term carries a factor of NADPH (`Q`). With `Q = 0` those terms vanish identically, so the
predicted rate is **flat across CO₂/Ru5P level** while the data fall ~2.3× (CO₂, 6A) / ~1.6× (Ru5P,
7A). All conditions in each figure collapse onto one trajectory.

This is **not** the fitted `Kd_CO2` (11 mM) or the unconstrained `Kd_Ru5P`: setting `Kd_CO2 → 10⁻⁹ M`
(10⁷× tighter) leaves the CO₂ = 2200 µM rate **bit-identical** to CO₂ = 0, because the term's
numerator (`Q·R·C`) is exactly zero when `Q = 0`. `Kd_CO2` is pinned only by the reverse Villet 1972
rows (all products present); `Kd_Ru5P` is unconstrained precisely because the forward corpus never
supplies the NADPH needed to activate its terms. NADPH inhibition (Weisz 5A/5B) *is* captured (NADPH
is the last-released product ⇒ free-E competitor), and the reverse Villet 1972 panels *do* separate by
CO₂ (there NADPH/Ru5P/CO₂ are co-substrates) — the law reproduces product inhibition exactly when the
partner products are present, and only then.

**Biochemical reading:** Weisz's single-product CO₂/Ru5P inhibition is direct evidence that PGD
product release is **not strictly ordered** — CO₂ and Ru5P inhibit while binding free/partial enzyme
forms (dead-end complexes, or random release), which the strictly-ordered `:full_re` law forbids by
construction. Addressing it is a **mechanism** change, deferred to follow-up phases: **Phase 4** adds
CO₂/Ru5P free-E dead-end terms (fittable single-product inhibition, +2 params); **Phase 5** relaxes to
random product release. This limitation does not affect the forward-flux deployment target (the
forward corpus does not constrain CO₂/Ru5P free-E binding regardless), but it is the honest boundary
of the ordered law.

## Verdict

- **Fits the corpus better than `cha_base`** — CV 8–11× lower with a tighter SE, fewer parameters.
  ✅
- **Passes Holten intrinsically** — data-identifies a sub-µM competitive Kd_NADPH (0.87–2.1 µM)
  and reproduces the ~28× dilution `cha_base` suffers, **without** the literature Ki_NADPH pin that
  `cha_base` requires. ✅ (This is the validation target the law exists for.)
- **Correct product-inhibition fingerprint** — NADPH competitive vs both substrates; Ru5P/CO₂
  conditional on later products. ✅
- **Cannot fit Weisz 6A–7B single-product CO₂/Ru5P inhibition** — structural (those assays have
  NADPH = 0, so the ordered CO₂/Ru5P terms vanish); not a parameter artifact. ⚠️ Deferred to
  Phase 4 (CO₂/Ru5P dead-ends) / Phase 5 (random release). Does not affect the forward-flux target.
- **α is the one soft direction** — the corpus does not support a distinct α. **Deploy α = 1**
  (single-site RE): ~4% loss cost, rescues Kd_Ru5P, lands Km_NADP in the human band. Keep α-free
  as a diagnostic only.
- **Km_PGA is not data-identified forward** (corpus lacks sub-Km [6PG]); it is set by the mode-2/3
  anchor (60/38 µM), exactly as for `cha_base`. Not a regression.

**Bottom line:** the fully-RE law is at least competitive with `cha_base` on every fit metric and
strictly better on generalization and on the NADPH-inhibition physics it was built for. The
recommended deployable form is **α = 1, modes 2/3** (Km_PGA-anchored), fiber-free.

## Deployment note (still gated — separate follow-up)

Deployment into `PentosePhosphatePathway.jl` remains **out of scope and gated**. The Phase-3 run
emits the deployable inputs: the fiber-free `:full_re` `micro_parameters.jl` deploy block (6 RE
constants per mode, no `koff`/`kon`; `K_NADP_EPGA = α·Kd_NADP` by detailed balance), the α verdict,
and the Holten / CV head-to-head above. Remaining gates before wiring into the coupled ODE:
(1) resolve **α = 1 vs α-free against the flux gate** (this evaluation recommends α = 1);
(2) the PGD Km_PGA anchor is a literature imposition, not a data identification — the coupled-flux
sensitivity to 38 vs 60 µM must be checked; (3) the over-reduction axis is entangled with the
G6PD_Vmax issue (spec §Scope). Those are scoped separately when the over-reduction axis is ready.

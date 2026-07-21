# Parsimony-first reframing of diagnostic language (docs/prose only)

## Motivation

James is shifting G6PD (and likely other enzymes') modeling philosophy away from strict
mechanistic-consensus derivation (fit every micro-constant, expect it to match
historical/literature assay values, treat disagreement as a problem to explain away) and
toward **parsimony + coupled-ODE-simulation accuracy as the primary model-selection
criterion**. If the simplest law that gives accurate simulated flux disagrees with a
historical assay value, that may mean the assay value isn't what the live parameter
should be — not that the simpler law is wrong.

The motivating case (2026-07-20 session): three G6PD dead-end-dropped variants were
fit and compared (`results/G6PD_deadend_variants_report.md`). One (`no_g6p_atp_deadend`)
is a clean parsimony win — forward constants unchanged, CV improves an order of
magnitude — despite the dropped term's literature Ki disagreeing with nothing in
particular; it's just an unconstrained nuisance dimension. The other two are NOT
recommended, because dropping them destabilizes `alpha`'s identifiability and CV,
independent of literature agreement either way.

Some of the package's existing language reads as if literature/mechanistic disagreement
is inherently a defect. This spec updates that language in three specific places. **No
computation, banner logic, CSV schema, or test changes** — this is prose only, confirmed
against `test/` (no test asserts on the exact wording being changed).

An adjacent idea — a built-in ablation-grid runner/leaderboard to automate the
multi-`run_all`-calls-plus-hand-diffing process used to produce the deadend-variants
report — was discussed and explicitly **descoped**; not pursued in this round.

## Change 1 — `unconstrained` classification: candidate-for-ablation, not defect-or-verdict

**Where:** `AGENTS.md`, wherever `unconstrained`/railing is discussed for macro
coordinates (the `Ki_NADPH` conflation section at minimum; any general description of
the `data_identified`/`unconstrained`/`literature_pinned`/`derived` classification).

**What changes:** state explicitly that an `unconstrained` macro coordinate means the
data don't pin it — that's a candidate worth testing by ablation (refit without it,
compare CV/identifiability/coupled flux against the full law), **not** proof either way
on its own. Explicitly warn against assuming `unconstrained`-in-isolation predicts
drop-safety: cite the 2026-07-20 finding that `Ki_ATP_EG` (unconstrained) was safe to
drop, while dropping `Ki_NADPH` (`data_identified`, not unconstrained) destabilized
`alpha`'s identifiability and CV — the two are not correlated 1:1.

Suggested language (adapt in place, don't necessarily quote verbatim):

> An `unconstrained` macro coordinate means the data do not pin this term — it is a
> candidate for ablation (refit without it, compare CV/identifiability/coupled flux),
> not evidence the term is safe to drop on its own. Whether dropping it costs anything
> can only be answered by actually dropping it and checking; classification status
> alone does not predict drop-safety across terms.

## Change 2 — `anchor_reverse=false` / "NOT DEPLOYABLE" banner: keep the hard flag, fix the stated reason

**Where:** `src/run.jl` (the `micro_parameters.jl` header comment block, ~line 393-397;
the `report.md` blockquote, ~line 426-429) and `AGENTS.md`'s description of the
`anchor_reverse` switch (~line 200-207).

**What changes:** wording only, not the banner's force or trigger condition. This case
is NOT a parsimony trade-off — `anchor_reverse=false` doesn't reduce parameter count, it
removes an identifiability anchor, leaving `Ki_NADPH` railed/structurally undetermined
at the *same* parameter count. That's a distinct failure mode from "simpler law, still
explains the data" and stays unconditionally non-deployable regardless of modeling
philosophy. Reword so the stated reason reads as an identifiability failure, not a
literature-disagreement complaint — e.g. swap "conflated diagnostic" framing for
language like "Ki_NADPH is structurally undetermined without this anchor — this is not
a simpler model, it is a non-identifiable one."

## Change 3 — Mode1/Mode2 framing: formalize mode1-as-deploy-default

**Where:** `AGENTS.md`'s "Modes" section (~line 184-224).

**What changes:** the section already documents an emerging convention for the
`:no_atp` law ("mode 1 is the deploy source of truth; mode 2/3 are diagnostics/
validation only"). Generalize this from a one-off convention to the standing rule for
all *future* law re-derivations:

- Mode 1 (free fit) is the deploy source of truth going forward.
- Mode 2 (and mode 3 for PGD) are literature-sensitivity diagnostics only — never an
  independent deploy path.
- The deploy gate is mode1's fit + coupled-flux validation (in `PPP_Experiments`,
  downstream of this package) — not literature agreement.
- Explicitly grandfather the currently-deployed full G6PD law (mode2-pinned
  `Ki_ATP`/`Ki_NADPH`) and PGD (mode2/3-pinned) — this does **not** retroactively
  un-deploy anything; they keep their existing status until themselves re-derived.

## Out of scope

- Any change to `cha_classify.jl`, `cha_fit.jl`, CSV column schemas, or banner
  *trigger* logic.
- The ablation-grid runner / AIC-BIC leaderboard / new CLI-API surface — descoped,
  not deferred to a specific future date; revisit only if it comes up again.
- Retroactively changing the deploy status of any currently-deployed law.

## Testing

None needed beyond confirming the existing test suite still passes unchanged (no test
asserts on the wording being changed, confirmed via
`grep -rn "NOT DEPLOYABLE\|conflated diagnostic\|CONFLATED DIAGNOSTIC" test/` returning
nothing). This is a documentation/comment-string change; standard `Pkg.test()` is
sufficient verification that nothing else broke.

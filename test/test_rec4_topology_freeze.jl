using FitRateEquation
using EnzymeRates
using Test

# FROZEN DECISION (panel review 2026-06-07): adding steady-state release steps does NOT
# free forward Ki_NADPH — the blocker is the shared E·NADPH symbol, not SS-step count.
# This test codifies WHY: in V2 the bare [NADPH] denominator term is governed by the SAME
# micro symbol that defines the reverse Km_NADPH (the NADPH-release step constant). The
# only de-conflation route is a DECOUPLED forward symbol (V3), not more SS steps.
@testset "REC-4: V2 [NADPH] coeff and reverse Km_NADPH share the NADPH-release symbol" begin
    v2 = FitRateEquation.consensus_variants(:PGD)[2].mech
    _, den = EnzymeRates._raw_symbolic_rate_polys(typeof(v2))
    metset = Set(EnzymeRates.metabolites(v2))
    # The bare-[NADPH]^1 denominator monomials (exactly one concentration symbol: NADPH).
    bare_nadph = [mono for mono in keys(den)
                  if Set(s for (s,_e) in mono if s in metset) == Set([:NADPH])]
    @test !isempty(bare_nadph)
    # Rate-constant symbols appearing in those monomials.
    ksyms = Set(s for mono in bare_nadph for (s,_e) in mono if !(s in metset))
    # The NADPH-release step (upstream canonicalizes it as binding, E + NADPH -> E(NADPH),
    # so NADPH lands on lhs_mets). Its rate constants are koff_NADPH_E / kon_NADPH_E.
    steps = _mechanism_steps(v2)
    rel = steps[findfirst(st -> :NADPH in st.rhs_mets || :NADPH in st.lhs_mets, steps)]
    @test rel.is_eq == false   # V2 release is steady-state (the DOF that still didn't help)
    # EXACT shared symbol: the bare-[NADPH] term carries the release step's REVERSE (rebind)
    # constant kon_NADPH_E — the NADPH-rebinds-free-E rate that also defines the reverse
    # Km_NADPH (NOT a coincidental substring match). We assert the SS-release rebind symbol is
    # present and the forward release rate koff_NADPH_E is NOT what carries the term, so the
    # bare-[NADPH] coeff is governed by the reverse-Km micro constant.
    @test :kon_NADPH_E in ksyms                            # rebind (reverse) release constant carries it
    @test :koff_NADPH_E ∉ ksyms                            # not the forward release rate
end

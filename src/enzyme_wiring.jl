# Per-enzyme wiring registry. Each enzyme contributes its consensus mechanism
# variants, literature pin values, pin table, and the Ki/Kd/Km alias maps. The
# generic machinery (mechanisms/pins/macro_collect/identifiability/run) reads the
# enzyme-specific DATA from here, dispatched by an `enzyme::Symbol` that `fit_consensus_equation`
# derives from `cfg.name`. Bare calls (tests) infer the enzyme from the mechanism's
# signature metabolite via `_enzyme_of`.

# A monomial key (e.g. [:NADP => 1]) labels a single-metabolite denominator term.
const MonoKey = Vector{Pair{Symbol,Int}}

struct EnzymeWiring
    name::Symbol                                   # :G6PD / :PGD
    signature::Symbol                              # metabolite unique to this enzyme
    variants::Vector{<:NamedTuple}                 # [(name::Symbol, mech), ...]
    lit_values::Dict{Symbol,Float64}               # macro name => log10(M)
    pin_table::Dict{Symbol,Vector{Symbol}}         # variant => [macro names pinned in Mode 2]
    ki_map::Dict{MonoKey,Symbol}                   # [X]^1 => Ki_X (inhibitors)
    kd_map::Dict{MonoKey,Symbol}                   # [X]^1 => Kd_X (substrates, binary)
    km_ratio::Dict{MonoKey,Tuple{Symbol,MonoKey}}  # substrate => (Km name, cosubstrate mono)
    substrate_pair::MonoKey                        # sorted substrate-pair monomial
    # macro name => (carrier_mono, cross_mono) for a noncompetitive dead-end Ki read as
    # the cross-term ratio g[carrier]/g[carrier·ligand] (e.g. G6PD Ki_NADPH on E·G6P).
    ki_ratio::Dict{Symbol,Tuple{MonoKey,MonoKey}}
    # macro name => (ligand, free-form) for a dead-end Ki read MICRO-DIRECT (1/K{idx}),
    # used when the dead-end's bare [ligand] term lumps with another denominator term.
    ki_micro_direct::Dict{Symbol,Tuple{Symbol,Symbol}}
    # variant => macro names that are structurally/empirically conflated with a reverse
    # constant (reported :conflated_reverse, not data_identified).
    conflated::Dict{Symbol,Set{Symbol}}
end

# Outer constructor: 9-positional call sites keep working; the three kwarg maps
# (ki_ratio, ki_micro_direct, conflated) default to empty, so any call site that omits
# them is unaffected. (G6PD passes `ki_ratio` for the Ki_NADPH dead-end cross term; see
# enzymes/g6pd.jl.)
EnzymeWiring(name, signature, variants, lit_values, pin_table, ki_map, kd_map, km_ratio,
             substrate_pair;
             ki_ratio        = Dict{Symbol,Tuple{MonoKey,MonoKey}}(),
             ki_micro_direct = Dict{Symbol,Tuple{Symbol,Symbol}}(),
             conflated       = Dict{Symbol,Set{Symbol}}()) =
    EnzymeWiring(name, signature, variants, lit_values, pin_table, ki_map, kd_map, km_ratio,
                 substrate_pair, ki_ratio, ki_micro_direct, conflated)

const ENZYMES = Dict{Symbol,EnzymeWiring}()
register_enzyme!(w::EnzymeWiring) = (ENZYMES[w.name] = w)
_wiring(enz::Symbol) = ENZYMES[enz]

# Identify which registered enzyme a mechanism belongs to, by signature metabolite.
function _enzyme_of(mech)
    mets = Set(EnzymeRates.metabolites(mech))
    for (name, w) in ENZYMES
        w.signature in mets && return name
    end
    error("no registered enzyme matches mechanism metabolites $(collect(mets))")
end

# Accessors (thin; keep call sites readable).
_lit_values(enz::Symbol) = _wiring(enz).lit_values
_pin_table(enz::Symbol, variant::Symbol) = get(_wiring(enz).pin_table, variant, Symbol[])

# Public, enzyme-aware. `consensus_variants()` defaults to G6PD for back-compat with
# the serial-baseline test and the task-generation test, which call it with no arg.
consensus_variants(enz::Symbol=:G6PD) = _wiring(enz).variants
v1_mechanism() = _wiring(:G6PD).variants[1].mech
v2_mechanism() = _wiring(:G6PD).variants[2].mech

# Back-compat: `pin_table(variant)` (exported) resolves against G6PD.
pin_table(variant::Symbol) = _pin_table(:G6PD, variant)

# Per-(enzyme, variant) fit-mode set: Mode 1 (free) + Mode 2 (literature-pinned). Mode 2 runs
# only when it pins something Mode 1 does not -- otherwise it is a redundant refit with identical
# pins (e.g. PGD :full_re carries no forward Ki_NADPH coord to pin). Resolved at the default
# anchor_reverse/extra/scale: the pins Mode 2 adds are forward Ki's those settings never touch.
# PGD's Mode 3 (the Km_PGA override) was retired 2026-10-21.
function modes_for(enzyme::Symbol, variant::Symbol=:_deploy)
    p1 = keys(ChaFit.resolve_cha_pins(enzyme, variant, :mode1))
    p2 = keys(ChaFit.resolve_cha_pins(enzyme, variant, :mode2))
    issubset(p2, p1) ? (:mode1,) : (:mode1, :mode2)
end

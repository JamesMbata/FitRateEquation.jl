# promotable.jl — declarative, bounded registry of investigable slow steps per enzyme.
# G6PD: NADPH catalytic release (RE-vs-SS variant; fiber coord koffQ).
# PGD: none (NADPH release is fast >800/s; silent-variant class foreclosed).
module Promotable

const _PROMOTABLE = Dict{Symbol,Vector{Symbol}}(
    :G6PD => [:nadph_release],
    :PGD  => Symbol[],
    :HK1  => Symbol[],
)
const _FIBER = Dict{Tuple{Symbol,Symbol},Symbol}(
    (:G6PD, :nadph_release) => :koffQ,
)

promotable_steps(enz::Symbol) = _PROMOTABLE[enz]
fiber_coord(enz::Symbol, step::Symbol) = _FIBER[(enz, step)]

end # module

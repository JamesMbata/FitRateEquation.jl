# ##########################################################################################
#                 Frozen copy of the Cha loss as it stood before the speed rewrite           #
# ##########################################################################################
# Verbatim from src/cha_fit.jl at ad773b6 (0.10.0), which rebuilt everything on every call and
# dispatched at run time on every row. Do NOT tidy, speed up or re-derive this file: it is the
# reference test_loss_speed.jl compares the live loss against with `==`. Keeping the reference as
# code, compared in the same process, makes that check independent of the machine; hardcoded loss
# values would drift by ULPs across CPUs (see test_byte_identity.jl).
module FrozenChaLoss

using FitRateEquation
using FitRateEquation: Dataset, nrows, _SIGN_PENALTY
using FitRateEquation.ChaLaws
using FitRateEquation.ChaFit: CHA_KDRU_DEFAULT, CHA_DEPLOY_RELEASE_RATE, CHA_ABS_RELEASE_RATE

function cha_haldane_kr(enzyme::Symbol, coords::AbstractDict; keq::Real, release_rate::Real,
                        kf::Real, release_eq::Real = _default_release_eq(enzyme, coords),
                        variant::Symbol = :_deploy)
    if enzyme === :PGD && variant === :full_re
        return coords[:Kd_NADPH] * coords[:Kd_Ru5P] * coords[:Kd_CO2] * kf /
               (coords[:Kd_NADP] * coords[:Kd_PGA] * coords[:alpha] * keq)
    end
    if enzyme === :G6PD
        return coords[:Kd_6PGLn] * release_eq * kf /
               (coords[:Kd_NADP] * coords[:Kd_G6P] * coords[:alpha] * keq)
    elseif enzyme === :PGD
        return coords[:Kd_CO2] * coords[:Km_NADPH_rev] * release_eq * kf /
               (coords[:Kd_NADP] * coords[:Kd_PGA] * coords[:alpha] * keq)
    else
        error("cha_haldane_kr: unknown enzyme $enzyme (expected :G6PD or :PGD)")
    end
end

function _default_release_eq(enzyme::Symbol, coords::AbstractDict)
    enzyme === :G6PD && return coords[:Km_NADPH_rev]
    enzyme === :PGD  && return CHA_KDRU_DEFAULT
    error("_default_release_eq: unknown enzyme $enzyme")
end

_default_release_rate(enzyme::Symbol) =
    enzyme in (:G6PD, :PGD) ? CHA_DEPLOY_RELEASE_RATE :
    error("_default_release_rate: unknown enzyme $enzyme")

function cha_macro_tuple(enzyme::Symbol, coords::AbstractDict; keq::Real,
                         kf::Real = 1.0, Et::Real = 1.0,
                         release_rate::Real = _default_release_rate(enzyme),
                         release_eq::Real = _default_release_eq(enzyme, coords),
                         kr::Union{Nothing,Real} = nothing,
                         variant::Symbol = :_deploy)
    if enzyme === :PGD && variant === :full_re
        krv = kr === nothing ?
            cha_haldane_kr(enzyme, coords; keq=keq, release_rate=release_rate, kf=kf,
                           release_eq=release_eq, variant=variant) : kr
        tup = (; Kd_NADP  = coords[:Kd_NADP],
                 Kd_PGA   = coords[:Kd_PGA],
                 alpha    = coords[:alpha],
                 Kd_NADPH = coords[:Kd_NADPH],
                 Kd_Ru5P  = coords[:Kd_Ru5P],
                 Kd_CO2   = coords[:Kd_CO2],
                 kf = kf, kr = krv, Et = Et, Keq = keq)
        for s in (:Ki_ATP, :Ki_ATP_EN, :Ki_NADPH)
            haskey(coords, s) && (tup = merge(tup, NamedTuple{(s,)}((coords[s],))))
        end
        return tup
    end
    krv = kr === nothing ?
        cha_haldane_kr(enzyme, coords; keq=keq, release_rate=release_rate, kf=kf,
                       release_eq=release_eq) : kr
    if enzyme === :G6PD
        koffQ = release_rate
        konQ  = koffQ / release_eq
        return (; Kd_NADP   = coords[:Kd_NADP],
                  Kd_G6P    = coords[:Kd_G6P],
                  Kd_6PGLn  = coords[:Kd_6PGLn],
                  alpha     = coords[:alpha],
                  Km_NADPH_rev = coords[:Km_NADPH_rev],
                  Ki_NADPH  = get(coords, :Ki_NADPH,  Inf),
                  Ki_ATP    = get(coords, :Ki_ATP,    Inf),
                  Ki_ATP_EG = get(coords, :Ki_ATP_EG, Inf),
                  koffQ = koffQ, konQ = konQ, kf = kf, kr = krv, Et = Et,
                  Keq = keq)
    elseif enzyme === :PGD
        koff = release_rate
        kon  = koff / release_eq
        return (; Kd_NADP   = coords[:Kd_NADP],
                  Kd_PGA    = coords[:Kd_PGA],
                  alpha     = coords[:alpha],
                  Kd_CO2    = coords[:Kd_CO2],
                  Km_NADPH_rev = coords[:Km_NADPH_rev],
                  Ki_NADPH  = coords[:Ki_NADPH],
                  Ki_ATP    = get(coords, :Ki_ATP,    Inf),
                  Ki_ATP_EN = get(coords, :Ki_ATP_EN, Inf),
                  koff = koff, kon = kon, kf = kf, kr = krv, Et = Et,
                  Keq = keq)
    else
        error("cha_macro_tuple: unknown enzyme $enzyme (expected :G6PD or :PGD)")
    end
end

_cha_field(cc, s) = hasproperty(cc, s) ? getproperty(cc, s) : 0.0
function _cha_row_kwargs(enzyme::Symbol, cc)
    if enzyme === :G6PD
        return (NADP = _cha_field(cc, :NADP), G6P = _cha_field(cc, :G6P),
                NADPH = _cha_field(cc, :NADPH), PGLn = _cha_field(cc, :PGLn),
                ATP = _cha_field(cc, :ATP))
    elseif enzyme === :PGD
        return (NADP = _cha_field(cc, :NADP), PGA = _cha_field(cc, :PGA),
                Ru5P = _cha_field(cc, :Ru5P), CO2 = _cha_field(cc, :CO2),
                NADPH = _cha_field(cc, :NADPH), ATP = _cha_field(cc, :ATP))
    else
        error("_cha_row_kwargs: unknown enzyme $enzyme")
    end
end

function _cha_row_logratios!(logratio, enzyme::Symbol, mech, d::Dataset, coords::AbstractDict;
                             keq::Union{Nothing,Real}=nothing, kf::Real = 1.0, Et::Real = 1.0,
                             release_rate::Real = _default_release_rate(enzyme),
                             release_eq::Real = _default_release_eq(enzyme, coords),
                             kr::Union{Nothing,Real} = nothing, variant::Symbol = :_deploy)
    cha_rate_enz = enzyme === :G6PD ? ChaLaws.cha_rate_G6PD :
                   enzyme === :PGD  ? (variant === :full_re ? ChaLaws.cha_rate_PGD_fullRE :
                                                              ChaLaws.cha_rate_PGD) :
                   error("_cha_row_logratios!: unknown enzyme $enzyme")
    penalty = 0.0
    groups = [findall(==(g), d.group) for g in unique(d.group)]
    for idx in groups
        keq_g = if keq === nothing
            ks = unique(d.keq[idx])
            length(ks) == 1 ||
                error("_cha_row_logratios!: keq not uniform within figure $(d.group[idx[1]]): $ks")
            ks[1]
        else
            keq
        end
        m = cha_macro_tuple(enzyme, coords; keq=keq_g, kf=kf, Et=Et,
                            release_rate=release_rate, release_eq=release_eq, kr=kr,
                            variant=variant)
        for i in idx
            vunit = cha_rate_enz(m; _cha_row_kwargs(enzyme, d.concs[i])...)
            eti = isnan(d.Et[i]) ? 1.0 : d.Et[i]
            v = eti * vunit
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

function cha_centered_logratio_loss(enzyme::Symbol, mech, d::Dataset,
                                    coords::AbstractDict; keq::Union{Nothing,Real}=nothing,
                                    kf::Real = 1.0, Et::Real = 1.0,
                                    release_rate::Real = _default_release_rate(enzyme),
                                    release_eq::Real = _default_release_eq(enzyme, coords),
                                    kr::Union{Nothing,Real} = nothing,
                                    variant::Symbol = :_deploy)
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

function cha_absolute_logratio_loss(enzyme::Symbol, mech, d::Dataset, coords::AbstractDict;
        keq::Union{Nothing,Real}=nothing, kf::Real=1.0,
        release_rate::Real = CHA_ABS_RELEASE_RATE,
        release_eq::Real = _default_release_eq(enzyme, coords),
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

# `release_rate` is not in the original signature: it mirrors cha_koffq_report's objective, which
# passes it straight through to the centered loss.
function _cha_loss_with_pins(enzyme, mech, d, u, coords_syms, pins;
                             keq::Union{Nothing,Real}=nothing, variant::Symbol=:_deploy,
                             scale::Symbol=:relative,
                             release_rate::Real=_default_release_rate(enzyme))
    if !isempty(pins)
        u = collect(u)
        for (idx, k) in enumerate(coords_syms)
            haskey(pins, k) && (u[idx] = pins[k])
        end
    end
    coords_dict = Dict(coords_syms .=> 10 .^ u)
    if scale === :absolute
        kcat = pop!(coords_dict, :kcat)
        return cha_absolute_logratio_loss(enzyme, mech, d, coords_dict;
                   keq=keq, kf=kcat, release_rate=CHA_ABS_RELEASE_RATE, variant=variant)
    end
    cha_centered_logratio_loss(enzyme, mech, d, coords_dict; keq=keq, variant=variant,
                               release_rate=release_rate)
end

end # module FrozenChaLoss

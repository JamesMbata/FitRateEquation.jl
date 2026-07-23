using FitRateEquation
using EnzymeRates
using Test

using FitRateEquation.ChaLaws

# A representative fully-RE macro tuple. kr is set per-test (Haldane or arbitrary).
_base_tuple(; kr) = (; Kd_NADP=1e-5, Kd_PGA=4e-5, alpha=1.4, Kd_NADPH=1e-6,
                       Kd_Ru5P=5e-5, Kd_CO2=1e-4, kf=1.0, kr=kr, Et=1.0, Keq=0.079)

@testset "fully-RE law: forward-only reduces to RE-random bireactant" begin
    m = _base_tuple(kr=0.0)                     # kr irrelevant when products are 0
    A, B = 7e-6, 3e-5
    gAB = A*B/(m.alpha*m.Kd_NADP*m.Kd_PGA)
    expected = m.Et*m.kf*gAB / (1 + A/m.Kd_NADP + B/m.Kd_PGA + gAB)
    got = cha_rate_PGD_fullRE(m; NADP=A, PGA=B)   # Ru5P=CO2=NADPH=0
    @test isapprox(got, expected; rtol=1e-12)
end

@testset "fully-RE law: Haldane v=0 at equilibrium" begin
    # kr = kf*Kq*Kr*Kc/(Keq*alpha*Ka*Kb); then v=0 iff Q*R*C = Keq*A*B.
    mt0 = _base_tuple(kr=0.0)
    kr  = mt0.kf*mt0.Kd_NADPH*mt0.Kd_Ru5P*mt0.Kd_CO2 /
          (mt0.Keq*mt0.alpha*mt0.Kd_NADP*mt0.Kd_PGA)
    m   = _base_tuple(kr=kr)
    A, B = 1e-5, 4e-5                            # A*B = 4e-10
    Q, R = 1e-6, 1e-3                            # Q*R = 1e-9
    Cc  = mt0.Keq*A*B / (Q*R)                    # => Q*R*Cc = Keq*A*B
    v   = cha_rate_PGD_fullRE(m; NADP=A, PGA=B, NADPH=Q, Ru5P=R, CO2=Cc)
    vfwd = cha_rate_PGD_fullRE(m; NADP=A, PGA=B) # characteristic scale
    @test abs(v) < 1e-10*abs(vfwd)
end

@testset "fully-RE law: alpha=1 => apparent Km_NADP == Kd_NADP" begin
    m = _base_tuple(kr=0.0)
    m = merge(m, (; alpha=1.0))
    Bsat = 1.0                                   # saturating PGA (>> Kd_PGA)
    Vmax_app = cha_rate_PGD_fullRE(m; NADP=1.0, PGA=Bsat)          # NADP saturating too
    half     = cha_rate_PGD_fullRE(m; NADP=m.Kd_NADP, PGA=Bsat)   # A = Kd_NADP
    @test isapprox(half, Vmax_app/2; rtol=1e-3)
end

using FitRateEquation.ChaInvert

# The registered fully-RE PGD mechanism (:RE_rate_eq = V1): random RE substrate binding,
# catalysis SS, CO2/Ru5P/NADPH releases all RE. Carries both ATP dead-ends, no NADPH dead-end.
function _pgd_re_mech()
    vs = FitRateEquation.consensus_variants(:PGD)
    vs[findfirst(v -> Symbol(v.name) === :RE_rate_eq, vs)].mech
end

@testset "exactness gate: cha_rate_PGD_fullRE == rate_equation (:RE_rate_eq, rtol 1e-10)" begin
    m    = _pgd_re_mech()
    mets = EnzymeRates.metabolites(m)           # (:NADP,:PGA,:CO2,:NADPH,:Ru5P,:ATP)
    grid = [
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=0.0,  NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=5e-5, NADPH=0.0,  ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=0.0,  Ru5P=0.0,  NADPH=5e-6, ATP=0.0),
        (; NADP=5e-6, PGA=40e-6, CO2=1e-4, Ru5P=5e-5, NADPH=5e-6, ATP=0.0),
        (; NADP=2e-5, PGA=80e-6, CO2=2e-4, Ru5P=1e-4, NADPH=8e-6, ATP=5e-4),
    ]
    for _ in 1:20
        free = free_params(m)
        logθ = -3 .+ 2 .* rand(length(free))
        keq  = 0.079
        mac  = cha_macro_readoffs_PGD_fullRE(m, logθ; keq=keq)
        # Characteristic forward rate: floor for the rel-error denominator.
        vsat = abs(EnzymeRates.rate_equation(m,
            NamedTuple{Tuple(mets)}(Tuple(s in (:NADP,:PGA) ? 1e-2 : 0.0 for s in mets)),
            build_params(m, logθ; keq=keq)))
        for conc in grid
            cc   = NamedTuple{Tuple(mets)}(Tuple(getfield(conc, s) for s in mets))
            vref = EnzymeRates.rate_equation(m, cc, build_params(m, logθ; keq=keq))
            vcha = cha_rate_PGD_fullRE(mac; conc...)
            @test isapprox(vcha, vref; rtol=1e-10, atol=1e-10 * vsat)
        end
    end
end

using FitRateEquation
using FitRateEquation: v2_mechanism
using EnzymeRates
using Test

using FitRateEquation.ChaLaws
using FitRateEquation.ChaInvert

@testset "read-offs match analytic apparent-Km formulas (corrected C=1+kf/koffQ)" begin
    m = v2_mechanism()
    for _ in 1:10
        free = free_params(m)
        logθ = -3 .+ 2 .* rand(length(free))
        keq  = 10.0
        mac  = cha_macro_readoffs_G6PD(m, logθ; keq=keq)
        C = 1 + mac.kf / mac.koffQ                       # koffQ alone (two-SS topology)
        @test isapprox(mac.Km_NADP_apparent, mac.alpha*mac.Kd_NADP/C; rtol=1e-9)
        @test isapprox(mac.Km_G6P_apparent,  mac.alpha*mac.Kd_G6P/C;  rtol=1e-9)
        @test isapprox(mac.Km_NADPH_rev, mac.koffQ/mac.konQ; rtol=1e-12)
    end
end

@testset "macro->micro->law: Kd & alpha/C invariant under koffQ sweep" begin
    # The fiber map holds the forward-identifiable macro constants (Kd's, Km_NADPH_rev, and
    # alpha/C => apparent Km's) invariant by construction, regardless of koffQ.
    m = v2_mechanism()
    free = free_params(m)
    logθ = -3 .+ 2 .* rand(length(free))
    keq  = 10.0
    mac  = cha_macro_readoffs_G6PD(m, logθ; keq=keq)
    C_ref = 1 + mac.kf / mac.koffQ
    alphaC_ref = mac.alpha / C_ref
    for koffQ in (1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3)
        fib = cha_micro_from_macro_G6PD(mac; koffQ=koffQ)
        @test isapprox(fib.Kd_NADP,      mac.Kd_NADP;      rtol=1e-12)
        @test isapprox(fib.Kd_G6P,       mac.Kd_G6P;       rtol=1e-12)
        @test isapprox(fib.Kd_6PGLn,     mac.Kd_6PGLn;     rtol=1e-12)
        @test isapprox(fib.Km_NADPH_rev, mac.Km_NADPH_rev; rtol=1e-12)
        @test isapprox(fib.Ki_ATP_EG,    mac.Ki_ATP_EG;    rtol=1e-12)
        C_fib = 1 + fib.kf / fib.koffQ
        @test isapprox(fib.alpha / C_fib, alphaC_ref; rtol=1e-10)
    end
end

@testset "koffQ is data-identified (NOT silent) on the corrected two-SS topology" begin
    # CORRECTED-TOPOLOGY DESIGN FINDING (cha_derive_g6pd.py Property 5): koffQ is silent
    # ONLY along the pure-Vmax gauge (uniform rescale of kf,kr,koffQ,konQ). The
    # `cha_micro_from_macro_G6PD` fiber holds the APPARENT macro constants (Kd's, alpha/C,
    # Km_NADPH_rev) invariant -- but it does NOT hold the Vmax-shape kf/alpha, so the
    # observable RATE moves as koffQ slides. The catalysis-reverse kr populating E_C makes
    # koffQ identifiable from PGLn-product-inhibition / both-product data. This testset
    # asserts the (verified) fact that the observable rate is NOT koffQ-invariant under the
    # macro-preserving fiber -- i.e. koffQ is data-identified, not a free gauge.
    m = v2_mechanism()
    free = free_params(m)
    logθ = -3 .+ 2 .* rand(length(free))
    keq  = 10.0
    mac  = cha_macro_readoffs_G6PD(m, logθ; keq=keq)

    pts = [(; NADP=5e-6, G6P=40e-6, NADPH=0.0,  PGLn=1e-4, ATP=0.0),   # PGLn inhibition
           (; NADP=5e-6, G6P=40e-6, NADPH=5e-6, PGLn=1e-4, ATP=0.0)]   # both products
    vref = [cha_rate_G6PD(mac; p...) for p in pts]

    moved = false
    for koffQ in (1e-2, 1e-1, 1.0, 1e1, 1e2)
        fib = cha_micro_from_macro_G6PD(mac; koffQ=koffQ)
        for (i, p) in enumerate(pts)
            isapprox(cha_rate_G6PD(fib; p...), vref[i]; rtol=1e-3) || (moved = true)
        end
    end
    @test moved      # rate changes with koffQ => koffQ is data-identified (not silent)
end

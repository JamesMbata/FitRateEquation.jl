# =========================================================================================
#           EnzymeRates name layer: the installed EnzymeRates replays the 0.5.0 oracle
# =========================================================================================
# FRE keeps EnzymeRates 0.5.0's composition-semantic parameter names (K_<ligand>_<form>,
# koff_/kon_<ligand>_<form>, k_<from>_to_<to>) as its own vocabulary; core/er_names.jl
# translates the installed EnzymeRates' names to them at the boundary. This replays
# fixtures/er050_oracle.toml (generated on 0.5.0 by fixtures/gen_er050_oracle.jl): every
# registered mechanism must show the same names, order, step view and rate values.
using FitRateEquation
using EnzymeRates
using TOML
using Test

const _ORACLE = TOML.parsefile(joinpath(@__DIR__, "fixtures", "er050_oracle.toml"))

@testset "EnzymeRates name layer replays the 0.5.0 oracle" begin
    @test _ORACLE["enzymerates_version"] == "0.5.0"
    for rec in _ORACLE["mechanisms"]
        enz, var = Symbol(rec["enzyme"]), Symbol(rec["variant"])
        m = only(filter(v -> v.name === var, FitRateEquation.consensus_variants(enz))).mech
        @testset "$(enz) $(var)" begin
            @test string.(collect(FitRateEquation.er_fitted_params(m))) == rec["fitted_params"]
            @test string.(collect(FitRateEquation.er_parameters(m))) == rec["parameters"]
            dep, indep = FitRateEquation.er_dependent_param_exprs(m)
            @test sort(string.(collect(keys(dep)))) == rec["dependent_params"]
            @test string.(collect(indep)) == rec["fitted_params"]
            @test string.(collect(EnzymeRates.metabolites(m))) == rec["metabolites"]
            @test string.(collect(FitRateEquation.er_substrates(m))) == rec["substrates"]
            steps = FitRateEquation._mechanism_steps(m)
            @test length(steps) == length(rec["steps"])
            for (st, o) in zip(steps, rec["steps"])
                @test string.(st.lhs) == o["lhs"] && string.(st.rhs) == o["rhs"]
                @test st.is_eq == o["is_eq"] && st.group == o["group"] &&
                      st.param_index == o["param_index"]
            end
            mets = Tuple(Symbol.(rec["metabolites"]))
            fp = Tuple(Symbol.(rec["fitted_params"]))
            maxrel = 0.0
            for pt in rec["points"]
                p = NamedTuple{(fp..., :Keq, :E_total)}(
                    (map(s -> pt["params"][string(s)], fp)..., pt["keq"], 1.0))
                c = NamedTuple{mets}(map(s -> pt["concs"][string(s)], mets))
                v = FitRateEquation.er_rate_equation(m, c, p)
                maxrel = max(maxrel, abs(v - pt["rate"]) / abs(pt["rate"]))
            end
            @test maxrel <= 1e-12
        end
    end
end

@testset "legacy_param_name translates the structural scheme" begin
    ln = FitRateEquation.legacy_param_name
    @test ln(:K_EG6P_to_E_G6P) === :K_G6P_E                       # RE binding (dissociation)
    @test ln(:K_EG6PNADP_to_ENADP_G6P) === :K_G6P_ENADP
    @test ln(:K_EATPinh_to_E_ATPinh) === :K_ATPinh_E               # dead-end inhibitor copy
    @test ln(:k_ENADPH_to_E_NADPH) === :koff_NADPH_E               # SS release (off)
    @test ln(:k_E_NADPH_to_ENADPH) === :kon_NADPH_E                # SS binding (on)
    @test ln(:k_EG6PNADP_to_ENADPHPGLn) === :k_EG6PNADP_to_ENADPHPGLn  # chemistry: unchanged
    @test ln(:Keq) === :Keq && ln(:E_total) === :E_total
    @test_throws ErrorException ln(:K_ES_to_EP)                    # RE isomerization: no 0.5.0 form
end

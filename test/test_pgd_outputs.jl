using FitRateEquation
using FitRateEquation: pgd_config
using CSV, DataFrames
using Test

@testset "PGD fit_consensus_equation outputs (Cha macro-coord, smoke)" begin
    outdir = mktempdir()
    res = FitRateEquation._fit_consensus(pgd_config(); outdir=outdir, n_restarts=2, maxiter=150, maxtime=120.0, seed=1)
    for f in ("macro_constants.csv","goodness_of_fit.csv","fit_corpus.csv",
              "identifiable_functions.csv","micro_parameters.jl","report.md")
        @test isfile(joinpath(outdir, f))
    end
    # Deploy-variant-only: PGD = 1 variant (:cha_base) × 2 modes (mode1/mode2) = 2 cells.
    @test length(res) == 2
    @test all(r -> r.variant === :cha_base, res)
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    @test Set(string.(mc.mode)) == Set(["mode1","mode2"])
    # Forward shape constants: Cha coords (Kd_*/alpha/Ki_*/Km_NADPH_rev) plus the DERIVED
    # apparent Km_PGA / Km_NADP readoffs. Ki_NADPH is a Cha coord on cha_base.
    names = Set(string.(mc.name))
    @test "Km_PGA" in names && "Kd_NADP" in names && "Kd_PGA" in names && "Ki_NADPH" in names
    # Classes are the Cha macro vocabulary (derived apparent Km's carry :derived).
    @test issubset(Set(string.(mc.class)),
        Set(["data_identified","unconstrained","literature_pinned","derived"]))
    # Mode 2 pins ONLY Ki_NADPH: Ki_ATP is never literature-pinned (the corpus has no ATP).
    @test !any((mc.name .== "Ki_ATP") .& (string.(mc.class) .== "literature_pinned"))
    # Km_PGA is the DERIVED apparent constant (alpha·Kd_PGA/C), never pulled or pinned.
    km_pga = mc[mc.name .== "Km_PGA", :]
    @test nrow(km_pga) == 2 && all(string.(km_pga.class) .== "derived")
    # Forward Ki_NADPH is a Cha coord: :literature_pinned (17 µM) in Mode 2, data-driven in
    # Mode 1, NEVER :conflated_reverse.
    @test !any(string.(mc.class) .== "conflated_reverse")
    ki_nadph_m2 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_nadph_m2) >= 1 && all(string.(ki_nadph_m2.class) .== "literature_pinned")
    @test all(isapprox.(ki_nadph_m2.value, 17e-6; rtol=1e-6))
    @test all(in.(string.(mc[mc.name .== "Ki_NADPH", :].class),
                  Ref(Set(["data_identified","unconstrained","literature_pinned"]))))
    # Report carries the PGD title, the single mode1↔mode2 agreement, and the note; the
    # retired Km_PGA-gap block (a mode1-vs-anchored-mode2 contrast) is gone.
    rep = read(joinpath(outdir, "report.md"), String)
    @test occursin("(PGD)", rep)
    @test occursin("### mode1 ↔ mode2", rep)
    @test !occursin("mode3", rep)
    @test !occursin("Km_PGA gap", rep)
    @test occursin("FROZEN: do not add SS steps", rep)
    # micro_parameters.jl carries the closed-form deploy block for the deploy variant.
    mp = read(joinpath(outdir, "micro_parameters.jl"), String)
    @test occursin("cha_base", mp)
    @test occursin("cha_deploy_micro", mp)
    @test !occursin("deploy unavailable", mp)
end

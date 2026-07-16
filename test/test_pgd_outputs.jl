using FitRateEquation
using CSV, DataFrames
using Test

@testset "PGD run_all outputs (Cha macro-coord, smoke)" begin
    outdir = mktempdir()
    res = run_all(pgd_config(); outdir=outdir, n_restarts=2, maxiter=150, maxtime=120.0, seed=1)
    for f in ("macro_constants.csv","goodness_of_fit.csv",
              "identifiable_functions.csv","micro_parameters.jl","report.md")
        @test isfile(joinpath(outdir, f))
    end
    # Deploy-variant-only: PGD = 1 variant (:cha_base) × 3 modes (mode1/mode2/mode3) = 3 cells.
    @test length(res) == 3
    @test all(r -> r.variant === :cha_base, res)
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    @test Set(string.(mc.mode)) == Set(["mode1","mode2","mode3"])
    # Forward shape constants: Cha coords (Kd_*/alpha/Ki_*/Km_NADPH_rev) plus the DERIVED
    # apparent Km_PGA / Km_NADP readoffs. Ki_NADPH is a Cha coord on cha_base.
    names = Set(string.(mc.name))
    @test "Km_PGA" in names && "Kd_NADP" in names && "Kd_PGA" in names && "Ki_NADPH" in names
    # Classes are the Cha macro vocabulary (derived apparent Km's carry :derived).
    @test issubset(Set(string.(mc.class)),
        Set(["data_identified","unconstrained","literature_pinned","derived"]))
    # Ki_ATP literature-pinned in mode2 (hard Cha coord pin).
    ki_atp_m2 = mc[(mc.name .== "Ki_ATP") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_atp_m2) >= 1 && all(string.(ki_atp_m2.class) .== "literature_pinned")
    # Km_PGA is the DERIVED apparent constant (alpha·Kd_PGA/C); its Mode-3 hard override and
    # Mode-2 soft anchor are realized through the anchors mechanism, not as a coord class, so
    # Km_PGA rows carry :derived (a readoff), never :override.
    km_pga = mc[mc.name .== "Km_PGA", :]
    @test nrow(km_pga) >= 1 && all(string.(km_pga.class) .== "derived")
    # Forward Ki_NADPH is a Cha coord: :literature_pinned (17 µM) in Mode 2/3, data-driven in
    # Mode 1, NEVER :conflated_reverse.
    @test !any(string.(mc.class) .== "conflated_reverse")
    ki_nadph_m2 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_nadph_m2) >= 1 && all(string.(ki_nadph_m2.class) .== "literature_pinned")
    @test all(isapprox.(ki_nadph_m2.value, 17e-6; rtol=1e-6))
    @test all(in.(string.(mc[mc.name .== "Ki_NADPH", :].class),
                  Ref(Set(["data_identified","unconstrained","literature_pinned"]))))
    # Report carries the PGD title, the 3-way agreement, the Km_PGA gap warning, and the note.
    rep = read(joinpath(outdir, "report.md"), String)
    @test occursin("(PGD)", rep)
    @test occursin("Km_PGA gap", rep)
    @test occursin("WARNING — Km_PGA is NOT data-identified", rep)
    @test occursin("FROZEN: do not add SS steps", rep)
    # micro_parameters.jl carries the closed-form deploy block for the deploy variant.
    mp = read(joinpath(outdir, "micro_parameters.jl"), String)
    @test occursin("cha_base", mp)
    @test occursin("cha_deploy_micro", mp)
end

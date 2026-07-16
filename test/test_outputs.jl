using FitRateEquation
using CSV, DataFrames
using Test


@testset "run_all outputs (Cha macro-coord)" begin
    outdir = mktempdir()
    res = run_all(g6pd_config(); outdir=outdir, n_restarts=2, maxiter=150, maxtime=5.0, seed=1)
    for f in ("macro_constants.csv","goodness_of_fit.csv",
              "identifiable_functions.csv","micro_parameters.jl","report.md")
        @test isfile(joinpath(outdir, f))
    end
    # Deploy-variant-only: G6PD = 1 variant (:SS_NADPH_release_rate_eq) × 2 modes = 2 cells.
    @test length(res) == 2
    @test all(r -> r.variant === :SS_NADPH_release_rate_eq, res)
    mc = CSV.read(joinpath(outdir, "macro_constants.csv"), DataFrame)
    @test :mode in propertynames(mc)
    @test Set(string.(mc.mode)) == Set(["mode1","mode2"])
    # Forward shape constants: the Cha coords (Kd_*/alpha/Ki_*) plus the DERIVED apparent
    # Km_G6P readoff. Km_G6P in-band (~10–200 µM → between 1e-6 and 1e-2 M).
    names = Set(string.(mc.name))
    @test "Km_G6P" in names && "Kd_G6P" in names && "Kd_NADP" in names && "Ki_ATP" in names
    km = mc[(mc.name .== "Km_G6P") .& (string.(mc.mode) .== "mode1"), :]
    @test nrow(km) >= 1
    @test all(1e-6 .< km.value .< 1e-2)
    # Classes are the Cha macro vocabulary (derived apparent Km's carry :derived).
    @test issubset(Set(string.(mc.class)),
        Set(["data_identified","unconstrained","literature_pinned","derived"]))
    # Ki_ATP literature-pinned in Mode 2 (a hard Cha coord pin).
    ki_atp_m2 = mc[(mc.name .== "Ki_ATP") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_atp_m2) >= 1 && all(string.(ki_atp_m2.class) .== "literature_pinned")
    # Forward Ki_NADPH is a Cha coord (E·G6P·NADPH dead-end cross term): :literature_pinned
    # (15 µM) in Mode 2, data-driven in Mode 1; NEVER :conflated_reverse.
    ki_nadph_m2 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode2"), :]
    @test nrow(ki_nadph_m2) >= 1 && all(string.(ki_nadph_m2.class) .== "literature_pinned")
    @test all(isapprox.(ki_nadph_m2.value, 15e-6; rtol=1e-6))
    ki_nadph_m1 = mc[(mc.name .== "Ki_NADPH") .& (string.(mc.mode) .== "mode1"), :]
    @test nrow(ki_nadph_m1) >= 1
    @test all(in.(string.(ki_nadph_m1.class), Ref(Set(["data_identified","unconstrained"]))))
    @test !any(string.(mc.class) .== "conflated_reverse")
    # report.md carries the 2-way (mode1<->mode2) agreement + the G6PD koffQ hybrid block.
    rep = read(joinpath(outdir,"report.md"), String)
    @test occursin("agreement across modes", rep)
    @test occursin("mode1 ↔ mode2", rep)
    @test occursin("koffQ hybrid", rep)
    @test occursin("data-identified koffQ", rep)
    # micro_parameters.jl carries the closed-form deploy block for the deploy variant.
    mp = read(joinpath(outdir, "micro_parameters.jl"), String)
    @test occursin("SS_NADPH_release_rate_eq", mp)
    @test occursin("cha_deploy_micro", mp)
end

@testset "provenance records deploy_keq" begin
    d    = FitRateEquation.load_dataset(g6pd_config())
    meta = (n_restarts=1, maxiter=100, maxtime=1.0, seed=1,
            n_rows=FitRateEquation.nrows(d))
    dir  = mktempdir()
    FitRateEquation._write_provenance(dir, d, meta; deploy_keq=13.655)
    txt = read(joinpath(dir, "provenance.toml"), String)
    @test occursin("deploy_keq", txt)
    @test occursin("13.655", txt)
end

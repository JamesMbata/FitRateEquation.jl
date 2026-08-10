using Test, FitRateEquation, CSV, DataFrames
using FitRateEquation: g6pd_config, read_corpus

# Minimal 1-row G6PD corpus with the canonical columns.
function _write_canonical_csv(path; drop=nothing, extra=false)
    cols = Dict(
        "[NADP] (uM)" => 100.0, "[G6P] (uM)" => 200.0, "[NADPH] (uM)" => 0.0,
        "[PGLn] (uM)" => 0.0, "[ATP] (uM)" => 0.0,
        "Rate_V" => 1.0, "Article" => "T", "Fig" => "1", "Apparent_Keq" => 13.0,
    )
    drop === nothing || delete!(cols, drop)
    extra && (cols["Notes"] = "annotation")
    CSV.write(path, DataFrame(cols))
    path
end

@testset "schema validation" begin
    dir = mktempdir()

    # conforming corpus loads
    ok = _write_canonical_csv(joinpath(dir, "ok.csv"))
    @test read_corpus(g6pd_config(; data_csv=ok)) isa DataFrame

    # extra annotation column is allowed (ignored)
    ex = _write_canonical_csv(joinpath(dir, "extra.csv"); extra=true)
    @test read_corpus(g6pd_config(; data_csv=ex)) isa DataFrame

    # missing required column errors, and the message names it + the expected header
    bad = _write_canonical_csv(joinpath(dir, "bad.csv"); drop="[G6P] (uM)")
    err = try read_corpus(g6pd_config(; data_csv=bad)); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("[G6P] (uM)", err.msg)
    @test occursin("Rate_V", err.msg)   # expected header is printed
end

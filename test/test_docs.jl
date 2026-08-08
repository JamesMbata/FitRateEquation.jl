using Test

@testset "docs sanity" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    agents = read(joinpath(@__DIR__, "..", "AGENTS.md"), String)
    @test occursin("run_g6pd", readme)
    @test occursin("Pkg.add", readme)
    @test occursin("AGENTS.md", readme)
    @test !occursin("ConsensusMacro", agents)
    @test !occursin("run from the repo root", agents)   # repo-root rule removed

    # 0.2.0 artifact + contract changes must be documented, not just implemented.
    @test occursin("fit_corpus.csv", readme)
    @test occursin("fit_corpus.csv", agents)
    @test occursin("read_fit_corpus", agents)    # the helper that replaced build_plot_df
    @test occursin("seven artifacts", agents)    # was six

    # The determinism section must describe the single structural gate as implemented.
    # Anchor on phrases unique to that section: a bare "structural" needle would pass on
    # unrelated prose elsewhere in the file and would keep passing if the section vanished.
    @test occursin("variant/mode/name", agents)
    @test occursin("fitted values are not compared", agents)
    @test !occursin("FITRATEEQ_BYTE_IDENTITY", agents)   # the retired opt-in
end

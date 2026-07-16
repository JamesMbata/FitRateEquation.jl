using Test

@testset "docs sanity" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    agents = read(joinpath(@__DIR__, "..", "AGENTS.md"), String)
    @test occursin("run_g6pd", readme)
    @test occursin("Pkg.add", readme)
    @test occursin("AGENTS.md", readme)
    @test !occursin("ConsensusMacro", agents)
    @test !occursin("run from the repo root", agents)   # repo-root rule removed
end

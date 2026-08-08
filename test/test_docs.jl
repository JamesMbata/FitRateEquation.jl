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

    # 0.3.0 contract change: corpus= is required on the exported write_outputs (breaking).
    # Anchored on wording unique to that sentence — bare "write_outputs" / "corpus" needles
    # were already satisfied by the Layout file map and many unrelated corpus mentions, so they
    # would have passed before the sentence existed and protected nothing. The anchor is
    # bold-independent (it carries no `**`), so reflowing the paragraph or dropping the bold
    # does not turn it red for a purely cosmetic edit. It is NOT markup-free: the code span in
    # `corpus=` is part of the needle, and removing it would be a red the wording did not earn.
    @test occursin("0.3.0", agents)
    @test occursin("requires `corpus=`", agents)
    @test occursin("breaking change", agents)

    # 0.3.0 also changed a run-dir ARTIFACT, not just the API: fit_corpus.csv now emits its
    # metabolite columns in sorted order, so a consumer reading it positionally breaks
    # silently. Same rule as the contract break above — artifact changes must be documented.
    @test occursin("sorted alphabetically", agents)
end

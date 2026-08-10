using Test

@testset "docs sanity" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    agents = read(joinpath(@__DIR__, "..", "AGENTS.md"), String)
    @test occursin("run_g6pd", readme)
    @test occursin("Pkg.add", readme)
    @test occursin("AGENTS.md", readme)
    @test !occursin("ConsensusMacro", agents)
    @test !occursin("run from the repo root", agents)   # repo-root rule removed

    # Unified entry point: fit_consensus_equation(:enzyme; …) is THE documented entry;
    # run_g6pd/run_pgd/run_hk1 survive only as thin aliases (Task 8 of the refactor).
    @test occursin("fit_consensus_equation", readme)
    @test occursin("fit_consensus_equation", agents)

    # The per-variant convenience runners were deleted — a variant is now a keyword on
    # the single entry point, never its own function. These names must not resurface as
    # live API in the user guide (README). AGENTS.md keeps a past-tense migration note,
    # so it is deliberately exempt from the absence check.
    @test !occursin("run_g6pd_noatp", readme)
    @test !occursin("run_pgd_fullre", readme)

    # Config builders (g6pd_config/pgd_config/hk1_config) are UN-EXPORTED internals; the
    # own-data path is fit_consensus_equation(:enzyme; data_csv=…). They must not appear
    # anywhere in the human user guide as an API the reader is told to call.
    @test !occursin("g6pd_config", readme)
    @test !occursin("pgd_config", readme)
    @test !occursin("hk1_config", readme)

    # Own-data contract: data_csv= entry + canonical (non-remappable) column schema.
    @test occursin("data_csv", readme)
    @test occursin("canonical", readme)
    @test occursin("canonical", agents)

    # CLI: --variant NAME replaces the removed g6pd-noatp subcommand, and --data now works
    # for any enzyme subcommand. Both README §8 and the AGENTS.md CLI block document it.
    @test occursin("--variant", readme)
    @test occursin("--variant", agents)

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

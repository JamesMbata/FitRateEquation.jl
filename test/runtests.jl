using Test

@testset "FitRateEquation" begin
    include("test_data_path.jl")
    include("test_data.jl")
    include("test_mechanisms.jl")
    include("test_cv.jl")
    include("test_promotable.jl")
    include("test_config_deploy_keq.jl")
    # Cha core + exactness anchors
    include("test_cha_laws.jl")
    include("test_cha_invert.jl")
    include("test_cha_fit.jl")
    include("test_cha_classify.jl")
    include("test_cha_deploy.jl")
    include("test_cha_koffq.jl")
    include("test_cha_silence.jl")
    include("test_cha_noatp.jl")
    include("test_rec4_topology_freeze.jl")
    # Run / outputs / modes / parallel determinism
    include("test_run_fit.jl")
    include("test_outputs.jl")
    include("test_mode_agreement.jl")
    include("test_modes_cascade.jl")
    include("test_parallel_equivalence.jl")
    include("test_runners.jl")
    # PGD
    include("test_cha_pgd_laws.jl")
    include("test_pgd_data.jl")
    include("test_pgd_mechanisms.jl")
    include("test_pgd_outputs.jl")
    include("test_pgd_macro_collect.jl")
    # HK1 (auto-skip while HK1_AVAILABLE is false)
    include("test_cha_hk1_laws.jl")
    include("test_hk1_deploy.jl")
    include("test_hk1_fit.jl")
    # Plotting (non-render assertions; render covered in Task 8)
    include("test_plot_consensus_fit.jl")
end

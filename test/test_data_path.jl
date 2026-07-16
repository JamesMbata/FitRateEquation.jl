using Test, FitRateEquation

@testset "data path resolves from any cwd" begin
    cfg = g6pd_config()
    @test isabspath(cfg.data_csv)
    @test isfile(cfg.data_csv)
    mktempdir() do d
        cd(d) do                      # a cwd with no repo-relative data dir
            ds = load_dataset(g6pd_config())
            @test nrows(ds) > 0
        end
    end
    @test isfile(pgd_config().data_csv)
end

using FitRateEquation
using FitRateEquation: g6pd_config
using FitRateEquation.ChaFit
using FitRateEquation.ChaLaws
using Test

@testset "absolute loss: Et is a linear prefactor (log-shift)" begin
    cfg = g6pd_config(data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    m = FitRateEquation.v2_mechanism()
    # planted coords (any physically valid point); kcat as :kf via scale-absolute coords (Task 4)
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    n = FitRateEquation.nrows(d)
    lr = fill(NaN, n)
    ChaFit._cha_row_logratios!(lr, :G6PD, m, d, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE, Et=1.0)   # Et=1 baseline
    # Doubling every row's Et must shift each finite log-ratio by exactly log(2).
    d2 = Dataset(d.concs, d.rate, d.group, d.keq, 2 .* d.Et)
    lr2 = fill(NaN, n)
    ChaFit._cha_row_logratios!(lr2, :G6PD, m, d2, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
    for i in 1:n
        (isfinite(lr[i]) && isfinite(lr2[i])) || continue
        @test lr2[i] - lr[i] ≈ log(2) atol=1e-10
    end
end

@testset "absolute loss: uncentered sum of squares" begin
    cfg = g6pd_config(data_csv=joinpath(@__DIR__, "fixtures", "g6pd_abs_mini.csv"))
    d = FitRateEquation.load_dataset(cfg)
    m = FitRateEquation.v2_mechanism()
    coords = Dict(:Kd_NADP=>5e-5, :Kd_G6P=>2e-4, :Kd_6PGLn=>2e-4, :alpha=>1.0,
                  :Ki_NADPH=>2e-5, :Ki_ATP=>1.5e-3, :Ki_ATP_EG=>3e-2, :Km_NADPH_rev=>3.9e-6)
    L = ChaFit.cha_absolute_logratio_loss(:G6PD, m, d, coords; keq=13.655, kf=178.0)
    lr = fill(NaN, FitRateEquation.nrows(d))
    pen, _ = ChaFit._cha_row_logratios!(lr, :G6PD, m, d, coords; keq=13.655, kf=178.0,
        release_rate=ChaFit.CHA_ABS_RELEASE_RATE)
    expect = (pen + sum(x -> x^2, filter(isfinite, lr))) / FitRateEquation.nrows(d)
    @test L ≈ expect
end

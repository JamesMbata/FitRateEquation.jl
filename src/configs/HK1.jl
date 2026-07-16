# HK1 config for the FitRateEquation extractor. data_csv resolves to the
# bundled corpus under the package's data/ directory via pkgdir, so it works
# from any cwd; pass the data_csv keyword to point at a different corpus.
# The runner and tests call hk1_config(). The Choe HK1 corpus is already in Molar, so
# every metabolite declares unit :M (the loader divides by 1e6 only for :uM). The
# phosphate law/mechanism symbol is :Pi but the CSV column is "Phosphate".
function hk1_config(; data_csv = joinpath(pkgdir(FitRateEquation), "data", "Choe_HK1_kinetic_data.csv"))
    (
        name = "HK1",
        data_csv = data_csv,
        rate_col = "Rate",
        article_col = "Article",
        fig_col = "Fig",
        keq_col = "Apparent_Keq",
        deploy_keq = 2700.0,               # unchanged; TODO recompute at pH 7.2 when HK1 is re-fit.
        metabolites = Dict(                # symbol => (csv column, unit)
            :Glucose => ("Glucose",   :M),
            :ATP     => ("ATP",       :M),
            :G6P     => ("G6P",       :M),
            :ADP     => ("ADP",       :M),
            :Pi      => ("Phosphate", :M),
        ),
    )
end

# PGD config for the FitRateEquation extractor. data_csv resolves to the
# bundled corpus under the package's data/ directory via pkgdir, so it works
# from any cwd; pass the data_csv keyword to point at a different corpus.
# The runner and tests call pgd_config().
function pgd_config(; data_csv = joinpath(pkgdir(FitRateEquation), "data", "PGD_EnzymeData_with_CO2.csv"))
    (
        name = "PGD",
        data_csv = data_csv,
        rate_col = "Rate_V",
        article_col = "Article",
        fig_col = "Fig",
        keq_col = "Apparent Keq (M)",
        deploy_keq = 0.17,                 # pH-flat apparent Keq at 37C (Villet 38C, CO2 aq;
                                           # matches deployed PGD_Keq=0.17). READOUT + DEPLOY only.
        metabolites = Dict(                # symbol => (csv column, unit)
            :NADP  => ("[NADP] (uM)",  :uM),
            :PGA   => ("[6PGA] (uM)",  :uM),
            :CO2   => ("[CO2] (uM)",   :uM),
            :Ru5P  => ("[Ru5P] (uM)",  :uM),
            :NADPH => ("[NADPH] (uM)", :uM),
            :ATP   => ("[ATP] (uM)",   :uM),   # zero-filled placeholder (Task 2)
        ),
    )
end

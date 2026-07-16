# PGD config for the FitRateEquation extractor. Paths are repo-root-relative (the
# data loader reads data_csv relative to cwd, so run from repo root). The
# runner and tests call pgd_config().
function pgd_config()
    (
        name = "PGD",
        data_csv = "fitting/PGD/PGD_EnzymeData_with_CO2.csv",
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

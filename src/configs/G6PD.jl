# G6PD config for the FitRateEquation extractor. Paths are relative to the
# PPP_Experiments repo root (the data loader reads data_csv relative to
# cwd, so run from repo root). Loader-relevant fields mirror the mechanism_id
# G6PD config; the runner and tests call g6pd_config().
function g6pd_config()
    (
        name = "G6PD",
        data_csv = "fitting/G6PD/rate_eq/G6PD_all_EnzymeData.csv",
        rate_col = "Rate_V",
        article_col = "Article",
        fig_col = "Fig",
        keq_col = "Apparent_Keq",
        deploy_keq = 13.655,               # apparent Keq at cellular pH 7.2 / 37C (eQuilibrator;
                                           # = corpus Mbata2026; matches deployed G6PD_Keq=13.7).
                                           # READOUT + DEPLOY only; the FIT uses per-figure d.keq.
        metabolites = Dict(                # symbol => (csv column, unit)
            :NADP  => ("[NADP] (uM)",  :uM),
            :G6P   => ("[G6P] (uM)",   :uM),
            :NADPH => ("[NADPH] (uM)", :uM),
            :PGLn  => ("[PGLn] (uM)",  :uM),
            :ATP   => ("[ATP] (uM)",   :uM),
        ),
    )
end

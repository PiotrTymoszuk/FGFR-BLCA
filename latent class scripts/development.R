# Development of the clustering solution by poLCA.
# A six-class solution seems to work the best

  insert_head()

# container ------

  lca_dev <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  lca_dev$variables <- lca_globals$variables

  lca_dev$data <- lca_globals$factor_data$genie %>%
    column_to_rownames('sample_id')

# parallelization --------

  insert_msg('Parallel backend')

  plan('multisession')

# tuning of the poLCA models ---------

  insert_msg('Tuninig of the poLCA models')

  lca_dev$tune_object <- tune_lca(formula = lca_globals$formula,
                                  data = lca_dev$data,
                                  class_range = 1:10,
                                  seed = 12345,
                                  cust_theme = globals$common_theme,
                                  maxiter = 3000,
                                  nrep = 5)

# Caching the results ---------

  insert_msg('Caching the results')

  lca_dev$data <- NULL

  lca_dev <- compact(lca_dev)

  save(lca_dev, file = './cache/lca_dev.RData')

# END ------

  plan('sequential')

  insert_tail()

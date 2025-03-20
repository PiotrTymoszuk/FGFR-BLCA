# SHAP variable importance for the Elastic Net and Random Forest models in
# the test portion of the TCGA data set.
# The models with the full set of predictors proposed by Kamoun et al.
#
# The background: a randomly chosen 30% of the TCGA gene expression data
# set.


  insert_head()

# container ------

  mibc_shap <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## matrix of explanatory factors,
  ## and responses in the test portions

  mibc_shap$data <- mibc_globals$data$tcga_test %>%
    column_to_rownames('sample_id')

  mibc_shap$x <- mibc_globals$x$tcga_test
  mibc_shap$y <- mibc_globals$y$tcga_test

  ## background matrix

  mibc_shap$bg_idx <- createDataPartition(mibc_shap$y, p = 0.3)[[1]]

  mibc_shap$bg_x <- mibc_shap$x[mibc_shap$bg_idx, ]

  ## best tunes

  mibc_shap$best_tunes <- mibc_train$tune_obj %>%
    map(~.x$best_tune)

# Models --------

  insert_msg('Models')

  set.seed(12345)

  mibc_shap$models$elnet <-
    glmnet(x = mibc_shap$x,
           y = mibc_shap$y,
           family = 'multinomial',
           alpha = 0.5,
           lambda = mibc_shap$best_tunes$elnet$lambda)

  mibc_shap$models$ranger <-
    ranger(x = mibc_shap$x,
           y = mibc_shap$y,
           num.trees = 500,
           mtry = mibc_shap$best_tunes$ranger$mtry,
           min.node.size = mibc_shap$best_tunes$ranger$min.node.size,
           splitrule = mibc_shap$best_tunes$ranger$splitrule,
           probability = TRUE)

# SHAP objects -------

  insert_msg('SHAP objects')

  mibc_shap$shap_obj$elnet <-
    kernelshap(mibc_shap$models$elnet,
               X = mibc_shap$x,
               bg_X = mibc_shap$bg_x,
               pred_fun = function(object, X, ...){

                 predict(object, X, type = 'response', ...)[, , 1]

               })

  mibc_shap$shap_obj$ranger <-
    kernelshap(mibc_shap$models$ranger,
               X = mibc_shap$x,
               bg_X = mibc_shap$bg_x)

# Caching ----------

  insert_msg('Caching')

  mibc_shap <- mibc_shap["shap_obj"]

  save(mibc_shap, file = './cache/mibc_shap.RData')

# END ------

  insert_tail()

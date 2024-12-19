# SHAP variable importance for the Elastic Net and Random Forest models in
# the TCGA data set.
#
# The background: a randomly chosen 30% of the TCGA gene expression data
# set.


  insert_head()

# container ------

  ml_shap <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## matrix of explanatory factors,
  ## and responses in the TCGA training cohort

  ml_shap$data <- ml_globals$data$tcga %>%
    column_to_rownames('sample_id')

  ml_shap$x <- ml_globals$x$tcga
  ml_shap$y <- ml_globals$y$tcga

  ## background matrix

  ml_shap$bg_idx <- createDataPartition(ml_shap$y, p = 0.3)[[1]]

  ml_shap$bg_x <- ml_shap$x[ml_shap$bg_idx, ]

  ## best tunes

  ml_shap$best_tunes <- ml_train$tune_obj %>%
    map(~.x$best_tune)

# Models --------

  insert_msg('Models')

  set.seed(12345)

  ml_shap$models$elnet <-
    glmnet(x = ml_shap$x,
           y = ml_shap$y,
           family = 'multinomial',
           alpha = 0.5,
           lambda = ml_shap$best_tunes$elnet$lambda)

  ml_shap$models$ranger <-
    ranger(x = ml_shap$x,
           y = ml_shap$y,
           num.trees = 500,
           mtry = ml_shap$best_tunes$ranger$mtry,
           min.node.size = ml_shap$best_tunes$ranger$min.node.size,
           splitrule = ml_shap$best_tunes$ranger$splitrule,
           probability = TRUE)

# SHAP objects -------

  insert_msg('SHAP objects')

  ml_shap$shap_obj$elnet <-
    kernelshap(ml_shap$models$elnet,
               X = ml_shap$x,
               bg_X = ml_shap$bg_x,
               pred_fun = function(object, X, ...){

                 predict(object, X, type = 'response', ...)[, , 1]

               })

  ml_shap$shap_obj$ranger <-
    kernelshap(ml_shap$models$ranger,
               X = ml_shap$x,
               bg_X = ml_shap$bg_x)

# Caching ----------

  insert_msg('Caching')

  ml_shap <- ml_shap["shap_obj"]

  save(ml_shap, file = './cache/ml_shap.RData')

# END ------

  insert_tail()

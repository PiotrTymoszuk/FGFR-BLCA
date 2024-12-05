# Predictions of consensus class assignment by the Elastic Net, Random
# Forest, and Ensemble models trained with FGF/FGFR/FGFBP expression data
# in the TCGA BLCA cohort.

  insert_head()

# container -------

  ml_pred <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## modeling matrices and the actual cluster assignment

  ml_pred$x <- ml_globals$x
  ml_pred$y <- ml_globals$y

  ## the optimized Elastic Net and Random Forest models

  ml_pred$models <- ml_train$tune_obj %>%
    map(~.x$model)

  ml_pred$models$ensemble <- ml_ens$tune_obj$model

  ## class levels

  ml_pred$class_levels <- ml_globals$class_levels

# Elastic Net predictions ---------

  insert_msg('Predictions')

  ## predicted class assignment

  ml_pred$predictions$elnet$class <- ml_pred$x %>%
    map(predict,
        object = ml_pred$models$elnet,
        type = 'class') %>%
    map(~.x[, 1]) %>%
    map(factor, ml_globals$class_levels)

  ## class assignment probabilities

  ml_pred$predictions$elnet$p_values <-  ml_pred$x %>%
    map(predict,
        object = ml_pred$models$elnet,
        type = 'response') %>%
    map(~.x[, , 1]) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id')

  ## the complete prediction tables

  ml_pred$predictions$elnet <-
    list(x = ml_pred$predictions$elnet$p_values,
         y = ml_pred$predictions$elnet$class,
         z = ml_pred$y) %>%
    pmap(function(x, y, z) x %>%
           mutate(pred = y,
                  obs = z))

# Random Forest predictions ---------

  insert_msg('Random Forest predictions')

  # class assignment p values

  ml_pred$predictions$ranger$p_values <- ml_pred$x %>%
    map(predict,
        object = ml_pred$models$ranger) %>%
    map(~.x[['predictions']]) %>%
    map(as.data.frame)

  ## class assignment by nodal voting

  for(i in names(ml_pred$predictions$ranger$p_values)) {

    ml_pred$predictions$ranger$class[[i]] <-
      1:nrow(ml_pred$predictions$ranger$p_values[[i]]) %>%
      map(~which(ml_pred$predictions$ranger$p_values[[i]][.x, ] ==
                   max(ml_pred$predictions$ranger$p_values[[i]][.x, ]))) %>%
      map_dbl(~.x[[1]])

    ml_pred$predictions$ranger$class[[i]] <-
      ml_pred$class_levels[ml_pred$predictions$ranger$class[[i]]]

  }

  ## the entire prediction frames

  ml_pred$predictions$ranger <-
    list(x = ml_pred$predictions$ranger$p_values,
         y = ml_pred$predictions$ranger$class,
         z = ml_pred$y,
         w = ml_pred$x) %>%
    pmap(function(x, y, z, w) x %>%
           mutate(pred = factor(y, ml_pred$class_levels),
                  obs = z,
                  sample_id = rownames(w)))

# modeling matrices for the Ensemble model -------

  insert_msg('modeling matrices for the Ensemble model')

  for(i in names(ml_pred$predictions)) {

    ml_pred$ens_x[[i]] <- ml_pred$predictions[[i]] %>%
      map(select, all_of(ml_pred$class_levels)) %>%
      map(~set_names(.x, paste(names(.x), i, sep = '_')))

  }

  ml_pred$ens_x <- ml_pred$ens_x %>%
    transpose %>%
    map(reduce, cbind) %>%
    map(as.matrix)

# Predictions for the Ensemble model --------

  insert_msg('Predictions for the ensemble model')

  ## predicted class assignment

  ml_pred$predictions$ensemble$class <- ml_pred$ens_x %>%
    map(predict,
        object = ml_pred$models$ensemble,
        type = 'class') %>%
    map(~.x[, 1]) %>%
    map(factor, ml_globals$class_levels)

  ## class assignment probabilities

  ml_pred$predictions$ensemble$p_values <-  ml_pred$ens_x %>%
    map(predict,
        object = ml_pred$models$ensemble,
        type = 'response') %>%
    map(~.x[, , 1]) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id')

  ## the complete prediction tables

  ml_pred$predictions$ensemble <-
    list(x = ml_pred$predictions$ensemble$p_values,
         y = ml_pred$predictions$ensemble$class,
         z = ml_pred$y) %>%
    pmap(function(x, y, z) x %>%
           mutate(pred = y,
                  obs = z))

# Brier scores for the model and random predictions --------

  insert_msg('Brier scores')

  ml_pred$predictions <- ml_pred$predictions %>%
    map(map, brier_squares) %>%
    map(map, brier_reference, n = 1000)

# END ------

  ml_pred <- ml_pred["predictions"]

  insert_tail()

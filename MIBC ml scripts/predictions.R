# Predictions of consensus class assignment by the Elastic Net, Random
# Forest, and Ensemble models trained with expression data proposed by
# Kamoun et al. in the TCGA BLCA cohort.

  insert_head()

# container -------

  mibc_pred <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## modeling matrices and the actual cluster assignment

  mibc_pred$x <- mibc_globals$x
  mibc_pred$y <- mibc_globals$y

  ## the optimized Elastic Net and Random Forest models

  mibc_pred$models <- mibc_train$tune_obj %>%
    map(~.x$model)

  mibc_pred$models$ensemble <- mibc_ens$tune_obj$model

  ## class levels

  mibc_pred$class_levels <- mibc_globals$class_levels

# Elastic Net predictions ---------

  insert_msg('Predictions')

  ## predicted class assignment

  mibc_pred$predictions$elnet$class <- mibc_pred$x %>%
    map(predict,
        object = mibc_pred$models$elnet,
        type = 'class') %>%
    map(~.x[, 1]) %>%
    map(factor, mibc_globals$class_levels)

  ## class assignment probabilities

  mibc_pred$predictions$elnet$p_values <-  mibc_pred$x %>%
    map(predict,
        object = mibc_pred$models$elnet,
        type = 'response') %>%
    map(~.x[, , 1]) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id')

  ## the complete prediction tables

  mibc_pred$predictions$elnet <-
    list(x = mibc_pred$predictions$elnet$p_values,
         y = mibc_pred$predictions$elnet$class,
         z = mibc_pred$y) %>%
    pmap(function(x, y, z) x %>%
           mutate(pred = y,
                  obs = z))

# Random Forest predictions ---------

  insert_msg('Random Forest predictions')

  # class assignment p values

  mibc_pred$predictions$ranger$p_values <- mibc_pred$x %>%
    map(predict,
        object = mibc_pred$models$ranger) %>%
    map(~.x[['predictions']]) %>%
    map(as.data.frame)

  ## class assignment by nodal voting

  for(i in names(mibc_pred$predictions$ranger$p_values)) {

    mibc_pred$predictions$ranger$class[[i]] <-
      1:nrow(mibc_pred$predictions$ranger$p_values[[i]]) %>%
      map(~which(mibc_pred$predictions$ranger$p_values[[i]][.x, ] ==
                   max(mibc_pred$predictions$ranger$p_values[[i]][.x, ]))) %>%
      map_dbl(~.x[[1]])

    mibc_pred$predictions$ranger$class[[i]] <-
      mibc_pred$class_levels[mibc_pred$predictions$ranger$class[[i]]]

  }

  ## the entire prediction frames

  mibc_pred$predictions$ranger <-
    list(x = mibc_pred$predictions$ranger$p_values,
         y = mibc_pred$predictions$ranger$class,
         z = mibc_pred$y,
         w = mibc_pred$x) %>%
    pmap(function(x, y, z, w) x %>%
           mutate(pred = factor(y, mibc_pred$class_levels),
                  obs = z,
                  sample_id = rownames(w)))

# modeling matrices for the Ensemble model -------

  insert_msg('modeling matrices for the Ensemble model')

  for(i in names(mibc_pred$predictions)) {

    mibc_pred$ens_x[[i]] <- mibc_pred$predictions[[i]] %>%
      map(select, all_of(mibc_pred$class_levels)) %>%
      map(~set_names(.x, paste(names(.x), i, sep = '_')))

  }

  mibc_pred$ens_x <- mibc_pred$ens_x %>%
    transpose %>%
    map(reduce, cbind) %>%
    map(as.matrix)

# Predictions for the Ensemble model --------

  insert_msg('Predictions for the ensemble model')

  ## predicted class assignment

  mibc_pred$predictions$ensemble$class <- mibc_pred$ens_x %>%
    map(predict,
        object = mibc_pred$models$ensemble,
        type = 'class') %>%
    map(~.x[, 1]) %>%
    map(factor, mibc_globals$class_levels)

  ## class assignment probabilities

  mibc_pred$predictions$ensemble$p_values <-  mibc_pred$ens_x %>%
    map(predict,
        object = mibc_pred$models$ensemble,
        type = 'response') %>%
    map(~.x[, , 1]) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id')

  ## the complete prediction tables

  mibc_pred$predictions$ensemble <-
    list(x = mibc_pred$predictions$ensemble$p_values,
         y = mibc_pred$predictions$ensemble$class,
         z = mibc_pred$y) %>%
    pmap(function(x, y, z) x %>%
           mutate(pred = y,
                  obs = z))

# Brier scores for the model and random predictions --------

  insert_msg('Brier scores')

  mibc_pred$predictions <- mibc_pred$predictions %>%
    map(map, brier_squares) %>%
    map(map, brier_reference, n = 1000)

# END ------

  mibc_pred <- mibc_pred["predictions"]

  insert_tail()

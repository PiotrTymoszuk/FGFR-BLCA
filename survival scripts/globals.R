# Data for survival analyses. The expression levels are ComBat-processed.

  insert_head()

# container -------

  surv_globals <- list()

# analysis data ------

  insert_msg('Analysis globals')

  ## variables of interest and their lexicon

  surv_globals$lexicon <-
    globals[c("ligands", "receptors", "binding_proteins")] %>%
    map(~tibble(variable = .x)) %>%
    compress(names_to = 'protein_type')

  surv_globals$data <- combat$expression[c("tcga", "imvigor", "bcan")] %>%
    map(select, sample_id, any_of(surv_globals$lexicon$variable))

  surv_globals$variables <- surv_globals$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  surv_globals$lexicon <- surv_globals$lexicon %>%
    filter(variable %in% surv_globals$variable)

  ## first and second order terms of log2 gene expression levels

  surv_globals$data <- surv_globals$data %>%
    map(select, sample_id, all_of(surv_globals$variables))

  for(i in names(surv_globals$data)) {

    surv_globals$data[[i]][paste0(surv_globals$variables, '_sec')] <-
      surv_globals$data[[i]][surv_globals$variables] %>%
      map_dfc(~.x^2)

  }

  surv_globals$variables <- surv_globals$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  ## survival information

  surv_globals$clinic <- list(tcga = tcga,
                              imvigor = imvigor,
                              bcan = bcan) %>%
    map(~.x$clinic[c('sample_id', 'os_days', 'death')])

  surv_globals$data <-
    map2(surv_globals$data,
         surv_globals$clinic,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate, os_days = ifelse(os_days <= 0,
                                 1,
                                 os_days + 1))

# Survival objects and model matrices ------

  insert_msg('Survival objects and x matrices')

  surv_globals$y <- surv_globals$data %>%
    map(~Surv(.x$os_days, .x$death))

  surv_globals$x <- surv_globals$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(surv_globals$variables)) %>%
    map(as.matrix)

# CV fold used for tuning ---------

  insert_msg('CV folds used for tuning')

  surv_globals$n_rep <- 200

  set.seed(1234)

  surv_globals$cv_folds <- 1:surv_globals$n_rep %>%
    map(function(x) createFolds(surv_globals$data$tcga$death,
                                list = FALSE,
                                returnTrain = TRUE)) %>%
    set_names(paste0('rep_', 1:surv_globals$n_rep))

# END -----

  surv_globals$clinic <- NULL

  surv_globals <- compact(surv_globals)

  insert_tail()

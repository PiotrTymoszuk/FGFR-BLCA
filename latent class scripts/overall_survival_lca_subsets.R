# Checking whether a combination of genetic information (Genetic Subsets) and
# transcriptomic classification (consensus clusters of MIBC) may improve
# prediction of overall survival.
#
# We're comparing three models:
#
# 1) a Cox model with the consensus MIBC subset as an independent variable
#
# 2) a Cox model with the genetic subset as an independent variable
#
# 3) a Cox model with the consensus and genetic subsets as independent variables
#
# A demonstrated both in the training data set and in 5-repeats 5-fold
# cross-validation, there are no differences in performance between the models
# 1 and 3, and, hence, no benefot from inclusion of the genetic information.


  insert_head()

# container -------

  lca_mos <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  lca_mos$data$transcriptome <- subtypes$assignment$tcga %>%
    select(sample_id, consensusClass)

  lca_mos$data$genome <- lca_pred$assignment$tcga

  lca_mos$data$clinic <- tcga$clinic %>%
    select(sample_id, os_days, death)

  lca_mos$data <- lca_mos$data %>%
    reduce(inner_join, by = 'sample_id') %>%
    filter(complete.cases(.))

# Cox models and their cross-validation --------

  insert_msg('Cox models and CV')

  lca_mos$models <-
    list(transcriptome = Surv(os_days, death) ~ consensusClass,
         genome = Surv(os_days, death) ~ clust_id,
         combi = Surv(os_days, death) ~ consensusClass + clust_id) %>%
    map(~call2('coxph',
               formula = .x,
               data = lca_mos$data,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map(as_coxex, lca_mos$data)

  lca_mos$cv_models <- lca_mos$models %>%
    map(cvCox, n_folds = 5, n_repeats = 5)

# Assumptions, fit stats, and inference, training data ------

  insert_msg('Assumptions, fit, and inference')

  ## assumptions: more or less hold by all models

  lca_mos$assumptions <- lca_mos$models %>%
    map(summary, 'assumptions')

  ## fit stats

  lca_mos$stats <- list(train = lca_mos$models,
                        cv = lca_mos$cv_models) %>%
    map(map, summary, 'fit')

  lca_mos$stats$cv <-   lca_mos$stats$cv %>%
    map(filter, summary_stat == 'mean') %>%
    map(select, -summary_stat)

  lca_mos$stats <- lca_mos$stats %>%
    transpose %>%
    map(~map2_dfr(.x, names(.x), ~mutate(.x, data_set = .y)))

  ## inference

  lca_mos$inference <- lca_mos$models %>%
    map(summary, 'inference')

# END -------

  lca_mos$data <- NULL
  lca_mos$models <- NULL
  lca_mos$cv_models <- NULL

  lca_mos <- compact(lca_mos)

  insert_tail()

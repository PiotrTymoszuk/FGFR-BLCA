# Evaluation of the Elastic Net Cox model of overall survival.
#
# For predictions of linear predictor scores in the TCGA, IMvigor, and BCAN
# cohorts, univariable Cox models of OS are constructed.
# Their performance is assessed by C-index, IBC, and R^2.

  insert_head()

# sontainer -------

  surv_eval <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## survival data sets and modeling matrices

  surv_eval$data <- surv_globals$data %>%
    map(select, sample_id, os_days, death)

  surv_eval$x <- surv_globals$x

  ## the GLMNET model

  surv_eval$glmnet_model <- surv_train$tune_obj$model

# Linear predictor scores -------

  insert_msg('LP scores for the training and validation cohorts')

  surv_eval$lp_scores <- surv_eval$x %>%
    map(predict, object = surv_eval$glmnet_model) %>%
    map(as.data.frame) %>%
    map(set_names, 'lp_score') %>%
    map(rownames_to_column, 'sample_id') %>%
    map2(., surv_eval$data,
         inner_join, by = 'sample_id') %>%
    map(as_tibble)

# Univariable Cox models of OS with the LP scores as explanatory variables -------

  insert_msg('Univariable Cox models')

  surv_eval$models <- surv_eval$lp_scores %>%
    map(~call2('coxph',
               formula = Surv(os_days, death) ~ lp_score,
               data = .x,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map2(., surv_eval$lp_scores, as_coxex)

# Model assumptions, fit stats, and inference -------

  insert_msg('Model assumptions, stats, inference')

  ## the PH assumption is met!

  surv_eval$assumptions <- surv_eval$models %>%
    map(summary, 'assumptions')

  ## fit stats

  surv_eval$stats <- surv_eval$models %>%
    map(summary, 'fit')

  ## inference

  surv_eval$inference <- surv_eval$models %>%
    map(summary, 'inference')

  surv_eval[c("stats", "inference")] <-
    surv_eval[c("stats", "inference")] %>%
    map(compress, names_to = 'cohort') %>%
    map(mutate,
        data_set = ifelse(cohort == 'tcga',
                          'training', 'test'),
        data_set = factor(data_set,
                          c('training', 'test')))

# Plot of global performance stats -------

  insert_msg('Plot of global performance stats')

  surv_eval$stat_plot <- surv_eval$stats %>%
    ggplot(aes(x = c_index,
               y = 1 - ibs_model,
               fill = data_set,
               color = data_set,
               size = raw_rsq)) +
    geom_hline(yintercept = 0.75,
               linetype = 'dashed') +
    geom_vline(xintercept = 0.5,
               linetype = 'dashed') +
    geom_point(shape = 21,
               color = 'black') +
    geom_text_repel(aes(label = globals$cohort_labs[cohort] %>%
                          paste(signif(raw_rsq, 2),
                                sep = '\nR\u00B2 = ')),
                    size = 2.5,
                    hjust = 0,
                    box.padding = 0.5,
                    show.legend = FALSE) +
    scale_fill_manual(values = c(training = 'indianred3',
                                 test = 'steelblue'),
                      name = 'data set') +
    scale_color_manual(values = c(training = 'indianred3',
                                 test = 'steelblue'),
                      name = 'data set') +
    scale_size_area(max_size = 5,
                    limits = c(0, max(surv_eval$stats$raw_rsq)),
                    name = expression(R^2)) +
    globals$common_theme +
    labs(title = 'Prediction of overall survival',
         subtitle = 'model performance',
         x = 'accuracy, C-index',
         y = 'calibration, 1 - IBS')

# Plots of Brier scores for unique time points --------

  insert_msg('Plots of Brier scores for unique time points')

  surv_eval$brier_plots <- surv_eval$models %>%
    map(surv_brier) %>%
    map(plot,
        cust_theme = globals$common_theme)

  surv_eval$brier_plots <-
    list(x = surv_eval$brier_plots,
         y = globals$cohort_labs[names(surv_eval$brier_plots)]) %>%
    pmap(function(x, y) x +
           scale_color_manual(values = c(reference = 'gray40',
                                         training = 'orangered3'),
                              labels = c(reference = 'NULL model',
                                         training = 'Cox model')))

# Model calibration: survival in tertiles of LP scores -------

  insert_msg('Survival in tertiles of LP scores')

  plan('multisession')

  surv_eval$calibrator_obj <- surv_eval$models %>%
    future_map(calibrate.coxex,
               n = 3,
               labels = c('low', 'int', 'high'),
               .options = furrr_options(seed = TRUE))

  plan('sequential')

  ## testing for differences in survival between the tertiles

  surv_eval$tertile_test <- surv_eval$calibrator_obj %>%
    map(~.x$surv_fit) %>%
    map(surv_pvalue, method = 'S1') %>%
    compress(names_to = 'cohort') %>%
    mutate(p_value = pval) %>%
    re_adjust(method = 'none') %>%
    mutate(raw_significance = significance) %>%
    re_adjust %>%
    mutate(plot_lab = paste0('raw: ', raw_significance,
                             '\nFDR: ', significance))

  ## legends for the KM plots with total numbers of observations and events

  surv_eval$tertile_legends <- surv_eval$calibrator_obj %>%
    map(~.x$strata_calibration) %>%
    map(function(dat) list(x = dat$label,
                           y = dat$n,
                           z = dat$km_events) %>%
          pmap_chr(function(x, y, z) paste0(x,
                                            '\ntotal: n = ', y,
                                            '\nevents: n = ', z)))

  ## Kaplan-Meier plot, displaying the event numbers
  ## in the plot legend

  surv_eval$tertile_km_plots <-
    list(x = surv_eval$calibrator_obj,
         title = globals$cohort_labs[names(surv_eval$calibrator_obj)]) %>%
    pmap(plot,
         cust_theme = globals$common_theme,
         show_cox = FALSE,
         xlab = 'overall survival, days') %>%
    map(tag2subtitle)

  surv_eval$tertile_km_plots <-
    list(x = surv_eval$tertile_km_plots,
         y = surv_eval$tertile_test$plot_lab,
         z = surv_eval$tertile_legends) %>%
    pmap(function(x, y, z) x +
           annotate('text',
                    label = y,
                    x = 0,
                    y = 0,
                    hjust = 0,
                    vjust = 0,
                    size = 2.5) +
           scale_color_manual(values = c('steelblue', 'gray30', 'firebrick'),
                              labels = z,
                              name = 'LP\nscore'))

# END -------

  surv_eval$data <- NULL
  surv_eval$x <- NULL
  surv_eval$glmnet_model <- NULL
  surv_eval$models <- NULL
  surv_eval$calibrator_obj <- NULL

  surv_eval <- compact(surv_eval)

  insert_tail()

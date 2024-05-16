# Univariable and multivariable modeling of overall survival as a function
# of genetic subset, age, sex, and pathological stage.
# Done for the GENIE and MSK cohorts


# container --------

  lca_cox <- list()

# analysis data --------

  insert_msg('Analysis data')

  ## analysis data: GENIE and TCGA

  lca_cox$data <- list(genie = genie, msk = msk) %>%
    map(~.x$clinic) %>%
    map(select,
        sample_id,
        os_days, death,
        age, sex,
        any_of('pt_stage'))

  lca_cox$data <- map2(lca_pred$assignment[names(lca_cox$data)],
                      lca_cox$data,
                      inner_join, by = 'sample_id')

  lca_cox$data <- lca_cox$data %>%
    map(~filter(.x, complete.cases(.x)))

  ## simplifying the stage: T3 and T4 together,
  ## age in decades

  lca_cox$data$msk <- lca_cox$data$msk %>%
    mutate(pt_stage = car::recode(pt_stage,
                                  "'T1' = 'T1/2';
                                  'T2' = 'T1/2';
                                  'T3' = 'T3/4';
                                  'T4' = 'T3/4'"),
           pt_stage = factor(pt_stage, c('T1/2', 'T3/4')))

  lca_cox$data <- lca_cox$data %>%
    map(mutate,
        age = age/10)

# Numbers of cases and events -------

  insert_msg('Number of cases and events')

  lca_cox$n_numbers <- lca_cox$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x$death))) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort_lab = globals$cohort_labs[cohort],
           n_caption = paste0('total: n = ', n_total,
                              ', events: n = ', n_events),
           facet_lab = paste(cohort_lab, n_caption, sep = '\n'))

# Cox models --------

  insert_msg('Cox models')

  ## univariable models and models abjusted for clinical confounders

  ## the GENIE cohort; age and sex as confounders

  lca_cox$models$genie[c('uni', 'full')] <-
    list(Surv(os_days, death) ~ clust_id,
         Surv(os_days, death) ~ age + sex + clust_id) %>%
    map(~call2(.fn = 'coxph',
               formula = .x,
               data = lca_cox$data$genie,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map(as_coxex, lca_cox$data$genie)

  ## the MSK cohort; age, sex, and pT stage as confounders

  lca_cox$models$msk[c('uni', 'full')] <-
    list(Surv(os_days, death) ~ clust_id,
         Surv(os_days, death) ~ age + sex + pt_stage + clust_id) %>%
    map(~call2(.fn = 'coxph',
               formula = .x,
               data = lca_cox$data$msk,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map(as_coxex, lca_cox$data$msk)

# Assumptions and fit stats ---------

  insert_msg('Proportional hazard assumption and fit stats')

  ## the proportional hazard assumption is significantly violated for
  ## the full model of survival in the MSK data set

  lca_cox$assumptions <- lca_cox$models %>%
    map(map, summary, 'assumptions')

  ## fit stats

  lca_cox$stats <- lca_cox$models %>%
    map(map, summary, 'fit') %>%
    map(compress, names_to = 'model_type') %>%
    compress(names_to = 'cohort')

# Scatter plot of the general model performance stats -------

  insert_msg('Scatter/bubble plots of the fit stats')

  lca_cox$stat_plot <- lca_cox$stats %>%
    ggplot(aes(x = c_index,
               y = 1 - ibs_model,
               size = raw_rsq,
               fill = model_type,
               color = model_type)) +
    geom_vline(xintercept = 0.5,
               linetype = 'dashed') +
    geom_hline(yintercept = 0.75,
               linetype = 'dashed') +
    geom_errorbarh(aes(xmin = lower_ci,
                       xmax = upper_ci),
                   height = 0,
                   linewidth = 0.5,
                   show.legend = FALSE) +
    geom_point(shape = 21,
               color = 'black') +
    facet_grid(. ~ cohort,
               labeller = as_labeller(set_names(lca_cox$n_numbers$facet_lab,
                                                lca_cox$n_numbers$cohort))) +
    scale_size_area(max_size = 5,
                    limits = c(0, 0.25),
                    name = expression('R'^2))  +
    scale_fill_manual(values = c(uni = 'steelblue',
                                 full = 'indianred3'),
                      labels = c(uni = 'genetic subsets',
                                 full = 'genetic subsets\n+ confounders'),
                      name = 'Model type') +
    scale_color_manual(values = c(uni = 'steelblue',
                                 full = 'indianred3'),
                      labels = c(uni = 'genetic subsets',
                                 full = 'genetic subsets\n+ confounders'),
                      name = 'Model type') +
    globals$common_theme +
    labs(title = 'Cox models of overall survival',
         x = 'C-index \u00B1 95% CI',
         y = '1 - IBS')

# LRT for the genetic subset term -------

  insert_msg('LRT')

  lca_cox$lrt <- lca_cox$models %>%
    map(~anova(as_coxph(.x[[1]]),
               as_coxph(.x[[2]]))) %>%
    map(as.data.frame) %>%
    map(~.x[2, ]) %>%
    map(set_names,
        c('log_lik', 'chisq', 'df', 'p_value')) %>%
    compress(names_to = 'cohort') %>%
    re_adjust(method = 'none') %>%
    mutate(plot_lab = paste0('\u03C7\u00B2(', df, ') = ',
                             signif(chisq, 2),
                             ', ', significance)) %>%
    as_tibble

# Inference -------

  insert_msg('Inference')

  lca_cox$inference <- lca_cox$models %>%
    map(map, summary, 'inference') %>%
    map(compress, names_to = 'model_type') %>%
    ## hazard ratios
    map(mutate,
        model_type = factor(model_type, c('uni', 'full')),
        hr = exp(estimate),
        hr_lower = exp(lower_ci),
        hr_upper = exp(upper_ci),
        hr_label = paste0(signif(hr, 2),
                          ' [', signif(hr_lower, 2),
                          ' to ', signif(hr_upper, 2), ']')) %>%
    ## plot labels
    map(mutate,
        axis_label = car::recode(variable,
                                 "'clust_id' = '';
                                 'pt_stage' = 'stage';
                                 'age' = 'age, per decade'"),
        axis_label = ifelse(level != '',
                            ifelse(axis_label == '',
                                   level,
                                   paste(axis_label, level, sep = ': ')),
                            axis_label),
        regulation = ifelse(p_value > 0.05, 'ns',
                            ifelse(estimate < 0, 'favorable',
                                   ifelse(estimate > 0, 'unfavorable', 'ns'))),
        regulation = factor(regulation,
                            c('favorable', 'unfavorable', 'ns')))

# Forest plots of HR and their confidence intervals ------

  insert_msg('Forest plots')

  lca_cox$forest_plots <-
    list(x = lca_cox$inference,
         y = globals$cohort_labs[names(lca_cox$inference)],
         z = lca_cox$n_numbers$n_caption,
         w = lca_cox$lrt$plot_lab) %>%
    pmap(function(x, y, z, v, w) x %>%
           ggplot(aes(x = hr,
                      y = reorder(axis_label, hr),
                      color = regulation)) +
           geom_vline(xintercept = 1,
                      linetype = 'dashed') +
           geom_errorbarh(aes(xmin = hr_lower,
                              xmax = hr_upper),
                          height = 0) +
           geom_point(shape = 16,
                      size = 2) +
           geom_text(aes(label = hr_label),
                     size = 2.5,
                     hjust = 0.15,
                     vjust = -1.2) +
           scale_color_manual(values = c(favorable = 'steelblue',
                                         unfavorable = 'firebrick',
                                         ns = 'gray60'),
                              name = 'Survival marker') +
           facet_grid(model_type ~ .,
                      scales = 'free',
                      space = 'free',
                      labeller = as_labeller(c('uni' = 'genetic subsets',
                                               'full' = 'genetic subsets\n+confounders'))) +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = paste(z, w, sep = '\n'),
                x = 'HR \u00B1 95% CI'))

# END -----

  insert_tail()

# Evaluation of performance of the ML models of consensus classes that
# use FGFR, FGF, and FGFBP genes as explanatory factors.

  insert_head()

# container --------

  ml_eval <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## predictions

  ml_eval$predictions <- ml_pred$predictions

  ## algorithm labels

  ml_eval$algo_labs <-
    c(elnet = 'Elastic Net',
      ranger = 'Random Forest',
      ensemble = 'Ensemble')

# Performance stats --------

  insert_msg('Performance stats')

  ## global performance stats: ROC stats, accuracy, kappa, and Brier scores

  ml_eval$stats <- ml_eval$predictions %>%
    map(map, ml_summary) %>%
    map(compress, names_to = 'cohort') %>%
    compress(names_to = 'algorithm') %>%
    mutate(data_set = ifelse(cohort == 'tcga',
                             'training', 'test'),
           data_set = factor(data_set, c('training', 'test')))

  ## dot plot

  ml_eval$stat_plot <- ml_eval$stats %>%
    ggplot(aes(x = Kappa,
               y = 2 - brier_score,
               size = AUC,
               fill = data_set,
               color = data_set)) +
    facet_grid(. ~ algorithm,
               labeller = as_labeller(ml_eval$algo_labs)) +
    geom_rect(aes(xmin = -Inf,
                  ymin = 2 - min(ml_eval$stats$brier_reference),
                  xmax = Inf,
                  ymax = 2 - max(ml_eval$stats$brier_reference)),
              fill = 'gray80',
              color = NA,
              alpha = 0.5,
              show.legend = FALSE) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_point(shape = 21,
               color = 'black') +
    geom_text_repel(aes(label = paste(globals$cohort_labs[cohort],
                                      signif(AUC, 2),
                                      sep = '\nAUC = ')),
                    size = 2.3,
                    hjust = 0,
                    box.padding = 0.6,
                    show.legend = FALSE) +
    scale_fill_manual(values = c(training = 'indianred3',
                                 test = 'steelblue'),
                      name = 'data set') +
    scale_color_manual(values = c(training = 'indianred3',
                                 test = 'steelblue'),
                      name = 'data set') +
    scale_size_area(max_size = 5,
                    limits = c(0.5, 1),
                    name = 'AUC') +
    guides(size = guide_legend(override.aes = list(fill = 'steelblue'))) +
    globals$common_theme +
    labs(title = 'Prediction of consensus classes',
         subtitle = 'model performance',
         x = "accuracy, Cohen's \u03BA",
         y = 'calibration, 2 - Brier score')

# Confusion matrices --------

  insert_msg('Confusion matrices')

  ## confusion matrices

  ml_eval$confusion_mtx <- ml_eval$predictions %>%
    map(map,
        ~table(obs = .x$obs,
               pred = .x$pred))

  ## heat map representations

  for(i in names(ml_eval$confusion_mtx)) {

    ml_eval$confusion_plots[[i]] <-
      list(mtx = ml_eval$confusion_mtx[[i]],
           plot_title = globals$cohort_labs[names(ml_eval$confusion_mtx[[i]])] %>%
             paste(ml_eval$algo_labs[[i]], sep = ', ')) %>%
      pmap(plot_confusion) %>%
      map(~.x +
            scale_y_discrete(labels = function(x) stri_replace(x,
                                                               fixed = '-',
                                                               replacement = '\n')) +
            scale_fill_gradient2(low = 'steelblue',
                                 high = 'firebrick',
                                 mid = 'white',
                                 limits = c(0, 45),
                                 midpoint = 22.5,
                                 name = '% of all samples'))

  }

# Plots of Brier squares --------

  insert_msg('Plots of brier squares')

  for(i in names(ml_eval$predictions)) {

    ## box plots

    ml_eval$brier_square_plots[[i]] <-
      map2(ml_eval$predictions[[i]] %>%
             map(mutate,
                 square_inv = 2 - square_dist),
           globals$cohort_labs[names(ml_eval$predictions[[i]])] %>%
             paste(ml_eval$algo_labs[[i]], sep = ', '),
           function(x, y) x %>%
             plot_variable(variable = 'square_inv',
                           split_factor = 'obs',
                           type = 'box',
                           cust_theme = globals$common_theme +
                             theme(legend.position = 'none'),
                           y_lab = '2 - Brier square distance',
                           x_lab = 'observed consensus class',
                           plot_subtitle = 'model calibration',
                           plot_title = y,
                           x_n_labs = TRUE))

    ## styling: colors and squares expected for a nonsense model

    ml_eval$brier_square_plots[[i]] <-
      list(x = ml_eval$brier_square_plots[[i]],
           y = ml_eval$stats %>%
             filter(algorithm == i) %>%
             .$brier_reference) %>%
      pmap(function(x, y) x +
             geom_hline(yintercept = y,
                        linetype = 'dashed') +
             scale_fill_manual(values = globals$sub_colors) +
             scale_y_continuous(limits = c(0, 2)))

  }

# Calibration curves: predicted vs observed probability for single classes -------

  insert_msg('Calibration curves')

  for(i in names(ml_eval$predictions)) {

    ml_eval$calibration_plots[[i]] <-
      list(prediction = ml_eval$predictions[[i]],
           plot_title = globals$cohort_labs[names(ml_eval$predictions[[i]])] %>%
             paste(ml_eval$algo_labs[[i]], sep = ', ')) %>%
      pmap(plot_calibration,
           discrete = FALSE,
           smooth = 'none',
           logistic.cal = TRUE,
           allowPerfectPredictions = TRUE)

  }

# END -------

  rm(i)

  ml_eval <-
    ml_eval[c("stats", "stat_plot",
              "confusion_mtx", "confusion_plots",
              "brier_square_plots", "calibration_plots")]

  insert_tail()

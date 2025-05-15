# Comparison of performance of the models of consensus classes of MIBC.
#
# The models that use the full set of predictors proposed by Kamoun et al.
# and the FGFR/FGF/FGFBP-related genes are compared.

  insert_head()

# container --------

  mibc_comp <- list()

# Analysis globals --------

  insert_msg('Analysis globals')

  ## labels and colors for the model types

  mibc_comp$model_labels <-
    c(fgfr_genes = 'FGFR/FGF/FGFBP',
      mibc_genes = 'Kamoun et al.')

  mibc_comp$model_colors <-
    c(fgfr_genes = 'aquamarine4',
      mibc_genes = 'bisque4')

  ## predictions

  mibc_comp$predictions <-
    list(fgfr_genes = ml_pred,
         mibc_genes = mibc_pred) %>%
    map(~.x$predictions) %>%
    map(map, map, as_tibble) %>%
    transpose %>%
    map(transpose) %>%
    map(map, compress, names_to = 'model')

  ## performance stats, ready to use plot labels
  ## with the cohort name

  mibc_comp$stats <-
    list(fgfr_genes = ml_eval,
         mibc_genes = mibc_eval) %>%
    map(~.x$stats) %>%
    compress(names_to = 'model') %>%
    mutate(plot_lab = ifelse(model == 'fgfr_genes',
                             globals$cohort_labs[cohort], NA),
           cohort = factor(cohort, names(mibc_comp$predictions[[1]])))

  ## calibration data (optional)

  #mibc_comp$calibration <-
   # list(fgfr_genes = ml_eval,
    #     mibc_genes = mibc_eval) %>%
    #map(~.x$calibration_plots) %>%
    #map(map, map, ~.x$data) %>%
    #transpose %>%
    #map(transpose) %>%
    #map(map, compress, names_to = 'model')

# Plots of the performance stats --------

  insert_msg('Plots of the performance stats')

  mibc_comp$stat_plot <- mibc_comp$stats %>%
    filter(cohort != 'tcga') %>%
    ggplot(aes(x = Kappa,
               y = 2 - brier_score,
               fill = model,
               color = model,
               shape = model)) +
    facet_grid(. ~ algorithm,
               labeller = as_labeller(globals$algo_labs)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    annotate('rect',
             xmin = -Inf,
             xmax = Inf,
             ymin = 2 - max(mibc_comp$stats$brier_reference),
             ymax = 2 - min(mibc_comp$stats$brier_reference),
             fill = 'gray60',
             alpha = 0.25) +
    geom_line(aes(group = cohort),
              color = 'gray40') +
    geom_point(size = 2,
               color = 'black') +
    geom_text_repel(aes(label = plot_lab),
                    size = 2.5,
                    show.legend = FALSE,
                    box.padding = 0.6,
                    force = 2) +
    scale_color_manual(values = mibc_comp$model_colors,
                       labels = mibc_comp$model_labels,
                       name = 'predictor genes') +
    scale_fill_manual(values = mibc_comp$model_colors,
                       labels = mibc_comp$model_labels,
                       name = 'predictor genes') +
    scale_shape_manual(values = c(fgfr_genes = 21,
                                  mibc_genes = 23),
                       labels = mibc_comp$model_labels,
                       name = 'predictor genes') +
    globals$common_theme +
    labs(title = 'Prediction of consensus classes',
         subtitle = 'model performance',
         x = "accuracy, Cohen's \u03BA",
         y = 'calibration, 2 - Brier score')

# Plots of Brier squares -------

  insert_msg('Plots of Brier squares')

  ## violin plots with point overlay, to make them compatible
  ## with the journal policy

  for(i in names(mibc_comp$predictions)) {

    mibc_comp$brier_square_plots[[i]] <-
      list(x = mibc_comp$predictions[[i]] %>%
             map(mutate,
                 obs = fct_recode(obs, `Stroma\nrich` = 'Stroma-rich')),
           y = mibc_comp$stats %>%
             #filter(algorithm == i) %>%
             blast(cohort),
           z = globals$cohort_labs[names(mibc_comp$predictions[[i]])] %>%
             paste(globals$algo_labs[i], sep = ', ')) %>%
      pmap(function(x, y, z) x %>%
             ggplot(aes(x = obs,
                        y = 2 - square_dist,
                        fill = model)) +
             geom_violin(alpha = 0.5,
                         position = position_dodge(0.9),
                         scale = "width") +
             geom_point(shape = 21,
                        position = position_jitterdodge(dodge.width = 0.9,
                                                        jitter.width = 0.1,
                                                        jitter.height = 0)) +
             geom_hline(yintercept = 2 - min(y$brier_reference),
                        linetype = 'dashed') +
             scale_fill_manual(values = mibc_comp$model_colors,
                               labels = mibc_comp$model_labels,
                               name = 'predictor genes') +
             scale_y_continuous(limits = c(0, 2)) +
             globals$common_theme +
             labs(title = z,
                  x = 'observed consensus class',
                  y = '2 - square distance to outcome'))

  }

# END ------

  insert_tail()

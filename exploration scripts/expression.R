# Expression of FGFR, FGF, and FGFBP genes and proteins.
# The later are investigated by IHC in the Human Protein Atlas.

  insert_head()

# container ------

  expl_expr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variables of interest and the corresponding expression levels

  expl_expr$lexicon <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    map(~tibble(variable = .x)) %>%
    compress(names_to = 'protein_type') %>%
    mutate(protein_type = factor(protein_type,
                                 c('receptors', 'ligands', 'binding_proteins')))

  expl_expr$variables <- expl_expr$lexicon$variable

  expl_expr$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    c(list(hpa = hpa$ihc_score)) %>%
    compact

  expl_expr$data <- expl_expr$data %>%
    map(select, any_of(expl_expr$variables))

  expl_expr$variables <- expl_expr$data %>%
    map(names)

  ## expression units

  expl_expr$txt_units <-
    c(tcga = 'log\u2082 gene count',
      imvigor = 'log\u2082 gene count',
      bcan = 'Z-score',
      hpa = 'IHC score')

  ## plotting data in a long format:
  ## for box plot panels and a heat map of Z-scores

  expl_expr$plot_data <-
    map2(expl_expr$data,
         expl_expr$variables,
         ~pivot_longer(.x,
                       cols = all_of(.y),
                       names_to = 'variable',
                       values_to = 'exp_level'))

  expl_expr$plot_data <- expl_expr$plot_data %>%
    map(left_join,
        expl_expr$lexicon,
        by = 'variable') %>%
    map(mutate,
        z_score = zScores(exp_level))

# Descriptive stats --------

  insert_msg('Descriptive stats')

  expl_expr$stats <- expl_expr$data %>%
    map2_dfr(., names(.),
             ~mutate(.x,
                     cohort = factor(.y, names(expl_expr$data)))) %>%
    fast_num_stats(variables = expl_expr$lexicon$variable,
                   split_fct = 'cohort') %>%
    map_dfc(~ifelse(stri_detect(.x, fixed = 'Inf'),
                    NA, .x))

  ## in the RNA seq data sets, all entries are complete

  expl_expr$stats[, c("tcga", "imvigor", "bcan")] <-
    expl_expr$stats[, c("tcga", "imvigor", "bcan")] %>%
    map_dfc(stri_replace, regex = '\\ncomplete.*', replacement = '')

  ## appending with the unit

  expl_expr$stats <-
    list(expl_expr$stats[1, ],
         expl_expr$txt_units %>%
           as.list %>%
           as_tibble %>%
           mutate(variable = 'Expression unit'),
         expl_expr$stats[-1, ]) %>%
    reduce(full_rbind)

  ## numbers of complete cases in the HPA data set

  expl_expr$hpa_complete <- expl_expr$data$hpa %>%
    colComplete

  expl_expr$hpa_complete <- expl_expr$hpa_complete %>%
    map2_chr(names(.), .,
             paste, sep = ', n = ') %>%
    set_names(names(expl_expr$hpa_complete))

# Box plot panels --------

  insert_msg('Box plot panels')

  expl_expr$plots <-
    list(x = expl_expr$plot_data,
         u = c('gray60', 'gray60', 'gray60', 'black'),
         t = c(0.5, 0.5, 0.5, 1),
         y = globals$cohort_labs[names(expl_expr$plot_data)],
         z = map_dbl(expl_expr$data, nrow),
         w = globals$cohort_axis_labs[names(expl_expr$plot_data)]) %>%
    pmap(function(x, u, t, y, z, w) x %>%
           ggplot(aes(x = exp_level,
                      y = reorder(variable,
                                  exp_level,
                                  FUN = function(x) quantile(x, 0.75, na.rm = TRUE)),
                      fill = protein_type)) +
           facet_grid(protein_type ~ .,
                      labeller = as_labeller(function(x) stri_replace(x,
                                                                      fixed = '_',
                                                                      replacement = '\n')),
                      scales = 'free',
                      space = 'free') +
           geom_boxplot(outlier.color = NA,
                        alpha = 0.35) +
           geom_point(shape = 16,
                      size = t,
                      alpha = 0.5,
                      position = position_jitter(width = 0.2,
                                                 height = 0.1),
                      color = u,
                      show.legend = FALSE) +
           scale_fill_manual(values = c('indianred3', 'plum4', 'aquamarine4'),
                             name = 'protein\ntype') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.title.x = element_markdown()) +
           labs(title = y,
                subtitle = paste('total: n =', z),
                x = w))

  ## additional styling

  expl_expr$plots[c("tcga", "imvigor", "bcan")] <-
    expl_expr$plots[c("tcga", "imvigor", "bcan")] %>%
    map(~.x + theme(axis.text.y = element_text(face = 'italic')))

  expl_expr$plots$hpa <- expl_expr$plots$hpa +
    scale_y_discrete(labels = expl_expr$hpa_complete)

# Heat map of expression Z-scores for the RNA seq data sets -------

  insert_msg('Heat map of expression Z-scores, RNA seq data sets')

  ## plotting data: expression medians and X axis labels with sample numbers

  expl_expr$hm_data <- expl_expr$plot_data[c("tcga", "imvigor", "bcan")] %>%
    map(group_by, variable, protein_type) %>%
    map(summarise,
        z_score = median(z_score, na.rm = TRUE)) %>%
    map(ungroup) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, names(expl_expr$plot_data)))

  expl_expr$hm_axis_labs <-
    map2(globals$cohort_labs[names(expl_expr$data)],
         map_dbl(expl_expr$data, nrow),
         paste, sep = '\nn = ') %>%
    set_names(names(expl_expr$data))

  ## heat map

  expl_expr$hm_plot <- expl_expr$hm_data %>%
    ggplot(aes(x = cohort,
               y = reorder(variable,
                           z_score,
                           FUN = function(x) median(x, na.rm = TRUE)),
               fill = z_score)) +
    facet_grid(protein_type ~ .,
               labeller = as_labeller(function(x) stri_replace(x,
                                                               fixed = '_',
                                                               replacement = '\n')),
               scales = 'free',
               space = 'free') +
    geom_tile() +
    scale_fill_gradient2(low = 'steelblue',
                         mid = 'black',
                         high = 'firebrick',
                         midpoint = 0,
                         #limits = c(-2, 2),
                         oob = scales::squish,
                         name = 'Z-score') +
    scale_x_discrete(labels = expl_expr$hm_axis_labs) +
    globals$common_theme +
    theme(axis.title = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(face = 'italic')) +
    labs(title = 'FGFR-, FGF-, and FGFBP-coding genes, RNA sequencing')

# END -------

  expl_expr$lexicon <- NULL
  expl_expr$hm_data <- NULL
  expl_expr$data <- NULL

  expl_expr <- compact(expl_expr)

  insert_tail()

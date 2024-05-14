# Functions for the project

# tools -------

  library(tidyverse)
  library(stringi)
  library(trafo)
  library(rlang)

  library(microViz)

  library(clustTools)
  library(jaccard)
  library(expm)

  library(ranger)

  library(cowplot)
  library(ggrepel)
  library(ggtext)

  set_rownames <- trafo::set_rownames

# data import ------

  annotate_raw_symbol <- function(data) {

    data %>%
      filter(!is.na(entrez_id),
             entrez_id != '') %>%
      mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                  column = 'SYMBOL',
                                  keys = entrez_id,
                                  keytype = 'ENTREZID')) %>%
      mutate(gene_symbol = unlist(gene_symbol)) %>%
      filter(gene_symbol != '') %>%
      filter(complete.cases(.))

  }

  extract_mutations <- function(data,
                                id_variable = 'Tumor_Sample_Barcode',
                                gene_id = 'Hugo_Symbol',
                                as_numeric = TRUE) {

    ## counts mutations per gene and sample

    data[[gene_id]] <- factor(data[[gene_id]])

    data_splits <- split(data[[gene_id]], data[[id_variable]])

    mut_counts <- data_splits %>%
      map(table)

    sample_ids <- names(mut_counts)

    mut_counts <- do.call('rbind', mut_counts)

    mut_counts <- ifelse(mut_counts > 0, 1, mut_counts)

    mut_counts <- as_tibble(mut_counts)

    if(!as_numeric) {

      mut_counts <- mut_counts %>%
        future_map_dfc(~car::recode(as.character(.x),
                                    "'0' = 'WT'; '1' = 'mutated'"),
                       .options = furrr_options(seed = TRUE)) %>%
        future_map_dfc(factor, c('WT', 'mutated'),
                       .options = furrr_options(seed = TRUE))

    }

    mut_counts %>%
      mutate(sample_id = sample_ids) %>%
      relocate(sample_id)

  }

  box_from_stats <- function(data,
                             width = 0.8,
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             y_lab = 'log<sub>2</sub> expression',
                             labeller = c('expression' = 'raw',
                                          'corr_expression' = 'COMBAT-corrected'),
                             scales = 'free') {

    ## plots a box plot given a data frame with the columns
    ## `median`, `perc25`, `perc75`, `perc025` and `perc975`
    ## storing the median expression, interquartile and 95 percentile ranges

    data <- data %>%
      mutate(center_pos = as.numeric(cohort))

    data %>%
      ggplot(aes(x = cohort,
                 y = median,
                 fill = cohort)) +
      facet_grid(type ~ .,
                 labeller = as_labeller(labeller),
                 scales = scales) +
      geom_errorbar(aes(ymin = perc025,
                        ymax = perc975),
                    width = 0) +
      geom_rect(aes(xmin = center_pos - 0.5 * width,
                    xmax = center_pos + 0.5 * width,
                    ymin = perc25,
                    ymax = perc75),
                color = 'black') +
      geom_segment(aes(x = center_pos - 0.5 * width,
                       xend = center_pos + 0.5 * width,
                       y = median,
                       yend = median),
                   color = 'black') +
      scale_fill_manual(values = globals$cohort_colors,
                        labels = globals$cohort_labs,
                        name = '') +
      scale_x_discrete(labels = globals$cohort_labs) +
      guides(fill = 'none') +
      globals$common_theme +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_markdown(),
            axis.text.x = element_text(hjust = 1,
                                       vjust = 0.5,
                                       angle = 90)) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           y = y_lab)

  }

# data transformation -------

  safely_mutate <- function(.data, ...) {

    ## enables for error-resistant change of data frame columns
    ## in case of any error, the genuine data frame is returned

    new_data <- try(mutate(.data, ...),
                    silent = TRUE)

    if(inherits(new_data, 'try-error')) return(.data)

    return(new_data)

  }

  safely_select <- function(.data, cols, fill = NULL) {

    ## selects columns with a character vector of column names
    ## in safe mode, i.e. without errors if the column does not exist
    ##
    ## if `fill` is specified, the non-existing columns are added and filled
    ## with the defined value

    stopifnot(is.character(cols))

    new_data <- select(.data, any_of(cols))

    if(is.null(fill)) return(new_data)

    absent_cols <- cols[!cols %in% names(.data)]

    for(i in absent_cols) {

      new_data[[i]] <- fill

    }

    new_data

  }

  minimum_shift <- function(df) {

    ## shifts values of variables of a numeric data frame by their minima

    min_values <- map_dbl(df, min, ns.rm = TRUE)

    map2_dfc(df, abs(min_values), `+`)

  }

# Result table formatting ------

  format_desc <- function(df) {

    ## formats the output of 'explore()' function

    levs <- names(df)

    df %>%
      reduce(left_join, by = 'variable') %>%
      set_names(c('variable', levs))

  }

  merge_stat_test <- function(stats, test, by = 'variable') {

    left_join(stats, test[c(by, 'significance', 'eff_size')], by = by)

  }

  format_res_tbl <- function(df,
                             dict = globals$clinic_lexicon,
                             remove_complete = FALSE, ...) {

    if(!is.null(dict)) {

      df <- df %>%
        exchange(variable = 'variable',
                 dict = dict, ...)

    }

    df <- df %>%
      map_dfc(stri_replace,
              fixed = 'Complete',
              replacement = 'complete') %>%
      map_dfc(stri_replace,
              regex = '^Mean.*\\nMedian\\s{1}=\\s{1}',
              replacement = '') %>%
      map_dfc(stri_replace,
              fixed = 'Range',
              replacement = 'range') %>%
      map_dfc(stri_replace,
              regex = '^no.*\\nyes:\\s{1}',
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^NA.*',
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^\\n.*',
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^\\n',
              replacement = '')

    if(!remove_complete) return(df)

    df %>%
      map_dfc(stri_replace,
              regex = '\\ncomplete.*',
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^complete.*',
              replacement = '')


  }

# Exploratory analysis -------

  bin2stats <- function(data,
                        split_factor = 'cohort') {

    ## frequency stats for the strata given a 0/1 data frame

    if(!is.factor(data[[split_factor]])) {

      data <- data %>%
        mutate(!!split_factor := factor(.data[[split_factor]]))

    }

    stats <- data[names(data) != 'sample_id'] %>%
      blast(all_of(split_factor), .skip = TRUE) %>%
      map(colSums) %>%
      map(compress,
          names_to = 'variable',
          values_to = 'n') %>%
      compress(names_to = split_factor)

    ## appending with the numbers of observations in the clusters
    ## and percentages

    totals <- count(data, .data[[split_factor]], name = 'n_total')

    stats <- left_join(stats, totals, by = split_factor)

    stats %>%
      mutate(percent = n/n_total * 100) %>%
      arrange(-percent)

  }

# Mutation type analysis -------

  find_domain <- function(position,
                          domain_tbl) {

    if(is.na(position)) return('not assigned')

    sel_vector1 <- position >= domain_tbl$start_position
    sel_vector2 <- position <= domain_tbl$end_position

    sel_vector <- sel_vector1 & sel_vector2

    if(sum(sel_vector) != 1) return('not assigned')

    domain_tbl[sel_vector, ][['description']][1]

  }

  append_domain_ <- function(mutation_data,
                             domain_tbl) {

    ## adds information on protein domain to the mutation data
    ## based on the protein position variable

    mutation_data$protein_domain <-
      mutation_data$protein_start_position %>%
      map_chr(find_domain, domain_tbl)

    mutation_data

  }

  append_domain <- function(mutation_data, domain_lst) {

    ## adds information on protein domain to the mutation data
    ## based on the protein position variable

    mutation_lst <- blast(mutation_data, gene_symbol)[names(domain_lst)]

    map2_dfr(mutation_lst,
             domain_lst,
             append_domain_)

  }

  plot_stack <- function(data,
                         x_var = 'percent',
                         y_var = 'gene_symbol',
                         pos_var = 'plot_pos',
                         fill_var = 'mutation_type',
                         y_limits = NULL,
                         plot_title = NULL,
                         plot_subtitle = NULL,
                         fill_title = 'Mutation type',
                         palette = NULL,
                         x_lab = '% of samples',
                         reverse_fill = TRUE) {

    # draws a stack plot with mutation or variant types

    if(reverse_fill) {

      data <- data %>%
        mutate(!!fill_var := factor(.data[[fill_var]], rev(levels(.data[[fill_var]]))))

    }

    if(is.null(y_limits)) {

      stack_plot <- data %>%
        ggplot(aes(x = .data[[x_var]],
                   y = reorder(.data[[y_var]], .data[[x_var]]),
                   fill = .data[[fill_var]]))

    } else {

      stack_plot <- data %>%
        ggplot(aes(x = .data[[x_var]],
                   y = .data[[y_var]],
                   fill = .data[[fill_var]])) +
        scale_y_discrete(limits = y_limits)

    }

    if(!is.null(palette)) {

      stack_plot <- stack_plot +
        scale_fill_manual(values = palette,
                          drop = FALSE)

    }

    stack_plot +
      geom_bar(position = position_stack(),
               color = 'black',
               stat = 'identity') +
      globals$common_theme +
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           fill = fill_title)

  }

  position_plot <- function(data,
                            position_var = 'protein_start_position',
                            split_by = c('mutation_type', 'variant_type'),
                            fill_var = NULL,
                            point_size = 2,
                            point_shape = 21,
                            point_alpha = 0.75,
                            point_hjitter = 0.1,
                            point_color = 'steelblue',
                            palette = NULL,
                            plot_title = NULL,
                            plot_subtitle = NULL,
                            x_lab = 'protein position',
                            y_lab = NULL,
                            palette_title = NULL,
                            reverse_split = TRUE,
                            protein_length = NULL, ...) {

    ## scatter plot: each point represents a single mutation
    ## split by the mutation or variant type

    ## metadata -----

    split_by <- match.arg(split_by[1],
                          c('mutation_type', 'variant_type'))

    if(reverse_split) {

      data <- data %>%
        mutate(!!split_by := factor(.data[[split_by]],
                                    rev(levels(.data[[split_by]]))))

    }

    if(is.null(plot_subtitle)) {

      plot_subtitle <- paste('total mutations: n =', nrow(data))

    }

    ## plotting -------

    if(is.null(fill_var)) {

      pos_plot <- data %>%
        ggplot(aes(x = .data[[position_var]],
                   y = .data[[split_by]]))

      if(point_shape == 21) {

        pos_plot <- pos_plot +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     position = position_jitter(height = point_hjitter),
                     fill = point_color,
                     color = 'black')

      } else {

        pos_plot <- pos_plot +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     position = position_jitter(height = point_hjitter),
                     color = point_color)

      }

    } else {

      if(point_shape == 21) {

        pos_plot <- data %>%
          ggplot(aes(x = .data[[position_var]],
                     y = .data[[split_by]],
                     fill = .data[[fill_var]])) +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     position = position_jitter(height = point_hjitter),
                     color = 'black') +
          scale_fill_manual(values = palette, ...)


      } else {

        pos_plot <- data %>%
          ggplot(aes(x = .data[[position_var]],
                     y = .data[[split_by]],
                     color = .data[[fill_var]])) +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     position = position_jitter(height = point_hjitter)) +
          scale_color_manual(values = palette, ...)

      }

    }

    if(!is.null(protein_length)) {

      pos_plot <- pos_plot + scale_x_continuous(limits = c(0, protein_length))

    }

    pos_plot +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab,
           fill = palette_title,
           color = palette_title)


  }

  plot_domains <- function(domain_tbl,
                           protein_name = 'FGFR3',
                           protein_length = globals$receptor_lengths["FGFR3"],
                           plot_title = html_italic('FGFR3'),
                           plot_subtitle = NULL,
                           x_lab = 'Protein position',
                           palette = globals$domain_colors,
                           fill_title = 'Protein domain') {

    ## draws a scheme of the domain structure

    domain_plot <- domain_tbl %>%
      ggplot(aes(xmin = start_position,
                 xmax = end_position,
                 y = protein_name,
                 ymin = 0.5,
                 ymax = 1.5,
                 fill = description))  +
      #annotate('segment',
       #        xmin = 0,
        #       xmax = protein_length,
         #      ymin = 1,
          #     ymax = 1) +
      geom_rect(color = 'black') +
      scale_x_continuous(limits = c(0, protein_length)) +
      scale_fill_manual(values = palette) +
      globals$common_theme +
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           fill = fill_title)

    domain_plot

  }

  append_domain_scheme <- function(position_plot,
                                   domain_tbl,
                                   protein_name = 'FGFR3',
                                   protein_length = globals$receptor_lengths["FGFR3"],
                                   palette = globals$domain_colors,
                                   fill_title = 'Protein domain',
                                   rm_position_legend = TRUE,
                                   rm_domain_legend = FALSE,
                                   rel_heights = c(1, 3),
                                   align = 'hv',
                                   cust_theme = globals$common_theme +
                                     theme(plot.title = element_markdown(),
                                           axis.title.x = element_blank()), ...) {

    ## adds a plot of the domain structure to a scatter plot of mutation
    ## positions

    domain_plot <-
      plot_domains(domain_tbl = domain_tbl,
                   protein_name = protein_name,
                   protein_length = protein_length,
                   plot_title = position_plot$labels$title,
                   plot_subtitle = position_plot$labels$subtitle,
                   x_lab = position_plot$labels$x,
                   palette = palette,
                   fill_title = fill_title) +
      cust_theme +
      theme(axis.title.y = element_blank())

    if(rm_position_legend) {

      position_plot <- position_plot + theme(legend.position = 'none')

    }

    if(rm_domain_legend) {

      domain_plot <- domain_plot + theme(legend.position = 'none')

    }

    plot_grid(domain_plot,
              position_plot +
                theme(plot.title = element_blank(),
                      plot.subtitle = element_blank()),
              nrow = 2,
              rel_heights = rel_heights,
              align = align,
              axis = 'tblr', ...)

  }

# Co-occurrence ------

  mds_genetics <- function(data,
                           distance_method = 'jaccard',
                           kdim = 2, ...) {

    ## performs multi-dimensional scaling of 0/1 coded genetic alteration data
    ## computes frequencies (percentages) of the alterations
    ## and a distance matrix as well

    ## distance matrix and MDS coordinates -------

    dist_mtx <- calculate_dist(data = data, method = distance_method)

    mds_tbl <- reduce_data(data = data,
                           distance_method = distance_method,
                           kdim = kdim,
                           red_fun = 'mds')

    mds_tbl <- mds_tbl %>%
      extract %>%
      relocate(observation)

    names(mds_tbl)[1] <- 'variable'

    ## frequencies ------

    freq_tbl <- data %>%
      rowSums %>%
      compress(names_to = 'variable',
               values_to = 'n') %>%
      mutate(n_total = ncol(data),
             percent = n/n_total * 100) %>%
      arrange(-percent)

    list(dist_mtx = dist_mtx,
         mds_tbl = mds_tbl,
         freq_tbl = freq_tbl)

  }

  plot_mds <- function(mds_output,
                       plot_title = NULL,
                       plot_subtitle = NULL,
                       palette = c('mutation' = unname(globals$mut_status_color["mutated"]),
                                   'amplification' = unname(globals$amp_status_color["amplified"]),
                                   'deletion' = unname(globals$del_status_color["deleted"])),
                       max_size = 6,
                       size_limits = c(0, 25),
                       point_size = 2,
                       point_alpha = 0.75,
                       point_hjittter = 0,
                       point_wjitter = 0,
                       txt_size = 2.75,
                       full_label = FALSE, ...) {

    ## draws a scatter plot of coordinates of MDS dimensions.
    ## Each point represents a single feature. If max_size is not NULL,
    ## point size codes for the feature frequency.

    ## plotting data and labels -------

    if(is.null(plot_subtitle)) {

      plot_subtitle <-
        paste('total samples: n =', mds_output$freq_tbl$n_total[[1]])

    }

    plot_data <- left_join(mds_output$mds_tbl,
                           mds_output$freq_tbl,
                           by = 'variable') %>%
      mutate(gene_symbol = stri_split_fixed(variable,
                                            pattern = '_',
                                            simplify = TRUE)[, 1],
             alteration = stri_split_fixed(variable,
                                           pattern = '_',
                                           simplify = TRUE)[, 2])

    ## plotting -------

    if(!is.null(max_size)) {

      sc_plot <- plot_data %>%
        ggplot(aes(x = comp_1,
                   y = comp_2,
                   fill = alteration,
                   color = alteration,
                   size = percent)) +
        geom_point(shape = 21,
                   color = 'black',
                   alpha = point_alpha) +
        scale_size_area(max_size = max_size,
                        limits = size_limits)

    } else {

      sc_plot <- plot_data %>%
        ggplot(aes(x = comp_1,
                   y = comp_2,
                   fill = alteration,
                   color = alteration)) +
        geom_point(shape = 21,
                   color = 'black',
                   size = point_size,
                   alpha = point_alpha)

    }

    if(full_label) {

      sc_plot <- sc_plot +
        geom_text_repel(aes(label = paste(gene_symbol, alteration)),
                        size = txt_size,
                        show.legend = FALSE, ...)

    } else {

      sc_plot <- sc_plot +
        geom_text_repel(aes(label = gene_symbol),
                        size = txt_size,
                        fontface = 'italic',
                        show.legend = FALSE, ...)

    }

    sc_plot +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = 'MDS1',
           y = 'MDS2',
           fill = '',
           size = '% of samples') +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette)

  }

  post_hoc_jaccard <- function(data, B = 1000, ...) {

    ## post-hoc bootstrap Jaccard test for overlap

    pairs <- combn(rownames(data), m = 2, simplify = FALSE)

    test <- pairs %>%
      future_map(~jaccard.test.bootstrap(x = data[.x[1], ],
                                         y = data[.x[2], ],
                                         verbose = FALSE,
                                         B = B, ...),
                 .options = furrr_options(seed = TRUE))

    map2_dfr(pairs, test,
             ~tibble(variable1 = .x[1],
                     variable2 = .x[2],
                     j = 1 - calculate_dist(data[.x, ], method = 'jaccard'),
                     p_value = .y$pvalue)) %>%
      re_adjust(method = 'BH') %>%
      mutate(pair_id = paste(variable1, variable2, sep = '|')) %>%
      relocate(pair_id)

  }

  plot_binary_overlap <- function(data,
                                  variable1,
                                  variable2,
                                  jaccard_test,
                                  scale_max = TRUE,
                                  title_suffix = 'GENIE BLCA',
                                  txt_size = 2.75) {

    ## heat map visualization of a contingency matrix

    ## plotting labels --------

    ax_labs <- list(x = variable1,
                    y = variable2) %>%
      map(stri_split_fixed, pattern = '_', simplify = TRUE) %>%
      map(~paste(html_italic(.x[, 1]), .x[, 2]))

    n_total <- nrow(data)

    jaccard_test <- jaccard_test %>%
      filter(.data[['variable1']] == .env[['variable1']],
             .data[['variable2']] == .env[['variable2']])

    plot_subtitle <- paste0('J = ', signif(jaccard_test$j, 2),
                            ', ', jaccard_test$significance)

    plot_title <- list(x = variable1,
                       y = variable2) %>%
      map(stri_split_fixed, pattern = '_', simplify = TRUE) %>%
      map_chr(~paste(html_italic(.x[, 1]), .x[, 2])) %>%
      paste(collapse = ':')

    if(is.null(title_suffix)) {

      plot_title <- paste(plot_title, title_suffix, sep = ', ')

    }

    ## plotting data -------

    count_tbl <- table(x_var = data[[variable1]],
                       y_var = data[[variable2]]) %>%
      as.data.frame %>%
      mutate(percent = Freq/n_total * 100,
             tile_lab = paste(signif(percent, 2), Freq, sep = '%\nn = '))

    ## heat map

    hm_plot <- count_tbl %>%
      ggplot(aes(x = factor(x_var),
                 y = factor(y_var),
                 fill = percent)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = tile_lab),
                size = txt_size) +
      globals$common_theme +
      scale_x_discrete(labels = c('0' = 'no', '1' = 'yes')) +
      scale_y_discrete(labels = c('0' = 'no', '1' = 'yes')) +
      theme(axis.title.x = element_markdown(),
            axis.title.y = element_markdown(),
            plot.title = element_markdown()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = ax_labs[['x']],
           y = ax_labs[['y']])

    if(!scale_max) {

      hm_plot <- hm_plot +
        scale_fill_gradient(low = 'white',
                            high = 'firebrick',
                            name = '% of samples')

    } else {

      plot_max <- count_tbl %>%
        filter(x_var == 1 | y_var == 1) %>%
        .$percent

      hm_plot <- hm_plot +
        scale_fill_gradient(low = 'white',
                            high = 'firebrick',
                            limits = c(0, max(plot_max)[1]),
                            name = '% of samples',
                            oob = scales::squish)

    }

    hm_plot

  }

# Random Forests ------

  tune_rf <- function(data,
                      formula,
                      tune_grid, ...) {

    ## tunes a Random Forest model by minimizing the OOB prediction error

    ## construction of the models ---------

    tune_grid <- tune_grid %>%
      mutate(model_id = paste0('rf_', 1:nrow(.)))

    tune_models <- tune_grid %>%
      select(-model_id) %>%
      pmap(ranger,
           formula = formula,
           data = data, ...) %>%
      set_names(tune_grid$model_id)

    ## OOB stats ----------

    oob_errors <- tune_models %>%
      map_dbl(~.x$prediction.error)

    tune_stats <- tune_grid %>%
      mutate(oob_error = oob_errors) %>%
      relocate(model_id, oob_error) %>%
      as_tibble

    best_tune <- tune_stats %>%
      filter(oob_error == min(oob_error, na.rm = TRUE))

    best_tune <- best_tune[1, ]

    tune_stats <- tune_stats %>%
      mutate(best = ifelse(model_id == best_tune$model_id,
                           'yes', 'no'))

    list(stats = tune_stats,
         best_tune = best_tune)

  }

  plot_ranger_tuning <- function(tune_stats,
                                 plot_title = NULL,
                                 plot_subtitle = NULL,
                                 x_lab = '# variables per random tree, mtry',
                                 y_lab = 'Classification error, OOB',
                                 split_color = c(gini = 'cornflowerblue',
                                                 extratrees = 'orangered3',
                                                 hellinger = 'darkolivegreen4',
                                                 variance = 'darkolivegreen4',
                                                 maxstat = 'steelblue'),
                                 split_labels = c(gini = 'Gini index',
                                                  extratrees = 'ExtraTrees',
                                                  hellinger = 'Hellinger',
                                                  variance = 'Variance',
                                                  maxstat = 'MaxStat'),
                                 split_title = 'Split rule',
                                 line_alpha = 0.5,
                                 point_alpha = 0.75,
                                 split_by_node_size = FALSE) {


    ## plots results of tuning of a Random Forest classifier or a regressor

    ## plot subtitle --------

    if(is.null(plot_subtitle)) {

      sub_tbl <- tune_stats %>%
        filter(best == 'yes')

      min_error <- paste('Minimal error =', signif(sub_tbl$oob_error, 2))

      sub_tbl <- sub_tbl %>%
        select(-model_id, -oob_error, -best)

      tune_params <- sub_tbl %>%
        as.list %>%
        map2_chr(names(.), .,
                 paste, sep = ' = ') %>%
        paste(collapse = ', ')

      plot_subtitle <- paste(min_error, tune_params, sep = '\n')

    }

    ## plotting -------

    tune_plot <- tune_stats %>%
      ggplot(aes(x = mtry,
                 y = oob_error,
                 fill = splitrule,
                 color = splitrule))

    if(split_by_node_size) {

      tune_plot <- tune_plot +
        facet_grid(. ~ factor(min.node.size),
                   scales = 'free',
                   space = 'free',
                   labeller = as_labeller(function(x) paste('min.node =', x)))

    }

    tune_plot <- tune_plot +
      geom_path(alpha = line_alpha) +
      geom_point(shape = 16,
                 size = 2,
                 alpha = point_alpha) +
      geom_smooth(se = FALSE,
                  show.legend = FALSE)  +
      scale_color_manual(values = split_color,
                         labels = split_labels,
                         name = split_title) +
      scale_fill_manual(values = split_color,
                        labels = split_labels,
                        name = split_title) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    tune_plot

  }

# END -------

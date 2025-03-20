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
  library(zoo)

  library(glmnet)
  library(DescTools)
  library(performanceEstimation)
  library(CalibrationCurves)
  library(ranger)

  library(cowplot)
  library(ggrepel)
  library(ggtext)

  library(xml2)
  library(curl)

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

# import of the HPA data -------

  download_hpa <- function(url, file, ...) {

    ## manages downloads from HPA

    idx <- try(download_xml(url = url, file = file, ...),
               silent = TRUE)

    if(inherits(idx, 'try-error')) message(paste('Failed for:', url))

  }

  extract_patients <- function(xml_obj) {

    ## extracts the nodes with patient/sample staining information

    as_list(xml_find_all(xml_obj, './/patient'))

  }

  cancer_search <- function(x, regex = '(U|u)rothe') {

    ## extracts elements of a branched list with HPA
    ## patient's data, that are matched by a regular expression
    ## applied to the sample/snomedParameters node

    snomeds <- x %>%
      map(~.x$sample$snomedParameters) %>%
      map(map, attr, 'tissueDescription') %>%
      map(unlist)

    lgl_vec <- snomeds %>%
      map(stri_detect, regex = regex) %>%
      map_lgl(any)

    x[lgl_vec]

  }

  cancer_format <- function(x) {

    ## coerces the branched list with patient's information and IHC
    ## staining scores into a handy table

    ## removal of patients with incomplete information,
    ## unique names for the list elements

    required_names <-
      c('sex', 'age', 'patientId', 'level', 'quantity', 'location', 'sample')

    complete_patients <- x %>%
      map(names) %>%
      map(~required_names %in% .x) %>%
      map_lgl(all)

    x <- x[complete_patients]

    ## stripping the most essential information

    x <- x %>%
      map(~.x[1:8]) %>%
      map(set_names,
          c('sex', 'age', 'patient_id',
            'staining', 'intensity', 'quantity', 'location',
            'sample'))

    patient_tbl <- x %>%
      map_dfr(as_tibble)

    ## plain data frame
    ## for each patient, we just need an entry
    ## with the snomed pathology description.
    ## It is an empty list with attributes.
    ## We also need information with image URLs stored in the non-empty
    ## elements

    patient_tbl[1:7] <- patient_tbl[1:7] %>%
      map_dfc(map_chr, ~.x[[1]])

    patient_tbl[[8]] <- patient_tbl[[8]] %>%
      map(~.x[1][[1]])

    empty_lst <- patient_tbl[[8]] %>%
      map_dbl(length)

    empty_lst <- empty_lst == 0

    ## extraction of image URLS

    image_tbl <- patient_tbl[!empty_lst, c('patient_id', 'sample')] %>%
      set_names(c('patient_id', 'image_url'))

    image_tbl[[2]] <- image_tbl[[2]] %>%
      map_chr(~unlist(.x[[1]]))

    ## appending the table with tissue descriptions

    patient_tbl <- patient_tbl[empty_lst, ]

    patient_tbl[[8]] <- patient_tbl[[8]] %>%
      map_chr(attr, 'tissueDescription')

    left_join(patient_tbl,
              image_tbl,
              by = 'patient_id')

  }

  process_hpa_entry <- function(xml_obj,
                                regex = '(U|u)rothe') {

    ## executes the complete extraction, cancer sample search, and
    ## formatting procedure for a HPA entry downloaded as a XML document

    xml_obj %>%
      extract_patients %>%
      cancer_search(regex = regex) %>%
      cancer_format

  }

  download_image <- function(url, file, ...) {

    ## fetches an image from HPA

    idx <- try(curl_download(url = url, destfile = file, ...),
               silent = TRUE)

    if(inherits(idx, 'try-error')) {

      message(paste('Failed for:', url))

      print(idx)

    }

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
        mutate(variable = ifelse(stri_detect(variable,
                                             regex = '^(Samples|Patients)'), #
                                 variable,
                                 exchange(variable,
                                          dict = dict, ...)))

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
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^yes:\\s{1}',
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
                            y_limits = NULL,
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

      pos_plot <- pos_plot +
        scale_x_continuous(limits = c(0, protein_length))

    }

    if(!is.null(y_limits)) {

      pos_plot <- pos_plot +
        scale_y_discrete(limits = y_limits)

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
                cust_theme +
                theme(plot.title = element_blank(),
                      plot.subtitle = element_blank(),
                      legend.position = 'none'),
              nrow = 2,
              rel_heights = rel_heights,
              align = align,
              axis = 'tblr', ...)

  }

  plot_hot_spot <- function(data,
                            color_variable = 'protein_domain',
                            color_palette = globals$domain_colors,
                            hot_cutoff = NULL,
                            plot_title = NULL,
                            plot_subtitle = NULL,
                            x_lab = 'Protein position',
                            y_lab = '% of cancer samples',
                            color_lab = 'protein\ndomain',
                            point_size = 2,
                            point_alpha = 1,
                            point_hjitter = 0.5,
                            point_shape = 21,
                            txt_size = 2.5,
                            txt_color = NULL, ...) {

    ## plot of percentages of mutations of particular protein positions

    ## plotting data and meta-data -------

    if(!is.null(hot_cutoff)) {

      data <- data %>%
        mutate(hot_label = ifelse(percent >= hot_cutoff,
                                  paste0(signif(percent, 2), '%\n',
                                         protein_change),
                                  NA),
               highlight = ifelse(percent >= hot_cutoff, 'yes', 'no'))

    }

    cohort_labeller <- data %>%
      filter(!duplicated(cohort))

    cohort_labeller <-
      paste(globals$cohort_labs[as.character(cohort_labeller$cohort)],
            cohort_labeller$n_total,
            sep = '\nn = ') %>%
      set_names(as.character(cohort_labeller$cohort))

    ## plotting ---------

    hot_plot <- data %>%
      ggplot(aes(x = protein_start_position,
                 y = percent,
                 color = .data[[color_variable]],
                 fill = .data[[color_variable]])) +
      geom_hline(yintercept = 0,
                 color = 'gray60',
                 linewidth = 0.25) +
      facet_grid(cohort ~ .,
                 labeller = as_labeller(cohort_labeller))

    if(is.null(hot_cutoff)) {

      if(point_shape %in% 21:25) {

        hot_plot <- hot_plot +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     color = 'black',
                     position = position_jitter(width = 0,
                                                height = point_hjitter))

      } else {

        hot_plot <- hot_plot +
          geom_point(shape = point_shape,
                     size = point_size,
                     alpha = point_alpha,
                     position = position_jitter(width = 0,
                                                height = point_hjitter))

      }

    } else {

      segment_data <- data %>%
        filter(highlight == 'yes')

      hot_plot <- hot_plot +
        geom_segment(aes(x = protein_start_position,
                         y = 0,
                         xend = protein_start_position,
                         yend = percent),
                     linewidth = 0.5) +
        scale_size_manual(values = c(no = 1,
                                     yes = point_size),
                          name = 'hotspot')

      if(point_shape %in% 21:25) {

        hot_plot <- hot_plot +
          geom_point(aes(size = highlight),
                     shape = point_shape,
                     alpha = point_alpha,
                     color = 'black',
                     position = position_jitter(width = 0,
                                                height = point_hjitter))

      } else {

        hot_plot <- hot_plot +
          geom_point(aes(size = highlight),
                     shape = point_shape,
                     alpha = point_alpha,
                     position = position_jitter(width = 0,
                                                height = point_hjitter))

      }

    }

    if(!is.null(hot_cutoff)) {

      if(is.null(txt_color)) {

        hot_plot <- hot_plot +
          geom_text_repel(aes(label = hot_label),
                          size = txt_size,
                          show.legend = FALSE, ...)

      } else {

        hot_plot <- hot_plot +
          geom_text_repel(aes(label = hot_label),
                          size = txt_size,
                          color = txt_color,
                          show.legend = FALSE, ...)

      }

    }

    hot_plot <- hot_plot +
      scale_color_manual(values = color_palette,
                         name = color_lab) +
      scale_fill_manual(values = color_palette,
                        name = color_lab) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    hot_plot

  }

  slide_window_density <- function(data,
                                   stat = c('mean', 'sum'),
                                   k = 5,
                                   from = 1,
                                   to = max(data$protein_start_position),
                                   padding = 2 * k,
                                   n_total = NULL) {

    ## takes a numeric vector of positions of a protein with mutations
    ## and counts mutations within a sliding window of width specified
    ## by k

    ## position data frame

    stat <- match.arg(stat[1], c('mean', 'sum'))

    event_counts <- count(data, protein_start_position)

    pos_df <-
      left_join(tibble(protein_start_position = from:to),
                event_counts,
                by = 'protein_start_position') %>%
      mutate(n = ifelse(is.na(n), 0, n))

    pos_df <- rbind(tibble(protein_start_position = rep(0, padding),
                           n = rep(0, padding)),
                    pos_df,
                    tibble(protein_start_position = rep(0, padding),
                           n = rep(0, padding)))

    if(stat == 'mean') {

      pos_df$density <- pos_df$n %>%
        rollmean(k = k, fill = NA)

    } else {

      pos_df$density <- pos_df$n %>%
        rollsum(k = k, fill = NA)

    }

    if(!is.null(n_total)) {

      pos_df$percent <- pos_df$density/n_total * 100

    }

    pos_df %>%
      filter(protein_start_position >= from,
             protein_start_position <= to)




  }


  plot_hot_density <- function(data,
                               plot_variable = c('percent', 'density'),
                               color_palette = globals$domain_colors,
                               color_lab = 'protein\ndomain',
                               plot_title = NULL,
                               plot_subtitle = NULL,
                               x_lab = 'Protein position',
                               y_lab = 'density, % of cancer samples',
                               point_size = 1,
                               line_color = 'gray80') {


    ## plots mutation density: raw or percent-scaled

    ## plotting data and metadata -------

    cohort_labeller <- data %>%
      filter(!duplicated(cohort))

    cohort_labeller <-
      paste(globals$cohort_labs[as.character(cohort_labeller$cohort)],
            cohort_labeller$n_total,
            sep = '\nn = ') %>%
      set_names(as.character(cohort_labeller$cohort))

    plot_variable <- match.arg(plot_variable[1], c('percent', 'density'))

    ## plot -------

    dens_plot <- data %>%
      ggplot(aes(x = protein_start_position,
                 y = .data[[plot_variable]])) +
      facet_grid(cohort ~ .,
                 labeller = as_labeller(cohort_labeller)) +
      geom_hline(yintercept = 0,
                 linewidth = 0.25,
                 color = 'gray60') +
      geom_path(color = line_color) +
      geom_point(aes(color = protein_domain),
                 shape = 16,
                 size = 1) +
      scale_color_manual(values = color_palette,
                         name = color_lab) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    dens_plot

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

# ML modeling ---------

  augment_multi <- function(data,
                            minorities = c('LumU', 'LumNS', 'Stroma-rich'),
                            k = 7,
                            perc.over = rep(4, 3),
                            perc.under = rep(1, 3)) {

    ## multi-class augmentation

    data_lst <- minorities %>%
      map(~mutate(data,
                  consensusClass = ifelse(consensusClass == .x, .x, 'other'),
                  consensusClass = factor(consensusClass, c('other', .x)))) %>%
      set_names(minorities)

    synth_lst <-
      list(data = data_lst,
           perc.over = perc.over,
           perc.under = perc.under) %>%
      pmap(smote,
           form = consensusClass ~ .) %>%
      map_dfr(filter,
              consensusClass %in% minorities)

    major_tbl <- data %>%
      filter(!consensusClass %in% minorities)

    if(nrow(major_tbl) != 0) synth_lst <- rbind(major_tbl, synth_lst)

    synth_lst %>%
      mutate(consensusClass = factor(consensusClass,
                                     levels(data$consensusClass)))

  }

  remove_missing <- function(data) {

    ## removes observations with missing expression data

    row_complete <- rowMissing(data[names(data) != 'consensusClass'])

    row_complete <- row_complete == 0

    data[row_complete, ]

  }

  tune_glmnet <- function(x, y, indexes, ...) {

    ## finds optimal lambda by repeated cross-validation

    ## tuning --------

    mod_lst <- indexes %>%
      future_map(function(fold) cv.glmnet(x = x,
                                          y = y,
                                          foldid = fold, ...),
                 .options = furrr_options(seed = TRUE))

    mod_stats <- mod_lst %>%
      map(~as_tibble(.x[c('lambda', 'cvm')])) %>%
      map(filter, cvm == min(cvm)) %>%
      map_dfr(~.x[1, ])

    best_tune <- mod_stats %>%
      filter(cvm == min(cvm))

    ## training ---------

    glmnet_model <- glmnet(x = x,
                           y = y,
                           lambda = best_tune$lambda, ...)

    ## output --------

    list(stats = mod_stats,
         best_tune = best_tune,
         model = glmnet_model)

  }

  tune_ranger <- function(x, y,
                          tune_grid,
                          tune_stat = c('prediction.error',
                                        'mse',
                                        'rsq'), ...) {

    ## tunes the Ranger RF model by minimizing the OOB prediction error
    ## or maximizing R^2

    ## tuning models and stats ---------

    tune_stat <-
      match.arg(tune_stat[1], c('prediction.error', 'mse', 'rsq'))

    tune_grid <- tune_grid %>%
      mutate(model_id = paste0('model_', 1:nrow(.)))

    tune_models <- tune_grid %>%
      select(- model_id) %>%
      future_pmap(ranger,
                  x = x,
                  y = y,
                  ...,
                  .options = furrr_options(seed = TRUE))

    tune_stats <- tune_models %>%
      map(~.x[c('prediction.error', 'r.squared', 'prediction.error')]) %>%
      map(set_names, c('mse_oob', 'rsq_oob', 'class_error_oob')) %>%
      map(compact) %>%
      map_dfr(as_tibble)

    tune_stats <- cbind(tune_grid, tune_stats) %>%
      as_tibble

    if(tune_stat == 'mse') {

      best_tune <- tune_stats %>%
        filter(mse_oob == min(mse_oob))

    } else if(tune_stat == 'mse') {

      best_tune <- tune_stats %>%
        filter(rsq_oob == max(rsq_oob))

    } else {

      best_tune <- tune_stats %>%
        filter(class_error_oob == min(class_error_oob))

    }

    best_tune <- best_tune[1, ]

    ## the best-performing model

    ranger_args <-
      best_tune[!names(best_tune) %in% c('mse_oob',
                                           'rsq_oob',
                                           'class_error_oob',
                                           'model_id')]

    if(tune_stat == 'prediction.error') {

      best_model <- call2(.fn = 'ranger',
                          x = x,
                          y = y,
                          probability = TRUE,
                          !!!ranger_args,
                          ...)

    } else {

      best_model <- call2(.fn = 'ranger',
                          x = x,
                          y = y,
                          !!!ranger_args,
                          ...) %>%
        eval

    }

    best_model <- eval(best_model)

    ## output -------

    list(stats = tune_stats,
         best_tune = best_tune,
         model = best_model)

  }

  oof_glmnet <- function(x_train,
                         x_test,
                         y_train,
                         y_test, ...) {

    ## retrieves OOF predictions for lists of X matrices and Y response vectors

    ## models ---------

    models <-
      list(x = x_train,
           y = y_train) %>%
      pmap(glmnet, ...)

    ## predictions -------

    oof_probs <-
      list(object = models,
           newx = x_test) %>%
      pmap(predict, type = 'response') %>%
      map(~.x[, , 1]) %>%
      map(as.data.frame) %>%
      map(~set_colnames(.x, paste0(colnames(.x), '_elnet')))

    oof_classes <-
      list(object = models,
           newx = x_test) %>%
      pmap(predict, type = 'class')

    oof_preds <-
      list(x = oof_probs,
           y = oof_classes,
           z = y_test) %>%
      pmap(function(x, y, z) x %>%
             mutate(obs = factor(y, levels(y_test[[1]])),
                    pred = z))

    oof_preds %>%
      map(rownames_to_column, 'sample_id') %>%
      compress(names_to = 'fold_id')

  }

  oof_ranger <- function(x_train,
                         x_test,
                         y_train,
                         y_test, ...) {

    ## computes OOF predictions for a ranger model

    ## models --------

    models <-
      list(x = x_train,
           y = y_train) %>%
      pmap(ranger,
           probability = TRUE, ...)

    class_levels <- levels(y_train[[1]])

    ## predictions --------

    oof_probs <-
      list(object = models,
           data = x_test) %>%
      pmap(predict) %>%
      map(~.x$predictions) %>%
      map(as.data.frame) %>%
      map(~set_colnames(.x, paste0(colnames(.x), '_ranger')))

    oof_classes <- list()

    for(i in seq_along(oof_probs)) {

      oof_classes[[i]] <- 1:nrow(oof_probs[[i]]) %>%
        map(~which(oof_probs[[i]][.x, ] == max(oof_probs[[i]][.x, ]))) %>%
        map_dbl(~.x[[1]])

      oof_classes[[i]] <-
        factor(class_levels[oof_classes[[i]]], class_levels)

    }

    oof_preds <-
      list(x = oof_probs,
           y = oof_classes,
           z = y_test,
           w = x_test) %>%
      pmap(function(x, y, z, w) x %>%
             mutate(obs = factor(y, levels(y_test[[1]])),
                    pred = z,
                    sample_id = rownames(w)))

    oof_preds %>%
      compress(names_to = 'fold_id')

  }

  brier_squares <- function(prediction) {

    ## calculates Brier squares, i.e. squared differences between the observed
    ## outcome and the class assignment probability

    obs_levs <- levels(prediction[['obs']])

    if(length(obs_levs) == 2) {

      pos_class <- obs_levs[2]

      prediction <- mutate(prediction,
                           .num_outcome = as.numeric(obs) - 1,
                           square_dist = (.num_outcome - .data[[pos_class]])^2)

      prediction <- prediction %>%
        select(-.num_outcome)

    } else {

      outcomes_num <-
        as.data.frame(Dummy(prediction[['obs']],
                            method = 'full',
                            levels = obs_levs))

      fitted_num <- prediction[obs_levs]

      sqd <- map2(outcomes_num, fitted_num,
                  function(x, y) (x - y)^2)

      sqd <- reduce(sqd, function(x, y) x + y)

      prediction <-
        mutate(prediction, square_dist = sqd)

    }

    prediction

  }

  brier_reference <- function(prediction, n = 100) {

    ## calculates reference Brier squares expected for a NULL model

    if(n == 1) {

      resamp_data <- prediction

      resamp_data[['obs']] <- resamp_data[['obs']] %>%
        sample(size = nrow(prediction), replace = FALSE)

      resamp_data <- brier_squares(resamp_data)

      prediction[['square_ref']] <- resamp_data[['square_dist']]

      return(prediction)

    } else {

      refs <- map(1:n,
                  function(x) brier_reference(prediction, n = 1)) %>%
        map(~.x$square_ref)

      refs <- reduce(refs, `+`)/n

      prediction[['square_ref']] <- refs

      return(prediction)

    }

  }

  caret_stats <- function(prediction) {

    ## a helper wrapper around multiClassSummary()

    ## re-leveling, because caret assumes the response of interest as the
    ## first level

    levs <- rev(levels(prediction[['obs']]))

    prediction[c("obs", "pred")] <- prediction[c("obs", "pred")] %>%
      map_dfc(factor, levs)

    caret_stats <- multiClassSummary(prediction, lev = levs)

    caret_stats <- c(caret_stats,
                     brier_score = mean(prediction[['square_dist']]),
                     brier_reference = mean(prediction[['square_ref']]))

    caret_names <- names(caret_stats)

    matrix(caret_stats, nrow = 1) %>%
      as.data.frame %>%
      set_names(caret_names) %>%
      mutate(Kappa = ifelse(Kappa < 0, 0, Kappa)) %>%
      as_tibble


  }

  ml_summary <- function(prediction) {

    ## computes numeric metrics of classification model performance

    ## predictions for the training, canonical CV data

    if(!'resample' %in% names(prediction)) {

      return(caret_stats(prediction))

    }

    if(!all(stri_detect(prediction$resample, fixed = 'Rep'))) {

      stats <- prediction %>%
        blast(resample) %>%
        map_dfr(caret_stats) %>%
        map_dfc(mean)

    } else {

      prediction <- prediction %>%
        mutate(fold_no = stri_split_fixed(resample,
                                          pattern = '.',
                                          simplify = TRUE)[, 1],
               rep_no = stri_split_fixed(resample,
                                         pattern = '.',
                                         simplify = TRUE)[, 2])

      stats <- prediction %>%
        blast(rep_no) %>%
        map(blast, fold_no) %>%
        map(map_dfr, caret_stats) %>%
        map_dfr(map_dfc, mean) %>%
        map_dfc(mean)

    }

    return(stats)

  }

  plot_confusion <- function(mtx,
                             plot_title = NULL,
                             x_lab = 'observed class',
                             y_lab = 'predicted class',
                             txt_size = 2.3) {

    ## heat map visualization of confusion matrices

    plot_tbl <- mtx %>%
      as.data.frame %>%
      mutate(count = Freq,
             percent = count/sum(count) * 100,
             plot_lab = paste0(signif(percent, 2), '%\n(n = ',
                               count, ')'))

    plot_tbl %>%
      ggplot(aes(x = obs,
                 y = pred,
                 fill = .data[['percent']])) +
      geom_tile(color = 'black') +
      geom_text(aes(label = plot_lab),
                size = txt_size) +
      scale_fill_gradient(low = 'white',
                          high = 'firebrick',
                          limits = c(0, 100),
                          name = '% of all samples') +
      globals$common_theme +
      labs(title = plot_title,
           x = x_lab,
           y = y_lab)

  }

  plot_calibration <- function(prediction,
                               discrete = TRUE,
                               smooth = FALSE,
                               logistic.cal = TRUE,
                               bins = NULL,
                               plot_title = NULL,
                               x_lab = 'mean predicted probability',
                               y_lab = 'observed fraction',
                               fill_title = 'consensus\nclass',
                               point_size = 1,
                               line_size = 0.5,
                               color_palette = globals$sub_colors, ...) {

    ## plots a calibration curve for a multinomial model.
    ## the general idea is taken from Python's sklearn calibration_curve()
    ## and R package calibrationCurves, for discrete and flexible
    ## calibration, respectively.
    ## Dots stand for arguments passed to valProbggplot

    ## plotting data ----------

    ## list of predictions for each of the classes

    classes <- levels(prediction$obs)

    pred_lst <- classes %>%
      map(~transmute(prediction,
                     p_value = .data[[.x]],
                     class = ifelse(obs == .x, 1, 0))) %>%
      set_names(classes)

    ## plotting: discrete bins ------------

    if(discrete) {

      ## binning: simple heuristics is used by
      ## default (Sturges). By we're trying to coerce it
      ## to have at least one positive class assignment in a bin

      if(is.null(bins)) {

        class_n <- table(prediction$obs)

        bins <- ceiling(1 + log2(class_n))

        bin_n <- list(x = pred_lst,
                      y = bins) %>%
          pmap(function(x, y) x %>%
                 mutate(p_strata = cut(p_value,
                                       breaks = y))) %>%
          map(group_by, p_strata) %>%
          map(mutate, n_bin = sum(class)) %>%
          map(~.x$n_bin)

        bins <-
          map2(bin_n, bins,
               function(x, y) if(any(x == 0)) y - 1 else y)

      }

      pred_lst <-
        list(x = pred_lst,
             y = bins) %>%
        pmap(function(x, y) x %>%
               mutate(p_strata = cut(p_value,
                                     breaks = y)))


      ## mean probabilities

      bin_means <- pred_lst %>%
        map(group_by, p_strata) %>%
        map(summarise,
            mean_predicted = mean(p_value, na.rm = TRUE),
            mean_observed = mean(class, na.rm = TRUE),
            n_bin = sum(class)) %>%
        map(ungroup)

      bin_means <- bin_means %>%
        compress(names_to = 'class')

      ## plot

      cal_plot <- bin_means %>%
        ggplot(aes(x = mean_predicted,
                   y = mean_observed,
                   color = class)) +
        geom_abline(slope = 1,
                    intercept = 0,
                    linetype = 'dashed')

      if(smooth) {

        cal_plot <- cal_plot +
          geom_smooth(..., se = FALSE)

      } else {

        cal_plot <- cal_plot +
          geom_path()

      }

      cal_plot <- cal_plot +
        geom_point(shape = 16,
                   size = point_size) +
        scale_color_manual(values = color_palette,
                           name = fill_title) +
        globals$common_theme +
        labs(title = plot_title,
             x = x_lab,
             y = y_lab)

      return(cal_plot)

    }

    ## plotting: non-discrete calibration curves -------------

    plot_tbl <-
      list(y = map(pred_lst, ~.x$class),
           p = map(pred_lst, ~.$p_value)) %>%
      pmap(valProbggplot,
           logistic.cal = logistic.cal, ...)

    if(logistic.cal) {

      plot_tbl <- plot_tbl %>%
        map(~.x$CalibrationCurves[['LogisticCalibration']])

    } else {

      plot_tbl <- plot_tbl %>%
        map(~.x$CalibrationCurves[['FlexibleCalibration']])

    }

    plot_tbl <- plot_tbl %>%
      compress(names_to = 'class') %>%
      mutate(class = factor(class, classes))

    cal_plot <- plot_tbl %>%
      ggplot(aes(x = x,
                 y = y,
                 color = class)) +
      geom_abline(slope = 1,
                  intercept = 0,
                  linetype = 'dashed') +
      geom_path() +
      scale_color_manual(values = color_palette,
                         name = fill_title) +
      globals$common_theme +
      labs(title = plot_title,
           x = x_lab,
           y = y_lab)

    return(cal_plot)

  }

  predict_loocv <- function(caret_model) {

    ## extracts the observed outcome,
    ## and predictions in the training data set and resamples
    ## from a caret model

    ## observed outcome and sample identifiers

    obs <- caret_model$trainingData$.outcome
    sample_ids <- rownames(caret_model$trainingData)

    ## predictions for the training data

    train_preds <- data.frame(sample_id = sample_ids,
                              obs = obs,
                              pred = unname(predict(caret_model)))

    ## predictions for the resamples

    cv_preds <- caret_model$pred %>%
      mutate(sample_id = sample_ids[as.numeric(rowIndex)]) %>%
      select(sample_id, obs, pred)


    list(training = train_preds,
         cv = cv_preds)

  }

  ml_summary_loocv <- function(pred_lst) {

    ## computes performance stats for regression models:
    ## training data set and OOF predictions in LOOCV

    caret_stats <- pred_lst %>%
      map(~postResample(.x$pred, .x$obs)) %>%
      map(matrix, nrow = 1) %>%
      map(set_colnames, c('rmse', 'caret_rsq', 'mae')) %>%
      map(as_tibble) %>%
      compress(names_to = 'data_set')

    ## correlation coefficients between the predicted and observed
    ## outcome

    corr_stats <- pred_lst %>%
      map(~tibble(pearson = cor(.x[['pred']], .x[['obs']]),
                  spearman = cor(.x[['pred']], .x[['obs']],
                                 method = 'spearman'))) %>%
      compress(names_to = 'data_set')

    ## pseudo-R-square

    variances <- pred_lst %>%
      map(~.x$obs) %>%
      map(~(.x - mean(.x))^2)

    mse <- pred_lst %>%
      map(~(.x$obs - .x$pred)^2)

    rsq <- map2(mse, variances,
                ~1 - sum(.x)/sum(.y)) %>%
      map_dbl(~ifelse(.x < 0, 0, .x))


    left_join(caret_stats,
              corr_stats,
              by = 'data_set') %>%
      mutate(rsq = rsq)

  }

  my_shap_bee <- function(shapviz_obj,
                          variables = NULL,
                          top_variables = 10,
                          title_suffix = NULL,
                          point_size = 2,
                          point_alpha = 0.75,
                          point_jitter = 0.1,
                          show_violin = TRUE,
                          violin_alpha = 0.25,
                          violin_fill = 'cornsilk',
                          panels = FALSE,
                          rel_widths = c(0.65, 0.35),
                          add_legend = TRUE,
                          rel_heights = c(0.85, 0.15),
                          n_breaks = NULL,
                          cust_theme = globals$common_theme) {

    ## a customized bee-swarm plot of SHAP values as a function of feature
    ## value (min/max) for multinomial classification models

    class_levels <- names(shapviz_obj)

    plot_titles <- class_levels

    if(!is.null(title_suffix)) {

      plot_titles <- paste(plot_titles, title_suffix, sep = ', ')

    }

    ## plotting data: min/max-normalized feature values and SHAP --------

    plot_data <- list()

    plot_data$ft <- get_feature_values(shapviz_obj) %>%
      map(min_max)

    plot_data$shap <- get_shap_values(shapviz_obj)

    plot_data$shap <-
      map2(plot_data$shap,
           plot_data$ft,
           ~set_rownames(.x, rownames(.y)))

    if(!is.null(variables)) {

      plot_data$shap <- plot_data$shap %>%
        map(~.x[, variables])

    }

    ft_names <- plot_data$ft %>%
      map(colnames)

    ## mean absolute values of SHAP per feature as a variable importance metric -------

    ## the correlation of SHAP and value will be used to classify the
    ## explanatory factors as variables that associate positively and negatively
    ## with the outcome

    shap_imp <- list()

    shap_imp$absmean <- plot_data$shap %>%
      map(abs) %>%
      map(colMeans) %>%
      map(compress,
          names_to = 'variable',
          values_to = 'meanabs_shap')

    shap_imp$corr <-
      map2(plot_data$ft,
           plot_data$shap,
           function(x, y) map2_dbl(as.data.frame(x),
                                   as.data.frame(y),
                                   cor,
                                   method = 'spearman')) %>%
      map(compress,
          names_to = 'variable',
          values_to = 'corr')

    shap_imp <- shap_imp %>%
      map(compress, names_to = 'class') %>%
      reduce(left_join, by = c('class', 'variable')) %>%
      mutate(sign = ifelse(corr < 0, 'negative',
                           ifelse(corr > 0, 'positive', 'ns')),
             sign = factor(sign, c('positive', 'negative', 'ns')),
             class = factor(class, class_levels)) %>%
      blast(class)

    ## top most influential variables to be shown in the plots ------

    ## per class, by mean(abs(shap))

    top_vars <- shap_imp %>%
      map(slice_max, meanabs_shap, n = top_variables) %>%
      map(~.x$variable)

    plot_data <- plot_data %>%
      transpose %>%
      map(map, as.data.frame) %>%
      map(map, rownames_to_column, 'sample_id') %>%
      map2(., top_vars,
           ~.x %>%
             map(select, sample_id, all_of(.y))) %>%
      map(map, as_tibble)

    shap_imp <-
      map2(shap_imp,
           top_vars,
           function(x, y) x %>%
             filter(variable %in% y) %>%
             mutate(variable = factor(variable, rev(y))))

    ## long format of the plotting data, variable order, association sign ---------

    plot_data <- plot_data %>%
      map(function(class) map2(class,
                               c('value', 'shap'),
                               ~pivot_longer(.x,
                                             cols = all_of(names(.x)[-1]),
                                             values_to = .y,
                                             names_to = 'variable'))) %>%
      map(reduce, left_join, by = c('sample_id', 'variable'))

    plot_data <-
      map2(plot_data, top_vars,
           ~mutate(.x,
                   variable = factor(variable, rev(.y))))

    plot_data <-
      map2(plot_data,
           map(shap_imp, ~.x[c('variable', 'sign')]),
           left_join, by = 'variable')

    ## common plotting ranges for the plots

    bee_range <- plot_data %>%
      map(~.x$shap) %>%
      reduce(c) %>%
      range

    imp_range <- shap_imp %>%
      map(~.x$meanabs_shap) %>%
      reduce(c) %>%
      range

    ## bee-swarm plots --------

    bee_plots <-
      list(x = plot_data,
           y = plot_titles) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = shap,
                        y = variable,
                        color = value)) +
             facet_grid(sign ~ .,
                        space = 'free',
                        scales = 'free') +
             labs(title = y,
                  x = 'SHAP value'))

    if(show_violin) {

      bee_plots <- bee_plots %>%
        map(~.x +
              geom_violin(scale = 'width',
                          fill = violin_fill,
                          alpha = violin_alpha,
                          color = 'black',
                          trim = FALSE))

    }

    bee_plots <- bee_plots %>%
      map(~.x +
            geom_vline(xintercept = 0, linetype = 'dashed') +
            geom_point(shape = 16,
                       size = point_size,
                       alpha = point_alpha,
                       position = position_jitter(width = 0,
                                                  height = point_jitter)) +
            cust_theme +
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_text(face = 'italic')) +
            scale_x_continuous(limits = bee_range) +
            scale_color_gradient2(low = 'steelblue',
                                  mid = 'gray80',
                                  high = 'firebrick',
                                  limits = c(0, 1),
                                  midpoint = 0.5,
                                  name = 'min/max scaled\nfeature value'))

    ## bar plots of mean absolute variable importance --------

    imp_plots <-
      list(x = shap_imp,
           y = plot_titles) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = meanabs_shap,
                        y = variable,
                        fill = sign)) +
             facet_grid(sign ~ .,
                        scales = 'free',
                        space = 'free') +
             geom_bar(color = 'white',
                      stat = 'identity') +
             geom_vline(xintercept = 0,
                        linetype = 'dashed') +
             scale_x_continuous(limits = c(0, imp_range[[2]]),
                                n.breaks = n_breaks) +
             scale_fill_manual(values = c(positive = 'firebrick',
                                          negative = 'steelblue',
                                          ns = 'gray60'),
                               name = 'association\nwith response') +
             cust_theme +
             theme(axis.title.y = element_blank(),
                   axis.text.y = element_text(face = 'italic')) +
             labs(title = y,
                  x = 'mean |SHAP value|'))

    ## optional: plot panels --------

    if(!panels) return(list(bee = bee_plots,
                            bar = imp_plots))

    panels <-
      map2(bee_plots,
           imp_plots,
           ~plot_grid(.x +
                        theme(legend.position = 'none',
                              strip.background.y = element_blank(),
                              strip.text.y = element_blank()),
                      .y +
                        labs(title = '') +
                        theme(legend.position = 'none',
                              axis.text.y = element_blank(),
                              strip.background.y = element_blank(),
                              strip.text.y = element_blank()),
                      ncol = 2,
                      align = 'h',
                      axis = 'tblr',
                      rel_widths = rel_widths))

    legend_grob <-
      get_legend(bee_plots[[1]] + theme(legend.position = 'bottom'))

    if(add_legend) {

      panels <- panels %>%
        map(~.x %>%
              plot_grid(legend_grob,
                        nrow = 2,
                        rel_heights = rel_heights))

      return(panels)

    }

    return(list(panels = panels,
                legend = legend_grob))

  }

# Representative HPA images -------

  choose_representative <- function(x,
                                    probs = c(0, 0.5, 1),
                                    break_ties = TRUE) {

    ## chooses indexes of a vector that correspond to values as close as possible
    ## to quantiles specified by probs argument

    ## naming the vector with its original numbers

    x <- set_names(x, 1:length(x))

    ## vectors of absolute differences; finding its minima
    ## and the corresponding values

    quants <- quantile(x, probs, na.rm = TRUE)

    diff_lst <- quants %>%
      map(~abs(x - .x))

    diff_min <- diff_lst %>%
      map_dbl(min)

    diff_selection <-
      map2(diff_lst, diff_min, ~.x[.x == .y]) %>%
      map(names) %>%
      map(as.numeric)

    if(!break_ties) return(diff_selection)

    ## tie breaking ------

    diff_single <- diff_selection %>%
      map_dbl(function(x) if(length(x) == 1) x else sample(x, 1))

    attempts <- 1

    while((length(unique(diff_single)) < length(diff_single)) & attempts < 10) {

      diff_single <- diff_selection %>%
        map_dbl(function(x) if(length(x) == 1) x else sample(x, 1))

      attempts <- attempts + 1

    }

    diff_single

  }

  stitch_hpa_images <- function(lexicon,
                                caption_type = c('none', 'full', 'short'),
                                txt_size = 7,
                                gene_size = 8,
                                y_pos = c(0, 0.9)) {

    ## stitches an image panel

    caption_type <- match.arg(caption_type[1],
                              c('none', 'full', 'short'))

    caps <-
      switch(caption_type,
             none = NULL,
             full = lexicon$caption,
             short = lexicon$caption_short)

    image_lst <-
      list(x = lexicon$image_path,
           y = paste('IHC:', lexicon$ihc_score),
           z = paste('ID:', lexicon$patient_id)) %>%
      pmap(function(x, y, z) ggdraw() +
             draw_image(x) +
             draw_text(y,
                       x = 0.05,
                       y = y_pos[2],
                       size = txt_size,
                       hjust = 0, vjust = 0) +
             draw_text(z,
                       x = 0.05,
                       y = y_pos[1],
                       size = txt_size,
                       hjust = 0, vjust = 0))

    plot_panel <- image_lst %>%
      plot_grid(plotlist = .,
                ncol = nrow(lexicon),
                labels = caps,
                label_size = 8) %>%
      plot_grid(ggdraw() +
                  draw_text(lexicon$gene_symbol[[1]],
                            size = gene_size,
                            angle = 90),
                .,
                ncol = 2,
                rel_widths = c(0.1, 0.9))

    plot_panel

  }

# Drug resistance analyses --------

  classify_moa <- function(x) {

    # creates a simplified classification of drugs based on
    # mechanism of action and molecular targets

    ## classification criteria

    class_lst <-
      list('DNA' = c('DNA', 'Genome', 'genome', 'replication', 'nucleo',
                     'topoisomerase', 'helicase', 'telomerase', 'pyrimidine',
                     'Nelarabine', 'platin', 'mitomycin', 'NAMPT', 'Ara-'),
           'cytoskeleton' = c('Cytoskeleton', 'cytoskeleton', 'microtubule',
                              'Microtubulin', 'kinesin', 'Kinesin', 'tubulin',
                              'Tubulin', 'Vincristine', 'vincristine',
                              'PACLITAX', 'paclitax',
                              'docetax'),
           'GF signaling' = c('RTK', 'kit', 'KIT', 'Abl', 'ABL',
                              'ABL', 'IGF1R', 'ephrin', 'IRAK', 'PAK',
                              'Ephrin', 'EPH', 'FGFR', 'PI3K', 'mTOR', 'MTOR',
                              'RAF', 'RAS', 'PKC', 'Bruton', 'FLT', 'flt', 'MET',
                              'STAT3', 'SRC', 'WNT', 'focal\\s{1}adhes',
                              'insulin', 'rho', 'exportin', 'BTK', 'FYN', 'LCK',
                              'TAK', 'NIK', 'PIM', 'AKT'),
           'MAPK' = c('MAPK', 'ERK', 'MEK', 'p38', 'JNK', 'Jnk'),
           'EGFR/ERBB' = c('EGFR', 'Egfr', 'ERBB', 'Erbb', 'HER'),
           'VEGFR/PDGFR' = c('VEGF', 'PDGFR', 'KDR'),
           'cell cycle' = c('cell\\s{1}cycle', 'Cell\\s{1}cycle',
                            'aurora', 'Aurora', 'CDK', 'PARP',
                            'Mitosis', 'mitosis', 'cyclin', 'Cyclin',
                            'polo-like', 'PLK', 'p53', 'TP53', 'AURK', 'MDM2',
                            'ADP-ribose', 'GADD34'),
           'cell death' = c('apoptosis', 'Apoptosis',
                            'Bcl', 'BCL', 'Bid', 'BID', 'Bax', 'BAX',
                            'SMAC', 'Smac', 'Caspase', 'caspase', 'DIABLO',
                            'Mcl', 'MCL', 'IAP', 'Iap'),
           'epigenetics' = c('chromatin', 'Chromatin', 'histone', 'Histone',
                             'HDAC', 'DNMT', 'methyltransf', 'acetylase', 'BRD',
                             'Bromo', 'bromo'),
           'protein' = c('Protein\\s{1}stability',
                         'protein\\s{1}stability', 'HSP', 'proteasom',
                         'protein\\s{1}translation', 'ubiquitin',
                         'eEF2K', 'farnesyl'),
           'hormones' = c('Hormone', 'hormone', 'aromatase', 'progesterone',
                          'estrogen', 'androgen', 'ESR', 'AR'),
           'metabolism' = c('metabolism', 'Metabolism', 'HMG-CoA', 'anhydrase',
                            'fatty', 'ATP\\s{1}synthase', 'hydroxyphenylpyruvate',
                            'Dihydrorotenone', 'CYP3A4')) %>%
      map(paste, collapse = '|')

    ## classification

    x$moa_short <- NA

    x <- x %>%
      mutate(targets_tmp = map_chr(targets, paste, collapse = ', '))

    for(i in seq_along(class_lst)) {

      x <- x %>%
        mutate(moa_short = ifelse(stri_detect(variable, regex = class_lst[[i]]),
                                  names(class_lst)[[i]], moa_short)) %>%
        mutate(moa_short = ifelse(stri_detect(drug_name, regex = class_lst[[i]]),
                                  names(class_lst)[[i]], moa_short)) %>%
        mutate(moa_short = ifelse(stri_detect(moa, regex = class_lst[[i]]),
                                  names(class_lst)[[i]], moa_short)) %>%
        mutate(moa_short = ifelse(stri_detect(targets_tmp, regex = class_lst[[i]]),
                                  names(class_lst)[[i]], moa_short))

    }

    x <- x %>%
      mutate(moa_short = ifelse(drug_name == 'Vincristine',
                                'cytoskeleton', moa_short),
             moa_short = ifelse(is.na(moa_short), 'other', moa_short),
             moa_short = factor(moa_short),
             moa_short = fct_relevel(moa_short,
                                     'other',
                                     after = Inf))

    return(x)

  }

  make_target_dict <- function(x) {

    ## generates a list of drug variable vectors: each element represents
    ## one molecular target

    ## target replacement: DNA and cytoskeleton --------

    x$targets <-
      map2(x$moa_short,
           x$targets,
           function(x, y) {

             if(stri_detect(x, fixed = 'cytoskel')) {

               return(c(y, 'cytoskeleton'))

             }

             return(y)

           })

    x$targets <-
      map2(x$moa_short,
           x$targets,
           function(x, y) {

             if(stri_detect(x, fixed = 'DNA')) {

               return(c(y, 'DNA'))

             }

             return(y)

           })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '(m|M)icrotub'))) {

          return(c(x, 'microtubule'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '(KDR|VEGF)'))) {

          return(c(x, 'VEGFR'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'DNA'))) {

          return(c(x, 'DNA'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'MEK(1|2)'))) {

          return(stri_replace(x,
                              fixed = 'MEK',
                              replacement = 'MAP2K'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '^c\\-'))) {

          return(stri_replace(x,
                              regex = '^c\\-',
                              replacement = ''))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'RAS'))) {

          return(c(x, 'RAS'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'mTOR'))) {

          return(stri_replace(x,
                              fixed = 'mTOR',
                              replacement = 'MTOR'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '^HDAC'))) {

          return(c(x, 'HDAC'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'Alkylating agent'))) {

          return(c(x, 'DNA alkylating agent'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '(Anti\\-meta|Antimeta)'))) {

          return(c(x, 'Antimetabolite'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'PARP'))) {

          return(c(x, 'PARP'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = 'PI3K|PI\\-3K'))) {

          return(c(x, 'PI3K'))

        }

        return(x)

      })

    x$targets <- x$targets %>%
      map(function(x) {

        if(all(is.na(x))) return(x)

        if(any(stri_detect(na.omit(x), regex = '^IAP'))) {

          return(stri_replace(x,
                              regex = '^IAP',
                              replacement = 'XIAP'))

        }

        return(x)

      })

    ## dictionary list --------

    all_targets <- na.omit(unique(unlist(x$targets)))

    dict <- all_targets %>%
      set_names(all_targets) %>%
      map(~safely(reglook)(x[, c('variable', 'targets')], .x)) %>%
      map(~.x$result)

    dict %>%
      compact %>%
      map(~.x$variable)

  }

# END -------

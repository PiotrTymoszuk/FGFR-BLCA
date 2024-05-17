# Co-occurrence of the top most frequent somatic mutations, amplifications
# and deletions.
#
# Co-occurrence is investigated by two-dimensional MDS of the Jaccard distance
# matrices. The analysis is done for the GENIE and TCGA cohort.

  insert_head()

# container -------

  lca_occur <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data -----

  insert_msg('Analysis globals')

  ## variables and lexicon

  lca_occur$variables <- lca_globals$variables

  lca_occur$lexicon <- lca_globals$lexicon

  lca_occur$data <- lca_globals$data %>%
    map(column_to_rownames, 'sample_id')

# MDS, frequencies, and distances -------

  insert_msg('MDS, frequences of alterations, and Jaccard distances')

  lca_occur$mds <- lca_occur$data %>%
    map(t) %>%
    map(mds_genetics,
        distance_method = 'jaccard')

# Scatter plots of MDS dimensions, frequency coded by point size -------

  insert_msg('Plots of MDS dimensions')

  lca_occur$mds_plots <-
    list(mds_output = lca_occur$mds,
         plot_title = globals$cohort_labs[names(lca_occur$mds)]) %>%
    pmap(plot_mds,
         size_limits = c(0, 75),
         min.segment.length = 0.25,
         max.overlaps = 20,
         txt_size = 2.3,
         force_pull = 0.5,
         force = 1)

# Stats for overlap of alterations ------

  insert_msg('Stats for overlap of alterations')

  for(i in names(lca_occur$data)) {

    lca_occur$stats[[i]] <- lca_occur$variables %>%
      map(~count_binary(lca_occur$data[[i]],
                        split_fct = .x)) %>%
      set_names(lca_occur$variables)

  }

# Pairwise Jaccard tests --------

  insert_msg('Pairwise Jaccard tests')

  lca_occur$test <- lca_occur$data %>%
    map(t) %>%
    map(post_hoc_jaccard,
        B = 10000)

  ## significant effects. common significant effects
  ## are shared by both cohorts

  lca_occur$significant <- lca_occur$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$pair_id)

  lca_occur$common_significant <- lca_occur$significant %>%
    reduce(intersect)

# Plots ---------

  insert_msg('Plots for single associations')

  ## the plots are generated for the significant variable pairs
  ## and selected pairs of interest

  lca_occur$plot_pairs <- lca_occur$common_significant %>%
    c('FGFR3_mutation|CCND1_amplification',
      'FGFR3_mutation|RB1_mutation') %>%
    stri_split_fixed(pattern = '|', simplify = FALSE) %>%
    set_names(c(lca_occur$common_significant,
                'FGFR3_mutation|CCND1_amplification',
                'FGFR3_mutation|RB1_mutation'))

  ## plots

  for(i in names(lca_occur$data)) {

    lca_occur$hm_plots[[i]] <- lca_occur$plot_pairs %>%
      map(~plot_binary_overlap(data = lca_occur$data[[i]],
                               variable1 = .x[[1]],
                               variable2 = .x[[2]],
                               jaccard_test = lca_occur$test[[i]],
                               title_suffix = globals$cohort_labs[i]))

  }

# caching the analysis results -------

  insert_msg('Caching the analysis results')

  lca_occur$data <- NULL

  lca_occur <- compact(lca_occur)

  save(lca_occur, file = './cache/lca_occur.RData')

# END -----

  plan('sequential')

  rm(i)

  insert_tail()

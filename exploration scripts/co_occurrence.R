# Co-occurrence of somatic mutations, amplifications and deletions of
# FGF/FGFR genes.
#
# Co-occurrence is investigated by two-dimensional MDS of the Jaccard distance
# matrices. The analysis is done for the GENIE and TCGA cohort.

  insert_head()

# container ------

  expl_occur <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data ------

  insert_msg('Analysis data')

  expl_occur$data <- list(genie = genie,
                          msk = msk,
                          tcga = tcga) %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    map(map,
        select,
        sample_id,
        any_of(reduce(globals[c("receptors", "ligands")], union))) %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map(~map2(., names(.),
              ~set_names(.x, paste(names(.x), .y, sep = '_')))) %>%
    map(map, rownames_to_column, 'sample_id') %>%
    map(reduce, inner_join, by = 'sample_id')

  ## common variables and variables coding for receptor mutations

  expl_occur$variables <- expl_occur$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  expl_occur$rec_variables <- paste0(globals$receptors, '_mutation')

  expl_occur$rec_variables <-
    expl_occur$rec_variables[expl_occur$rec_variables %in% expl_occur$variables]

  expl_occur$data <- expl_occur$data %>%
    map(select, sample_id, expl_occur$variables) %>%
    map(column_to_rownames, 'sample_id')

# MDS, frequencies, and distances -------

  insert_msg('MDS, frequences of alterations, and Jaccard distances')

  expl_occur$mds <- expl_occur$data %>%
    map(t) %>%
    map(mds_genetics,
        distance_method = 'jaccard')

# Scatter plots of MDS dimensions, frequency coded by point size -------

  insert_msg('Plots of MDS dimensions')

  expl_occur$mds_plots <-
    list(mds_output = expl_occur$mds,
         plot_title = globals$cohort_labs[names(expl_occur$mds)]) %>%
    pmap(plot_mds,
         size_limits = c(0, 25),
         min.segment.length = 0.25,
         max.overlaps = 25,
         txt_size = 2.3,
         force_pull = 0.25,
         force = 1)

# Stats for overlap of alterations ------

  insert_msg('Stats for overlap of alterations')

  for(i in names(expl_occur$data)) {

    expl_occur$stats[[i]] <- expl_occur$variables %>%
      map(~bin2stats(expl_occur$data[[i]],
                     split_factor = .x)) %>%
      set_names(expl_occur$variables )

  }

# Pairwise Jaccard tests --------

  insert_msg('Pairwise Jaccard tests')

  expl_occur$test <- expl_occur$data %>%
    map(t) %>%
    map(post_hoc_jaccard,
        B = 10000)

  ## significant effects. common significant effects
  ## are shared by both cohorts

  expl_occur$significant <- expl_occur$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$pair_id)

  expl_occur$common_significant <- expl_occur$significant %>%
    reduce(intersect)

# Plots ---------

  insert_msg('Plots for single associations')

  ## the plots are generated for the significant variable pairs
  ## and the FGFR receptors

  expl_occur$plot_variables <- expl_occur$significant %>%
    map(stri_split_fixed, pattern = '|', simplify = FALSE) %>%
    unlist %>%
    unname %>%
    unique %>%
    c(expl_occur$rec_variables)

  expl_occur$plot_pairs <- expl_occur$plot_variables %>%
    combn(m = 2, simplify = FALSE)

  expl_occur$plot_pairs <- expl_occur$plot_pairs %>%
    set_names(map_chr(expl_occur$plot_pairs, paste, collapse = '|'))

  ## plots

  for(i in names(expl_occur$data)) {

    expl_occur$hm_plots[[i]] <- expl_occur$plot_pairs %>%
      map(~plot_binary_overlap(data = expl_occur$data[[i]],
                               variable1 = .x[[1]],
                               variable2 = .x[[2]],
                               jaccard_test = expl_occur$test[[i]],
                               title_suffix = globals$cohort_labs[i]))

  }

# caching the analysis results -------

  insert_msg('Caching the analysis results')

  expl_occur$data <- NULL

  expl_occur <- compact(expl_occur)

  save(expl_occur, file = './cache/expl_occur.RData')

# END -----

  plan('sequential')

  rm(i)

  insert_tail()

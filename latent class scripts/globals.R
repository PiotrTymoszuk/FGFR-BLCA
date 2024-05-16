# Analysis globals:
#
# ready to use data frames with 0/1 coded top most frequent genetic features
# present in at least 5% of the samples
#
# feature lexicons
#
# Definition of the training and test subset of the GENIE cohort.

  insert_head()

# container -------

  lca_globals <- list()

# analysis data -------

  insert_msg('Data frames')

  ## variables: the most frequent genetic features
  ## and MTAP deletion as a hallmark of more general 9p21 chromosome
  ## region deletion

  lca_globals$variables <- expl_genet$top_alterations

  lca_globals$variables$deletion <-
    c(lca_globals$variables$deletion, 'MTAP')

  ## data

  lca_globals$data <- list(genie = genie, msk = msk, tcga = tcga) %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')])

  for(i in names(lca_globals$data)) {

    lca_globals$data[[i]] <-
      map2(lca_globals$data[[i]],
           lca_globals$variables,
           ~select(.x, sample_id, all_of(.y))) %>%
      map2(., names(.),
           ~set_names(.x,
                      c('sample_id', paste(names(.x)[-1], .y, sep = '_')))) %>%
      reduce(inner_join, by = 'sample_id') %>%
      filter(complete.cases(.))

  }

  lca_globals$variables <- lca_globals$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  ## data sets in a factor-coded form used in modeling
  ## they contain the hallmark genetic alterations identified by co-occurrence
  ## analysis: RB1 mutations, TP53 mutations, FGFR3 mutations,
  ## 9p21del, and 11q13amp

  lca_globals$factor_data <- lca_globals$data

  for(i in names(lca_globals$factor_data)) {

    lca_globals$factor_data[[i]][['amp11q13']] <-
      lca_globals$factor_data[[i]][, c("FGF3_amplification",
                                       "FGF4_amplification",
                                       "FGF19_amplification",
                                       "CCND1_amplification")] %>%
      reduce(`+`)

    lca_globals$factor_data[[i]][['del9p21']] <-
      lca_globals$factor_data[[i]][, c("CDKN2A_deletion",
                                       "CDKN2B_deletion",
                                       "MTAP_deletion")] %>%
      reduce(`+`)


  }

  lca_globals$modeling_variables <-
    lca_globals$variables[!lca_globals$variables %in% c(c("FGF3_amplification",
                                                         "FGF4_amplification",
                                                         "FGF19_amplification",
                                                         "CCND1_amplification",
                                                         "CDKN2A_deletion",
                                                         "CDKN2B_deletion",
                                                         "MTAP_deletion"))] %>%
    c('amp11q13', 'del9p21')

  lca_globals$factor_data <- lca_globals$factor_data %>%
    map(select, sample_id, all_of(lca_globals$modeling_variables))

  for(i in names(lca_globals$factor_data)) {

    lca_globals$factor_data[[i]][lca_globals$modeling_variables] <-
      lca_globals$factor_data[[i]][lca_globals$modeling_variables] %>%
      map_dfc(~ifelse(.x > 0, 'yes', 'no')) %>%
      map_dfc(factor, c('no', 'yes'))

  }

  ## modeling formula

  lca_globals$formula <-
    paste0('cbind(',
           paste(lca_globals$modeling_variables, collapse = ', '),
           ') ~ 1') %>%
    as.formula

# variable lexicon ------

  insert_msg('Variable lexicon')

  lca_globals$lexicon <-
    tibble(variable = lca_globals$variables) %>%
    mutate(gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           label = paste(gene_symbol, alteration),
           label_html = paste(html_italic(gene_symbol), alteration))

# Training and test subset of the GENIE data set --------

  insert_msg('Test and training subset of the GENIE cohort')

  set.seed(12345)

  lca_globals$genie_train_indexes <- lca_globals$data$genie$TP53_mutation %>%
    createDataPartition %>%
    unlist %>%
    unname

# END ------

  insert_tail()

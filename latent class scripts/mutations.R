# Enrichment and differential alteration analysis for somatic mutations in
# Genetic Clusters.
#
# The analysis is done in the GENIE, MSK, and TCGA cohort by weighted
# permutation testing with perTest/perich. The background, i.e. sample weights
# for the permutations are obtained by counting all mutated genes per sample.
#
# Significant enrichment in a cluster is assumed for enrichment score >= 1.2
# and FDR-corrected p value < 0.05.

  insert_head()

# container -------

  lca_mut <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## analysis data and background vectors with counts of all mutated genes

  lca_mut$data <- list(genie = genie,
                       msk = msk,
                       tcga = tcga) %>%
    map(~.x$mutation) %>%
    map(column_to_rownames, 'sample_id')

  lca_mut$backgrounds <- lca_mut$data

  ## cluster assignment vectors

  lca_mut$f_vectors <- lca_pred$assignment %>%
    map(~set_names(.x$clust_id, .x$sample_id))

  lca_mut$f_vectors <-
    map2(lca_mut$f_vectors,
         lca_mut$data,
         ~.x[rownames(.y)])

# Frequencies and mutations of interest present in at least 1% samples -------

  insert_msg('Frequencies of mutations, mutations of interest')

  ## overall mutation frequencies,
  ## mutations of interest shared by all data sets

  lca_mut$frequencies <- lca_mut$data %>%
    map(count_binary)

  lca_mut$variables <- lca_mut$frequencies %>%
    map(filter, percent >= 1) %>%
    map(~.x$variable) %>%
    reduce(intersect)

  lca_mut$data <- lca_mut$data %>%
    map(~.x[lca_mut$variables])

  ## frequencies of mutations in genetic clusters

  lca_mut$stats <-
    map2(lca_mut$data,
         lca_mut$f_vectors,
         ~mutate(.x, clust_id = .y)) %>%
    map(count_binary,
        variables = lca_mut$variables,
        split_fct = 'clust_id')

# Permutation tests -------

  insert_msg('Permutation tests')

  plan('multisession')

  lca_mut$test <-
    list(x = lca_mut$data,
         f = lca_mut$f_vectors,
         background = lca_mut$backgrounds,
         n_iter = c(1000, 10000, 100000)) %>%
    pmap(perTest,
         alternative = 'greater',
         as_data_frame = TRUE,
         compress = TRUE,
         adj_method = 'BH',
         .parallel = TRUE)

  plan('sequential')

# Formatting the testing results, significant effects --------

  insert_msg('Formatting the testing results, significant effects')

  lca_mut$test <- lca_mut$test %>%
    map(re_adjust, 'global_p_value') %>%
    map(mutate, global_significance = significance) %>%
    map(re_adjust) %>%
    map(mutate,
        clust_id = factor(strata, levels(lca_mut$f_vectors[[1]])),
        global_eff_size = paste('V =', signif(cramer_v, 2)),
        regulation = ifelse(p_adjusted < 0.05 & ES >= 1.2,
                            'enriched', 'ns'),
        regulation = factor(regulation, c('enriched', 'ns'))) %>%
    map(as_tibble)

  ## significant alterations: overall and in single clusters

  lca_mut$significant <- lca_mut$test %>%
    map(filter, regulation == 'enriched') %>%
    map(~.x$variable) %>%
    map(unique)

  lca_mut$enriched <- lca_mut$test %>%
    map(filter, regulation == 'enriched') %>%
    map(blast, strata) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## common significant alterations shared by at least two cohorts

  lca_mut$common_significant <- lca_mut$significant %>%
    shared_features(m = 2) %>%
    as.character

  lca_mut$common_enriched <- lca_mut$enriched %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Caching -------

  insert_msg('Caching')

  lca_mut$data <- NULL
  lca_mut$backgrounds <- NULL
  lca_mut$f_vectors <- NULL

  lca_mut <- compact(lca_mut)

  save(lca_mut, file = './cache/lca_mut.RData')

# END -------

  insert_tail()

# Enrichment and differential alteration analysis for deletions in
# Genetic Clusters.
#
# The analysis is done in the GENIE, MSK, and TCGA cohort by weighted
# permutation testing with perTest/perich. The background, i.e. sample weights
# for the permutations are obtained by counting all deleted genes per sample.
#
# Significant enrichment in a cluster is assumed for enrichment score >= 1.2
# and FDR-corrected p value < 0.05.

  insert_head()

# container -------

  lca_del <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## analysis data and background vectors with counts of all mutated genes

  lca_del$data <- list(genie = genie,
                       msk = msk,
                       tcga = tcga) %>%
    map(~.x$deletion) %>%
    map(column_to_rownames, 'sample_id')

  lca_del$backgrounds <- lca_del$data

  ## cluster assignment vectors

  lca_del$f_vectors <- lca_pred$assignment %>%
    map(~set_names(.x$clust_id, .x$sample_id))

  lca_del$f_vectors <-
    map2(lca_del$f_vectors,
         lca_del$data,
         ~.x[rownames(.y)])

# Frequencies and deletions of interest present in at least 1% samples -------

  insert_msg('Frequencies of mutations, mutations of interest')

  ## overall deletion frequencies,
  ## deletions of interest shared by all data sets

  lca_del$frequencies <- lca_del$data %>%
    map(count_binary)

  lca_del$variables <- lca_del$frequencies %>%
    map(filter, percent >= 1) %>%
    map(~.x$variable) %>%
    reduce(intersect)

  lca_del$data <- lca_del$data %>%
    map(~.x[lca_del$variables])

  ## frequencies of deletions in genetic clusters

  lca_del$stats <-
    map2(lca_del$data,
         lca_del$f_vectors,
         ~mutate(.x, clust_id = .y)) %>%
    map(count_binary,
        variables = lca_del$variables,
        split_fct = 'clust_id')

# Permutation tests -------

  insert_msg('Permutation tests')

  plan('multisession')

  lca_del$test <-
    list(x = lca_del$data,
         f = lca_del$f_vectors,
         background = lca_del$backgrounds,
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

  lca_del$test <- lca_del$test %>%
    map(re_adjust, 'global_p_value') %>%
    map(mutate, global_significance = significance) %>%
    map(re_adjust) %>%
    map(mutate,
        clust_id = factor(strata, levels(lca_del$f_vectors[[1]])),
        global_eff_size = paste('V =', signif(cramer_v, 2)),
        regulation = ifelse(p_adjusted < 0.05 & ES >= 1.2,
                            'enriched', 'ns'),
        regulation = factor(regulation, c('enriched', 'ns'))) %>%
    map(as_tibble)

  ## significant alterations: overall and in single clusters

  lca_del$significant <- lca_del$test %>%
    map(filter, regulation == 'enriched') %>%
    map(~.x$variable) %>%
    map(unique)

  lca_del$enriched <- lca_del$test %>%
    map(filter, regulation == 'enriched') %>%
    map(blast, strata) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## common significant alterations shared by at least two cohorts

  lca_del$common_significant <- lca_del$significant %>%
    shared_features(m = 2) %>%
    as.character

  lca_del$common_enriched <- lca_del$enriched %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Caching -------

  insert_msg('Caching')

  lca_del$data <- NULL
  lca_del$backgrounds <- NULL
  lca_del$f_vectors <- NULL

  lca_del <- compact(lca_del)

  save(lca_del, file = './cache/lca_del.RData')

# END -------

  insert_tail()

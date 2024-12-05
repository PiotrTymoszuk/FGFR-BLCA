# Enrichment and differential alteration analysis for gene amplifications in
# Genetic Clusters.
#
# The analysis is done in the GENIE, MSK, and TCGA cohort by weighted
# permutation testing with perTest/perich. The background, i.e. sample weights
# for the permutations are obtained by counting all amplified genes per sample.
#
# Significant enrichment in a cluster is assumed for enrichment score >= 1.2
# and FDR-corrected p value < 0.05.

  insert_head()

# container -------

  lca_amp <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## analysis data and background matrix with counts of all amplified genes

  lca_amp$data <- list(genie = genie,
                       msk = msk,
                       tcga = tcga) %>%
    map(~.x$amplification) %>%
    map(column_to_rownames, 'sample_id')

  lca_amp$backgrounds <- lca_amp$data

  ## cluster assignment vectors

  lca_amp$f_vectors <- lca_pred$assignment %>%
    map(~set_names(.x$clust_id, .x$sample_id))

  lca_amp$f_vectors <-
    map2(lca_amp$f_vectors,
         lca_amp$data,
         ~.x[rownames(.y)])

# Frequencies and amplifications of interest present in at least 1% samples -------

  insert_msg('Frequencies of amplifications, amplifications of interest')

  ## overall amplification frequencies,
  ## mutations of interest shared by all data sets

  lca_amp$frequencies <- lca_amp$data %>%
    map(count_binary)

  lca_amp$variables <- lca_amp$frequencies %>%
    map(filter, percent >= 1) %>%
    map(~.x$variable) %>%
    reduce(intersect)

  lca_amp$data <- lca_amp$data %>%
    map(~.x[lca_amp$variables])

  ## frequencies of mutations in genetic clusters

  lca_amp$stats <-
    map2(lca_amp$data,
         lca_amp$f_vectors,
         ~mutate(.x, clust_id = .y)) %>%
    map(count_binary,
        variables = lca_amp$variables,
        split_fct = 'clust_id')

# Permutation tests -------

  insert_msg('Permutation tests')

  plan('multisession')

  lca_amp$test <-
    list(x = lca_amp$data,
         f = lca_amp$f_vectors,
         background = lca_amp$backgrounds,
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

  lca_amp$test <- lca_amp$test %>%
    map(re_adjust, 'global_p_value') %>%
    map(mutate, global_significance = significance) %>%
    map(re_adjust) %>%
    map(mutate,
        clust_id = factor(strata, levels(lca_amp$f_vectors[[1]])),
        global_eff_size = paste('V =', signif(cramer_v, 2)),
        regulation = ifelse(p_adjusted < 0.05 & ES >= 1.2,
                            'enriched', 'ns'),
        regulation = factor(regulation, c('enriched', 'ns'))) %>%
    map(as_tibble)

  ## significant alterations: overall and in single clusters

  lca_amp$significant <- lca_amp$test %>%
    map(filter, regulation == 'enriched') %>%
    map(~.x$variable) %>%
    map(unique)

  lca_amp$enriched <- lca_amp$test %>%
    map(filter, regulation == 'enriched') %>%
    map(blast, strata) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## common significant alterations shared by at least two cohorts

  lca_amp$common_significant <- lca_amp$significant %>%
    shared_features(m = 2) %>%
    as.character

  lca_amp$common_enriched <- lca_amp$enriched %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Caching -------

  insert_msg('Caching')

  lca_amp$data <- NULL
  lca_amp$backgrounds <- NULL
  lca_amp$f_vectors <- NULL

  lca_amp <- compact(lca_amp)

  save(lca_amp, file = './cache/lca_amp.RData')

# END -------

  insert_tail()

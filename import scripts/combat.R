# ComBat adjustment of gene expression in the TCGA and Imvigor cohorts

  insert_head()

# container ------

  combat <- list()

# raw expression data sets ------

  insert_msg('Raw expression data sets')

  ## log2-transformed, restricted to the common genes

  combat$raw_expression <- list(tcga = tcga, imvigor = imvigor) %>%
    map(~.x$expression)

  combat$genes <- combat$raw_expression %>%
    map(names) %>%
    map(~.x[!.x %in% c('sample_id', 'patient_id', 'tissue')]) %>%
    reduce(intersect)

  combat$raw_expression <- combat$raw_expression %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, combat$genes) %>%
    map(t)

# Batch effect adjustment --------

  insert_msg('Batch effect adjustment')

  combat$expression <- pre_process(train = combat$raw_expression$tcga,
                                   test = combat$raw_expression$imvigor)

  combat$expression <- combat$expression[c("train", "test")] %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble) %>%
    set_names(names(combat$raw_expression))

# Cachinig the results -------

  insert_msg('Caching')

  combat <- combat[c("genes", "expression")]

  save(combat, file = './data/combat.RData')

# END -----

  insert_tail()

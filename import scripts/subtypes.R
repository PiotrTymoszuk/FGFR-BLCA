# Assignment of the TCGA BLCA and IMvigor samples to the consensus molecular
# subtypes based on the ComBat-adjusted whole-genome expression.

  insert_head()

# container -------

  subtypes <- list()

# analysis data -----

  insert_msg('Analysis data')

  subtypes$data <- combat$expression %>%
    map(column_to_rownames, 'sample_id') %>%
    map(t)

# Assignment to the subtypes -------

  insert_msg('Assignment to the subtypes')

  subtypes$assignment <- subtypes$data %>%
    map(getConsensusClass,
        gene_id = 'hgnc_symbol') %>%
    map(rownames_to_column, 'sample_id') %>%
    map(mutate,
        consensusClass = factor(consensusClass ,
                                c('LumP', 'LumU', 'LumNS',
                                  'Stroma-rich', 'Ba/Sq',
                                  'NE-like'))) %>%
    map(as_tibble)

# END -------

  subtypes$data <- NULL

  subtypes <- compact(subtypes)

  insert_tail()

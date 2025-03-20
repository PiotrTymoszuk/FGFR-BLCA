# Import and clearing of the baseline expression data for the PRISM project
# cell lines. The data are available as log2(x + 1) transformed TPM.
# This is basically the same data set as in case of CTRP2.

  insert_head()

# container -----

  prism_exp <- list()

# reading the expression data -------

  insert_msg('Reading the expression data')

  prism_exp$expression <-
    read_csv('./data/PRISM/Expression_Public_23Q4.csv')

  names(prism_exp$expression)[1] <- 'sample_id'

# Annotation -------

  insert_msg('Annotation')

  prism_exp$annotation <-
    tibble(probe_id = names(prism_exp$expression)[-1]) %>%
    mutate(entrez_id = mapIds(org.Hs.eg.db,
                              keys = probe_id,
                              keytype = 'SYMBOL',
                              column = 'ENTREZID'),
           entrez_id = ifelse(!is.na(entrez_id),
                              entrez_id,
                              mapIds(org.Hs.eg.db,
                                     keys = probe_id,
                                     keytype = 'ALIAS',
                                     column = 'ENTREZID'))) %>%
    filter(complete.cases(.)) %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    filter(complete.cases(.))

  ## re-naming the expression data set: official symbols

  prism_exp$expression <- prism_exp$expression %>%
    select(sample_id, all_of(prism_exp$annotation$probe_id)) %>%
    set_names(c('sample_id', prism_exp$annotation$gene_symbol))

# Removal of duplicated features --------

  insert_msg('Unique genes only')

  prism_exp$unique_genes <- names(prism_exp$expression)[-1] %>%
    unique

  prism_exp$expression <- prism_exp$expression %>%
    select(sample_id, all_of(prism_exp$unique_genes))

# Caching the results -------

  insert_msg('Caching')

  prism_exp$unique_genes <- NULL

  prism_exp <- compact(prism_exp)

  save(prism_exp, file = './data/PRISM/prism_exp.RData')

# END ------

  insert_tail()

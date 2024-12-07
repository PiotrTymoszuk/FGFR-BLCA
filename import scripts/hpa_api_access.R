# Access to XML files with patient and staining information in the Human Protein
# Atlas repository

  insert_head()

# container -------

  hpa_access <- list()

# genes of interest -------

  insert_msg('Genes of interest')

  hpa_access$lexicon <-
    tibble(gene_symbol = reduce(globals[c("receptors",
                                          "ligands",
                                          "binding_proteins")],
                                union)) %>%
    mutate(ensembl_id = mapIds(org.Hs.eg.db,
                               keys = gene_symbol,
                               keytype = 'SYMBOL',
                               column = 'ENSEMBL')) %>%
    filter(complete.cases(.))

  ## specifying the URLs and file names

  hpa_access$lexicon <- hpa_access$lexicon %>%
    mutate(url = paste0('https://www.proteinatlas.org/', ensembl_id, '.xml'),
           file = paste0('./data/HPA/', ensembl_id, '.xml'))

# Downloads ---------

  insert_msg('Downloads')

  hpa_access$lexicon[, c("url", "file")] %>%
    pwalk(download_hpa)

# END ------

  insert_tail()

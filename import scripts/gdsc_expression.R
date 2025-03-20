# GDSC baseline cell line gene expression data: RMA normalized 
# and log2-transformed

  insert_head()
  
# container -------
  
  gdsc_exp <- list()
  
# Reading the expression data -------
  
  insert_msg('Reading the gene expression data')
  
  gdsc_exp$expression <- 
    read_tsv('./data/GDSC/Cell_line_RMA_proc_basalExp.txt') %>% 
    mutate(probe_id = paste0('probe_', 1:nrow(.)))
  
# Probe/gene annotation ------
  
  insert_msg('Probe/gene annotation')
  
  gdsc_exp$annotation <- gdsc_exp$expression %>% 
    transmute(probe_id = probe_id, 
              gene_symbol = GENE_SYMBOLS) %>% 
    filter(complete.cases(.))
  
  ## Entrez ID: some of 'gene symbols' are in fact aliases
  
  gdsc_exp$annotation <- gdsc_exp$annotation %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              keytype = 'SYMBOL', 
                              column = 'ENTREZID')) %>% 
    mutate(entrez_id = ifelse(!is.na(entrez_id), entrez_id, 
                              mapIds(org.Hs.eg.db, 
                                     keys = gene_symbol, 
                                     keytype = 'ALIAS', 
                                     column = 'ENTREZID'))) %>% 
    filter(complete.cases(.))
  
  ## official gene symbols
  
  gdsc_exp$annotation <- gdsc_exp$annotation %>% 
    mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                                keys = entrez_id, 
                                keytype = 'ENTREZID', 
                                column = 'SYMBOL')) %>% 
    filter(complete.cases(.))
  
# Clearing of the expression data set -----
  
  insert_msg('Clearing of the expression data set')
  
  gdsc_exp$expression<- gdsc_exp$expression %>% 
    select(-GENE_SYMBOLS, -GENE_title) %>% 
    column_to_rownames('probe_id') %>% 
    integrate_expression(annotation = gdsc_exp$annotation)
  
  ## sample ID: official COSMIC identifier
  
  gdsc_exp$expression <- gdsc_exp$expression %>% 
    mutate(sample_id = stri_replace(sample_id, 
                                    regex = '^DATA.',
                                    replacement = 'cosmic_'))

# Caching the cleared data sets -------
  
  insert_msg('Caching the cleared data sets')
  
  save(gdsc_exp, file = './data/GDSC/gdsc_exp.RData')
  
# END ------
  
  insert_tail()
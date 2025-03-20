# Import and clearing of the baseline expression data for the CTRP2 project 
# cell lines. The data are available as log2(x + 1) transformed TPM.

  insert_head()
  
# container -----
  
  ctrp2_exp <- list()
  
# reading the expression data -------
  
  insert_msg('Reading the expression data')
  
  ctrp2_exp$expression <- 
    read_csv('./data/CTRP2/Expression_Public_23Q4.csv')
  
  names(ctrp2_exp$expression)[1] <- 'sample_id'
  
# Annotation -------
  
  insert_msg('Annotation')
  
  ctrp2_exp$annotation <- 
    tibble(probe_id = names(ctrp2_exp$expression)[-1]) %>% 
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
  
  ctrp2_exp$expression <- ctrp2_exp$expression %>% 
    select(sample_id, all_of(ctrp2_exp$annotation$probe_id)) %>% 
    set_names(c('sample_id', ctrp2_exp$annotation$gene_symbol))
  
# Removal of duplicated features --------
  
  insert_msg('Unique genes only')
  
  ctrp2_exp$unique_genes <- names(ctrp2_exp$expression)[-1] %>% 
    unique
  
  ctrp2_exp$expression <- ctrp2_exp$expression %>% 
    select(sample_id, all_of(ctrp2_exp$unique_genes))
  
# Caching the results -------
  
  insert_msg('Caching')
  
  ctrp2_exp$unique_genes <- NULL
  
  ctrp2_exp <- compact(ctrp2_exp)

  save(ctrp2_exp, file = './data/CTRP2/ctrp2_exp.RData')
  
# END ------
  
  insert_tail()
# Computation of ssGSEA scores of gene signatures associated with FGFR
# signaling.
# The signatures are extracted from MSig database version 7.5.1

  insert_head()

# container ------

  sig <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# reading the signature database, signatures of interest ------

  insert_msg('Reading the signature database')

  sig$db <- load_dbsig('./data/msig/msigdb.v7.5.1.symbols.gmt')

  sig$lexicon <- sig$db %>%
    filter(stri_detect(sign_name, regex = '^(REACTOME|PID|KEGG|WP)'),
           stri_detect(sign_name, regex = 'FGF|FGFR'))

# computing the ssGSEA scores for the TCGA, Imvigor, and BCAN cohorts ------

  insert_msg('Computing the ssGSEA scores')

  ## whole transcriptome expression

  sig$expression <- list(tcga = tcga, imvigor = imvigor, bcan = bcan) %>%
    map(~.x$expression) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(~.x[!names(.x) %in% c('patient_id', 'tissue_type')])

  sig$scores <- sig$expression %>%
    future_map(calculate.dbsig,
               x = sig$lexicon,
               .options = furrr_options(seed = TRUE))

# Edits of the lexicon -------

  insert_msg('Edits of the lexicon')

  sig$sign_names <- sig$scores %>%
    map(names) %>%
    reduce(intersect)

  sig$lexicon <- sig$lexicon %>%
    filter(sign_name %in% sig$sign_names) %>%
    mutate(variable = sign_name,
           label = stri_replace_all(variable, fixed = '_', replacement = ' '))

# sample IDs for the score data frames -------

  insert_msg('Setting the sample IDs')

  sig$scores <-
    map2(sig$scores,
         sig$expression,
         ~mutate(.x, sample_id = rownames(.y))) %>%
    map(relocate, sample_id)

# Caching ------

  insert_msg('Caching')

  sig <- sig[c("lexicon", "scores")]

  save(sig, file = './data/sig.RData')

# END ------

  plan('sequential')

  insert_tail()

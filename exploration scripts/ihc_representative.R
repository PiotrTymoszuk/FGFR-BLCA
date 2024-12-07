# Representative IHC staining images.
# We're showing for each of them a sample with the
# lowest expression, median, and 75% percentile.

  insert_head()

# container -------

  expl_ihc <- list()

# analysis globals ---------

  insert_msg('Genes of interest and their image paths')

  expl_ihc$genes <- hpa$data$gene_symbol %>%
    unique

  expl_ihc$lexicon <- hpa$data %>%
    filter(gene_symbol %in% expl_ihc$genes) %>%
    mutate(gene_symbol = factor(gene_symbol, expl_ihc$genes),
           caption = paste0(sex, ', ', age, 'y\n',
                            grade, ' grade'),
           caption_short = paste0(sex, ', ', age, 'y')) %>%
    select(patient_id, gene_symbol, ihc_score,
           image_path, caption, caption_short)

  ## selection of the samples: work on that!!!

  set.seed(123)

  expl_ihc$lexicon <- expl_ihc$lexicon %>%
    blast(gene_symbol)

  expl_ihc$representative_idx <- expl_ihc$lexicon %>%
    map(~.x$ihc_score) %>%
    map(choose_representative,
        probs = c(0, 0.6, 1))

  expl_ihc$lexicon <-
    map2(expl_ihc$lexicon,
         expl_ihc$representative_idx,
         ~.x[.y, ])

# staining panels saved as calls to spare memory -------

  insert_msg('Staining panels')

  expl_ihc$panels <- expl_ihc$lexicon %>%
    map(~call2(.fn = 'stitch_hpa_images',
               lexicon = .x,
               y_pos = c(0.12, 0.88)))

# END -------

  insert_tail()

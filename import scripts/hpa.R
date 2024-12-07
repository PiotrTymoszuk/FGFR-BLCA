# Processing of XML entries with patient's and sample staining information for
# the FGFR, FGF, and FGFBP proteins. Resorting to locally saved XML files.

  insert_head()

# container --------

  hpa <- list()

# available XML files --------

  insert_msg('Available XML files and gene lexicon')

  ## genes and the expected file paths

  hpa$annotation <-
    tibble(gene_symbol = reduce(globals[c("receptors",
                                          "ligands",
                                          "binding_proteins")],
                                union)) %>%
    mutate(ensembl_id = mapIds(org.Hs.eg.db,
                               keys = gene_symbol,
                               keytype = 'SYMBOL',
                               column = 'ENSEMBL'),
           entrez_id =  mapIds(org.Hs.eg.db,
                               keys = ensembl_id,
                               keytype = 'ENSEMBL',
                               column = 'ENTREZID')) %>%
    filter(complete.cases(.))

  hpa$annotation <- hpa$annotation %>%
    mutate(file = paste0('./data/HPA/', ensembl_id, '.xml'),
           exists = file.exists(file))

  if(sum(hpa$annotation$exists) < 10) {

    insert_msg('Too few files found, accessing HPA')

    source_all('./import scripts/hpa_api_access.R',
               message = TRUE, crash = TRUE)

  }

  hpa$annotation <- hpa$annotation %>%
    filter(exists)

# Extraction of the patient and staining data ------

  insert_msg('Extraction')

  hpa$raw_data <- hpa$annotation$file %>%
    map(read_xml) %>%
    set_names(hpa$annotation$gene_symbol)

  hpa$raw_data <- hpa$raw_data %>%
    future_map(safely(process_hpa_entry),
               .options = furrr_options(seed = TRUE)) %>%
    map(~.x$result) %>%
    compact

# Formatting ---------

  insert_msg('Formatting')

  hpa$data <- hpa$raw_data %>%
    map(mutate,
        sex = factor(tolower(sex), c('female', 'male')),
        age = as.numeric(age),
        staining = factor(tolower(staining),
                          c('not detected', 'low', 'medium', 'high')),
        intensity = factor(tolower(intensity),
                           c('negative', 'weak', 'moderate', 'strong')),
        quantity = factor(tolower(quantity),
                          c('none', '<25%', '75%-25%', '>75%')),
        grade = stri_extract(tolower(sample), regex = '(high|low)\\s{1}grade'),
        grade = stri_extract(grade, regex = '(high|low)'),
        grade = factor(grade, c('low', 'high')),
        image_name = stri_extract(image_url, regex = '\\d+/.*$'),
        image_name = stri_replace_all(image_name,
                                      fixed = '/',
                                      replacement = '__'),
        image_path = paste0('./data/HPA/images/', image_name))

  ## averaging the patients, IHC intensity score is a product of
  ## intensity and quantity

  hpa$data <- hpa$data %>%
    map(blast, patient_id) %>%
    map(map,
        mutate,
        staining = factor(levels(staining)[mean(as.numeric(staining))],
                          levels(staining)),
        intensity = factor(levels(intensity)[mean(as.numeric(intensity))],
                           levels(intensity)),
        quantity = factor(levels(quantity)[mean(as.numeric(quantity))],
                          levels(quantity))) %>%
    map(map_dfr, ~.x[1, ]) %>%
    compress(names_to = 'gene_symbol') %>%
    mutate(ihc_score = as.numeric(intensity) * as.numeric(quantity))

# Separate data frames for staining, intensity, quantity and IHC scores ------

  insert_msg('Analysis data set for IHC quantity variables')

  for(i in c('staining', 'intensity', 'quantity', 'ihc_score')) {

    hpa[[i]] <- hpa$data %>%
      select(patient_id, gene_symbol, all_of(i)) %>%
      pivot_wider(names_from = 'gene_symbol',
                  values_from = all_of(i))

  }

# A data frame with clinical information -------

  insert_msg('A data frame with clinical information')

  hpa$clinic <-
    hpa$data[, c("patient_id", "sex", "age", "grade", "sample")] %>%
    filter(!duplicated(patient_id))

# Downloading IHC images ------

  insert_msg('Downloading IHC images')

  hpa$data <- hpa$data %>%
    mutate(image_exists = file.exists(image_path))

  if(sum(hpa$data$image_exists) < 10) {

    insert_msg('Too few images found, accessing HPA')

    source_all('./import scripts/hpa_images.R',
               message = TRUE, crash = TRUE)

  }

# Caching -------

  insert_msg('Caching')

  save(hpa, file = './data/hpa.RData')

# END -------

  insert_tail()

 download_image(hpa$data$image_url[[1]], hpa$data$image_path[[1]])


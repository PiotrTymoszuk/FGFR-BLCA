# Import of the Genie public data set for bladder and urothelial cancer

  library(cbioportalR)

  library(tidyverse)

  #Sys.setenv(CBIOPORTAL_TOKEN = '')

  set_cbioportal_db(db = 'genie.cbioportal.org')

  query <- read_tsv('./data/GENIE/genie_public_clinical_data.tsv')

  clinical <- get_clinical_by_patient(study_id = 'genie_public',
                                      patient_id = query$`Patient ID`)

  sample <- get_clinical_by_sample(study_id = 'genie_public',
                                   sample_id = query$`Sample ID`)

  mutations <- get_mutations_by_sample(study_id = 'genie_public',
                                       sample_id = query$`Sample ID`)

  cna <- get_cna_by_sample(study_id = 'genie_public',
                           sample_id = query$`Sample ID`)

  ## saving the results

  genie_raw <- list(patient = clinical,
                    sample = sample,
                    mutations = mutations,
                    cna = cna)

  save(genie_raw, file = './data/GENIE/genie_raw.RData')

# END ------

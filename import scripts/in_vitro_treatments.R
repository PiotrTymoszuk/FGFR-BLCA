# Import of in vitro treatment data:
# urothelial cancer cells, WT and FGFR3 mutants, challenged with the pan-FGFR
# inhibitor Erdafitnib.
# The readouts are measurements of confluence of the cell culture on day four
# following treatment with 0, 1, 4, 7, and 10 µM erdafitinib.
#
# The data are read from a Graph Pad Prism

  insert_head()

# container ------

  iv <- list()

# reading the content of the .prism file -------

  insert_msg("Reading the prism file")

  ## there's only one data table in th file;
  ## the "test_1/2/3" columns are empty

  iv$raw <-
    read_prism("./data/in vitro treatments/Erdafitinib FGFR report.prism")[[1]]

  iv$raw <- iv$raw %>%
    select(-starts_with("test"))

# wrangling: long data format -------

  insert_msg("Wrangling")

  iv$data <- iv$raw %>%
    pivot_longer(cols = all_of(names(iv$raw)[-1]),
                 names_to = "condition",
                 values_to = "confluence")

  names(iv$data)[1] <- "concentration"

  ## information about the replicate, cell line name, and genotype,
  ## HTML labels for plots

  iv$data <- iv$data %>%
    mutate(cell_name = stri_split_fixed(condition,
                                        pattern = " ",
                                        simplify = TRUE)[, 1],
           genotype = stri_extract(condition, regex = "\\(.*\\)"),
           genotype = stri_extract(genotype, regex = "WT|(\\w{1}\\d{3}\\w{1})"),
           genotype = fct_relevel(factor(genotype), "WT", "S249C"),
           FGFR3_mutation = ifelse(genotype == "WT",
                                   "WT", "mutated"),
           FGFR3_mutation = factor(FGFR3_mutation, c("WT", "mutated")),
           replicate = stri_extract(condition, regex = "_\\d{1}$"),
           replicate = stri_extract(replicate, regex = "\\d{1}"),
           replicate = factor(replicate),
           html_label = paste0(cell_name, " (<em>FGFR3</em><sup>",
                               genotype, "</sup>)"))

# END ------

  insert_tail()

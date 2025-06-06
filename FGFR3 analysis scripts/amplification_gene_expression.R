# Effects of FGF3/4/19 amplification on expression of FGF and FGFR genes.
# The analysis is done for the TCGA cohort.

  insert_head()

# container -------

  fgfr_ampl <- list()

# analysis data and variables -------

  insert_msg('Analysis data')

  ## expression data

  fgfr_ampl$variables <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    reduce(union)

  fgfr_ampl$data <- tcga$expression %>%
    select(sample_id, CCND1, any_of(fgfr_ampl$variables))

  fgfr_ampl$variables <-
    fgfr_ampl$variables[fgfr_ampl$variables %in% names(fgfr_ampl$data)]

  fgfr_ampl$variables <- c("CCND1", fgfr_ampl$variables)

  ## genetics data

  fgfr_ampl$genetic <-
    tcga$amplification[c("sample_id", "FGF3", "FGF4", "FGF19")]

  fgfr_ampl$genetic$amp11q13 <-
    fgfr_ampl$genetic[c("FGF3", "FGF4", "FGF19")] %>%
    reduce(`+`)

  fgfr_ampl$genetic <- fgfr_ampl$genetic %>%
    mutate(amp11q13 = ifelse(amp11q13 > 0, 'amplified', 'WT'),
           amp11q13 = factor(amp11q13, c('WT', 'amplified'))) %>%
    select(sample_id, amp11q13)

  fgfr_ampl$data <- inner_join(fgfr_ampl$genetic,
                               fgfr_ampl$data,
                               by = 'sample_id')

# Descriptive stats --------

  insert_msg('Descriptive stats')

  fgfr_ampl$stats <- fgfr_ampl$data %>%
    fast_num_stats(variables = fgfr_ampl$variables,
                   split_fct = 'amp11q13')

# Testing for differences between WT and amplified cancers -------

  insert_msg('Testing for differences')

  fgfr_ampl$test <-
    f_wilcox_test(fgfr_ampl$data[fgfr_ampl$variables],
                  f = fgfr_ampl$data$amp11q13,
                  exact = FALSE,
                  as_data_frame = TRUE,
                  adj_method = 'BH',
                  safely = TRUE) %>%
    re_adjust %>%
    p_formatter(text = TRUE) %>%
    mutate(eff_size = paste('r =', signif(biserial_r, 2)),
           plot_cap = paste(eff_size, significance, sep = ', ')) %>%
    as_tibble

  ## significant differences

  fgfr_ampl$significant <- fgfr_ampl$test %>%
    filter(p_adjusted < 0.05) %>%
    .$variable

# Box plots --------

  insert_msg('Box plots')

  fgfr_ampl$plots <-
    list(variable = fgfr_ampl$test$variable,
         plot_title = fgfr_ampl$test$variable %>%
           html_italic %>%
           paste(globals$cohort_labs["tcga"], sep = ', '),
         plot_subtitle = fgfr_ampl$test$plot_cap) %>%
    pmap(plot_variable,
         fgfr_ampl$data,
         split_factor = 'amp11q13',
         type = 'box',
         cust_theme = globals$common_theme +
           theme(plot.title = element_markdown(),
                 axis.title.x = element_markdown()),
         x_n_labs = TRUE,
         y_lab = expression('log'[2] * ' expression'),
         x_lab = html_italic('chromosome 11q13')) %>%
    map(~.x + scale_fill_manual(values = globals$amp_status_color)) %>%
    set_names(fgfr_ampl$test$variable)

# Result table -------

  insert_msg('Result table')

  fgfr_ampl$result_tbl <-
    merge_stat_test(fgfr_ampl$stats,
                    fgfr_ampl$test) %>%
    format_res_tbl(dict = NULL,
                   remove_complete = TRUE) %>%
    set_names(c('Variable',
                levels(fgfr_ampl$data$amp11q13),
                'Significance', 'Effect size'))

# END ------

  fgfr_ampl$data <- NULL

  fgfr_ampl <- compact(fgfr_ampl)

  insert_tail()

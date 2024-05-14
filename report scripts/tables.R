# Tables for the analysis report

  insert_head()

# container ------

  rep_tabs <- list()

# characteristic of the cohorts -------

  insert_msg('Cohort characteristic')

  rep_tabs$cohorts <- expl_cohorts$result_tbl %>%
    set_names(c('Variable',
                globals$cohort_labs[c("genie", "tcga", "imvigor")],
                'Significance',
                'Effect size')) %>%
    mdtable(label = 'cohort_characteristic',
            ref_name = 'cohorts',
            caption = paste('Characteristic of the investigated cohorts.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))

# Frequency of genetic alterations -------

  insert_msg('Frequency of genetic features')

  rep_tabs$genetics <- expl_genet$result_tbl %>%
    compress(names_to = 'Alteration type') %>%
    relocate(`Alteration type`) %>%
    mdtable(label = 'genetic_features',
            ref_name = 'genetics',
            caption = paste('Frequency of the most common somatic mutations',
                            'and copy number variants in the',
                            'investigated cohorts.',
                            'Frequencies of alterations present in at least 5%',
                            'of cancer samples are displayed.',
                            'Statistical significance for differences between',
                            'the cohorts was assessed by chi-squared test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Frequency of genetic alterations of the FGF- and FGFR-coding genes -------

  insert_msg('Frequency of alterations in FGF and FGFR genes')

  rep_tabs$fgf_genetics <- expl_genet$fgf_result_tbl %>%
    compress(names_to = 'Alteration type') %>%
    relocate(`Alteration type`) %>%
    mdtable(label = 'genetic_fgf_fgfr_features',
            ref_name = 'fgf_genetics',
            caption = paste('Frequency of somatic mutations and copy number',
                            'variants of FGF- and FGFR-coding genes in the',
                            'investigated cohorts.',
                            'Statistical significance for differences between',
                            'the cohorts was assessed by chi-squared test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Classification of somatic mutations of the FGFR genes ------

  insert_msg('Classification of mutations')

  rep_tabs$fgfr_class <- expl_fgfr$result_tbl

  for(i in names(rep_tabs$fgfr_class)) {

    names(rep_tabs$fgfr_class[[i]])[3] <-
      'Mutation classification'

  }

  rep_tabs$fgfr_class <- rep_tabs$fgfr_class %>%
    set_names(c('Mutation type', 'Variant type', 'Protein domain')) %>%
    compress(names_to = 'Classification scheme') %>%
    relocate(`Classification scheme`)

  rep_tabs$fgfr_class <- rep_tabs$fgfr_class %>%
    mdtable(label = 'fgfr_mutation_classification',
            ref_name = 'fgfr_class',
            caption = paste('Frequency of mutations in the FGFR1/2/3/4 genes',
                            'in the GENIE BLCA and TCGA BLCA cohorts',
                            'classified by mutation type, variant type, and',
                            'protein domain.',
                            'The table is available as a supplementary',
                            'Excel file.'))


# Clinical characteristic of samples with and without FGFR3 mutations ------

  insert_msg('Clinical characteristic, FGFR3 WT and FGFR3 mutated samples')

  rep_tabs$fgfr3_clinic <- fgfr_clinic$result_tbl

  names(rep_tabs$fgfr3_clinic)[c(3, 4)] <-
    c('FGFR3 wild-type', 'FGFR3 mutated')

  rep_tabs$fgfr3_clinic <- rep_tabs$fgfr3_clinic %>%
    mdtable(label = 'fgfr3_clinical_features',
            ref_name = 'fgfr3_clinic',
            caption = paste('Demographic, clinical, and pathological',
                            'characteristic of cancers with WT and mutated',
                            'FGFR3.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))

# Differential gene expression, FGFR3 WT and FGFR3 mutated -------

  insert_msg('Differential gene expresssion, FGFR3 WT and mutated')

  rep_tabs$fgfr3_dge <- expl_dge$result_tbl$FGFR3_mut %>%
    mdtable(label = 'fgfr3_differential_gene_expression',
            ref_name = 'fgfr3_dge',
            caption = paste('Differential expression of FGF- and FGFR-coding',
                            'genes in urothelial cancers with and without',
                            'FGFR3 mutations.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance was determined by',
                            'two-tailed T test with',
                            "Cohen's d effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Genetic alterations in the consensus molecular subtypes ----------

  insert_msg('Gnetics of the consensus classes')

  rep_tabs$cons_genetics <- sub_genet$result_tbl %>%
    mdtable(label = 'consensus_subsets_genetics',
            ref_name = 'cons_genetics',
            caption = paste('Frequency of somatic mutations and copy number',
                            'variants of FGF- and FGFR-coding genes',
                            'in the consensus molecular classes of',
                            'urothelial cancers.',
                            'The frequencies are presented as percentages of the',
                            'consensus class.',
                            'Statistical significance of differences between',
                            'the consensus classes was determined by',
                            "chi-squared test with Cramer's V effect size",
                            'statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Differential gene expression in the consensus subsets -------

  insert_msg('Differential gene expression in the consensus subtypes')

  rep_tabs$cons_dge <- sub_dge$result_tbl %>%
    mdtable(label = 'consensus_subsets_expression',
            ref_name = 'cons_dge',
            caption = paste('Expression of FGF- and FGFR-coding genes in',
                            'the consensus molecular classes of urothelial',
                            'cancers.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance of differences between',
                            'the molecular classes was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Genetics of the genetic subsets --------

  insert_msg('Genetics of the genetic subsets')

  rep_tabs$lca_genetics <- lca_eval$result_tbl %>%
    mdtable(label = 'genetic_subsets_genetics',
            ref_name = 'lca_genetics',
            caption = paste('Frequency of the most frequent somatic mutations',
                            'and copy number alterations in the genetic subsets',
                            'of urothelial cancers.',
                            'The frequencies are presented as percentages of',
                            'samples with genetic alterations in the genetic',
                            'subsets.',
                            'Differences between the genetic subsets were',
                            'assessed by chi-squared test with',
                            "Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Differential gene expression in the genetic subsets -------

  insert_msg('Genetic subsets: differential gene expression')

  rep_tabs$lca_dge <- lca_dge$result_tbl %>%
    mdtable(label = 'genetic_subsets_gene_expression',
            ref_name = 'lca_dge',
            caption = paste('Expression of FGF- and FGFR-coding genes',
                            'in the genetic subsets of urothelial cancers of',
                            'the TCGA BLCA cohort.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance of differences between',
                            'the genetic subsets was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Clinical characteristic of the genetic subsets ------

  insert_msg('Clinical characteristic of the genetic subsets')

  rep_tabs$lca_clinic <- lca_clinic$result_tbl %>%
    mdtable(label = 'genetic_subsets_clinics',
            ref_name = 'lca_clinic',
            caption = paste('Demographic, clinical, and pathological',
                            'characteristic of the genetic subsets of',
                            'urothelial cancers.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))

# Saving the tables on the disc ------

  insert_msg('Saving the tables on the disc')

  rep_tabs %>%
    save_excel(path = './report/tables.xlsx')

# END -----

  insert_tail()

# R/deseq_functions.R
# Functions for building and using DESeq2 Data Set

source('config.R')
source('R/packages.R')

# Filter low-expression genes
# Inputs: Raw count matrix, count threshold, minimum samples meeting threshold
# Outputs: Filtered matrix, summary stats,
#         matrix of samples with detection per gene.
filter_low_counts = function(count_matrix, min_count = config$min_count,
                             min_samples = config$min_samples) {
  genes_detected = rowSums(count_matrix > 0)
  genes_above    = rowSums(count_matrix >= min_count)
  keep           = genes_above >= min_samples
  
  list(
    counts   = count_matrix[keep, ],
    summary  = data.frame(
      total_genes    = nrow(count_matrix),
      detected       = sum(genes_detected > 0),
      passing_filter = sum(keep),
      removed        = sum(!keep)
    ),
    detection_per_gene = genes_detected
  )
}

# Build DESeq2 data set with size factor & variance-stabilised transformation
# Inputs: Filtered count matrix, sample metadata, design formula
# Outputs: DESeq data set and variance stabilised data
build_deseq = function(count_matrix, metadata, design = ~ group) {
  dds = suppressMessages(DESeqDataSetFromMatrix(count_matrix, metadata, design))
  dds = suppressMessages(estimateSizeFactors(dds))
  vsd = suppressMessages(vst(dds, blind = TRUE))
  
  # Print size factors with sample names
  sf = sizeFactors(dds)
  cat('Size factors:\n')
  for (s in names(sf)) {
    cat(sprintf('  %-12s %.4f\n', s, sf[s]))
  }
  
  list(
    dds = dds,
    vsd = vsd
  )
}

# Drop samples from data and build DESeq2 data set without them
# Inputs: Filtered count matrix, sample metadata, samples to drop, design formula
# Outputs: DESeq data set and variance stabilised data
drop_samples_deseq = function(count_matrix, metadata, drop, design = ~ group) {
  keep_samples = setdiff(colnames(count_matrix), drop)
  build_deseq(count_matrix[, keep_samples],
              droplevels(metadata[keep_samples, ]),
              design)
}

# Extract a PCA object of the top n most variable genes from DESeq2 data set
# Input: VSD, number of top genes to include
# Output: PCA object built using the n most variable genes in the data set
pca_top_n = function(vsd, n = config$pca_n_top) {
  rv = rowVars(assay(vsd))
  top_n = order(rv, decreasing = TRUE)[1:n]
  pca_object = prcomp(t(assay(vsd)[top_n,]))
  
  pca_object
}

# Generate differential expression analysis of a single contrast
# Inputs: DDS, contrast, alpha, independent filter status
# For pairwise contrast: contrast = vector c('factor, numerator, denominator)
# For coefficient contrast: contrast = coefficient name from resultsNames(dds)
# A fold change threshold was not included as an input due to small sample sizes
# Output: Data frame of results sorted by padj
extract_contrast = function(dds, contrast, alpha = 0.05,
                            independent_filter = TRUE) {
  # Determine contrast type
  is_pairwise = length(contrast) == 3
  
  if (is_pairwise) {
    contrast_label = sprintf('%s vs %s', contrast[2], contrast [3])
  } else {
    # Validate coefficient name
    valid_coeff = resultsNames(dds)
    if (!contrast %in% valid_coeff) {
      stop(sprintf('Coefficient %s not found. Available: %s', contrast,
                   paste(valid_coeff, collapse = ', ')))
    } else {contrast_label = contrast}
  }

  if (is_pairwise) {
    # Raw results for hypothesis testing and p-values
    res_raw    = results(dds, contrast = contrast, alpha = alpha,
                         independentFiltering = independent_filter)
    # Shrunk log fold change using ashr to estimate FC with small sample sizes
    res_shrunk = lfcShrink(dds, contrast = contrast, type = 'ashr', quiet = TRUE)
  } else {
    res_raw    = results(dds, name = contrast, alpha = alpha,
                         independentFiltering = independent_filter)
    res_shrunk = lfcShrink(dds, coef = contrast, type = 'ashr', quiet = TRUE)
  }

  # Assemble output
  output_df = data.frame(
    gene = rownames(res_raw),
    baseMean = res_raw$baseMean,
    log2FC_raw = res_raw$log2FoldChange,
    log2FC_shrunk = res_shrunk$log2FoldChange,
    lfcSE = res_raw$lfcSE,
    lfcSE_shrunk = res_shrunk$lfcSE,
    stat = res_raw$stat,
    pvalue = res_raw$pvalue,
    padj = res_raw$padj,
    contrast = contrast_label,
    direction = ifelse(
      is.na(res_raw$padj), 'NS',
      ifelse(res_raw$padj >= alpha, 'NS',
             ifelse(res_shrunk$log2FoldChange > 0, 'up', 'down'))
    )
  )
  
  # Order by adjusted p value
  output_df = output_df[order(output_df$padj, na.last = TRUE), ]
  output_df
}

# Run DEG analysis of multiple defined contrasts
# Inputs: DDS, contrast list, alpha, independent filter status
# Outputs: List with: - named list of data frames
#                     - Single combined data frame
#                     - Summary table
run_all_contrasts = function(dds, contrasts, alpha = 0.05,
                             independent_filter = TRUE) {
  
  # Extract contrasts and provide summary statistics to log
  per_contrast = list()
  
  for (name in names(contrasts)) {
    cat(sprintf('Extracting contrast: %s\n', name))
    per_contrast[[name]] = extract_contrast(dds, contrasts[[name]], alpha = alpha,
                                            independent_filter = independent_filter)
    
    sig = per_contrast[[name]] %>%
      filter(!is.na(padj) & padj < alpha)
    cat(sprintf('  Significant (padj < %.2f): %d up, %d down, %d total\n',
                alpha,
                sum(sig$direction == 'up'),
                sum(sig$direction == 'down'),
                nrow(sig)))
  }
  
  # Combine all significant results into one dataframe
  combined = bind_rows(per_contrast, .id = 'contrast')
  
  # Generate summary table
  summary_df = combined %>%
    filter(!is.na(padj)) %>%
    group_by(contrast, direction) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  
  # Output list
  list(
    per_contrast = per_contrast,
    combined = combined,
    summary = summary_df
  )
}

# Likelihood ratio test for overall group effect
# Use to identify DEGs, but do not use p-values from this test for significance
# testing with our sample sizes
# Inputs: DDS, alpha, reduced model formula 
# Reduced defaults to intercept only (~ 1) for pairwise comparisons
# Output: Data frame of LRT test results sorted by padj
likelihood_ratio = function(dds, alpha = 0.05, reduced = ~1) {
  dds_lrt = suppressMessages(DESeq(dds, test = 'LRT', reduced = reduced))
  res = results(dds_lrt, alpha = alpha)
  
  output_df = data.frame(
    gene = rownames(res),
    baseMean = res$baseMean,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj
  )
  
  output_df = output_df[order(output_df$padj, na.last = TRUE), ]
  cat(sprintf('LRT: %d genes significant at padj < %.2f (of %d tested)\n',
              sum(output_df$padj < alpha, na.rm = TRUE),
              alpha, nrow(output_df)))
  output_df
}

# Classify genes by their sensitivity to dsEHP_1's inclusion in the data
# - robust            = sig in both primary & sensitivity datasets
# - dsEHP_1_dependent = sig in sensitivity only
# - dsEHP_1_masked    = sig in primary only
# - not_significant
#
# Inputs: primary data, sensitivity data, alpha
# Outputs: df of genes with ds_EHP_1 sensitivity
classify_sensitivity = function(primary_df, sensitivity_df,
                                alpha = config$alpha) {
  # Build dfs of genes significant in either primary or sensitivity data
  sig_primary = primary_df %>%
    filter(!is.na(padj) & padj < alpha) %>%
    pull(gene)
  
  sig_sensitivity = sensitivity_df %>%
    filter(!is.na(padj) & padj < alpha) %>%
    pull(gene)
  
  # df of all genes in case filtering was different between datasets
  all_genes = union(primary_df$gene, sensitivity_df$gene)
  
  # Merge p-values and effect sizes into gene table for evaluation
  primary_info = primary_df %>%
    select(gene, log2FC_primary = log2FC_shrunk, padj_primary = padj)
  sensitivity_info = sensitivity_df %>%
    select(gene, log2FC_sensitivity = log2FC_shrunk, padj_sensitivity = padj)
  
  classification = data.frame(gene = all_genes) %>%
    left_join(primary_info, by = 'gene') %>%
    left_join(sensitivity_info, by = 'gene') %>%
    mutate(sig_in_primary = gene %in% sig_primary,
           sig_in_sensitivity = gene %in% sig_sensitivity,
           concordant = sign(log2FC_primary) == sign(log2FC_sensitivity),
           class = case_when(
             sig_in_primary & sig_in_sensitivity & concordant ~ 'robust',
             sig_in_primary & sig_in_sensitivity & !concordant ~ 'discordant',
             !sig_in_primary & sig_in_sensitivity ~ 'dsEHP_1_dependent',
             sig_in_primary & !sig_in_sensitivity ~ 'dsEHP_1_masked',
             TRUE ~ 'NS'
           ),
           class = factor(class, levels = c('robust', 'discordant',
                                            'dsEHP_1_dependent',
                                            'dsEHP_1_masked', 'NS'))) %>%
    arrange(class, padj_primary)
  
  classification
}


# Perform sensitivity analysis on a series of pairwise contrasts
# Inputs: primary data, sensitivity data, output path, list of contrasts, alpha
# Output: Per-contrast reports & overall summary
classify_sensitivity_mult = function(primary, sensitivity,
                                     comp_dir,
                                     pairwise_contrasts = config$pairwise_contrasts,
                                     alpha = config$alpha) {
  sens_results = list()
  
  for (name in names(pairwise_contrasts)) {
    comp = classify_sensitivity(
      primary$deg_results$per_contrast[[name]],
      sensitivity$deg_results$per_contrast[[name]],
      alpha = alpha
    )
    
    if ('discordant' %in% comp$class) {
      cat(sprintf('Discordant gene(s) found in %s: inspect manually\n', name))
    } 
    
    sens_results[[name]] = comp
    
    suppressMessages(write_csv(comp %>% filter(class != 'NS'),
                               file.path(comp_dir,
                                         sprintf('sensitivity_%s.csv', name))))
  }
  
  # Summary table
  comp_summary = bind_rows(
    lapply(names(sens_results), function(name) {
      sens_results[[name]] %>%
        filter(class != 'NS') %>%
        dplyr::count(class) %>%
        mutate(contrast = name)
    })
  ) %>% pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
    relocate(contrast)
  
  suppressMessages(write_csv(comp_summary,
                             file.path(comp_dir, 'sensitivity_summary.csv')))
  
  list(results = sens_results, comp_summary = comp_summary)
}
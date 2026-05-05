# R/pca_diagnostics.R
# Functions for extracting PCA axes to ID genes driving outliers

source('R/packages.R')

# Extract loadings from a PC, annotated with gene descriptions
# Inputs: PCA object, annotation table, PC choice
# Output: df of top loading genes by absolute loading
extract_pc_loadings = function(pca_obj, annotation_df, pc = 1) {
  
  loadings = data.frame(
    gene_id = rownames(pca_obj$rotation),
    loading = pca_obj$rotation[, pc]) %>%
    mutate(abs_loading = abs(loading),
           direction = if_else(loading > 0, 'positive', 'negative')) %>%
    arrange(desc(abs_loading))
  
  # Annotate gene info
  loadings = loadings %>%
    left_join(annotation_df %>% select(gene_id, symbol, description, annotated),
              by = 'gene_id') %>%
    mutate(label = case_when(
      annotated ~ description,
      !is.na(symbol) & symbol != gene_id ~ symbol,
      TRUE ~ gene_id
    ))
  
  # Variance explained by this PC
  var_explained = pca_obj$sdev^2 / sum(pca_obj$sdev^2)
  pct_var = round(100 * var_explained[pc], 1)
  
  # Attach values to loadings as metadata
  attr(loadings, 'pc') = pc
  attr(loadings, 'pct_var') = pct_var
  
  loadings
}

# Summarise loading attributes & functional themes using annotation terms
# Input: Loadings data frame, number of top genes to analyse
# Output: Summaries of loading printed to console
summarise_loading = function(loadings, n_top = 50) {
  pc_label = paste0('PC', attr(loadings, 'pc'))
  pct_var = attr(loadings, 'pct_var')
  cat(sprintf('\n%s accounts for %.1f%% of variance\n', pc_label, pct_var))
  
  top_loadings = head(loadings, n_top)
  described = filter(top_loadings, annotated)
  
  if (nrow(described) == 0) {
    cat(sprintf('\nNo annotated genes in %s top %d genes\n', pc_label, n_top))
    return(invisible(NULL))
  } 
  
  cat(sprintf('\nTop %d genes by absolute loading (Pos: %d, Neg: %d)\n',
              nrow(top_loadings),
              sum(top_loadings$direction == 'positive'),
              sum(top_loadings$direction == 'negative')))
  
  cat(sprintf(
    '\n%d of the top %d loaded genes in %s have functional annotation\n',
    nrow(described), nrow(top_loadings), pc_label))
  
  for (dir in c('positive', 'negative')) {
    subset = filter(described, direction == dir)
    if (nrow(subset) == 0) next
    
    cat(sprintf('\n%d annotated %s loadings in %s top %d. Top %d:\n',
                nrow(subset), dir, pc_label, n_top, min(10, nrow(subset))))
    
    for (i in seq_len(min(10, nrow(subset)))) {
      cat(sprintf('\n%s, (loading weight = %.4f): %s\n',
                  subset$gene_id[i], subset$loading[i], subset$description[i]))
    }
  }
  
  top_loadings
}

# Compare top loaded genes of a single sample to the rest of its group
# Inputs: VSD matrix (from assay(vsd)), metadata, loadings, sample name, number of genes
# Outputs: df of sample vs group mean for top loading genes
compare_sample_to_group = function(vsd_matrix, metadata, loadings,
                                   sample_name, n_top = 50) {
  pc_label = paste0('PC', attr(loadings, 'pc'))
  
  group = as.character(metadata[sample_name, 'group'])
  group_samples = rownames(metadata)[metadata$group == group]
  other_samples = setdiff(group_samples, sample_name)
  
  top_genes = head(loadings$gene_id, n_top)
  top_genes = top_genes[top_genes %in% rownames(vsd_matrix)]
  
  comparison = data.frame(
    gene_id = top_genes,
    sample_expression = vsd_matrix[top_genes, sample_name],
    group_mean = rowMeans(vsd_matrix[top_genes, other_samples, drop = FALSE])
  ) %>%
    mutate(
      log2_diff = sample_expression - group_mean,
      direction = if_else(log2_diff > 0, 'higher_in_sample', 'lower_in_sample')
    ) %>%
    left_join(loadings %>% select(gene_id, loading, symbol, description, label),
              by = 'gene_id') %>%
    arrange(desc(abs(log2_diff)))
  
  cat(sprintf('\n%s vs rest of %s group (top %d %s-loading genes):\n',
              sample_name, group, n_top, pc_label))
  cat(sprintf('\nHigher in %s: %d\nLower in %s: %d\n',
              sample_name, sum(comparison$direction == 'higher_in_sample'),
              sample_name, sum(comparison$direction == 'lower_in_sample')))
  cat(sprintf('\nMean log2 difference: %.2f\n', mean(abs(comparison$log2_diff))))

  comparison
}
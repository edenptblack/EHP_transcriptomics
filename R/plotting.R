# R/plotting.R
# Reusable variables and functions for the creation of plots

source('config.R')
source('R/packages.R')
source('R/annotation_functions.R')

# Aesthetics
group_colours = config$group_colours
group_labels = config$group_labels
direction_colours = config$direction_colours
group_order = config$group_order

# Library size bar plot
# Inputs: raw count matrix, sample metadata, title
# Output: Plot
plot_library_sizes = function(count_matrix,
                              metadata,
                              title = 'Library sizes (total mapped reads)') {
  lib_sizes = colSums(count_matrix)
  lib_df = data.frame(
    sample  = factor(names(lib_sizes), levels = names(lib_sizes)),
    millions = lib_sizes / 1e6,
    group   = metadata[names(lib_sizes), "group"]
  ) %>%
    mutate(group = fct_relevel(group, group_order)) %>%
    arrange(group) %>%
    mutate(sample = fct_inorder(sample))
  
  ggplot(lib_df, aes(x = sample, y = millions, fill = group)) +
    geom_col() +
    scale_fill_manual(values = group_colours, labels = group_labels) +
    labs(title = title, x = NULL, y = "Reads (millions)", fill = "Treatment") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    geom_hline(yintercept = mean(lib_sizes) / 1e6,
               linetype = "dashed", colour = "grey40")
}


# Gene detection distribution bar plot
# Input: Vector of number of samples containing > 0 count for each gene,
#       Number of samples in data set, title
# Output: Plot
plot_detection_distribution = function(genes_detected,
                                       n_samples,
                                       title = 'Gene detection across samples') {
  detect_df = data.frame(
    n_samples_detected = factor(genes_detected, levels = 0:n_samples)
  )
  
  ggplot(detect_df, aes(x = n_samples_detected)) +
    geom_bar() +
    labs(title = title,
         x = 'Number of samples with count > 0',
         y = 'Number of genes') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       breaks = seq(0, 10000, by = 1000)) +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}


# log2 count density distribution plot
# Inputs: Filtered count data, metadata, title
# Output: Plot
plot_count_density = function(filtered_count, metadata,
                              title = 'Count distribution per sample (filtered genes, log2)') {
  log_counts = log2(filtered_count + 1)
  log_melt = melt(log_counts)
  colnames(log_melt) = c('gene', 'sample', 'log2count')
  log_melt$group = metadata[as.character(log_melt$sample), 'group']
  log_melt$replicate = sub('^.+_', '', log_melt$sample)
  
  ggplot(log_melt, aes(x = log2count, 
                       colour = group,
                       linetype = replicate)) +
    geom_density() +
    scale_colour_manual(values = group_colours) +
    labs(title = title,
         x = 'log2(count + 1)',
         y = 'Density') +
    theme_bw(base_size = 13) +
    theme(legend.position = 'right')
}


# PCA plot using DESeq2 VST object
# Inputs: Variance-stabilised data, title
# Output: Plot
plot_pca = function(vsd, title = 'PCA of VST-transformed counts') {
  pca_data = plotPCA(vsd, intgroup = 'group', returnData = TRUE)
  percent_var = round(100 * attr(pca_data, 'percentVar'))
  
  ggplot(pca_data, aes(x = PC1, y = PC2, colour = group, label = name)) +
    geom_point() +
    geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = group_colours, labels = group_labels) +
    labs(title = title,
         x = sprintf('PC1 (%d%% variance)', percent_var[1]),
         y = sprintf('PC2 (%d%% variance)', percent_var[2])) +
    theme_bw(base_size = 13) +
    theme(legend.position = 'right')
}


# Scree plot
# Input: PCA object
# Output: Plot
plot_scree = function(pca_obj, title = 'Scree plot (500 most variable genes)') {
  var_explained = pca_obj$sdev^2 / sum(pca_obj$sdev^2) * 100
  scree_df = data.frame(PC = 1:length(var_explained),
                        variance = var_explained)
  
  ggplot(scree_df, aes(x = PC, y = variance)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(0, nrow(scree_df), by = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = 'Principal Component', y = '% Variance Explained') +
    theme_bw()
}


# Grid of PCA plots covering lower principal components
# Inputs: PCA object, metadata, number of principal components to include, title
# Output: Grid of PC pairwise scatter plots & density plots
plot_pca_grid = function(pca_obj, metadata, n = 3,
                         title = 'Grid of PCA pairwise combinations') {
  pc_df = data.frame(
    pca_obj$x[, 1:n],
    group = metadata[rownames(pca_obj$x), 'group']
  )
  
  ggpairs(pc_df, columns = 1:n, aes(colour = group),
          upper = list(continuous = 'blank'),
          diag = list(continuous = 'densityDiag')) +
    labs(title = title) +
    scale_colour_manual(values = group_colours) +
    scale_fill_manual(values = group_colours)
}


# Pearson correlation heatmap
# Inputs: VSD, metadata, title, file destination
# Output: Heatmap saved to file
plot_correlation_heatmap = function(vsd, metadata, filename,
                                    title = 'Sample-to-sample Pearson correlation (VST)') {
  correlation_matrix = cor(assay(vsd), method = 'pearson')
  
  annotation_columns = data.frame(
    Treatment = metadata[colnames(assay(vsd)), 'group'],
    row.names = colnames(assay(vsd))
  )
  
  pheatmap(
    correlation_matrix,
    color = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),
    breaks = seq(min(correlation_matrix) - 0.01, 1, length.out = 101),
    annotation_col = annotation_columns,
    annotation_colors = list(Treatment = group_colours),
    display_numbers = TRUE,
    number_format = '%.3f',
    fontsize_number = 7,
    fontsize = 11,
    main = title,
    filename = filename,
    width = 8,
    height = 8
  )
}


# Heatmap of a pre-selected gene matrix (Z-scored VST)
# Inputs: Gene expression matrix (genes x samples), metadata,
#         gene_labels named vector (or NULL), file name, title
# Output: Heatmap saved to file
plot_heatmap = function(gene_matrix, metadata, filename, title,
                        gene_labels = NULL) {
  gene_matrix_scaled = t(scale(t(gene_matrix)))
  
  # Replace rownames with readable labels where available
  if (!is.null(gene_labels)) {
    rownames(gene_matrix_scaled) = make.unique(
      resolve_labels(rownames(gene_matrix_scaled), gene_labels)
    )
  }
  
  annotation_columns = data.frame(
    Treatment = metadata[colnames(gene_matrix_scaled), 'group'],
    row.names = colnames(gene_matrix_scaled)
  )
  
  # Truncate labels to reasonable size
  rownames(gene_matrix_scaled) = strtrim(rownames(gene_matrix_scaled), 50)
  
  pheatmap(
    gene_matrix_scaled,
    color = colorRampPalette(c('blue', 'white', 'red'))(100),
    breaks = seq(-3, 3, length.out = 101),
    annotation_col = annotation_columns,
    annotation_colors = list(Treatment = group_colours),
    show_rownames = nrow(gene_matrix_scaled) <= 75,
    fontsize_row = 7,
    clustering_method = 'ward.D2',
    fontsize = 11,
    main = title,
    filename = filename,
    width = 8,
    height = 12
  )
}


# P-value histogram
# Input: Data frame from extract_contrast, title
# Output: Plot
plot_pvalue_histogram = function(deg_results_df, title = NULL) {
  if (is.null(title)) title = sprintf('Raw p-value distribution: %s',
                                      unique(deg_results_df$contrast))
  
  plot_df = deg_results_df %>% filter(!is.na(pvalue))
  
  ggplot(plot_df, aes(x = pvalue)) +
    geom_histogram(bins = 50) +
    geom_hline(yintercept = nrow(plot_df) * 0.02,
               linetype = 'dashed', colour = 'red') +
    labs(title = title, x = 'Raw p-value', y = 'Frequency') +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw(base_size = 13)
}

# Volcano plot
# Input: extract_contrast df, alpha, log2FC threshold, number of genes to label,
#       gene_labels named vector (or NULL), title (defaults to contrast name)
# Output: Plot
plot_volcano = function(deg_results_df, alpha = 0.05, lfc_line = 1, n_label = 20,
                        gene_labels = NULL, title = NULL) {
  if (is.null(title)) title = sprintf('Volcano plot: %s',
                                      unique(deg_results_df$contrast))
  
  plot_df = deg_results_df %>%
    filter(!is.na(padj)) %>%
    mutate(neg_log10p = -log10(pvalue),
           direction = factor(direction, levels = c('up', 'down', 'NS')))
  
  # Identify top genes for labeling
  top_genes = plot_df %>%
    filter(direction != 'NS') %>%
    slice_min(padj, n = n_label) %>%
    mutate(display_label = resolve_labels(gene, gene_labels))
  
  # Gene numbers for subtitle
  n_up = sum(plot_df$direction == 'up')
  n_down = sum(plot_df$direction == 'down')
  subtitle = sprintf('FDR < %.2f: %d up, %d down', alpha, n_up, n_down)
  
  ggplot(plot_df, aes(x = log2FC_shrunk, y = neg_log10p, colour = direction)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = direction_colours) +
    geom_vline(xintercept = c(-lfc_line, lfc_line),
               linetype = 'dashed', colour = 'grey') +
    geom_hline(yintercept = -log10(alpha),
               linetype = 'dashed', colour = 'grey') +
    geom_text_repel(data = top_genes, aes(label = display_label)) +
    labs(title = title, subtitle = subtitle,
         x = 'log2 fold change (shrunk)', y = '-log10 p-value') +
    theme_bw(base_size = 13)
}


# MA (log ratio (M) vs mean average (A)) plot
# Inputs: extract_contrast df, alpha, n_label,
#        gene_labels named vector (or NULL), title (defaults to contrast)
# Output: Plot
plot_ma = function(deg_results_df, alpha = 0.05, n_label = 15,
                   gene_labels = NULL, title = NULL) {
  if (is.null(title)) title = sprintf('MA plot: %s',
                                      unique(deg_results_df$contrast))
  
  plot_df = deg_results_df %>%
    filter(!is.na(padj)) %>%
    mutate(neg_log10p = -log10(pvalue),
           direction = factor(direction, levels = c('up', 'down', 'NS')))
  
  # Identify top genes for labeling
  top_genes = plot_df %>%
    filter(direction != 'NS') %>%
    slice_min(padj, n = n_label) %>%
    mutate(display_label = resolve_labels(gene, gene_labels))
  
  ggplot(plot_df, aes(x = log10(baseMean), y = log2FC_shrunk,
                      colour = direction)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = direction_colours) +
    geom_hline(yintercept = 0, colour = 'black') +
    geom_text_repel(data = top_genes, aes(label = display_label)) +
    labs(title = title,
         x = 'log10 mean normalised count', y = 'log2 fold change (shrunk)') +
    theme_bw(base_size = 13)
}


# Summary bar chart of DEGs across contrasts
# Inputs: summary df of all contrasts, alpha, title
# Output: Plot
plot_deg_summary = function(summary_df, alpha = 0.05,
                            title = 'Differentially expressed genes per contrast') {
  
  # Plot dataframe with number of DEGs per contrast separated by direction
  count_df = summary_df %>%
    filter(!is.na(padj) & padj < alpha) %>%
    mutate(direction = factor(direction, levels = c('up', 'down'))) %>%
    group_by(contrast, direction) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(directional_n = ifelse(direction == 'down', -n, n))
  
  ggplot(count_df, aes(x = contrast, y = directional_n, fill = direction)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = 'black') +
    geom_text(aes(label = n)) +
    scale_fill_manual(values = direction_colours[c('up', 'down')],
                      labels = c('Upregulated', 'Downregulated')) +
    labs(title = title, subtitle = sprintf('FDR < %.2f', alpha),
         y = 'Number of DEGs') +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}


# Boxplots of normalised counts for a single gene across groups
# Modified from earlier versions to highlight a sample of interest
# Inputs: Gene ID, DDS, metadata, title
# Output: Plot
plot_gene_counts = function(gene_id, dds, metadata,
                            highlight_sample = 'dsEHP_1',
                            gene_labels = NULL,
                            title = NULL) {
  if (is.null(title)) {
    title = resolve_labels(gene_id, gene_labels)
  }
  
  norm_counts = counts(dds, normalized = TRUE)
  if (!gene_id %in% rownames(norm_counts)) {
    cat(sprintf('Gene %s not found in DDS\n', gene_id))
    return(invisible(NULL))
  }
  
  plot_df = data.frame(
    sample = colnames(norm_counts),
    count = norm_counts[gene_id, ],
    group = metadata[colnames(norm_counts), 'group']
  ) %>%
    mutate(group = fct_relevel(group, config$group_order))
  
  # Highlight if the sample exists in this dataset
  has_highlight = highlight_sample %in% plot_df$sample
  
  if (has_highlight) {
    plot_df$is_highlight = plot_df$sample == highlight_sample
  }
  
  p = ggplot(plot_df, aes(x = group, y = count, colour = group))
  
  if (has_highlight) {
    p = p + geom_jitter(aes(shape = is_highlight),
                        width = 0.15) +
      scale_shape_manual(values = c('FALSE' = 16, 'TRUE' = 17),
                         guide = 'none')
  } else {
    p = p + geom_jitter(width = 0.15)
  }
  
  p + stat_summary(fun = mean, geom = 'crossbar', colour = 'black',
                   width = 0.4, linewidth = 0.3, fatten = 1.5) +
    stat_summary(fun.data = mean_se, geom = 'errorbar', colour = 'black',
                 width = 0.2, linewidth = 0.3) +
    scale_colour_manual(values = config$group_colours,
                        labels = config$group_labels) +
    labs(title = title, x = NULL, y = 'Normalised count') +
    theme_bw(base_size = 13)
}

# Concordance plot of log2FC values for sensitivity analysis
# Inputs: Sensitivity data, gene_labels named vector (or NULL), title
# Outputs: ggplot
plot_concordance = function(sens_df, title, gene_labels = NULL) {
  robust = sens_df %>% filter(class == 'robust')
  
  if (nrow(robust) < 3) {
    cat(sprintf('Contrast %s has too few robust genes\n',
                unique(sens_df$contrast)))
    return(invisible(NULL))
  }
  
  robust = robust %>%
    mutate(display_label = resolve_labels(gene, gene_labels))
  
  ggplot(robust, aes(x = log2FC_primary, y = log2FC_sensitivity)) +
    geom_point() +
    geom_text_repel(aes(label = display_label), size = 3, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'red') +
    labs(title = title, x = 'log2FC (primary)', y = 'log2FC (sensitivity)') +
    theme_bw(base_size = 13)
}

# Dot plot of GSEA enrichment
# Input: fgsea results df, n pathways to show, alpha, title
# Output: ggplot
plot_fgsea_dotplot = function(fgsea_res, n_show = 20, alpha = config$alpha,
                              title = 'GSEA enrichment') {
  sig = fgsea_res %>%
    filter(padj < alpha) %>%
    arrange(padj) %>%
    slice_head(n = n_show) %>%
    mutate(
      pathway = str_wrap(pathway, width = 40),
      pathway = fct_reorder(pathway, NES)
    )
  
  if (nrow(sig) == 0) {
    cat('No significant pathways to plot\n')
    return(invisible(NULL))
  }
  
  ggplot(sig, aes(x = NES, y = pathway, size = size, colour = padj)) +
    geom_point() +
    scale_colour_gradient(low = 'red', high = 'blue', name = 'padj') +
    scale_size_continuous(name = 'Gene set size', range = c(2, 8)) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey') +
    labs(title = title, x = 'Normalised Enrichment Score', y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
}

# Pathway-level heatmap of member genes (Z-scored VST expression)
# Inputs: Vector of pathway genes, vst matrix (genes x samples), metadata,
# filename, title, gene labels, annotations
# Output: Heatmap saved to file
plot_pathway_heatmap = function(gene_ids, vst_matrix, metadata, filename, title,
                                gene_labels = NULL, row_annotations = NULL) {
  
  present = gene_ids[gene_ids %in% rownames(vst_matrix)]
  
  if (length(present) < 2) {
    cat(sprintf('Too few genes (%d) for pathway heatmap: %s\n',
                length(present), title))
    return(invisible(NULL))
  }
  
  gene_matrix = vst_matrix[present, ]
  gene_matrix_scaled = t(scale(t(gene_matrix)))
  
  if (!is.null(gene_labels)) {
    rownames(gene_matrix_scaled) = make.unique(
      resolve_labels(rownames(gene_matrix_scaled), gene_labels)
    )
  }
  
  col_anno = data.frame(
    Treatment = metadata[colnames(gene_matrix_scaled), 'group'],
    row.names = colnames(gene_matrix_scaled)
  )
  
  row_anno = NULL
  anno_colors = list(Treatment = config$group_colours)
  
  if (!is.null(row_annotations)) {
    row_anno = row_annotations[present, , drop = FALSE]
    if (!is.null(gene_labels)) {
      rownames(row_anno) = make.unique(
        resolve_labels(present, gene_labels)
      )
    }
  }
  
  height = max(10, 5 + 0.4 * length(present))
  
  pheatmap(
    gene_matrix_scaled,
    color = colorRampPalette(c('blue', 'white', 'red'))(100),
    breaks = seq(-3, 3, length.out = 101),
    annotation_col = col_anno,
    annotation_row = row_anno,
    annotation_colors = anno_colors,
    show_rownames = TRUE,
    fontsize_row = 10,
    clustering_method = 'ward.D2',
    fontsize = 12,
    main = title,
    filename = filename,
    width = 12,
    height = height
  )
}
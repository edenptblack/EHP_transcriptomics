# R/de_pipeline.R
# Reusable DE pipeline applicable to all contrasts not using spore count
# Runs DESeq2 fitting, contrast extraction, LRT, plotting, and saving

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/plotting.R')
source('R/deseq_functions.R')
source('R/annotation_functions.R')

res_dir = 'results/03_differential_expression'

# Run complete DE analysis pipeline
#
# Inputs: dds; metadata;
# contrasts = 3-length vector for pairwise (c('factor', 'numerator', 'denominator'))
# string from resultsNames(dds) for coefficient;
# label (character string to ID run); alpha;
# gene_labels named vector (or NULL) for annotated plot labels;
# annotation_df (or NULL) for enriching exported CSV tables;
# lrt_reduced (formula for model - ~ 1 for pairwise,
# ~ vaccinated + challenged for multi-factor); vector of gene counts for heatmaps;
# number of top DEGs to plot per contrast; directory for processed data;
# dpi, height, and width for plots
#
# Outputs to disk: results/{label}/tables; results/{label}/plots; data/processed
# Returns: list - dds (fitted), deg_results, lrt_results, vsd
run_de_pipeline = function(dds, metadata, contrasts, label,
                           alpha = config$alpha,
                           gene_labels = NULL,
                           annotation_df = NULL,
                           lrt_reduced = ~ 1,
                           n_heatmap = c(50, 200),
                           n_gene_plots = 10,
                           proc_dir = config$proc_dir,
                           plot_dpi = config$plot_dpi,
                           plot_width = config$plot_width,
                           plot_height = config$plot_height) {
  # Set up directories
  plot_dir  = file.path(res_dir, label, 'plots')
  table_dir = file.path(res_dir, label, 'tables')
  gene_plot_dir = file.path(plot_dir, 'gene_counts')
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gene_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Fit DESeq2 model
  dds = suppressMessages(DESeq(dds))
  cat(sprintf('\nDispersion estimation complete. Design: %s\n',
              paste(as.character(design(dds)), collapse = ' ')))
  cat(sprintf('\nCoefficients available: %s\n',
              paste(resultsNames(dds), collapse = ', ')))
  
  # Dispersion plot
  png(file.path(plot_dir, '01_dispersion_estimates.png'),
      width = plot_width, height = plot_height, units = 'in', res = plot_dpi)
  plotDispEsts(dds, main = sprintf('Dispersion estimates (%s)', label))
  dev.off()
  
  
  # Extract and save pairwise contrasts
  deg_results = run_all_contrasts(dds, contrasts, alpha = alpha)
  
  # Annotate results for CSV export if annotation table provided
  if (!is.null(annotation_df)) {
    for (name in names(deg_results$per_contrast)) {
      deg_results$per_contrast[[name]] = annotate_results(
        deg_results$per_contrast[[name]], annotation_df
      )
    }
    deg_results$combined = annotate_results(deg_results$combined, annotation_df)
  }
  
  for (name in names(deg_results$per_contrast)) {
    suppressMessages(write_csv(deg_results$per_contrast[[name]],
                               file.path(table_dir, sprintf('DE_%s.csv', name))))
  }
  
  suppressMessages(write_csv(deg_results$combined,
                             file.path(table_dir, 'DE_all_contrasts_combined.csv')))
  
  suppressMessages(write_csv(deg_results$summary,
                             file.path(table_dir, 'DE_summary_counts.csv')))
  
  cat(sprintf('\nDifferentially expressed gene summary for %s\n', label))
  print(deg_results$summary)
  
  # Likelihood ratio test (effect across any groups)
  # Don't use p-value for analysis with pairwise due to small sample sizes
  # Especially dsPTP2
  lrt_results = likelihood_ratio(dds, alpha = alpha, reduced = lrt_reduced)
  suppressMessages(write_csv(lrt_results,
                             file.path(table_dir, 'DE_LRT.csv')))
  
  # P-value histograms for contrast diagnostics
  for (name in names(deg_results$per_contrast)) {
    p = plot_pvalue_histogram(deg_results$per_contrast[[name]])
    quiet_ggsave(file.path(plot_dir, sprintf('02_pvalue_hist_%s.png', name)),
                 p, width = plot_width, height = plot_height, dpi = plot_dpi)
  }
  
  # Volcano plots
  for (name in names(deg_results$per_contrast)) {
    p = plot_volcano(deg_results$per_contrast[[name]], alpha = alpha,
                     gene_labels = gene_labels)
    quiet_ggsave(file.path(plot_dir, sprintf('03_volcano_%s.png', name)),
                 p, width = plot_width, height = plot_height, dpi = plot_dpi)
  }
  
  # MA plots
  for (name in names(deg_results$per_contrast)) {
    p = plot_ma(deg_results$per_contrast[[name]], alpha = alpha,
                gene_labels = gene_labels)
    quiet_ggsave(file.path(plot_dir, sprintf('04_MA_%s.png', name)),
                 p, width = plot_width, height = plot_height, dpi = plot_dpi)
  }
  
  # Summary bar charts
  p_summary = plot_deg_summary(deg_results$combined, alpha = alpha)
  quiet_ggsave(file.path(plot_dir, '05_DEG_summary_barplot.png'),
               p_summary, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # Recompute VSD for trimmed data. Not blinded as model is already fit.
  vsd = suppressMessages(vst(dds, blind = FALSE))
  
  # Heatmaps of most significant DEGs
  for (n in n_heatmap) {
    top_genes = deg_results$combined %>%
      filter(!is.na(padj) & padj < alpha) %>%
      group_by(gene) %>%
      summarise(best_padj = min(padj), .groups = 'drop') %>%
      arrange(best_padj) %>%
      slice_head(n = n) %>%
      pull(gene)
    
    plot_heatmap(assay(vsd)[top_genes, ], metadata,
                 filename = file.path(plot_dir,
                                      sprintf('06_DEG_heatmap_top%d.png', n)),
                 title = sprintf('Heatmap of top %d DEGs (Z-scored VST)', n),
                 gene_labels = gene_labels)
  }
  
  # Per-gene count boxplots for top hits in key contrasts
  for (name in names(deg_results$per_contrast)) {
    res = deg_results$per_contrast[[name]]
    top_hits = res %>%
      filter(!is.na(padj) & padj < alpha) %>%
      slice_head(n = n_gene_plots)
    
    if (nrow(top_hits) == 0) {
      cat(sprintf('No significantly different genes for %s\n', name))
      next
    }
    
    for (i in seq_len(nrow(top_hits))) {
      gene_id = top_hits$gene[i]
      gene_display = resolve_labels(gene_id, gene_labels)
      p = plot_gene_counts(gene_id, dds, metadata,
                           title = sprintf('%s (padj = %.2f), %s',
                                           gene_display,
                                           top_hits$padj[i], name))
      quiet_ggsave(file.path(gene_plot_dir, sprintf('%s_rank_%d_%s.png',
                                                    name, i, gene_id)),
                   p, width = plot_width, height = plot_height, dpi = plot_dpi)
    }
    cat(sprintf('\nPlotted %d genes for %s\n', nrow(top_hits), name))
  }
  
  saveRDS(dds, file.path(proc_dir, sprintf('dds_%s_fitted.rds', label)))
  saveRDS(deg_results, file.path(proc_dir, sprintf('deg_results_%s.rds', label)))
  saveRDS(lrt_results, file.path(proc_dir, sprintf('lrt_results_%s.rds', label)))
  saveRDS(vsd, file.path(proc_dir, sprintf('vsd_%s.rds', label)))
  
  cat(sprintf('\nPipeline complete: %s\n', label))
  cat(sprintf('Tables: %s\n', table_dir))
  cat(sprintf('Plots:  %s\n', plot_dir))
  
  # Return objects for use by calling script
  invisible(list(
    dds = dds,
    deg_results = deg_results,
    lrt_results = lrt_results,
    vsd = vsd
  ))
}
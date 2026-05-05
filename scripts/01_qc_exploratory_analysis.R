# scripts/01_qc_exploratory_analysis.R
# QC and exploratory analysis of RNA-seq feature counts

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/metadata.R')
source('R/deseq_functions.R')
source('R/annotation_functions.R')
source('R/plotting.R')
source('R/pca_diagnostics.R')

main = function() {
  log_con = start_log('01_qc_exploratory_analysis')
  on.exit(end_log(log_con))
  cat('Configuration:\n')
  str(config)
  
  # Configs
  plot_dpi = config$plot_dpi
  plot_width = config$plot_width
  plot_height = config$plot_height
  
  res_dir = 'results/01_qc_exploratory_analysis'
  plot_dir = 'results/01_qc_exploratory_analysis/plots'
  proc_dir = config$proc_dir
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  detection_per_gene = readRDS('data/processed/detection_per_gene.rds')
  count_filtered = readRDS('data/processed/count_matrix_filtered.rds')
  coldata = readRDS('data/processed/coldata_full.rds')
  annotation_df = readRDS('data/processed/annotation_table.rds')
  
  # Build gene labels once for all plotting
  gene_labels = build_gene_labels(annotation_df)
  
  # Gene detection distribution
  p_detect = plot_detection_distribution(detection_per_gene, 12)
  ggsave(file.path(plot_dir, '00_gene_detection_distribution.png'),
         p_detect, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # Count density by sample
  p_density = plot_count_density(count_filtered, coldata)
  ggsave(file.path(plot_dir, '01_count_density.png'),
         p_density, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # Build DESeqDataSet & Variance stabilised data
  dds_list = build_deseq(count_filtered, coldata)
  
  # Check for flat mean:SD ratio in VSD
  png(file.path(plot_dir, '02_meanSD_vst.png'),
      width = plot_width, height = plot_height, units = 'in', res = plot_dpi)
  suppressMessages(meanSdPlot(assay(dds_list$vsd), ranks = FALSE))
  dev.off()
  
  # PCA of filtered counts
  p_pca = plot_pca(dds_list$vsd)
  ggsave(file.path(plot_dir, '03_pca.png'),
         p_pca, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # dsPTP2_1 is a clear outlier, accounting for all of PC1 (>50% total variance)
  # Potential tissue contamination. Analyse PC1
  cat('\n dsPTP2_1 outlier analysis\n')
  
  pca_obj_full = pca_top_n(dds_list$vsd)
  full_pca_loadings = extract_pc_loadings(pca_obj_full, annotation_df, pc = 1)
  top_genes_pc1_all_samples = summarise_loading(full_pca_loadings)
  dsPTP2_1_vs_group = compare_sample_to_group(assay(dds_list$vsd), coldata,
                                              full_pca_loadings, 'dsPTP2_1')
  
  write_csv(top_genes_pc1_all_samples,
            file.path(res_dir, 'all_samples_pc1_genes.csv'))
  write_csv(dsPTP2_1_vs_group, file.path(res_dir, 'dsPTP2_1_vs_group.csv'))
  
  # Rerun without dsPTP2_1
  coldata_drop_1 = build_metadata(
    colnames(subset(count_filtered, select = -dsPTP2_1)))
  
  dds_drop_outlier = drop_samples_deseq(count_filtered, coldata, 'dsPTP2_1')
  
  png(file.path(plot_dir, '04_meanSD_vst_drop_dsPTP2_1.png'),
      width = plot_width, height = plot_height, units = 'in', res = plot_dpi)
  suppressMessages(meanSdPlot(assay(dds_drop_outlier$vsd), ranks = FALSE))
  dev.off()
  
  p_pca_drop_dsPTP2_1 = plot_pca(dds_drop_outlier$vsd,
                                 title = 'PCA of VST-transformed counts (drop dsPTP2_1)')
  ggsave(file.path(plot_dir, '05_pca_drop_dsPTP2_1.png'),
         p_pca_drop_dsPTP2_1, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # dsEHP_1 is a smaller outlier. It is not clustering with EHP, making a
  # breakthrough infection explanation less likely. Spore count data shows
  # dsEHP_1 has ~50x lower spore count than other dsEHP_1 samples 
  # - PC1 as immune axis? Analyse PC1
  cat('\ndsEHP_1 outlier analysis\n')
  
  pca_obj_drop_1 = pca_top_n(dds_drop_outlier$vsd)
  drop_1_pca_loadings = extract_pc_loadings(pca_obj_drop_1, annotation_df, pc = 1)
  top_genes_pc1_drop_1 = summarise_loading(drop_1_pca_loadings)
  dsEHP_1_vs_group = compare_sample_to_group(assay(dds_drop_outlier$vsd),
                                             coldata_drop_1,
                                             drop_1_pca_loadings, 'dsEHP_1')
  
  write_csv(top_genes_pc1_drop_1,
            file.path(res_dir, 'drop_dsPTP2_1_pc1_genes.csv'))
  write_csv(dsEHP_1_vs_group, file.path(res_dir, 'dsEHP_1_vs_group.csv'))
  
  # Repeated without either for further exploration
  coldata_drop_2 = build_metadata(
    colnames(subset(count_filtered, select = -c(dsPTP2_1, dsEHP_1))))
  
  dds_drop_2_outliers = drop_samples_deseq(count_filtered,
                                           coldata, drop = c('dsPTP2_1', 'dsEHP_1'))
  
  png(file.path(plot_dir, '06_meanSD_vst_drop_outliers.png'),
      width = plot_width, height = plot_height, units = 'in', res = plot_dpi)
  suppressMessages(meanSdPlot(assay(dds_drop_2_outliers$vsd), ranks = FALSE))
  dev.off()
  
  p_pca_drop_outliers = plot_pca(dds_drop_2_outliers$vsd,
                                 title = 'PCA of VST-transformed counts (drop dsPTP2_1 & dsEHP_1)')
  ggsave(file.path(plot_dir, '07_pca_drop_outliers.png'),
         p_pca_drop_outliers, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # Extract 500 most variable genes for further analysis
  pca_obj_drop_2 = pca_top_n(dds_drop_2_outliers$vsd)
  p_scree = plot_scree(pca_obj_drop_2)
  ggsave(file.path(plot_dir, '08_scree_plot.png'),
         p_scree, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # PCA grid of 3 PCs accounting for >10% of variance
  p_pca_grid = plot_pca_grid(pca_obj_drop_2, coldata)
  ggsave(file.path(plot_dir, '09_pca_grid.png'),
         p_pca_grid, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  # Examine how all samples behave in pearson's heatmap
  plot_correlation_heatmap(dds_list$vsd,
                           coldata,
                           filename = file.path(plot_dir, '10_pearson_heatmap_all_samples.png'))
  
  # Heatmaps of most variable genes
  # Gene selection by variance, done here rather than inside the plotting function
  rv_drop_2 = rowVars(assay(dds_drop_2_outliers$vsd))
  
  for (n in c(50, 500)) {
    top_idx = order(rv_drop_2, decreasing = TRUE)[1:min(n, length(rv_drop_2))]
    plot_heatmap(assay(dds_drop_2_outliers$vsd)[top_idx, ],
                 coldata,
                 filename = file.path(plot_dir,
                                      sprintf('%d_top_%d_gene_heatmap.png',
                                              ifelse(n == 500, 11, 12), n)),
                 title = sprintf('Heatmap of %d most variable genes (Z-scored VST)', n),
                 gene_labels = gene_labels)
  }
  
  # Scree, PCA grid & Heatmaps including dsEHP_1
  pca_obj_dsEHP_1 = pca_top_n(dds_drop_outlier$vsd)
  p_scree_dsEHP_1 = plot_scree(pca_obj_dsEHP_1)
  ggsave(file.path(plot_dir, '13_scree_plot_dsEHP_1.png'),
         p_scree_dsEHP_1, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  p_pca_grid_dsEHP_1 = plot_pca_grid(pca_obj_dsEHP_1, coldata_drop_1)
  ggsave(file.path(plot_dir, '14_pca_grid_dsEHP_1.png'),
         p_pca_grid_dsEHP_1, width = plot_width, height = plot_height, dpi = plot_dpi)
  
  plot_correlation_heatmap(dds_drop_outlier$vsd,
                           coldata_drop_1,
                           filename = file.path(plot_dir, '15_pearson_heatmap_dsEHP_1.png'))
  
  # Heatmaps including dsEHP_1
  rv_drop_1 = rowVars(assay(dds_drop_outlier$vsd))
  
  for (n in c(50, 500)) {
    top_idx = order(rv_drop_1, decreasing = TRUE)[1:min(n, length(rv_drop_1))]
    plot_heatmap(assay(dds_drop_outlier$vsd)[top_idx, ],
                 coldata,
                 filename = file.path(plot_dir,
                                      sprintf('%d_top_%d_gene_heatmap_dsEHP_1.png',
                                              ifelse(n == 500, 16, 17), n)),
                 title = sprintf('Heatmap of %d most variable genes (Z-scored VST, with dsEHP_1)', n),
                 gene_labels = gene_labels)
  }
  
  # Save processed data for downstream use
  saveRDS(dds_list$dds, file.path(proc_dir, 'dds_full.rds'))
  saveRDS(dds_drop_outlier$dds, file.path(proc_dir, 'dds_drop_dsPTP2_1.rds'))
  saveRDS(dds_drop_2_outliers$dds, file.path(proc_dir, 'dds_drop_both.rds'))
  saveRDS(coldata_drop_2, file.path(proc_dir, 'coldata_drop_2.rds'))
  saveRDS(coldata_drop_1, file.path(proc_dir, 'coldata_drop_dsPTP2_1.rds'))
  
  
  cat(paste('\nProcessed data saved to', proc_dir))
}

main()
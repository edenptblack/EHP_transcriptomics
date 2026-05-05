# scripts/00_data_preparation.R
# Filtering and annotation of raw data

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/metadata.R')
source('R/plotting.R')
source('R/deseq_functions.R')
source('R/annotation_functions.R')

main = function() {
  log_con = start_log('00_data_preparation')
  on.exit(end_log(log_con))
  cat('Configuration:\n')
  str(config)
  
  # Constants & variables
  FT_URL  = paste0('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/',
                   config$ASSEMBLY_TAG, '/',
                   config$ASSEMBLY_TAG, '_feature_table.txt.gz')
  
  ref_dir = config$ref_dir
  proc_dir = config$proc_dir
  res_dir = 'results/00_data_preparation'
  dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir,  recursive = TRUE, showWarnings = FALSE)
  
  # Load and prepare data
  counts_raw = read_csv(file.path(config$raw_dir, 'EHP_featurecounts.csv'))
  count_matrix = as.matrix(counts_raw[, -1])
  rownames(count_matrix) = counts_raw$Geneid
  
  coldata = build_metadata(colnames(count_matrix), 'NaCl')
  
  # Library sizes
  lib_sizes = colSums(count_matrix)
  cat('Library sizes (millions)\n')
  
  for (i in names(lib_sizes)) {
    cat(sprintf('  %-12s %.1f M\n', i, lib_sizes[i] / 1e6))
  }
  cat(sprintf('Range: %.1f M - %.1f M\nFold difference: %.1fx\n',
              min(lib_sizes/1e6),
              max(lib_sizes/1e6),
              max(lib_sizes)/min(lib_sizes)))
  
  p_lib = plot_library_sizes(count_matrix, coldata)
  ggsave(file.path(res_dir, 'library_sizes.png'),
         p_lib, width = 8, height = 5, dpi = 300)
  
  # Filter low-expression genes
  filtered_genes = filter_low_counts(count_matrix)
  print(filtered_genes$summary)
  count_filtered = filtered_genes$counts
  
  # Download feature table
  ft_dest = file.path(ref_dir, basename(FT_URL))
  
  if (!file.exists(ft_dest)) {
    dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    cat('Downloading feature table...\n')
    download.file(FT_URL, ft_dest, mode = 'wb', quiet = TRUE)
  } else {
    cat('Using cached feature table\n')
  }
  
  ft = parse_feature_table(ft_dest)
  
  # Check assembly matches using LOC-prefixed genes from filtered data
  filtered_loc = grep('^LOC', rownames(count_filtered), value = TRUE)
  match_rate = mean(filtered_loc %in% ft$gene_id)
  
  cat(sprintf('Assembly check: %d / %d LOC genes matched (%.1f%%)\n',
              sum(filtered_loc %in% ft$gene_id),
              length(filtered_loc),
              100 * match_rate))
  
  if (match_rate < 0.95) {
    stop(sprintf('Only %.1f%% of LOC gene IDs match assembly %s. Check reference\n',
                 100 * match_rate, config$ASSEMBLY_TAG))
  }
  
  # Download gene2go, build & summarise annotation table
  g2g_dest = file.path(ref_dir, basename(config$G2G_URL))
  
  if(!file.exists(g2g_dest)) {
    cat('Downloading gene2go (may take several minutes)...\n')
    download.file(config$G2G_URL, g2g_dest, mode = 'wb', quiet = TRUE)
  } else {
    cat('Using cached gene2go\n')
  }
  
  annotation_df = build_annotation_table(ft_dest, g2g_dest, config$TAXID)
  
  summarise_annotation_coverage(annotation_df, rownames(count_filtered))
  
  # Save processed objects
  saveRDS(count_matrix, file.path(proc_dir, 'count_matrix_raw.rds'))
  saveRDS(count_filtered, file.path(proc_dir, 'count_matrix_filtered.rds'))
  saveRDS(coldata, file.path(proc_dir, 'coldata_full.rds'))
  saveRDS(annotation_df, file.path(proc_dir, 'annotation_table.rds'))
  saveRDS(filtered_genes$detection_per_gene,
          file.path(proc_dir, 'detection_per_gene.rds'))
  
  cat(sprintf('\nProcessed data saved to %s/\n', proc_dir))
}

main()

# scripts/02_build_subsets
# Build all data subsets required for differential expression analyses
# For those already built, apply structured naming convention

# Subsets to build:
# primary     - 10 samples by group (drop dsPTP2_1 & dsEHP_1)
# sensitivity - 11 samples by group (drop dsPTP2_1 only)
# multifactor - 10 samples by vaccinated + challenged + vaccinated:challenged
# multifactor sensitivity - Multifactor + dsEHP_1

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/metadata.R')
source('R/deseq_functions.R')

main = function() {
  log_con = start_log('02_build_subsets')
  on.exit(end_log(log_con))
  cat('Configuration:\n')
  str(config)
  
  proc_dir = config$proc_dir
  
  # Load data
  count_raw = readRDS(file.path(proc_dir, 'count_matrix_raw.rds'))
  count_filtered = readRDS(file.path(proc_dir, 'count_matrix_filtered.rds'))
  
  # Rename primary dataset (import from 01_qc_exploratory analysis)
  dds_primary = readRDS(file.path(proc_dir, 'dds_drop_both.rds'))
  coldata_primary = readRDS(file.path(proc_dir, 'coldata_drop_2.rds'))
  cat('\nGroup sizes dds_primary:\n')
  print(table(coldata_primary$group))
  
  # Rename dataset for dsEHP_1 sensitivity analysis (import from 01_qc)
  dds_sensitivity = readRDS(file.path(proc_dir, 'dds_drop_dsPTP2_1.rds'))
  coldata_sensitivity = readRDS(file.path(proc_dir, 'coldata_drop_dsPTP2_1.rds'))
  cat('\nGroup sizes dds_sensitivity:\n')
  print(table(coldata_sensitivity$group))
  
  # Build dataset for multifactor analysis
  coldata_multifactor = coldata_primary
  cat('\nMultifactor metadata check:\n')
  print(coldata_multifactor[, c('group', 'vaccinated', 'challenged')])
  
  dds_multifactor = drop_samples_deseq(count_filtered, coldata_multifactor,
                                       c('dsPTP2_1', 'dsEHP_1'),
                                       ~ vaccinated + challenged + vaccinated:challenged)
  dds_multifactor = dds_multifactor$dds
  
  # Sensitivity multifactor (which interaction terms are dsEHP_1 dependent?)
  dds_multifactor_sens = drop_samples_deseq(count_filtered, coldata_sensitivity,
                                            'dsPTP2_1',
                                            ~ vaccinated + challenged + vaccinated:challenged)
  dds_multifactor_sens = dds_multifactor_sens$dds
  
  # Save for downstream use
  datasets = list(
    primary      = list(dds = dds_primary,       coldata = coldata_primary),
    primary_sens = list(dds = dds_sensitivity,   coldata = coldata_sensitivity),
    multifactor  = list(dds = dds_multifactor,   coldata = coldata_multifactor),
    multifactor_sens = list(dds = dds_multifactor_sens, coldata = coldata_sensitivity)
  )
  
  for (name in names(datasets)) {
    saveRDS(datasets[[name]]$dds,
            file.path(proc_dir, sprintf('dds_%s.rds', name)))
    saveRDS(datasets[[name]]$coldata,
            file.path(proc_dir, sprintf('coldata_%s.rds', name)))
  }
  
  # Summarise datasets
  cat('\nDataset summary:\n')
  summary_df = data.frame(
    samples = sapply(datasets, function(x) ncol(x$dds)),
    design = sapply(datasets, function(x) deparse(design(x$dds)))
  )
  print(summary_df)
  
  cat(sprintf('\nAll datasets saved to %s\n', proc_dir))
}

main()
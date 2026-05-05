# scripts/03_differential_expression.R
# Runs DE pipeline on primary, sensitivty, and multifactor datasets
# Performs sensitivity analysis to dsEHP_1 by comparing primary vs sensitivity

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/metadata.R')
source('R/deseq_functions.R')
source('R/annotation_functions.R')
source('R/plotting.R')
source('R/de_pipeline.R')

main = function() {
  log_con = start_log('03_differential_expression')
  on.exit(end_log(log_con))
  cat('Configuration:\n')
  str(config)
  
  proc_dir = config$proc_dir
  alpha = config$alpha
  res_dir = 'results/03_differential_expression'
  sens_dir = file.path(res_dir, 'sensitivity_analysis')
  comp_dir = file.path(sens_dir, 'sensitivity')
  plot_dir = file.path(sens_dir, 'plots')
  pairwise_contrasts = config$pairwise_contrasts
  
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load annotation and build gene labels
  annotation_df = readRDS(file.path(proc_dir, 'annotation_table.rds'))
  gene_labels = build_gene_labels(annotation_df)
  
  # Primary dataset (no dsPTP2_1 or dsEHP_1) analysis
  cat('\nPrimary analysis:\n')
  
  dds_primary     = readRDS(file.path(proc_dir, 'dds_primary.rds'))
  coldata_primary = readRDS(file.path(proc_dir, 'coldata_primary.rds'))
  
  primary = run_de_pipeline(
    dds     = dds_primary,
    metadata = coldata_primary,
    contrasts = pairwise_contrasts,
    label   = 'primary_de',
    gene_labels = gene_labels,
    annotation_df = annotation_df
  )
  
  # Sensitivity dataset (no dsPTP2_1) analysis
  cat('\nSensivity analysis: \n')
  
  dds_sensitivity     = readRDS(file.path(proc_dir, 'dds_primary_sens.rds'))
  coldata_sensitivity = readRDS(file.path(proc_dir, 'coldata_primary_sens.rds'))
  
  sensitivity = run_de_pipeline(
    dds     = dds_sensitivity,
    metadata = coldata_sensitivity,
    contrasts = pairwise_contrasts,
    label   = 'sensitivity_de',
    gene_labels = gene_labels,
    annotation_df = annotation_df
  )
  
  # Multi-factor interaction model
  # vaccinated + challenged + vaccinated:challenged
  cat('Multifactor analysis:\n')
  
  dds_multifactor     = readRDS(file.path(proc_dir, 'dds_multifactor.rds'))
  coldata_multifactor = readRDS(file.path(proc_dir, 'coldata_multifactor.rds'))
  
  multifactor_contrasts = list(
    vaccination_effect = 'vaccinated_yes_vs_no',
    challenge_effect   = 'challenged_yes_vs_no',
    interaction        = 'vaccinatedyes.challengedyes'
  )
  
  multifactor = run_de_pipeline(
    dds = dds_multifactor,
    metadata = coldata_multifactor,
    contrasts = multifactor_contrasts,
    label = 'multifactor_de',
    gene_labels = gene_labels,
    annotation_df = annotation_df,
    lrt_reduced = ~ vaccinated + challenged
  )
  
  cat('\nMultifactor sensitivity (with dsEHP_1):\n')
  
  dds_multifactor_sens     = readRDS(file.path(proc_dir, 'dds_multifactor_sens.rds'))
  coldata_multifactor_sens = readRDS(file.path(proc_dir, 'coldata_multifactor_sens.rds'))
  
  multifactor_sens = run_de_pipeline(
    dds = dds_multifactor_sens,
    metadata = coldata_multifactor_sens,
    contrasts = multifactor_contrasts,
    label = 'multifactor_de_sensitivity',
    gene_labels = gene_labels,
    annotation_df = annotation_df,
    lrt_reduced = ~ vaccinated + challenged
  )
  
  # Sensitivity comparison
  cat('Sensitivity comparison:\n')
  
  sensitivity_all = classify_sensitivity_mult(primary, sensitivity,
                                              comp_dir = comp_dir)
  print(sensitivity_all$comp_summary)
  sens_res = sensitivity_all$results
  
  # Venn Diagrams of sensitivity analysis
  for (name in names(sens_res)) {
    sig_primary = sens_res[[name]] %>% 
      filter(sig_in_primary) %>% pull(gene)
    sig_sens = sens_res[[name]] %>% 
      filter(sig_in_sensitivity) %>% pull(gene)
    
    p_venn = ggVennDiagram(list(Primary = sig_primary, Sensitivity = sig_sens)) +
      labs(title = sprintf('DEG overlap: %s', name))
    quiet_ggsave(file.path(plot_dir, sprintf('venn_%s.png', name)),
                 p_venn, width = config$plot_width, height = config$plot_height,
                 dpi = config$plot_dpi)
  }
  
  # Concordance plots of fold changes in robust genes
  for (name in names(sens_res)) {
    p_conc = plot_concordance(sens_res[[name]],
                              title = sprintf(
                                'Robust DEG effect size concordance: %s', name),
                              gene_labels = gene_labels)
    
    if (!is.null(p_conc)) {
      quiet_ggsave(file.path(plot_dir, sprintf('concordance_%s.png', name)),
                   p_conc, width = config$plot_width, height = config$plot_height,
                   dpi = config$plot_dpi)
    }
  }
  
  # Interaction term sensitivity comparison
  cat('\nInteraction term sensitivity comparison:\n')
  
  mf_sens_comp = classify_sensitivity(
    multifactor$deg_results$per_contrast$interaction,
    multifactor_sens$deg_results$per_contrast$interaction
  )
  
  mf_class_counts = mf_sens_comp %>%
    filter(class != 'NS') %>%
    dplyr::count(class)
  cat('Interaction sensitivity classification:\n')
  print(as.data.frame(mf_class_counts))
  
  suppressMessages(write_csv(mf_sens_comp %>% filter(class != 'NS'),
                             file.path(comp_dir, 'sensitivity_interaction.csv')))
  
  # Venn diagram for interaction
  sig_mf_with_dsEHP1 = mf_sens_comp %>% filter(sig_in_sensitivity) %>% pull(gene)
  sig_mf_without_dsEHP1 = mf_sens_comp %>% filter(sig_in_primary) %>% pull(gene)
  
  p_venn_mf = ggVennDiagram(list('With dsEHP_1' = sig_mf_with_dsEHP1,
                                 'Without dsEHP_1' = sig_mf_without_dsEHP1)) +
    labs(title = 'Interaction DEG overlap: dsEHP_1 sensitivity')
  quiet_ggsave(file.path(plot_dir, 'venn_interaction_sensitivity.png'),
               p_venn_mf, width = config$plot_width, height = config$plot_height,
               dpi = config$plot_dpi)
  
  
  # Save sensitivity results
  saveRDS(sens_res, file.path(proc_dir, 'sensitivity_comparison.rds'))
  saveRDS(mf_sens_comp, file.path(proc_dir, 'sensitivity_interaction.rds'))
  
  cat(sprintf('Results saved to %s\n', res_dir))
}

main()
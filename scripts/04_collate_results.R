# scripts/04_collate_results.R
# Build master gene table merging all analyses, flag priority genes for HHPred
# Prepare GSEA ranked lists & ORA gene sets


source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/annotation_functions.R')

main = function() {
  log_con = start_log('04_collate_results')
  on.exit(end_log(log_con))
  
  proc_dir = config$proc_dir
  alpha    = config$alpha
  res_dir  = 'results/04_collate_results'
  gsea_dir = file.path(res_dir, 'gsea_ranks')
  ora_dir  = file.path(res_dir, 'ora_sets')
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load all results
  annotation_df = readRDS(file.path(proc_dir, 'annotation_table.rds'))
  
  primary_results = readRDS(file.path(proc_dir, 'deg_results_primary_de.rds'))
  multifactor_results = readRDS(file.path(proc_dir,
                                          'deg_results_multifactor_de.rds'))
  lrt_results = readRDS(file.path(proc_dir, 'lrt_results_primary_de.rds'))
  pairwise_sensitivity = readRDS(file.path(proc_dir,
                                           'sensitivity_comparison.rds'))
  interaction_sensitivity = readRDS(file.path(proc_dir,
                                              'sensitivity_interaction.rds'))
  
  cat('All result objects loaded\n')
  
  #Build base table with annotation
  master = primary_results$per_contrast[[1]] %>%
    select(gene, baseMean) %>%
    left_join(annotation_df, by = c('gene' = 'gene_id'))
  
  cat(sprintf('Base table: %d genes\n', nrow(master)))
  cat(sprintf('Annotated: %d (%.1f%%)\n',
              sum(master$annotated, na.rm = TRUE),
              100 * mean(master$annotated, na.rm = TRUE)))
  
  # Merge priority pairwise contrasts (primary analysis)
  priority_pairwise = c('EHP_vs_NaCl', 'dsEHP_vs_EHP',
                        'dsPTP2_vs_NaCl', 'dsEHP_vs_dsPTP2', 'dsEHP_vs_NaCl')
  
  cat(sprintf('\nMerging %d pairwise contrasts\n', length(priority_pairwise)))
  
  for (name in priority_pairwise) {
    res = primary_results$per_contrast[[name]] %>%
      select(gene, log2FC_shrunk, lfcSE_shrunk, padj, direction, stat)
    names(res)[-1] = paste0(c('lfc_', 'lfcSE_', 'padj_', 'dir_', 'stat_'), name)
    
    master = left_join(master, res, by = 'gene')
  }

  # Merge LRT results
  lrt_info = lrt_results %>%
    select(gene, padj_LRT = padj, stat_LRT = stat)
  
  master = left_join(master, lrt_info, by = 'gene')

  # Merge multifactor interaction term
  
  # Only the interaction term — the main effects duplicate the pairwise
  # contrasts (vaccination_effect = dsPTP2 vs NaCl, challenge_effect =
  # EHP vs NaCl, confirmed by identical DEG counts).
  
  int_res = multifactor_results$per_contrast$interaction %>%
    select(gene, log2FC_shrunk, lfcSE_shrunk, padj, direction, stat)
  names(int_res)[-1] = paste0(c('lfc_', 'lfcSE_', 'padj_', 'dir_', 'stat_'),
                              'interaction')
  
  master = left_join(master, int_res, by = 'gene')
  cat('Merged interaction term\n')
  
  # Merge pairwise sensitivity classification (dsEHP vs EHP)
  if ('dsEHP_vs_EHP' %in% names(pairwise_sensitivity)) {
    pw_sens = pairwise_sensitivity$dsEHP_vs_EHP %>%
      select(gene,
             sens_class_pairwise = class,
             lfc_sens_pairwise = log2FC_sensitivity,
             padj_sens_pairwise = padj_sensitivity)
    
    master = left_join(master, pw_sens, by = 'gene')
    cat('Merged pairwise dsEHP vs EHP sensitivity classification\n')
  }
  
  # 6. Merge interaction sensitivity classification
  
  # In classify_sensitivity, first arg = "primary" (without dsEHP_1),
  # second arg = "sensitivity" (with dsEHP_1).
  # So: sig_in_sensitivity = sig with dsEHP_1 included
  #     sig_in_primary = sig without dsEHP_1
  #     dsEHP_1_dependent = sig only when dsEHP_1 included
  #     dsEHP_1_masked = sig only when dsEHP_1 excluded
  
  int_sens = interaction_sensitivity %>%
    select(gene,
           sens_class_interaction = class,
           lfc_interaction_with_dsEHP1 = log2FC_sensitivity,
           padj_interaction_with_dsEHP1 = padj_sensitivity)
  
  master = left_join(master, int_sens, by = 'gene')
  cat('Merged interaction sensitivity classification\n')
  
  # Significance flags
  master = master %>%
    mutate(
      sig_EHP_vs_NaCl    = !is.na(padj_EHP_vs_NaCl) & padj_EHP_vs_NaCl < alpha,
      sig_dsEHP_vs_EHP   = !is.na(padj_dsEHP_vs_EHP) & padj_dsEHP_vs_EHP < alpha,
      sig_dsPTP2_vs_NaCl = !is.na(padj_dsPTP2_vs_NaCl) & padj_dsPTP2_vs_NaCl < alpha,
      sig_interaction    = !is.na(padj_interaction) & padj_interaction < alpha,
      sig_dsEHP_vs_NaCl  = !is.na(padj_dsEHP_vs_NaCl) & padj_dsEHP_vs_NaCl < alpha,
      sig_LRT            = !is.na(padj_LRT) & padj_LRT < alpha,
      n_priority_sig = sig_EHP_vs_NaCl + sig_dsEHP_vs_EHP +
        sig_dsPTP2_vs_NaCl + sig_interaction,
      # LRT excluded from sig_any — unreliable p-values at this sample size
      # Kept as separate flag for manual investigation
      sig_any = n_priority_sig > 0
    )
  
  cat('\nSignificance summary:\n')
  cat(sprintf('EHP vs NaCl: %d\n', sum(master$sig_EHP_vs_NaCl)))
  cat(sprintf('dsEHP vs EHP: %d\n', sum(master$sig_dsEHP_vs_EHP)))
  cat(sprintf('dsPTP2 vs NaCl: %d\n', sum(master$sig_dsPTP2_vs_NaCl)))
  cat(sprintf('Interaction: %d\n', sum(master$sig_interaction)))
  cat(sprintf('LRT (separate): %d\n', sum(master$sig_LRT)))
  cat(sprintf('Any (excl LRT): %d\n', sum(master$sig_any)))
  
  # Functional category classification
  immune_keywords = c(
    # Pattern recognition receptors
    'lectin', 'galectin', 'toll', 'peptidoglycan', 'beta.*glucan',
    'lipopolysaccharide', 'LPS', 'PGRP', 'LGBP', 'BGBP',
    'scavenger receptor', 'mannose receptor',
    # Signalling cascades
    'MyD88', 'TRAF', 'Dorsal', 'Relish', 'Cactus', 'IKK',
    'IRAK', 'pellino', 'NF.kB', 'IMD', 'JNK', 'JAK',
    'Spätzle', 'Spatzle', 'IRF', 'interferon regulatory factor',
    'STING', 'FoxO', 'MAPK', 'Wnt',
    'TGF.*beta', 'Basigin',
    # proPO cascade
    'phenoloxidase', 'prophenoloxidase', 'proPO', 'clip.*serine',
    'serine protease', 'serine proteinase', 'serpin',
    # Antimicrobial peptides and effectors
    'antimicrobial', 'anti.lipopolysaccharide', 'crustin', 'penaeidin',
    'lysozyme', 'defensin',
    # Oxidative and nitrosative burst
    'NADPH oxidase', 'NOX', 'superoxide dismutase',
    'peroxidase', 'catalase', 'glutathione.*transferase',
    'thioredoxin', 'peroxiredoxin', 'dual oxidase', 'DUOX',
    'nitric oxide synthase', '\\biNOS\\b',
    # Protease inhibitors
    'alpha.2.macroglobulin', 'kazal', 'kunitz', 'WAP',
    'cystatin', 'pacifastin',
    # Phagocytosis, clotting, complement
    'integrin', 'transglutaminase', 'clotting', 'hemocyanin',
    'complement', 'thioester',
    # Autophagy
    'autophagy', 'Beclin', '\\bATG\\d', '\\bmTOR\\b', '\\bTOR\\b',
    # Apoptosis
    'caspase', 'apoptosis', 'Bcl', 'IAP', 'inhibitor of apoptosis',
    '\\bp53\\b',
    # Cytokines and signalling molecules
    'interferon', 'cytokine', 'astakine', 'Vago',
    # Other immune-associated
    'heat shock', 'HSP',
    'venom', 'tachylectin', 'fibrinogen', 'ficolin',
    'Down syndrome.*cell', 'DSCAM',
    'matrix metalloproteinase', 'MMP',
    'cathepsin', 'chitinase',
    'NLRP', 'NACHT',
    'calreticulin', 'granulin',
    'Krüppel', 'Kruppel',
    'peritrophin',
    'JAK.STAT|\\bSTAT'
  )
  
  rnai_keywords = c(
    # siRNA pathway (somatic antiviral — core)
    'Dicer', 'Argonaute', '\\bAGO\\b',
    'R2D2', 'Loquacious', 'LOQS',
    # dsRNA uptake
    'SID.1', 'Sid.1',
    # miRNA pathway (somatic, upstream processing)
    'Drosha', 'Pasha', 'DGCR8', 'Ars2',
    # piRNA pathway
    'PIWI', 'piRNA', 'Tudor',
    # General RNAi terms
    'RNA.dependent RNA polymerase', 'RdRp',
    'RISC', 'RNA interference'
  )
  
  immune_pattern = paste(immune_keywords, collapse = '|')
  rnai_pattern   = paste(rnai_keywords, collapse = '|')
  
  master = master %>%
    mutate(
      functional_category = case_when(
        grepl(rnai_pattern, description, ignore.case = TRUE) ~ 'RNAi',
        grepl(immune_pattern, description, ignore.case = TRUE) ~ 'immune',
        annotated ~ 'other_annotated',
        TRUE ~ 'uncharacterised'
      )
    )
  
  cat('\nFunctional category breakdown:\n')
  print(table(master$functional_category))
  
  cat('\nImmune/RNAi genes significant per contrast:\n')
  genes_of_interest = master %>%
    filter(functional_category %in% c('immune', 'RNAi'))
  cat(sprintf('Total immune/RNAi genes: %d\n', nrow(genes_of_interest)))
  cat(sprintf('Sig in EHP vs NaCl: %d\n',
              sum(genes_of_interest$sig_EHP_vs_NaCl)))
  cat(sprintf('Sig in interaction: %d\n',
              sum(genes_of_interest$sig_interaction)))
  cat(sprintf('Sig in dsEHP vs EHP: %d\n',
              sum(genes_of_interest$sig_dsEHP_vs_EHP)))
  cat(sprintf('Sig in dsPTP2 vs NaCl: %d\n',
              sum(genes_of_interest$sig_dsPTP2_vs_NaCl)))
  cat(sprintf('Sig in dsEHP vs NaCl: %d\n',
              sum(genes_of_interest$sig_dsEHP_vs_NaCl)))
  
  # Priority scoring
  master = master %>%
    mutate(priority = case_when(
      # Tier 1: Immune/RNAi genes sig in BOTH interaction AND dsEHP vs EHP
      functional_category %in% c('immune', 'RNAi') &
        sig_interaction & sig_dsEHP_vs_EHP ~ 1L,
      # Tier 2: Immune/RNAi genes sig in EITHER interaction or dsEHP vs EHP
      functional_category %in% c('immune', 'RNAi') &
        (sig_interaction | sig_dsEHP_vs_EHP) ~ 2L,
      # Tier 3: Other annotated genes with convergent vaccine evidence
      annotated & sig_interaction & sig_dsEHP_vs_EHP ~ 3L,
      # Tier 4: Immune/RNAi genes sig in EHP vs NaCl or dsPTP2 vs NaCl
      functional_category %in% c('immune', 'RNAi') &
        (sig_EHP_vs_NaCl | sig_dsPTP2_vs_NaCl) ~ 4L,
      # Tier 5: Uncharacterised genes with convergent vaccine evidence
      functional_category == 'uncharacterised' &
        sig_interaction & sig_dsEHP_vs_EHP ~ 5L,
      # Tier 6: Uncharacterised genes sig in EITHER interaction or dsEHP vs EHP
      functional_category == 'uncharacterised' &
        (sig_interaction | sig_dsEHP_vs_EHP) ~ 6L,
      # Tier 7: Anything else significant in a priority contrast
      sig_any ~ 7L,
      # Tier 8: Not significant anywhere
      TRUE ~ 8L
    )) %>%
    arrange(priority, padj_interaction, padj_dsEHP_vs_EHP)
  
  cat('\nPriority tier distribution:\n')
  for (t in 1:8) {
    cat(sprintf('Tier %d: %d genes\n', t, sum(master$priority == t)))
  }
  
  # Priority subsets
  priority_genes = master %>% filter(priority < 8)
  
  unannotated_priority = master %>% filter(priority %in% c(5, 6))
  
  cat(sprintf('\n%d priority genes (tiers 1-7)\n', nrow(priority_genes)))
  cat(sprintf('\nHHPred candidates (tiers 5-6): %d genes\n',
              nrow(unannotated_priority)))
  cat(sprintf('With protein accession: %d\n',
              sum(!is.na(unannotated_priority$protein_accession))))
  cat(sprintf('Tier 5 (convergent): %d\n',
              sum(unannotated_priority$priority == 5)))
  cat(sprintf('Tier 6 (single analysis): %d\n',
              sum(unannotated_priority$priority == 6)))
  
  # Interaction sensitivity summary
  cat('\nInteraction sensitivity breakdown:\n')
  int_sens_summary = master %>%
    filter(!is.na(sens_class_interaction) & sens_class_interaction != 'NS') %>%
    dplyr::count(sens_class_interaction)
  print(as.data.frame(int_sens_summary))
  
  cat('\nInteraction sensitivity by priority tier (tiers 1-6):\n')
  int_sens_by_tier = master %>%
    filter(priority <= 6, !is.na(sens_class_interaction)) %>%
    dplyr::count(priority, sens_class_interaction) %>%
    pivot_wider(names_from = sens_class_interaction,
                values_from = n, values_fill = 0)
  print(as.data.frame(int_sens_by_tier))
  
  # GSEA ranked gene lists
  cat('\nPreparing GSEA ranked lists...\n')
  
  gsea_contrasts = config$gsea_stat_cols
  
  for (name in names(gsea_contrasts)) {
    stat_col = gsea_contrasts[[name]]
    
    ranked = master %>%
      select(gene, stat = all_of(stat_col)) %>%
      filter(!is.na(stat)) %>%
      arrange(desc(stat))
    
    write_tsv(ranked, file.path(gsea_dir, sprintf('gsea_ranks_%s.rnk', name)),
              col_names = FALSE)
    
    cat(sprintf('%s: %d genes ranked\n', name, nrow(ranked)))
  }
  
  # ORA gene sets
  cat('\nPreparing ORA gene sets...\n')
  
  writeLines(master$gene, file.path(ora_dir, 'background_all_tested.txt'))
  
  ora_contrasts = list(
    EHP_vs_NaCl    = list(padj = 'padj_EHP_vs_NaCl',    dir = 'dir_EHP_vs_NaCl'),
    dsEHP_vs_EHP   = list(padj = 'padj_dsEHP_vs_EHP',   dir = 'dir_dsEHP_vs_EHP'),
    dsPTP2_vs_NaCl = list(padj = 'padj_dsPTP2_vs_NaCl', dir = 'dir_dsPTP2_vs_NaCl'),
    interaction    = list(padj = 'padj_interaction',     dir = 'dir_interaction')
  )
  
  for (name in names(ora_contrasts)) {
    pcol = ora_contrasts[[name]]$padj
    dcol = ora_contrasts[[name]]$dir
    
    sig = master %>% filter(!is.na(.data[[pcol]]) & .data[[pcol]] < alpha)
    
    if (nrow(sig) == 0) {
      cat(sprintf('%s: no significant genes\n', name))
      next
    }
    
    writeLines(sig$gene, file.path(ora_dir, sprintf('ora_%s_all.txt', name)))
    
    sig_up = sig %>% filter(.data[[dcol]] == 'up')
    if (nrow(sig_up) > 0) {
      writeLines(sig_up$gene,
                 file.path(ora_dir, sprintf('ora_%s_up.txt', name)))
    }
    sig_down = sig %>% filter(.data[[dcol]] == 'down')
    if (nrow(sig_down) > 0) {
      writeLines(sig_down$gene,
                 file.path(ora_dir, sprintf('ora_%s_down.txt', name)))
    }
    
    cat(sprintf('%s: %d total (%d up, %d down)\n',
                name, nrow(sig), nrow(sig_up), nrow(sig_down)))
  }
  
  # dsEHP_1-dependent genes — separate exploratory sets from each sensitivity
  if (!is.null(master$sens_class_pairwise)) {
    dep_pw = master %>%
      filter(sens_class_pairwise == 'dsEHP_1_dependent') %>%
      pull(gene)
    if (length(dep_pw) > 0) {
      writeLines(dep_pw,
                 file.path(ora_dir, 'ora_dsEHP1_dependent_pairwise.txt'))
      cat(sprintf('dsEHP_1-dependent (pairwise): %d genes\n', length(dep_pw)))
    }
  }
  
  if (!is.null(master$sens_class_interaction)) {
    dep_int = master %>%
      filter(sens_class_interaction == 'dsEHP_1_dependent') %>%
      pull(gene)
    if (length(dep_int) > 0) {
      writeLines(dep_int,
                 file.path(ora_dir, 'ora_dsEHP1_dependent_interaction.txt'))
      cat(sprintf('dsEHP_1-dependent (interaction): %d genes\n',
                  length(dep_int)))
    }
  }
  
  # Robust dsEHP vs EHP — highest-confidence pairwise vaccine set
  if (!is.null(master$sens_class_pairwise)) {
    robust_pw = master %>%
      filter(sens_class_pairwise == 'robust', sig_dsEHP_vs_EHP) %>%
      pull(gene)
    if (length(robust_pw) > 0) {
      writeLines(robust_pw,
                 file.path(ora_dir, 'ora_dsEHP_vs_EHP_robust.txt'))
      cat(sprintf('dsEHP vs EHP robust: %d genes\n', length(robust_pw)))
    }
  }
  
  # KEGG/Entrez ID mapping
  kegg_map = annotation_df %>%
    filter(gene_id %in% master$gene, !is.na(ncbi_gene_id)) %>%
    select(gene_id, ncbi_gene_id) %>%
    distinct()
  
  suppressMessages(write_csv(kegg_map,
                             file.path(res_dir, 'gene_to_entrez_map.csv')))
  
  mapping_rate = nrow(kegg_map) / nrow(master) * 100
  cat(sprintf('\nKEGG mapping: %d of %d genes have Entrez IDs (%.1f%%)\n',
              nrow(kegg_map), nrow(master), mapping_rate))
  
  if (mapping_rate < 50) {
    cat('WARNING: Low KEGG mapping rate. Enrichment may be underpowered.\n')
  }
  
  # dsEHP_1-dependent gene table
  dependent_genes = master %>%
    filter(sens_class_interaction =='dsEHP_1_dependent' | 
             sens_class_pairwise == 'dsEHP_1_dependent') %>%
    arrange(sens_class_interaction, sens_class_pairwise, padj_dsEHP_vs_EHP)
  
  suppressMessages(write_csv(dependent_genes,
                             file.path(res_dir, 'dsEHP_dependent_genes.csv')))
  cat(sprintf('\ndsEHP_1_dependent genes: %d genes exported\n',
              nrow(dependent_genes)))
  
  # Save outputs
  saveRDS(master, file.path(proc_dir, 'master_gene_table.rds'))
  suppressMessages(write_csv(master,
                             file.path(res_dir, 'master_gene_table.csv')))
  suppressMessages(write_csv(priority_genes,
                             file.path(res_dir, 'priority_genes.csv')))
  suppressMessages(write_csv(unannotated_priority,
                             file.path(res_dir, 'unannotated_priority.csv')))
  
  # Summary
  cat('\n========== Collation summary ==========\n')
  cat(sprintf('Total genes: %d\n', nrow(master)))
  cat(sprintf('Significant (any, excl LRT): %d\n', sum(master$sig_any)))
  cat(sprintf('Priority genes (T1-7): %d\n', nrow(priority_genes)))
  cat(sprintf('T1 (immune + convergent): %d\n', sum(master$priority == 1)))
  cat(sprintf('T2 (immune + vaccine): %d\n', sum(master$priority == 2)))
  cat(sprintf('T3 (annotated + converg): %d\n', sum(master$priority == 3)))
  cat(sprintf('T4 (immune + baseline): %d\n', sum(master$priority == 4)))
  cat(sprintf('T5 (unchara + converg): %d\n', sum(master$priority == 5)))
  cat(sprintf('T6 (unchara + vaccine): %d\n', sum(master$priority == 6)))
  cat(sprintf('T7 (other significant): %d\n', sum(master$priority == 7)))
  cat(sprintf('HHPred candidates (T5-6): %d\n', nrow(unannotated_priority)))
  cat(sprintf('GSEA ranked lists: %s\n', gsea_dir))
  cat(sprintf('ORA gene sets: %s\n', ora_dir))
  cat(sprintf('Outputs: %s\n', res_dir))
}

main()
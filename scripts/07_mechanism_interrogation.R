# NOTE: DUE TO TIME CONSTRAINTS, THIS SCRIPT WAS WRITTEN BY AI TO MY SPECIFICATIONS 
# WITH HUMAN EDITING
#
# scripts/07_mechanism_interrogation.R
# Deep interrogation of interaction mechanism classes and dsEHP_1 biology
#
# Section 1: Non-additive gene class characterisation (emergent, vaccine-modified,
#            dsRNA-primed)
# Section 2: dsEHP_1-dependent genes within each mechanism class —
#            what distinguishes the extreme responder?
# Section 3: RNAi machinery expression profiles across all four groups —
#            does the silencing apparatus persist through challenge?

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/annotation_functions.R')
source('R/plotting.R')

main = function() {
  log_con = start_log('07_mechanism_interrogation')
  on.exit(end_log(log_con))
  
  proc_dir = config$proc_dir
  alpha    = config$alpha
  res_dir  = 'results/07_mechanism_interrogation'
  table_dir = file.path(res_dir, 'tables')
  plot_dir  = file.path(res_dir, 'plots')
  gene_dir  = file.path(plot_dir, 'gene_plots')
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_dpi    = config$plot_dpi
  plot_width  = config$plot_width
  plot_height = config$plot_height
  
  # ── Load data ──────────────────────────────────────────────────────────
  master = readRDS(file.path(proc_dir, 'master_gene_table.rds'))
  annotation_df = readRDS(file.path(proc_dir, 'annotation_table.rds'))
  gene_labels = build_gene_labels(annotation_df)
  
  # Multifactor sensitivity model objects (includes dsEHP_1 for visualising)
  dds_mf = readRDS(file.path(proc_dir, 'dds_multifactor_de_sensitivity_fitted.rds'))
  vsd_mf = readRDS(file.path(proc_dir, 'vsd_multifactor_de_sensitivity.rds'))
  coldata_mf = readRDS(file.path(proc_dir, 'coldata_multifactor_sens.rds'))
  vst_mf = assay(vsd_mf)
  
  # Mechanism classification from script 06
  interaction_crossref = read_csv(file.path(
    'results/06_candidate_pathways/tables',
    'interaction_mechanism_crossref.csv'), show_col_types = FALSE)
  
  cat('All data loaded\n')
  cat(sprintf('Mechanism classification: %d genes\n',
              nrow(interaction_crossref)))
  
  # ════════════════════════════════════════════════════════════════════════
  # SECTION 1: Non-additive gene class characterisation
  # ════════════════════════════════════════════════════════════════════════
  cat('\n══ Section 1: Non-additive gene class characterisation ══\n')
  
  # The three non-additive classes: genes where the combined effect of
  # vaccination + challenge deviates from additivity in ways that can't
  # be explained by simple pathway saturation
  
  non_additive_classes = c(
    'co-directional (sub-additive)',
    'co-directional (synergistic)',
    'opposing main effects',
    'interaction-only (emergent)',
    'infection-responsive, vaccine-amplified',
    'infection-responsive, vaccine-attenuated',
    'dsRNA-primed (synergistic)',
    'dsRNA-primed (sub-additive)'
  )
  
  for (mech_class in non_additive_classes) {
    # Pull full gene details from master via the crossref
    class_genes = interaction_crossref %>%
      filter(dsRNA_response == mech_class) %>%
      arrange(padj_interaction)
    
    class_gene_ids = class_genes$gene
    
    # Short label for filenames
    short_label = gsub('[^a-zA-Z0-9]', '_', mech_class)
    short_label = gsub('_+', '_', short_label)
    
    cat(sprintf('\n── %s: %d genes ──\n', mech_class, length(class_gene_ids)))
    
    # Functional composition
    func_summary = class_genes %>%
      dplyr::count(functional_category) %>%
      arrange(desc(n))
    cat('Functional composition:\n')
    print(as.data.frame(func_summary))
    
    # Direction summary
    cat(sprintf('Direction: %d up, %d down (interaction LFC)\n',
                sum(class_genes$lfc_interaction > 0),
                sum(class_genes$lfc_interaction < 0)))
    cat(sprintf('Median |interaction LFC|: %.2f\n',
                median(abs(class_genes$lfc_interaction))))
    
    # Sensitivity summary
    sens_counts = class_genes %>%
      dplyr::count(sens_class_interaction) %>%
      arrange(desc(n))
    cat('dsEHP_1 sensitivity:\n')
    print(as.data.frame(sens_counts))
    
    # Log all annotated genes (these are small enough to print in full)
    annotated_genes = class_genes %>%
      filter(annotated == TRUE)
    cat(sprintf('\n%d annotated genes (full list):\n', nrow(annotated_genes)))
    for (i in seq_len(nrow(annotated_genes))) {
      sens = ifelse(is.na(annotated_genes$sens_class_interaction[i]), '?',
                    annotated_genes$sens_class_interaction[i])
      cat(sprintf('  %s: %s\n    int=%+.2f±%.2f  dsPTP2=%+.2f±%.2f  EHP=%+.2f±%.2f\n    dsEHP-NaCl=%+.2f±%.2f  dsEHP-EHP=%+.2f±%.2f  [%s]\n'
                  , annotated_genes$gene[i]
                  , substr(annotated_genes$description[i], 1, 55)
                  , annotated_genes$lfc_interaction[i]
                  , annotated_genes$lfcSE_interaction[i]
                  , annotated_genes$lfc_dsPTP2_vs_NaCl[i]
                  , annotated_genes$lfcSE_dsPTP2_vs_NaCl[i]
                  , annotated_genes$lfc_EHP_vs_NaCl[i]
                  , annotated_genes$lfcSE_EHP_vs_NaCl[i]
                  , annotated_genes$lfc_dsEHP_vs_NaCl[i]
                  , annotated_genes$lfcSE_dsEHP_vs_NaCl[i]
                  , annotated_genes$lfc_dsEHP_vs_EHP[i]
                  , annotated_genes$lfcSE_dsEHP_vs_EHP[i]
                  , sens))
    }
    
    # Log uncharacterised genes (gene ID + LFCs for HHPred cross-reference)
    unchar = class_genes %>%
      filter(functional_category == 'uncharacterised')
    if (nrow(unchar) > 0) {
      cat(sprintf('\n%d uncharacterised genes:\n', nrow(unchar)))
      for (i in seq_len(nrow(unchar))) {
        sens = ifelse(is.na(unchar$sens_class_interaction[i]), '?',
                      unchar$sens_class_interaction[i])
        has_acc = ifelse(!is.na(unchar$protein_accession[i]),
                         'has accession', 'no accession')
        cat(sprintf('%s: int=%+.2f±%.2f  dsEHP-NaCl=%+.2f±%.2f  dsEHP-EHP=%+.2f±%.2f  [%s, %s]\n',
                    unchar$gene[i],
                    unchar$lfc_interaction[i],
                    unchar$lfcSE_interaction[i],
                    unchar$lfc_dsEHP_vs_NaCl[i],
                    unchar$lfcSE_dsEHP_vs_NaCl[i],
                    unchar$lfc_dsEHP_vs_EHP[i],
                    unchar$lfcSE_dsEHP_vs_EHP[i],
                    sens, has_acc))
      }
    }
    
    # Export full class table for manual inspection
    suppressMessages(write_csv(
      class_genes,
      file.path(table_dir, sprintf('mechanism_%s.csv', short_label))))
    
    # Heatmap (multifactor VST, includes dsEHP_1)
    if (length(class_gene_ids) >= 3) {
      row_anno = master %>%
        filter(gene %in% class_gene_ids) %>%
        transmute(
          gene,
          Interaction = ifelse(sig_interaction, 'sig', 'NS'),
          `EHP vs NaCl` = ifelse(sig_EHP_vs_NaCl, 'sig', 'NS'),
          `dsPTP2 vs NaCl` = ifelse(sig_dsPTP2_vs_NaCl, 'sig', 'NS'),
          dsEHP_1_sens = ifelse(
            is.na(sens_class_interaction), 'NS',
            sens_class_interaction)
        ) %>%
        column_to_rownames('gene') %>%
        as.data.frame()
      
      plot_pathway_heatmap(
        gene_ids = class_gene_ids,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        row_annotations = row_anno,
        filename = file.path(plot_dir,
                             sprintf('heatmap_%s.png', short_label)),
        title = sprintf('%s (%d genes)', mech_class, length(class_gene_ids))
      )
      cat(sprintf('Heatmap saved: %s\n', short_label))
    }
  }
  
  # Save full class gene table
  non_additive_genes = interaction_crossref %>%
    filter(dsRNA_response %in% non_additive_classes)
  suppressMessages(write_csv(
    non_additive_genes,
    file.path(table_dir, 'non_additive_genes.csv')))
  
  # ════════════════════════════════════════════════════════════════════════
  # SECTION 2: dsEHP_1-dependent genes within mechanism classes
  # ════════════════════════════════════════════════════════════════════════
  cat('\n\n══ Section 2: dsEHP_1-dependent genes by mechanism ══\n')
  
  # dsEHP_1 has near-complete parasite clearance (35 copies/µl vs 1400-2200
  # for groupmates). Genes that are dsEHP_1-dependent in the interaction
  # are genes where this extreme responder drives the signal.
  # This is an n=1 observation — hypothesis-generating, not conclusive.
  
  for (mech_class in non_additive_classes) {
    class_genes = interaction_crossref %>%
      filter(dsRNA_response == mech_class) %>%
      arrange(padj_interaction)
    
    
    dep = class_genes %>%
      filter(sens_class_interaction == 'dsEHP_1_dependent')
    robust = class_genes %>%
      filter(sens_class_interaction == 'robust')
    
    short_label = gsub('[^a-zA-Z0-9]', '_', mech_class)
    short_label = gsub('_+', '_', short_label)
    
    cat(sprintf('\n── %s ──\n', mech_class))
    cat(sprintf('  Total: %d | Robust: %d | dsEHP_1-dependent: %d\n',
                nrow(class_genes), nrow(robust), nrow(dep)))
    
    if (nrow(dep) == 0) {
      cat('  No dsEHP_1-dependent genes in this class\n')
      next
    }
    
    # What fraction is dsEHP_1-dependent?
    dep_frac = nrow(dep) / nrow(class_genes)
    cat(sprintf('  dsEHP_1-dependent fraction: %.1f%%\n', 100 * dep_frac))
    
    # Functional composition of dependent genes
    dep_func = dep %>%
      dplyr::count(functional_category) %>%
      arrange(desc(n))
    cat('  Functional composition (dsEHP_1-dependent):\n')
    for (i in seq_len(nrow(dep_func))) {
      cat(sprintf('    %s: %d\n', dep_func$functional_category[i],
                  dep_func$n[i]))
    }
    
    # Log annotated dependent genes with full LFC context
    dep_annotated = dep %>%
      filter(annotated == TRUE)
    if (nrow(dep_annotated) > 0) {
      cat('  Annotated dsEHP_1-dependent genes:\n')
      for (i in seq_len(nrow(dep_annotated))) {
        cat(sprintf('%s: %s\n int=%+.2f±%.2f  dsPTP2=%+.2f±%.2f  EHP=%+.2f±%.2f/n   dsEHP-NaCl=%+.2f±%.2f  dsEHP-EHP=%+.2f±%.2f\n'
                    , dep_annotated$gene[i]
                    , substr(dep_annotated$description[i], 1, 55)
                    , dep_annotated$lfc_interaction[i]
                    , dep_annotated$lfcSE_interaction[i]
                    , dep_annotated$lfc_dsPTP2_vs_NaCl[i]
                    , dep_annotated$lfcSE_dsPTP2_vs_NaCl[i]
                    , dep_annotated$lfc_EHP_vs_NaCl[i]
                    , dep_annotated$lfcSE_EHP_vs_NaCl[i]
                    , dep_annotated$lfc_dsEHP_vs_NaCl[i]
                    , dep_annotated$lfcSE_dsEHP_vs_NaCl[i]
                    , dep_annotated$lfc_dsEHP_vs_EHP[i]
                    , dep_annotated$lfcSE_dsEHP_vs_EHP[i]))
      }
    }
    
    # Log uncharacterised dependent genes
    dep_unchar = dep %>%
      filter(functional_category == 'uncharacterised')
    if (nrow(dep_unchar) > 0) {
      cat(sprintf('%d uncharacterised dsEHP_1-dependent genes:\n',
                  nrow(dep_unchar)))
      for (i in seq_len(nrow(dep_unchar))) {
        cat(sprintf('%s: int=%+.2f±%.2f  dsEHP-NaCl=%+.2f±%.2f  dsEHP-EHP=%+.2f±%.2f\n',
                    dep_unchar$gene[i],
                    dep_unchar$lfc_interaction[i],
                    dep_unchar$lfcSE_interaction[i],
                    dep_unchar$lfc_dsEHP_vs_NaCl[i],
                    dep_unchar$lfcSE_dsEHP_vs_NaCl[i],
                    dep_unchar$lfc_dsEHP_vs_EHP[i],
                    dep_unchar$lfcSE_dsEHP_vs_EHP[i]))
      }
    }
    
    # Compare direction of interaction LFC: dependent vs robust
    if (nrow(robust) > 0 && nrow(dep) > 0) {
      cat('\n  Direction comparison (interaction LFC):\n')
      cat(sprintf('    Robust - up: %d, down: %d\n',
                  sum(robust$lfc_interaction > 0),
                  sum(robust$lfc_interaction < 0)))
      cat(sprintf('    dsEHP_1-dep - up: %d, down: %d\n',
                  sum(dep$lfc_interaction > 0),
                  sum(dep$lfc_interaction < 0)))
      cat(sprintf('    Robust median |LFC|: %.2f\n',
                  median(abs(robust$lfc_interaction))))
      cat(sprintf('    dsEHP_1-dep median |LFC|: %.2f\n',
                  median(abs(dep$lfc_interaction))))
    }
    
    # Gene-level plots for dsEHP_1-dependent genes (capped at 15)
    plot_genes = dep %>%
      arrange(padj_interaction) %>%
      slice_head(n = 15)
    
    for (i in seq_len(nrow(plot_genes))) {
      gene_id = plot_genes$gene[i]
      gene_display = resolve_labels(gene_id, gene_labels)
      desc = ifelse(is.na(plot_genes$description[i]), 'uncharacterised',
                    substr(plot_genes$description[i], 1, 40))
      
      p = plot_gene_counts(
        gene_id, dds_mf, coldata_mf,
        highlight_sample = 'dsEHP_1',
        gene_labels = gene_labels,
        title = sprintf('dsEHP_1-dep [%s]\n%s (int padj=%.2g)',
                        mech_class, gene_display,
                        plot_genes$padj_interaction[i])
      )
      
      if (!is.null(p)) {
        quiet_ggsave(file.path(gene_dir,
                               sprintf('dsEHP1dep_%s_%02d_%s.png',
                                       short_label, i, gene_id)),
                     p, width = plot_width, height = plot_height,
                     dpi = plot_dpi)
      }
    }
    
    # Heatmap: dsEHP_1-dependent genes, with dsEHP_1 visible
    if (nrow(dep) >= 3) {
      plot_pathway_heatmap(
        gene_ids = dep$gene,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        filename = file.path(plot_dir,
                             sprintf('heatmap_dsEHP1dep_%s.png',
                                     short_label)),
        title = sprintf('dsEHP_1-dependent: %s (%d genes)',
                        mech_class, nrow(dep))
      )
    }
  }
  
  # dsEHP vs EHP exclusive genes
  # Investigating additive effects of vaccination which persist through challenge
  # Not AI - added post-hoc to fill gap in analysis
  cat('\ndsEHP vs EHP exlusive analysis\n')
  
  # Overlap categories & counts
  both_sig = sum(master$sig_dsEHP_vs_EHP & master$sig_interaction)
  pw_only = sum(master$sig_dsEHP_vs_EHP & !master$sig_interaction)
  int_only = sum(!master$sig_dsEHP_vs_EHP & master$sig_interaction)
  
  cat(sprintf('dsEHP vs EHP only: %d\n', pw_only))
  cat(sprintf('Interaction only: %d\n', int_only))
  cat(sprintf('Both: %d\n', both_sig))
  
  # Extract pairwise-only genes
  pw_exclusive = master %>%
    filter(sig_dsEHP_vs_EHP & !sig_interaction) %>%
    arrange(padj_dsEHP_vs_EHP)
  
  # Print direction and functional category
  pw_excl_functions = pw_exclusive %>%
    dplyr::count(functional_category) %>%
    arrange(desc(n))
  cat('Functional composition:\n')
  print(as.data.frame(pw_excl_functions))
  
  cat(sprintf('Direction: %d up, %d down\n',
              sum(pw_exclusive$dir_dsEHP_vs_EHP == 'up'),
              sum(pw_exclusive$dir_dsEHP_vs_EHP == 'down')))
  
  # Check that significance and direction match dsPTP2 vs NaCl
  pw_excl_also_dsPTP2 = sum(pw_exclusive$sig_dsPTP2_vs_NaCl)
  cat(sprintf('\nAlso significant in dsPTP2 vs NaCl: %d of %d (%.1f%%)\n',
              pw_excl_also_dsPTP2, nrow(pw_exclusive),
              100 * pw_excl_also_dsPTP2 / nrow(pw_exclusive)))
  
  if (pw_excl_also_dsPTP2 > 0) {
    concordant = pw_exclusive %>%
      filter(sig_dsPTP2_vs_NaCl) %>%
      mutate(concordant = sign(lfc_dsEHP_vs_EHP) == sign(lfc_dsPTP2_vs_NaCl))
    cat(sprintf('Direction concordant with dsPTP2: %d of %d\n',
                sum(concordant$concordant), nrow(concordant)))
  }
  
  # Analyse immune/RNAi genes in pairwise-only set
  pw_excl_immune = filter(pw_exclusive, functional_category %in% c('immune', 'RNAi'))
  
  cat(sprintf('\n%d immune/RNAi genes in dsEHP vs EHP-exclusive set:\n',
              nrow(pw_excl_immune)))
  
  if (nrow(pw_excl_immune) > 0) {
    for (i in seq_len(nrow(pw_excl_immune))) {
      also_ptp2 = ifelse(pw_excl_immune$sig_dsPTP2_vs_NaCl[i], 'sig', 'NS')
      also_ehp = ifelse(pw_excl_immune$sig_EHP_vs_NaCl[i], 'sig', 'NS')
      cat(sprintf('%s: %s\n  dsEHP-EHP=%+.2f±%.2f | dsPTP2=%+.2f±%.2f |\n  EHP=%+.2f±%.2f | dsEHP-NaCl=%+.2f±%.2f\n',
                  pw_excl_immune$gene[i],
                  substr(pw_excl_immune$description[i], 1, 50),
                  pw_excl_immune$lfc_dsEHP_vs_EHP[i],
                  pw_excl_immune$lfcSE_dsEHP_vs_EHP[i],
                  pw_excl_immune$lfc_dsPTP2_vs_NaCl[i],
                  pw_excl_immune$lfcSE_dsPTP2_vs_NaCl[i],
                  pw_excl_immune$lfc_EHP_vs_NaCl[i],
                  pw_excl_immune$lfcSE_EHP_vs_NaCl[i],
                  pw_excl_immune$lfc_dsEHP_vs_NaCl[i],
                  pw_excl_immune$lfcSE_dsEHP_vs_NaCl[i]))
    }
  }
  
  # Sensitivity classification of pairwise-only genes
  if ('sens_class_pairwise' %in% names(pw_exclusive)) {
    pw_excl_sens = pw_exclusive %>%
      filter(!is.na(sens_class_pairwise) & sens_class_pairwise != 'NS') %>%
      dplyr::count(sens_class_pairwise)
    cat('\ndsEHP vs EHP sensitivity (pairwise model):\n')
    print(as.data.frame(pw_excl_sens))
  }
  
  suppressMessages(write_csv(
    pw_exclusive,
    file.path(table_dir, 'dsEHP_vs_EHP_exclusive_genes.csv')))
  
  cat(sprintf('\nExported to %s\n',
              file.path(table_dir, 'dsEHP_vs_EHP_exclusive_genes.csv')))
  
  # ════════════════════════════════════════════════════════════════════════
  # SECTION 3: RNAi machinery expression across all four groups
  # ════════════════════════════════════════════════════════════════════════
  cat('\n\n══ Section 3: RNAi machinery expression profiles ══\n')
  
  # Search for all RNAi-related genes (broader than the immune keyword list,
  # which was tightened to core components — here we want the full picture)
  rnai_pattern = paste(c(
    'Dicer', 'Argonaute', 'AGO', 'R2D2', 'SID',
    'Drosha', 'Pasha', 'DGCR8',
    'RNA.dependent RNA polymerase', 'RdRp',
    'RISC', 'RNA interference',
    'TRBP', 'eIF6', 'exportin.5',  # RNAi accessory components
    'Staufen', 'Tudor'              # dsRNA-binding / piRNA-adjacent
  ), collapse = '|')
  
  rnai_genes = master %>%
    filter(grepl(rnai_pattern, description, ignore.case = TRUE)) %>%
    select(gene, description, baseMean,
           lfc_dsPTP2_vs_NaCl, lfcSE_dsPTP2_vs_NaCl, padj_dsPTP2_vs_NaCl,
           lfc_EHP_vs_NaCl, lfcSE_EHP_vs_NaCl, padj_EHP_vs_NaCl,
           lfc_interaction, lfcSE_interaction, padj_interaction,
           lfc_dsEHP_vs_NaCl, lfcSE_dsEHP_vs_NaCl, padj_dsEHP_vs_NaCl,
           lfc_dsEHP_vs_EHP, lfcSE_dsEHP_vs_EHP, padj_dsEHP_vs_EHP,
           sig_dsPTP2_vs_NaCl, sig_EHP_vs_NaCl,
           sig_interaction, sig_dsEHP_vs_EHP)
  
  cat(sprintf('\n%d RNAi or RNAi-related genes found\n', nrow(rnai_genes)))
  
  if (nrow(rnai_genes) > 0) {
    # Compute group mean expression
    rnai_genes = rnai_genes %>%
      mutate(
        # Does dsRNA upregulation persist through challenge?
        dsRNA_persists = abs(lfc_dsEHP_vs_EHP) > 0.5 &
          sign(lfc_dsEHP_vs_EHP) == sign(lfc_dsPTP2_vs_NaCl)
      )
    
    # Log each gene with its expression across contrasts
    cat('\nRNAi gene expression profile:\n')
    cat(sprintf('%-16s %-45s %8s %8s %8s %8s %8s %8s\n',
                'Gene', 'Description', 'dsPTP2', 'EHP', 'Int',
                'dsE-NaCl', 'dsE-E', 'Persist'))
    
    for (i in seq_len(nrow(rnai_genes))) {
      persist_flag = ifelse(rnai_genes$dsRNA_persists[i], 'YES', 'no')
      
      cat(sprintf('%-16s %-45s %+.2f±%.2f%s %+.2f±%.2f%s %+.2f±%.2f%s %+.2f±%.2f %+.2f±%.2f  %s\n',
                  rnai_genes$gene[i],
                  substr(rnai_genes$description[i], 1, 45),
                  rnai_genes$lfc_dsPTP2_vs_NaCl[i],
                  rnai_genes$lfcSE_dsPTP2_vs_NaCl[i],
                  ifelse(rnai_genes$sig_dsPTP2_vs_NaCl[i], '*', ' '),
                  rnai_genes$lfc_EHP_vs_NaCl[i],
                  rnai_genes$lfcSE_EHP_vs_NaCl[i],
                  ifelse(rnai_genes$sig_EHP_vs_NaCl[i], '*', ' '),
                  rnai_genes$lfc_interaction[i],
                  rnai_genes$lfcSE_interaction[i],
                  ifelse(rnai_genes$sig_interaction[i], '*', ' '),
                  rnai_genes$lfc_dsEHP_vs_NaCl[i],
                  rnai_genes$lfcSE_dsEHP_vs_NaCl[i],
                  rnai_genes$lfc_dsEHP_vs_EHP[i],
                  rnai_genes$lfcSE_dsEHP_vs_EHP[i],
                  persist_flag))
    }
    }
    
    # Summary: does the RNAi machinery persist through challenge?
    sig_dsrna = rnai_genes %>% filter(sig_dsPTP2_vs_NaCl)
    if (nrow(sig_dsrna) > 0) {
      cat(sprintf('\n%d RNAi genes significantly respond to dsRNA alone\n',
                  nrow(sig_dsrna)))
      cat(sprintf('  Of these, %d retain >0.5 LFC in same direction '
                  , sum(sig_dsrna$dsRNA_persists)))
      cat('during challenge (dsEHP-EHP)\n')
      cat(sprintf('  Median |dsPTP2 LFC|: %.2f (median SE: %.2f)\n',
                  median(abs(sig_dsrna$lfc_dsPTP2_vs_NaCl)),
                  median(sig_dsrna$lfcSE_dsPTP2_vs_NaCl)))
      cat(sprintf('  Median |dsEHP-NaCl|: %.2f (median SE: %.2f)\n',
                  median(abs(sig_dsrna$lfc_dsEHP_vs_NaCl)),
                  median(sig_dsrna$lfcSE_dsEHP_vs_NaCl)))
      cat(sprintf('  Median |dsEHP-EHP|: %.2f (median SE: %.2f)\n',
                  median(abs(sig_dsrna$lfc_dsEHP_vs_EHP)),
                  median(sig_dsrna$lfcSE_dsEHP_vs_EHP)))
      cat(sprintf('  Median survival fraction: %.2f\n',
                  median(abs(sig_dsrna$lfc_dsEHP_vs_EHP) /
                           abs(sig_dsrna$lfc_dsPTP2_vs_NaCl))))
    }
    
    # Save table
    suppressMessages(write_csv(
      rnai_genes,
      file.path(table_dir, 'rnai_machinery_expression.csv')))
    
    # Heatmap of RNAi genes
    rnai_in_vst = rnai_genes$gene[rnai_genes$gene %in% rownames(vst_mf)]
    if (length(rnai_in_vst) >= 3) {
      plot_pathway_heatmap(
        gene_ids = rnai_in_vst,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        filename = file.path(plot_dir, 'heatmap_rnai_machinery.png'),
        title = sprintf('RNAi machinery (%d genes)', length(rnai_in_vst))
      )
    }
    
    # Gene-level plots for dsRNA-responsive RNAi genes
    rnai_to_plot = rnai_genes %>%
      filter(sig_dsPTP2_vs_NaCl | sig_interaction) %>%
      arrange(padj_dsPTP2_vs_NaCl)
    
    for (i in seq_len(nrow(rnai_to_plot))) {
      gene_id = rnai_to_plot$gene[i]
      gene_display = resolve_labels(gene_id, gene_labels)
      persist = ifelse(rnai_to_plot$dsRNA_persists[i],
                       'persists', 'collapses')
      
      p = plot_gene_counts(
        gene_id, dds_mf, coldata_mf,
        highlight_sample = 'dsEHP_1',
        gene_labels = gene_labels,
        title = sprintf('RNAi: %s\ndsPTP2=%.2f, dsEHP-EHP=%.2f (%s)',
                        gene_display,
                        rnai_to_plot$lfc_dsPTP2_vs_NaCl[i],
                        rnai_to_plot$lfc_dsEHP_vs_EHP[i],
                        persist)
      )
      
      if (!is.null(p)) {
        quiet_ggsave(file.path(gene_dir,
                               sprintf('rnai_%02d_%s.png', i, gene_id)),
                     p, width = plot_width, height = plot_height,
                     dpi = plot_dpi)
      }
    }
  
  # ════════════════════════════════════════════════════════════════════════
  # SUMMARY
  # ════════════════════════════════════════════════════════════════════════
  cat('\n\n══ Summary ══\n')
  
  for (mech_class in non_additive_classes) {
    class_data = interaction_crossref %>%
      filter(dsRNA_response == mech_class)
    n_total = nrow(class_data)
    n_robust = sum(class_data$sens_class_interaction == 'robust',
                   na.rm = TRUE)
    n_dep = sum(class_data$sens_class_interaction == 'dsEHP_1_dependent',
                na.rm = TRUE)
    n_immune = sum(class_data$functional_category %in% c('immune', 'RNAi'))
    
    cat(sprintf('  %s: %d total (%d robust, %d dsEHP_1-dep, %d immune)\n',
                mech_class, n_total, n_robust, n_dep, n_immune))
  }
  
  cat(sprintf('\nRNAi genes responsive to dsRNA: %d\n',
              sum(rnai_genes$sig_dsPTP2_vs_NaCl)))
  cat(sprintf('RNAi genes with persistent effect: %d\n',
              sum(rnai_genes$sig_dsPTP2_vs_NaCl & rnai_genes$dsRNA_persists)))
  
  cat(sprintf('\nResults: %s\n', res_dir))
}

main()
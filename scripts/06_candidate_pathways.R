# scripts/06_candidate_pathways.R
# Thorough exploration of enriched pathways and immune genes
# Pathways enriched in KEGG & GO and immune pathways
# Search for dsRNA-triggered, challenge-triggered, and interaction-triggered

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/annotation_functions.R')
source('R/plotting.R')
source('R/enrichment_functions.R')

main = function() {
  log_con = start_log('06_candidate_pathways')
  on.exit(end_log(log_con))
  
  proc_dir   = config$proc_dir
  alpha      = config$alpha
  res_dir    = 'results/06_candidate_pathways'
  heatmap_dir = file.path(res_dir, 'pathway_heatmaps')
  gene_dir   = file.path(res_dir, 'gene_plots')
  table_dir  = file.path(res_dir, 'tables')
  dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_dpi    = config$plot_dpi
  plot_width  = config$plot_width
  plot_height = config$plot_height
  
  # Load data
  master = readRDS(file.path(proc_dir, 'master_gene_table.rds'))
  enrichment = readRDS(file.path(proc_dir, 'enrichment_results.rds'))
  annotation_df = readRDS(file.path(proc_dir, 'annotation_table.rds'))
  gene_labels = build_gene_labels(annotation_df)
  
  vsd_primary = readRDS(file.path(proc_dir, 'vsd_primary_de.rds'))
  coldata_primary = readRDS(file.path(proc_dir, 'coldata_primary.rds'))
  vst_primary = assay(vsd_primary)
  
  dds_mf = readRDS(file.path(proc_dir, 'dds_multifactor_de_sensitivity_fitted.rds'))
  vsd_mf = readRDS(file.path(proc_dir, 'vsd_multifactor_de_sensitivity.rds'))
  coldata_mf = readRDS(file.path(proc_dir, 'coldata_multifactor_sens.rds'))
  vst_mf = assay(vsd_mf)
  
  # Identify significant pathways
  cat('\nIdentifying pathways for investigation\n')
  
  # KEGG pathways
  kegg_results = enrichment$gsea[grep(
    '_KEGG$', names(enrichment$gsea), value = TRUE)]
  
  sig_kegg = NULL
  if (length(kegg_results) > 0) {
    sig_kegg = bind_rows(
      lapply(names(kegg_results), function(rname) {
        kegg_results[[rname]] %>%
          filter(padj < alpha) %>%
          mutate(contrast = sub('_KEGG$', '', rname),
                 source = 'KEGG')
      })
    )
  }
  
  if (is.null(sig_kegg) || nrow(sig_kegg) == 0) {
    cat('No significant KEGG pathways from GSEA\n')
    kegg_pathway_list = list()
    kegg_summary = NULL
  } else {
    kegg_summary = sig_kegg %>%
      group_by(pathway) %>%
      summarise(
        n_contrasts = n_distinct(contrast),
        mean_abs_NES = mean(abs(NES)),
        best_padj = min(padj),
        contrasts = paste(contrast, collapse = ', '),
        .groups = 'drop'
      ) %>%
      arrange(desc(n_contrasts), best_padj)
    
    cat(sprintf('%d unique significant KEGG pathways\n', nrow(kegg_summary)))
    print(as.data.frame(head(kegg_summary, 10)))
    
    suppressMessages(write_csv(kegg_summary,
                               file.path(table_dir,
                                         'kegg_pathway_investigation_list.csv')))
    
    kegg_pathway_list = enrichment$gene_sets$kegg[
      names(enrichment$gene_sets$kegg) %in% kegg_summary$pathway
    ]
  }
  
  # GO pathways
  go_results = enrichment$gsea[grep('_GO_', names(enrichment$gsea), value = TRUE)]
  
  sig_go = NULL
  go_summary = NULL
  if (length(go_results) > 0) {
    sig_go = bind_rows(
      lapply(names(go_results), function(rname) {
        go_results[[rname]] %>%
          filter(padj < alpha) %>%
          mutate(
            contrast = sub('_GO_.*$', '', rname),
            ontology = sub('^.*_GO_', '', rname),
            source = 'GO'
          )
      })
    )
    
    if (nrow(sig_go) > 0) {
      go_summary = sig_go %>%
        group_by(pathway, ontology) %>%
        summarise(
          n_contrasts = n_distinct(contrast),
          mean_abs_NES = mean(abs(NES)),
          best_padj = min(padj),
          contrasts = paste(contrast, collapse = ', '),
          .groups = 'drop'
        ) %>%
        arrange(ontology, desc(n_contrasts), best_padj)
      
      cat(sprintf('\n%d unique significant GO terms across all contrasts\n',
                  nrow(go_summary)))
      
      for (ont in unique(go_summary$ontology)) {
        ont_n = sum(go_summary$ontology == ont)
        cat(sprintf('  %s: %d terms\n', ont, ont_n))
      }
      
      suppressMessages(write_csv(go_summary,
                                 file.path(table_dir,
                                           'go_term_investigation_list.csv')))
    } else {
      cat('No significant GO terms from GSEA\n')
    }
  } else {
    cat('No GO GSEA results found\n')
  }
  
  # Flatten GO pathways (ontology separated -> all) for downstream use
  go_pathway_list = list()
  
  if (!is.null(go_summary) && nrow(go_summary) > 0) {
    for (ont in unique(go_summary$ontology)) {
      ont_sets = enrichment$gene_sets$go[[ont]]
      if (is.null(ont_sets)) next
      prefixed_sets = ont_sets
      names(prefixed_sets) = paste0(ont, ' | ', names(prefixed_sets))
      
      # Match against significant pathways
      ont_sig = go_summary %>% filter(ontology == ont)
      matching = prefixed_sets[names(prefixed_sets) %in% ont_sig$pathway]
      go_pathway_list = c(go_pathway_list, matching)
    }
    cat(sprintf('%d GO gene sets matched for heatmaps\n',
                length(go_pathway_list)))
  }
  
  # Pathway heatmaps
  cat('\nBuilding plots of significant pathways\n')
  
  # Row annotation builder function for heatmaps
  build_row_anno = function(gene_ids) {
    anno = master %>%
      filter(gene %in% gene_ids) %>%
      transmute(
        gene,
        Interaction = ifelse(sig_interaction, 'sig', 'NS'),
        `EHP vs NaCl` = ifelse(sig_EHP_vs_NaCl, 'sig', 'NS'),
        `dsEHP vs EHP` = ifelse(sig_dsEHP_vs_EHP, 'sig', 'NS'),
        `dsPTP2 vs NaCl` = ifelse(sig_dsPTP2_vs_NaCl, 'sig', 'NS')
      ) %>%
      column_to_rownames('gene') %>%
      as.data.frame()
    
    anno
  }
  
  # Function to convert pathway name to usable file name
  clean_pathway_name = function(pw_name, max_len = 80) {
    clean = gsub('[^A-Za-z0-9_]', '_', pw_name)
    clean = gsub('_+', '_', clean)
    substr(clean, 1, max_len)
  }
  
  # KEGG pathway heatmaps (multifactor VST)
  if (length(kegg_pathway_list) > 0) {
    cat(sprintf('\nPlotting %d KEGG pathway heatmaps\n',
                length(kegg_pathway_list)))
    
    for (pw_name in names(kegg_pathway_list)) {
      pw_genes = kegg_pathway_list[[pw_name]]
      clean_name = clean_pathway_name(pw_name)
      row_anno = build_row_anno(pw_genes)
      
      plot_pathway_heatmap(
        gene_ids = pw_genes,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        row_annotations = row_anno,
        filename = file.path(heatmap_dir,
                             sprintf('kegg_%s.png', clean_name)),
        title = pw_name
      )
      
      cat(sprintf('  KEGG: %s (%d genes)\n', pw_name,
                  sum(pw_genes %in% rownames(vst_mf))))
    }
  }
  
  # GO pathway heatmaps (limited to 20 max)
  if (length(go_pathway_list) > 0) {
    go_to_plot = names(go_pathway_list)
    if (length(go_to_plot) > 15) {
      # Prefer biological process terms
      bp_terms = go_to_plot[grepl('^Process', go_to_plot)]
      other_terms = setdiff(go_to_plot, bp_terms)
      go_to_plot = c(head(bp_terms, 15), head(other_terms, 5))
      cat(sprintf('Limiting GO heatmaps to %d (from %d significant)\n',
                  length(go_to_plot), length(go_pathway_list)))
    }
    
    cat(sprintf('\nPlotting %d GO pathway heatmaps\n', length(go_to_plot)))
    
    for (pw_name in go_to_plot) {
      pw_genes = go_pathway_list[[pw_name]]
      clean_name = clean_pathway_name(pw_name)
      row_anno = build_row_anno(pw_genes)
      
      plot_pathway_heatmap(
        gene_ids = pw_genes,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        row_annotations = row_anno,
        filename = file.path(heatmap_dir,
                             sprintf('go_%s.png', clean_name)),
        title = pw_name
      )
      
      cat(sprintf('  GO: %s (%d genes)\n', pw_name,
                  sum(pw_genes %in% rownames(vst_mf))))
    }
  }
  
  # Per-gene boxplots of multifactor DE for priority genes
  cat('\nBuilding gene-level plots\n')
  
  # Tier 1-3 genes from DE collation
  top_genes = master %>%
    filter(priority <= 3) %>%
    arrange(priority, padj_interaction)
  
  cat(sprintf('Plotting %d tier 1-3 genes\n', nrow(top_genes)))
  
  for (i in seq_len(nrow(top_genes))) {
    gene_id = top_genes$gene[i]
    tier = top_genes$priority[i]
    gene_display = resolve_labels(gene_id, gene_labels)
    
    padj_int = top_genes$padj_interaction[i]
    padj_pw  = top_genes$padj_dsEHP_vs_EHP[i]
    title = sprintf('T%d: %s\nInteraction padj=%.2g, Pairwise padj=%.2g',
                    tier, gene_display, padj_int, padj_pw)
    
    p = plot_gene_counts(
      gene_id, dds_mf, coldata_mf,
      highlight_sample = 'dsEHP_1',
      gene_labels = gene_labels,
      title = title
    )
    
    if (!is.null(p)) {
      quiet_ggsave(file.path(gene_dir,
                             sprintf('T%d_%02d_%s.png', tier, i, gene_id)),
                   p, width = plot_width, height = plot_height, dpi = plot_dpi)
    }
  }
  
  # Leading edge genes from interaction GSEA (KEGG)
  plotted_genes = top_genes$gene
  
  if ('interaction_KEGG' %in% names(enrichment$gsea)) {
    int_gsea = enrichment$gsea$interaction_KEGG
    sig_int = int_gsea %>% filter(padj < alpha)
    
    if (nrow(sig_int) > 0) {
      leading_edge_genes = unique(unlist(sig_int$leadingEdge))
      new_le_genes = setdiff(leading_edge_genes, plotted_genes)
      
      le_to_plot = master %>%
        filter(gene %in% new_le_genes) %>%
        arrange(desc(abs(stat_interaction))) %>%
        slice_head(n = 20)
      
      cat(sprintf('Plotting %d KEGG leading-edge genes not in tier 1-3\n',
                  nrow(le_to_plot)))
      
      for (i in seq_len(nrow(le_to_plot))) {
        gene_id = le_to_plot$gene[i]
        gene_display = resolve_labels(gene_id, gene_labels)
        
        p = plot_gene_counts(
          gene_id, dds_mf, coldata_mf,
          highlight_sample = 'dsEHP_1',
          gene_labels = gene_labels,
          title = sprintf('KEGG leading edge: %s', gene_display)
        )
        
        if (!is.null(p)) {
          quiet_ggsave(file.path(gene_dir,
                                 sprintf('LE_KEGG_%02d_%s.png', i, gene_id)),
                       p, width = plot_width, height = plot_height,
                       dpi = plot_dpi)
        }
      }
      
      plotted_genes = c(plotted_genes, le_to_plot$gene)
    }
  }
  
  # Leading edge genes from interaction GO process GSEA
  if ('interaction_GO_Process' %in% names(enrichment$gsea)) {
    go_int = enrichment$gsea[['interaction_GO_Process']]
    sig_go_int = go_int %>% filter(padj < alpha)
    
    if (nrow(sig_go_int) > 0) {
      go_le_genes = unique(unlist(sig_go_int$leadingEdge))
      new_go_le = setdiff(go_le_genes, plotted_genes)
      
      go_le_to_plot = master %>%
        filter(gene %in% new_go_le) %>%
        arrange(desc(abs(stat_interaction))) %>%
        slice_head(n = 15)
      
      cat(sprintf('Plotting %d GO Process leading-edge genes not yet plotted\n',
                  nrow(go_le_to_plot)))
      
      for (i in seq_len(nrow(go_le_to_plot))) {
        gene_id = go_le_to_plot$gene[i]
        gene_display = resolve_labels(gene_id, gene_labels)
        
        p = plot_gene_counts(
          gene_id, dds_mf, coldata_mf,
          highlight_sample = 'dsEHP_1',
          gene_labels = gene_labels,
          title = sprintf('GO Process leading edge: %s', gene_display)
        )
        
        if (!is.null(p)) {
          quiet_ggsave(file.path(gene_dir,
                                 sprintf('LE_GO_%02d_%s.png', i, gene_id)),
                       p, width = plot_width, height = plot_height,
                       dpi = plot_dpi)
        }
      }
      
      plotted_genes = c(plotted_genes, go_le_to_plot$gene)
    }
  }
  
  # Investigate immune pathway components (regardless of significance)
  cat('Immune gene search')
  
  candidate_terms = list(
    # Toll signalling pathway
    toll_pathway = 'Toll|MyD88|Dorsal|Cactus|Spätzle|Spatzle|pellino',
    # IMD signalling pathway
    imd_pathway = 'IMD|Relish|IKK|DREDD|FADD',
    # JAK-STAT signalling pathway
    jak_stat = 'JAK|\\bSTAT\\b|DOME|SOCS',
    # proPO / melanisation cascade
    propo_cascade = 'phenoloxidase|prophenoloxidase|proPO|clip.*serine|serpin',
    # Antimicrobial peptides
    amps = 'anti.lipopolysaccharide|crustin|penaeidin|lysozyme|defensin',
    # RNAi / siRNA machinery (somatic pathway only)
    rnai_machinery = 'Dicer|Argonaute|AGO|R2D2|SID|Drosha|Pasha',
    # Pattern recognition receptors
    prrs = 'lectin|galectin|LGBP|BGBP|PGRP|fibrinogen|DSCAM|ficolin',
    # Oxidative burst / antioxidant defence
    oxidative_burst = 'NADPH oxidase|superoxide dismutase|peroxidase|DUOX'
  )
  
  candidate_results = list()
  
  for (pathway_name in names(candidate_terms)) {
    pattern = candidate_terms[[pathway_name]]
    
    matches = master %>%
      filter(grepl(pattern, description, ignore.case = TRUE)) %>%
      select(gene, description, functional_category, priority,
             sig_interaction, sig_dsEHP_vs_EHP, sig_EHP_vs_NaCl,
             sig_dsPTP2_vs_NaCl,
             lfc_interaction, padj_interaction,
             lfc_dsEHP_vs_EHP, padj_dsEHP_vs_EHP,
             lfc_EHP_vs_NaCl, padj_EHP_vs_NaCl,
             lfc_dsPTP2_vs_NaCl, padj_dsPTP2_vs_NaCl)
    
    candidate_results[[pathway_name]] = matches
    
    n_sig_int = sum(matches$sig_interaction)
    n_sig_ehp = sum(matches$sig_EHP_vs_NaCl)
    n_sig_vac = sum(matches$sig_dsPTP2_vs_NaCl)
    
    cat(sprintf('\n  %s: %d genes found (%d sig interaction, %d sig EHP, '
                , pathway_name, nrow(matches), n_sig_int, n_sig_ehp))
    cat(sprintf('%d sig dsPTP2)\n', n_sig_vac))
    
    # Log individual genes if the set is small enough to be informative
    if (nrow(matches) > 0 && nrow(matches) <= 25) {
      for (j in seq_len(nrow(matches))) {
        sig_flags = paste0(
          ifelse(matches$sig_interaction[j], 'INT', ''),
          ifelse(matches$sig_dsEHP_vs_EHP[j], ' PW', ''),
          ifelse(matches$sig_EHP_vs_NaCl[j], ' EHP', ''),
          ifelse(matches$sig_dsPTP2_vs_NaCl[j], ' PTP2', ''))
        sig_flags = trimws(sig_flags)
        if (sig_flags == '') sig_flags = 'NS'
        
        cat(sprintf('    %s: %s [%s]\n',
                    matches$gene[j],
                    substr(matches$description[j], 1, 50),
                    sig_flags))
      }
    }
    
    # Per-gene count plots for candidate genes (dsEHP_1 highlighted)
    # Only plot genes not already plotted in section 3
    new_candidates = setdiff(matches$gene, plotted_genes)
    if (length(new_candidates) > 0) {
      for (j in seq_along(new_candidates)) {
        gene_id = new_candidates[j]
        gene_display = resolve_labels(gene_id, gene_labels)
        
        p = plot_gene_counts(
          gene_id, dds_mf, coldata_mf,
          highlight_sample = 'dsEHP_1',
          gene_labels = gene_labels,
          title = sprintf('%s: %s', pathway_name, gene_display)
        )
        
        if (!is.null(p)) {
          quiet_ggsave(file.path(gene_dir,
                                 sprintf('candidate_%s_%02d_%s.png',
                                         pathway_name, j, gene_id)),
                       p, width = plot_width, height = plot_height,
                       dpi = plot_dpi)
        }
      }
      plotted_genes = c(plotted_genes, new_candidates)
    }
    
    # Pathway-level heatmap if enough candidate genes in this group
    if (nrow(matches) >= 3) {
      row_anno = build_row_anno(matches$gene)
      
      plot_pathway_heatmap(
        gene_ids = matches$gene,
        vst_matrix = vst_mf,
        metadata = coldata_mf,
        gene_labels = gene_labels,
        row_annotations = row_anno,
        filename = file.path(heatmap_dir,
                             sprintf('candidate_%s.png', pathway_name)),
        title = sprintf('Candidate genes: %s', pathway_name)
      )
    }
  }
  
  # Save combined candidate gene table
  candidate_all = bind_rows(candidate_results, .id = 'pathway_group')
  suppressMessages(write_csv(candidate_all,
                             file.path(table_dir, 'candidate_genes.csv')))
  
  cat(sprintf('\n%d total candidate gene hits across %d pathway groups\n',
              nrow(candidate_all), length(candidate_results)))
  
  # Interaction/pairwise mechanism classification
  cat('\nInteraction/dsPTP2/EHP mechanism classification\n')
  
  interaction_genes = master %>%
    filter(sig_interaction | sens_class_interaction == 'dsEHP_1_dependent') %>%
    mutate(
      # Implied pairwise dsEHP vs EHP from the factorial model:
      # dsEHP = intercept + vaccination + challenge + interaction
      # EHP   = intercept + challenge
      # → dsEHP - EHP = vaccination + interaction
      # Retained to compare multifactor and pairwise model
      implied_dsEHP_vs_EHP = lfc_interaction + lfc_dsPTP2_vs_NaCl,
      
      dsRNA_response = case_when(
        # No dsRNA main effect
        
        # Interaction-only: neither main effect significant
        # Gene is silent in both dsPTP2 vs NaCl and EHP vs NaCl,
        # but the combination produces a significant effect.
        # Genuinely emergent — cannot be explained by either stimulus alone
        !sig_dsPTP2_vs_NaCl & !sig_EHP_vs_NaCl ~
          'interaction-only (emergent)',
        
        # Infection-responsive: no dsRNA main effect, but infection changes
        # the gene AND vaccination significantly alters that response.
        # Use sign(lfc_EHP_vs_NaCl) as reference direction (lfc_dsPTP2
        # is near zero so cannot serve as anchor).
        #
        # Amplified: β_int same sign as β_chal → vaccination pushes the
        # infection response further; |dsEHP-EHP| > |EHP alone|
        !sig_dsPTP2_vs_NaCl & sig_EHP_vs_NaCl &
          sign(lfc_interaction) == sign(lfc_EHP_vs_NaCl) ~
          'infection-responsive, vaccine-amplified',
        
        # Attenuated: β_int opposite sign to β_chal → vaccination dampens
        # or partially reverses the infection response; |dsEHP-EHP| < |EHP alone|
        !sig_dsPTP2_vs_NaCl & sig_EHP_vs_NaCl &
          sign(lfc_interaction) != sign(lfc_EHP_vs_NaCl) ~
          'infection-responsive, vaccine-attenuated',
        
        # dsRNA main effect present
        
        # Both stimuli push gene in same direction.
        # Interaction significant → combined effect deviates from additivity.
        # Check whether interaction amplifies or attenuates the shared direction:
        #   sign(interaction) == sign(main effects) → synergistic (super-additive)
        #   sign(interaction) != sign(main effects) → sub-additive (saturating)
        
        # Sub-additive: combined effect is LESS than sum of individual effects
        # The two stimuli partially saturate or are redundant
        sig_dsPTP2_vs_NaCl & sig_EHP_vs_NaCl &
          sign(lfc_EHP_vs_NaCl) == sign(lfc_dsPTP2_vs_NaCl) &
          sign(lfc_interaction) != sign(lfc_EHP_vs_NaCl) ~
          'co-directional (sub-additive)',
        
        # Synergistic: combined effect EXCEEDS sum of individual effects
        # Vaccination amplifies the infection response (or vice versa)
        sig_dsPTP2_vs_NaCl & sig_EHP_vs_NaCl &
          sign(lfc_EHP_vs_NaCl) == sign(lfc_dsPTP2_vs_NaCl) &
          sign(lfc_interaction) == sign(lfc_EHP_vs_NaCl) ~
          'co-directional (synergistic)',
        
        # dsRNA and infection push gene in opposite directions
        # Interaction significant → the net outcome deviates from the
        # tug-of-war prediction. Does not imply which stimulus "wins"
        sig_dsPTP2_vs_NaCl & sig_EHP_vs_NaCl &
          sign(lfc_EHP_vs_NaCl) != sign(lfc_dsPTP2_vs_NaCl) ~
          'opposing main effects',
        
        # dsRNA changes gene but infection alone does not.
        # Yet interaction is significant — challenge alters the dsRNA-primed
        # state in a way infection alone never produces.
        # Use sign(lfc_dsPTP2_vs_NaCl) as the reference direction
        #
        # Synergistic: β_int same sign as β_vacc → challenge amplifies
        # the primed response; |implied| > |dsPTP2 alone|
        sig_dsPTP2_vs_NaCl & !sig_EHP_vs_NaCl &
          sign(lfc_interaction) == sign(lfc_dsPTP2_vs_NaCl) ~
          'dsRNA-primed (synergistic)',
        
        # Sub-additive: β_int opposite sign to β_vacc → challenge attenuates
        # or partially reverses the primed state; |implied| < |dsPTP2 alone|
        sig_dsPTP2_vs_NaCl & !sig_EHP_vs_NaCl &
          sign(lfc_interaction) != sign(lfc_dsPTP2_vs_NaCl) ~
          'dsRNA-primed (sub-additive)'
      )
    )
  
  # Summary counts
  dsrna_summary = interaction_genes %>%
    dplyr::count(dsRNA_response)
  cat('\nInteraction DEGs by mechanism:\n')
  print(as.data.frame(dsrna_summary))
  
  # By functional category
  dsrna_by_category = interaction_genes %>%
    dplyr::count(dsRNA_response, functional_category) %>%
    pivot_wider(names_from = functional_category,
                values_from = n, values_fill = 0)
  cat('\nBy functional category:\n')
  print(as.data.frame(dsrna_by_category))
  
  # By dsEHP_1 sensitivity
  dsrna_by_sens = interaction_genes %>%
    filter(!is.na(sens_class_interaction)) %>%
    dplyr::count(dsRNA_response, sens_class_interaction) %>%
    pivot_wider(names_from = sens_class_interaction,
                values_from = n, values_fill = 0)
  cat('\nBy dsEHP_1 sensitivity:\n')
  print(as.data.frame(dsrna_by_sens))
  
  # Residual vaccine effect: does any dsRNA-primed advantage survive?
  # For sub-additive genes: if near zero → effects fully saturate
  # For synergistic genes: should exceed either main effect alone
  residual_by_class = interaction_genes %>%
    group_by(dsRNA_response) %>%
    summarise(
      n = n(),
      median_abs_interaction = median(abs(lfc_interaction)),
      median_abs_dsPTP2 = median(abs(lfc_dsPTP2_vs_NaCl)),
      median_abs_EHP = median(abs(lfc_EHP_vs_NaCl), na.rm = TRUE),
      median_abs_dsEHP_vs_EHP = median(abs(lfc_dsEHP_vs_EHP)),
      .groups = 'drop'
    )
  cat('\nResidual vaccine effect by mechanism:\n')
  print(as.data.frame(residual_by_class))
  
  # Survival fraction (dsRNA-responsive groups only)
  survival = interaction_genes %>%
    filter(sig_dsPTP2_vs_NaCl) %>%
    group_by(dsRNA_response) %>%
    summarise(
      n = n(),
      median_survival = median(
        abs(lfc_dsEHP_vs_EHP) / abs(lfc_dsPTP2_vs_NaCl)),
      .groups = 'drop'
    )
  cat('\nSurvival fraction (dsRNA-responsive groups only):\n')
  cat('(fraction of dsRNA-alone effect retained when challenged)\n')
  print(as.data.frame(survival))
  
  suppressMessages(write_csv(interaction_genes,
                             file.path(table_dir,
                                       'interaction_mechanism_crossref.csv')))
  
  # Flag immune/RNAi genes in each mechanism class
  for (mech in c('interaction-only (emergent)',
                 'infection-responsive, vaccine-amplified',
                 'infection-responsive, vaccine-attenuated',
                 'co-directional (synergistic)',
                 'co-directional (sub-additive)',
                 'opposing main effects',
                 'dsRNA-primed (synergistic)',
                 'dsRNA-primed (sub-additive)')) {
    immune_in_class = interaction_genes %>%
      filter(dsRNA_response == mech,
             functional_category %in% c('immune', 'RNAi'))
    
    if (nrow(immune_in_class) > 0) {
      cat(sprintf('\n%d immune/RNAi genes [%s]:\n',
                  nrow(immune_in_class), mech))
      for (i in seq_len(nrow(immune_in_class))) {
        cat(sprintf('  %s: %s\n    int=%.2f±%.2f, dsPTP2=%.2f±%.2f, EHP=%.2f±%.2f,\ndsEHP-NaCl=%.2f±%.2f, dsEHP-EHP=%.2f±%.2f\n'
                    , immune_in_class$gene[i]
                    , substr(immune_in_class$description[i], 1, 50)
                    , immune_in_class$lfc_interaction[i]
                    , immune_in_class$lfcSE_interaction[i]
                    , immune_in_class$lfc_dsPTP2_vs_NaCl[i]
                    , immune_in_class$lfcSE_dsPTP2_vs_NaCl[i]
                    , immune_in_class$lfc_EHP_vs_NaCl[i]
                    , immune_in_class$lfcSE_EHP_vs_NaCl[i]
                    , immune_in_class$lfc_dsEHP_vs_NaCl[i]
                    , immune_in_class$lfcSE_dsEHP_vs_NaCl[i]
                    , immune_in_class$lfc_dsEHP_vs_EHP[i]
                    , immune_in_class$lfcSE_dsEHP_vs_EHP[i]))
      }
    }
  }
  
  # Summary table of most important genes
  cat('Building summary table of important genes')
  
  mechanism_map = interaction_genes %>%
    select(gene, dsRNA_response)
  
  narrative = master %>%
    filter(priority <= 4) %>%
    left_join(mechanism_map, by = 'gene') %>%
    select(gene, description, symbol, functional_category, priority,
           # Interaction term (primary vaccine contrast)
           lfc_interaction, lfcSE_interaction, padj_interaction,
           # Pairwise dsEHP vs EHP
           lfc_dsEHP_vs_EHP, lfcSE_dsEHP_vs_EHP, padj_dsEHP_vs_EHP,
           # dsPTP2 vs NaCl (dsRNA-alone baseline)
           lfc_dsPTP2_vs_NaCl, lfcSE_dsPTP2_vs_NaCl, padj_dsPTP2_vs_NaCl,
           # EHP vs NaCl (infection response)
           lfc_EHP_vs_NaCl, lfcSE_EHP_vs_NaCl, padj_EHP_vs_NaCl,
           # dsEHP vs NaCl (combined vaccinated + challenged vs baseline)
           lfc_dsEHP_vs_NaCl, lfcSE_dsEHP_vs_NaCl, padj_dsEHP_vs_NaCl,
           # Sensitivity to dsEHP_1
           sens_class_pairwise, sens_class_interaction,
           # Mechanism classification (from section 5)
           dsRNA_response) %>%
    arrange(priority, padj_interaction)
  
  suppressMessages(write_csv(narrative,
                             file.path(table_dir, 'narrative_summary.csv')))
  
  cat(sprintf('%d genes in narrative summary (tiers 1-4)\n', nrow(narrative)))
  cat(sprintf('Tier 1 (immune + convergent): %d\n',
              sum(narrative$priority == 1)))
  cat(sprintf('Tier 2 (immune + vaccine): %d\n',
              sum(narrative$priority == 2)))
  cat(sprintf('Tier 3 (annotated + convergent): %d\n',
              sum(narrative$priority == 3)))
  cat(sprintf('Tier 4 (immune + baseline): %d\n',
              sum(narrative$priority == 4)))
  
  # Sensitivity composition of narrative genes
  if ('sens_class_interaction' %in% names(narrative)) {
    narrative_sens = narrative %>%
      filter(!is.na(sens_class_interaction)) %>%
      dplyr::count(priority, sens_class_interaction) %>%
      pivot_wider(names_from = sens_class_interaction,
                  values_from = n, values_fill = 0)
    cat('\nNarrative genes by dsEHP_1 sensitivity:\n')
    print(as.data.frame(narrative_sens))
  }
  
  # Log summary
  cat('\nCandidate pathway analysis complete\n')
  cat(sprintf('Pathway heatmaps: %s\n', heatmap_dir))
  cat(sprintf('Gene plots: %s\n', gene_dir))
  cat(sprintf('Tables: %s\n', table_dir))
  cat(sprintf('\nTotal genes plotted: %d\n', length(unique(plotted_genes))))
  cat(sprintf('KEGG pathways investigated: %d\n', length(kegg_pathway_list)))
  cat(sprintf('GO terms investigated: %d\n', length(go_pathway_list)))
  cat(sprintf('Candidate pathway groups: %d\n', length(candidate_terms)))
}

main()
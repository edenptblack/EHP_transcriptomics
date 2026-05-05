# scripts/05_pathway_enrichment.R
# KEGG & GO pathway enrichment via gene set enrichment analysis (GSEA, primary) &
# overrepresentation analysis (ORA, supplementary)

source('config.R')
source('R/logging.R')
source('R/packages.R')
source('R/annotation_functions.R')
source('R/enrichment_functions.R')
source('R/plotting.R')

main = function() {
  
  log_con = start_log('05_pathway_enrichment')
  on.exit(end_log(log_con))
  
  proc_dir   = config$proc_dir
  alpha      = config$alpha
  res_dir    = 'results/05_pathway_enrichment'
  table_dir  = file.path(res_dir, 'tables')
  plot_dir   = file.path(res_dir, 'plots')
  collation  = 'results/04_collate_results'
  gsea_dir   = file.path(collation, 'gsea_ranks')
  ora_dir    = file.path(collation, 'ora_sets')
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_dpi    = config$plot_dpi
  plot_width  = config$plot_width
  plot_height = config$plot_height
  
  # Load data
  master = readRDS(file.path(proc_dir, 'master_gene_table.rds'))
  annotation_df = readRDS(file.path(proc_dir, 'annotation_table.rds'))
  
  gene_universe = master$gene
  cat(sprintf('Gene universe: %d genes\n', length(gene_universe)))
  
  # Build gene sets
  cat('\n Building gene sets...\n')
  
  kegg_sets = build_kegg_gene_sets(annotation_df, gene_universe, label = 'primary')
  go_sets = build_go_gene_sets(annotation_df, gene_universe)
  
  custom_sets = list('CUSTOM: p.van immune & RNAi genes' = master %>%
                       filter(functional_category %in% c('immune', 'RNAi')) %>%
                       pull(gene))
  
  cat(sprintf('\nGene sets: %d KEGG, %d GO (%s), %d custom\n',
              length(kegg_sets),
              sum(sapply(go_sets, length)),
              paste(sapply(go_sets, length), collapse = '/'),
              length(custom_sets)))
  
  # GSEA
  cat('\nRunning GSEA...\n')
  
  gsea_contrasts = config$gsea_contrasts
  gsea_results = list()
  
  for (name in gsea_contrasts) {
    rnk_file = file.path(gsea_dir, sprintf('gsea_ranks_%s.rnk', name))
    if (!file.exists(rnk_file)) {
      cat(sprintf('WARNING: %s not found, skipping\n', rnk_file))
      next
    }
    
    ranked = read_tsv(rnk_file, col_names = c('gene', 'stat'),
                      show_col_types = FALSE)
    cat(sprintf('\n%s (%d genes ranked)\n', name, nrow(ranked)))
    
    # KEGG
    if (length(kegg_sets) > 0) {
      kegg_res = run_fgsea(ranked, kegg_sets)
      gsea_results[[paste0(name, '_KEGG')]] = kegg_res
      summarise_fgsea(kegg_res, alpha = alpha, label = paste(name, 'KEGG'))
      
      suppressMessages(write_csv(
        kegg_res %>% select(-leadingEdge),
        file.path(table_dir, sprintf('gsea_%s_KEGG.csv', name))))
    }
    
    # GO per ontology
    for (ont in names(go_sets)) {
      
      ont_sets = go_sets[[ont]]
      if (length(ont_sets) == 0) next
      
      # Prefix with full ontology name
      names(ont_sets) = paste0(ont, ' | ', names(ont_sets))
      
      go_res = run_fgsea(ranked, ont_sets)
      
      result_key = paste0(name, '_GO_', ont)
      gsea_results[[result_key]] = go_res
      
      summarise_fgsea(go_res, alpha = alpha,
                      label = paste(name, 'GO', ont))
      
      suppressMessages(write_csv(
        go_res %>% select(-leadingEdge),
        file.path(table_dir,
                  sprintf('gsea_%s_GO_%s.csv', name, ont))
      ))
    }
    
    # Custom set
    custom_res = run_fgsea(ranked, custom_sets, min_size = 3)
    gsea_results[[paste0(name, '_custom')]] = custom_res
    summarise_fgsea(custom_res, alpha = 1.0,
                    label = paste(name, 'custom immune'))
    
    suppressMessages(write_csv(
      custom_res %>% select(-leadingEdge),
      file.path(table_dir, sprintf('gsea_%s_custom.csv', name))))
  }
  
  # GSEA dot plots
  for (result_name in names(gsea_results)) {
    res = gsea_results[[result_name]]
    
    p = plot_fgsea_dotplot(res, n_show = 20, alpha = alpha,
                           title = sprintf('GSEA: %s', result_name))
    
    if (!is.null(p)) {
      quiet_ggsave(file.path(plot_dir,
                             sprintf('gsea_dotplot_%s.png', result_name)),
                   p, width = 10, height = plot_height + 2, dpi = plot_dpi)
    }
  }
  
  # Enrichment curves for top 3 pathways in KEGG & GO biological process enrichment
  # Limited to primary contrast (multifactor interaction term)
  enrich_curve_dir = file.path(plot_dir, 'enrichment_curves')
  dir.create(enrich_curve_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Contrasts to plot enrichment curves for
  curve_contrasts = c('interaction', 'dsPTP2_vs_NaCl')
  
  for (cname in curve_contrasts) {
    rnk_file = file.path(gsea_dir, sprintf('gsea_ranks_%s.rnk', cname))
    if (!file.exists(rnk_file)) next
    
    ranked = read_tsv(rnk_file, col_names = c('gene', 'stat'),
                      show_col_types = FALSE)
    ranks = setNames(ranked$stat, ranked$gene)
    ranks = ranks[!is.na(ranks)]
    ranks = ranks[!duplicated(names(ranks))]
    
    # KEGG curves
    kegg_key = paste0(cname, '_KEGG')
    if (kegg_key %in% names(gsea_results)) {
      top_kegg = gsea_results[[kegg_key]] %>%
        filter(padj < alpha) %>%
        arrange(padj) %>%
        slice_head(n = 3)
      
      if (nrow(top_kegg) > 0) {
        for (i in seq_len(nrow(top_kegg))) {
          pw = top_kegg$pathway[i]
          if (!pw %in% names(kegg_sets)) next
          
          p = fgsea::plotEnrichment(kegg_sets[[pw]], ranks) +
            labs(title = sprintf('%s\n%s (padj=%.2g, NES=%.2f)',
                                 cname, pw, top_kegg$padj[i], top_kegg$NES[i]))
          
          clean_pw = gsub('[^A-Za-z0-9_]', '_', pw)
          clean_pw = gsub('_+', '_', clean_pw)
          clean_pw = substr(clean_pw, 1, 60)
          
          quiet_ggsave(file.path(enrich_curve_dir,
                                 sprintf('enrichment_%s_KEGG_%02d_%s.png',
                                         cname, i, clean_pw)),
                       p, width = 8, height = 5, dpi = plot_dpi)
        }
      }
    }
    
    # GO BP curves
    process_key = paste0(cname, '_GO_Process')
    
    if (process_key %in% names(gsea_results)) {
      process_sets = go_sets[['Process']]
      if (length(process_sets) > 0) {
        names(process_sets) = paste0('Process | ', names(process_sets))
      }
      
      top_bp = gsea_results[[process_key]] %>%
        filter(padj < alpha) %>%
        arrange(padj) %>%
        slice_head(n = 3)
      
      if (nrow(top_bp) > 0) {
        for (i in seq_len(nrow(top_bp))) {
          pw = top_bp$pathway[i]
          if (!pw %in% names(process_sets)) next
          
          p = fgsea::plotEnrichment(process_sets[[pw]], ranks) +
            labs(title = sprintf('%s\n%s (padj=%.2g, NES=%.2f)',
                                 cname, pw, top_bp$padj[i], top_bp$NES[i]))
          
          clean_pw = gsub('[^A-Za-z0-9_]', '_', pw)
          clean_pw = gsub('_+', '_', clean_pw)
          clean_pw = substr(clean_pw, 1, 60)
          
          quiet_ggsave(file.path(enrich_curve_dir,
                                 sprintf('enrichment_%s_GO_Process_%02d_%s.png',
                                         cname, i, clean_pw)),
                       p, width = 8, height = 5, dpi = plot_dpi)
      }
    }
    }
  }
  
  
  
  # ORA
  cat('\nRunning ORA...\n')
  
  background = readLines(file.path(ora_dir, 'background_all_tested.txt'))
  
  # Define all ORA gene sets to test, with expected filenames from script 04
  ora_file_map = list(
    # Primary contrasts — all/up/down splits
    EHP_vs_NaCl_all = 'ora_EHP_vs_NaCl_all.txt',
    EHP_vs_NaCl_up = 'ora_EHP_vs_NaCl_up.txt',
    EHP_vs_NaCl_down = 'ora_EHP_vs_NaCl_down.txt',
    interaction_all = 'ora_interaction_all.txt',
    interaction_up = 'ora_interaction_up.txt',
    interaction_down = 'ora_interaction_down.txt',
    dsPTP2_vs_NaCl_all = 'ora_dsPTP2_vs_NaCl_all.txt',
    dsPTP2_vs_NaCl_up = 'ora_dsPTP2_vs_NaCl_up.txt',
    dsPTP2_vs_NaCl_down = 'ora_dsPTP2_vs_NaCl_down.txt',
    dsEHP_vs_EHP_all = 'ora_dsEHP_vs_EHP_all.txt',
    dsEHP_vs_EHP_up = 'ora_dsEHP_vs_EHP_up.txt',
    dsEHP_vs_EHP_down = 'ora_dsEHP_vs_EHP_down.txt',
    # Sensitivity-derived exploratory sets
    dsEHP1_dep_pairwise = 'ora_dsEHP1_dependent_pairwise.txt',
    dsEHP1_dep_interaction = 'ora_dsEHP1_dependent_interaction.txt',
    dsEHP_vs_EHP_robust = 'ora_dsEHP_vs_EHP_robust.txt'
  )
  
  ora_results = list()
  
  for (name in names(ora_file_map)) {
    ora_file = file.path(ora_dir, ora_file_map[[name]])
    if (!file.exists(ora_file)) {
      cat(sprintf('%s: file not found, skipping\n', name))
      next
    }
    
    sig_genes = readLines(ora_file)
    cat(sprintf('\n  %s: %d genes\n', name, length(sig_genes)))
    
    if (length(sig_genes) < 3) {
      cat('Too few genes for ORA, skipping\n')
      next
    }
    
    # KEGG
    if (length(kegg_sets) > 0) {
      kegg_ora = run_ora(sig_genes, background, kegg_sets, alpha = alpha)
      ora_key = paste0(name, '_KEGG')
      
      if (!is.null(kegg_ora)) {
        ora_results[[ora_key]] = kegg_ora
        summarise_ora(kegg_ora, alpha = alpha,
                      label = paste(name, 'KEGG'))
        
        res_df = as.data.frame(kegg_ora)
        if (nrow(res_df) > 0) {
          suppressMessages(write_csv(
            res_df,
            file.path(table_dir, sprintf('ora_%s_KEGG.csv', name))))
        }
      }
    }
    
    # GO biological processes
    # Only for larger gene sets where ORA has reasonable power
    # (interaction_all, dsPTP2_vs_NaCl_all, EHP_vs_NaCl_all)
    if (length(sig_genes) >= 50 && length(go_sets[['Process']]) > 0) {
      
      process_sets = go_sets[['Process']]
      names(process_sets) = paste0('Process | ', names(process_sets))
      
      go_ora = run_ora(sig_genes, background, process_sets, alpha = alpha)
      ora_key_go = paste0(name, '_GO_Process')
      
      if (!is.null(go_ora)) {
        ora_results[[ora_key_go]] = go_ora
        summarise_ora(go_ora, alpha = alpha,
                      label = paste(name, 'GO Process'))
        
        res_df = as.data.frame(go_ora)
        if (nrow(res_df) > 0) {
          suppressMessages(write_csv(
            res_df,
            file.path(table_dir, sprintf('ora_%s_GO_Process.csv', name))))
        }
      }
    }
  }
  
  # ORA dot plots
  for (ora_name in names(ora_results)) {
    ora_res = ora_results[[ora_name]]
    if (is.null(ora_res)) next
    
    res_df = as.data.frame(ora_res)
    if (nrow(res_df) == 0 || sum(res_df$p.adjust < alpha) == 0) next
    
    p = dotplot(ora_res, showCategory = 15,
                title = sprintf('ORA: %s', ora_name))
    quiet_ggsave(file.path(plot_dir,
                           sprintf('ora_dotplot_%s.png', ora_name)),
                 p, width = 10, height = plot_height + 2,
                 dpi = plot_dpi)
  }
  
  # Cross-contrast enrichment summary
  cat('\nCross-contrast pathway enrichment summary\n')
  
  # Sig KEGG pathways across all GSEA contrasts
  kegg_gsea_names = grep('_KEGG$', names(gsea_results), value = TRUE)
  
  cross_contrast = NULL
  
  if (length(kegg_gsea_names) > 0) {
    cross_contrast = bind_rows(
      lapply(kegg_gsea_names, function(rname) {
        gsea_results[[rname]] %>%
          filter(padj < alpha) %>%
          mutate(contrast = sub('_KEGG$', '', rname)) %>%
          select(contrast, pathway, NES, padj, size)
      })
    )
    
    if (nrow(cross_contrast) > 0) {
      pathway_counts = cross_contrast %>%
        dplyr::count(pathway, name = 'n_contrasts') %>%
        arrange(desc(n_contrasts))
      
      cat('\nKEGG pathways significant in multiple contrasts:\n')
      multi = pathway_counts %>% filter(n_contrasts > 1)
      if (nrow(multi) > 0) {
        for (i in seq_len(nrow(multi))) {
          cat(sprintf('  %s (%d contrasts)\n',
                      multi$pathway[i], multi$n_contrasts[i]))
        }
      } else {
        cat('  None\n')
      }
      
      suppressMessages(write_csv(cross_contrast,
                                 file.path(table_dir,
                                           'gsea_cross_contrast_KEGG.csv')))
      suppressMessages(write_csv(pathway_counts,
                                 file.path(table_dir,
                                           'pathway_contrast_counts_KEGG.csv')))
    }
  }
  
  # Sig GO biological process pathways across all GSEA contrasts
  bp_gsea_names = grep('_GO_Process$', names(gsea_results), value = TRUE)
  
  cross_contrast_bp = NULL
  
  if (length(bp_gsea_names) > 0) {
    cross_contrast_bp = bind_rows(
      lapply(bp_gsea_names, function(rname) {
        gsea_results[[rname]] %>%
          filter(padj < alpha) %>%
          mutate(contrast = sub('_GO_Process$', '', rname)) %>%
          select(contrast, pathway, NES, padj, size)
      })
    )
    
    if (nrow(cross_contrast_bp) > 0) {
      bp_pathway_counts = cross_contrast_bp %>%
        dplyr::count(pathway, name = 'n_contrasts') %>%
        arrange(desc(n_contrasts))
      
      cat('\nGO BP terms significant in multiple contrasts:\n')
      multi_bp = bp_pathway_counts %>% filter(n_contrasts > 1)
      if (nrow(multi_bp) > 0) {
        for (i in seq_len(min(10, nrow(multi_bp)))) {
          cat(sprintf('  %s (%d contrasts)\n',
                      multi_bp$pathway[i], multi_bp$n_contrasts[i]))
        }
        if (nrow(multi_bp) > 10) {
          cat(sprintf('  ... and %d more\n', nrow(multi_bp) - 10))
        }
      } else {
        cat('  None\n')
      }
      
      suppressMessages(write_csv(cross_contrast_bp,
                                 file.path(table_dir,
                                           'gsea_cross_contrast_GO_Process.csv')))
      suppressMessages(write_csv(bp_pathway_counts,
                                 file.path(table_dir,
                                           'pathway_contrast_counts_GO_Process.csv')))
    }
  }
  
  # GSEA/ORA agreement check for main contrast (interaction term, KEGG)
  # GSEA-only indicates modest shifts across pathway
  # ORA-only indicates a few strongly shifted individual genes
  cat('\nGSEA vs ORA agreement (multifactor interaction, KEGG)')
  
  gsea_int_key = 'interaction_KEGG'
  ora_int_key = 'interaction_all_KEGG'
  
  if (gsea_int_key %in% names(gsea_results) &&
      ora_int_key %in% names(ora_results)) {
    
    gsea_sig = gsea_results[[gsea_int_key]] %>%
      filter(padj < alpha) %>%
      pull(pathway)
    ora_sig = as.data.frame(ora_results[[ora_int_key]]) %>%
      filter(p.adjust < alpha) %>%
      pull(ID)
    
    both = intersect(gsea_sig, ora_sig)
    gsea_only = setdiff(gsea_sig, ora_sig)
    ora_only = setdiff(ora_sig, gsea_sig)
    
    cat(sprintf('Both GSEA & ORA: %d\n', length(both)))
    if (length(both) > 0) {
      for (pw in both) cat(sprintf('  %s\n', pw))
    }
    cat(sprintf('GSEA only: %d\n', length(gsea_only)))
    cat(sprintf('ORA only: %d\n', length(ora_only)))
    
    agreement_df = bind_rows(
      if (length(both) > 0) data.frame(pathway = both, status = 'both'),
      if (length(gsea_only) > 0) data.frame(pathway = gsea_only,
                                            status = 'GSEA_only'),
      if (length(ora_only) > 0) data.frame(pathway = ora_only,
                                           status = 'ORA_only')
    )
    
    if (nrow(agreement_df) > 0) {
      suppressMessages(write_csv(agreement_df,
                                 file.path(table_dir,
                                           'gsea_ora_agreement_interaction.csv')))
    }
  }
  
  # Save results
  all_enrichment = list(
    gsea = gsea_results,
    ora = ora_results,
    gene_sets = list(kegg = kegg_sets,
                     go = go_sets,
                     custom = custom_sets),
    cross_contrast_kegg = cross_contrast,
    cross_contrast_bp = cross_contrast_bp
  )
  
  saveRDS(all_enrichment, file.path(proc_dir, 'enrichment_results.rds'))
  
  cat(sprintf('\nEnrichment complete\n'))
  cat(sprintf('Tables: %s\n', table_dir))
  cat(sprintf('Plots:  %s\n', plot_dir))
}

main()
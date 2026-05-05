# NB: Due to time constraints, this script was written by AI with minimal
# human editing
#
# R/export_hhpred_candidates.R
# Export a pre-populated CSV for the HHPred tracker.
# Inputs: master_gene_table.rds (Tiers 5-6) + interaction_mechanism_crossref.csv
# Output:  results/hhpred_candidates.csv

source('config.R')
source('R/logging.R')
source('R/packages.R')

main = function() {
  log_con = start_log('export_hhpred_candidates')
  on.exit(end_log(log_con))
  
  proc_dir    = config$proc_dir
  crossref_file  = 'results/06_candidate_pathways/tables/interaction_mechanism_crossref.csv'
  out_file       = 'results/hhpred_candidates.csv'
  
  # ── Load ──────────────────────────────────────────────────────────────────
  master = readRDS(file.path(proc_dir, 'master_gene_table.rds'))
  
  crossref = read_csv(crossref_file, show_col_types = FALSE) %>%
    select(gene, dsRNA_response)
  
  cat(sprintf('Master table: %d genes\n', nrow(master)))
  cat(sprintf('Mechanism crossref: %d genes\n', nrow(crossref)))
  
  # ── Filter to HHPred candidates (Tiers 5-6) ───────────────────────────────
  candidates = master %>%
    filter(priority %in% c(5L, 6L)) %>%
    select(
      gene_id            = gene,
      symbol,
      description,
      protein_accession,
      priority,
      sig_interaction,
      sig_dsEHP_vs_EHP,
      sensitivity        = sens_class_interaction,
      lfc_interaction,
      padj_interaction,
      lfc_dsEHP_vs_EHP,
      padj_dsEHP_vs_EHP,
      baseMean
    ) %>%
    left_join(crossref, by = c('gene_id' = 'gene')) %>%
    rename(mechanism_class = dsRNA_response) %>%
    # Genes not in crossref are sig_dsEHP_vs_EHP only (not interaction-sig)
    mutate(mechanism_class = replace_na(mechanism_class, 'not interaction-significant')) %>%
    # Add empty columns matching tracker layout — fill manually after HHPred
    mutate(
      status            = 'Not started',
      hhpred_top_hit    = NA_character_,
      hit_database      = NA_character_,
      probability_pct   = NA_real_,
      e_value           = NA_character_,
      query_coverage_pct= NA_real_,
      pfam_domain       = NA_character_,
      proposed_function = NA_character_,
      reclassify        = NA_character_,
      thesis_relevance  = NA_character_,
      discussion_priority = NA_character_,
      notes             = NA_character_
    ) %>%
    arrange(priority, padj_interaction, padj_dsEHP_vs_EHP)
  
  n_t5 = sum(candidates$priority == 5L)
  n_t6 = sum(candidates$priority == 6L)
  n_acc = sum(!is.na(candidates$protein_accession))
  n_mech = sum(candidates$mechanism_class != 'not interaction-significant')
  
  cat(sprintf('\nHHPred candidates: %d total\n', nrow(candidates)))
  cat(sprintf('  Tier 5 (convergent):      %d\n', n_t5))
  cat(sprintf('  Tier 6 (single-contrast): %d\n', n_t6))
  cat(sprintf('  With protein accession:   %d\n', n_acc))
  cat(sprintf('  With mechanism class:     %d\n', n_mech))
  
  cat('\nMechanism class distribution:\n')
  print(as.data.frame(dplyr::count(candidates, priority, mechanism_class)))
  
  cat('\nSensitivity distribution:\n')
  print(as.data.frame(dplyr::count(candidates, priority, sensitivity)))
  
  # ── Write ─────────────────────────────────────────────────────────────────
  suppressMessages(write_csv(candidates, out_file))
  cat(sprintf('\nWritten: %s\n', out_file))
}

main()
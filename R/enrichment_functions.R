# R/enrichment_functions.R
# Functions for KEGG/GO pathway enrichment analysis

source('R/packages.R')
source('config.R')

# Build gene sets for KEGG - download KEGG mappings and convert to list for fgsea
# Input: df with gene_id & ncbi_id columns; gene universe (vector of genes to use);
# label for saving dataset (in case save multiple); organism (pvm); directory
# Output: Names list of gene id vectors, (pathway = gene IDs)
build_kegg_gene_sets = function(annotation_df, gene_universe, label,
                                organism = 'pvm', ref_dir = config$ref_dir) {
  kegg_cache = file.path(ref_dir, sprintf('kegg_gene_sets_%s.rds',
                                          label))
  
  if (file.exists(kegg_cache)) {
    cat('Loading cached KEGG gene sets from', kegg_cache, '\n')
    return(readRDS(kegg_cache))
  }
  
  cat('\nDownloading KEGG pathway data for', organism, '...\n')
  
  kegg_data = clusterProfiler::download_KEGG(organism)
  
  # Map KEGG Entrez IDs to gene IDs using annotation table
  pathway_genes = kegg_data$KEGGPATHID2EXTID %>%
    inner_join(annotation_df %>%
                 filter(gene_id %in% gene_universe) %>%
                 select(gene_id, ncbi_gene_id) %>%
                 distinct(ncbi_gene_id, .keep_all = TRUE),
               by = c('to' = 'ncbi_gene_id'))
  
  # Create list of pathways each containing vector of corresponding genes
  # (split for one-to-many relationship)
  kegg_pathways = split(pathway_genes$gene_id, pathway_genes$from)
  
  # Attach pathway names (setNames for one-to-one relationship)
  name_map = setNames(kegg_data$KEGGPATHID2NAME$to,
                      kegg_data$KEGGPATHID2NAME$from)
  
  names(kegg_pathways) = ifelse(
    names(kegg_pathways) %in% names(name_map),
    paste0(names(kegg_pathways), ': ', name_map[names(kegg_pathways)]),
    names(kegg_pathways)
  )
  
  cat(sprintf('Built %d KEGG gene sets (median size: %d genes; range: %d-%d)\n',
      length(kegg_pathways), median(sapply(kegg_pathways, length)),
      min(sapply(kegg_pathways, length)), max(sapply(kegg_pathways, length))))
  
  saveRDS(kegg_pathways, kegg_cache)
  cat('Saved to', kegg_cache, '\n')
  
  invisible(kegg_pathways)
}

# Build GO gene sets from annotation table
# Inputs: annotation_df, gene universe
# Output: Named list of gene ID vectors (1 per GO term)
build_go_gene_sets = function(annotation_df, gene_universe) {
  
  # Parse GO terms into a long format table
  go_long = annotation_df %>%
    filter(gene_id %in% gene_universe, !is.na(go_id)) %>%
    select(gene_id, go_id, go_categories, go_terms) %>%
    separate_rows(go_id, go_categories, go_terms, sep = '; ') %>%
    filter(go_id != '')
  
  if (nrow(go_long) == 0) {
    cat('WARNING: No GO annotations found in gene universe\n')
    return(list(Process = list(), Function = list(), Component = list()))
  }
  
  name_df = go_long %>%
    distinct(go_id, go_terms) %>%
    filter(go_terms != '')
    
  name_map = setNames(name_df$go_terms, name_df$go_id)
  
  # Build a list of gene sets for each GO ontology
  # Biological processes (primary analysis), molecular functions, cell component
  ontology_list = list()
  
  for (ont in unique(go_long$go_categories)) {
    ont_long = filter(go_long, go_categories == ont)
    gene_sets = split(ont_long$gene_id, ont_long$go_id)
    
    names(gene_sets) = ifelse(
      names(gene_sets) %in% names(name_map),
            paste0(names(gene_sets), ': ', name_map[names(gene_sets)]),
            names(gene_sets)
    )
    
    ontology_list[[ont]] = gene_sets
    cat(sprintf('%s: %d gene sets from %d genes\n', ont, length(gene_sets),
                length(unique(ont_long$gene_id))))
  }
  
  cat(sprintf('Built %d GO gene sets total\n', sum(sapply(ontology_list, length))))
  
  invisible(ontology_list)
}

# Gene Set Enrichment Analysis
# Run FGSEA on a ranked gene list against a collection of gene sets
# Inputs: ranked_df (cols: gene, stat); gene sets list; min/max size filters, seed
# Output: FGSEA results df
run_fgsea = function(ranked_df, gene_sets, min_size = 15,
                     max_size = 500, seed = 50) {
  
  # Turn ranked df into a named vector
  ranks = setNames(ranked_df$stat, ranked_df$gene)
  ranks = ranks[!is.na(ranks)]
  ranks = ranks[!duplicated(names(ranks))]
  
  # Perform gene set enrichment analysis
  fgsea_res = fgsea(pathways = gene_sets,
                    stats = ranks,
                    minSize = min_size,
                    maxSize = max_size)
  
  # Sort by significance and add column in CSV-friendly format
  fgsea_res = fgsea_res %>%
    arrange(padj) %>%
    mutate(leadingEdge_str = sapply(leadingEdge, paste, collapse = ', '))
  
  fgsea_res
}

# Summarise FGSEA results for logging
# Inputs: fgsea results df, alpha, label
# Output: Summary to console, invisible sig pathways df
summarise_fgsea = function(fgsea_res, alpha = config$alpha, label) {
  sig = filter(fgsea_res, padj < alpha)
  
  sig_up = filter(sig, NES > 0)
  sig_down = filter(sig, NES < 0)
  
  cat(sprintf('%s: %d / %d pathways significant (padj < %.2f)\n',
              label, nrow(sig), nrow(fgsea_res), alpha))
  cat(sprintf('Enriched (NES > 0): %d\n', nrow(sig_up)))
  cat(sprintf('Depleted (NES < 0): %d\n', nrow(sig_down)))
  
  if (nrow(sig) > 0) {
    top = head(sig, 5)
    cat(sprintf('Top %d:\n', nrow(top)))
    for (i in seq_len(nrow(top))) {
      cat(sprintf('%s (NES=%.2f, padj=%.3g, size=%d)\n',
                  top$pathway[i], top$NES[i], top$padj[i], top$size[i]))
    }
  }
  
  invisible(sig)
}

# Over-representation analysis
# Runs ORA (relies on significant genes - secondary to GSEA at our sample size)
# Inputs: - sig_genes (vector of sig gene IDs)
# - background (vector of all tested gene IDs)
# - gene_sets (named list of gene sets from KEGG or 1 GO ontology)
# - alpha
# Output: enrichResult object
run_ora = function(sig_genes, background, gene_sets, alpha = config$alpha) {
  
  # clusterProfiler needs term2gene & term2name dfs for enricher
  term2gene = data.frame(
    # Repeat each term for the length of the gene vector in that term
    term = rep(names(gene_sets), sapply(gene_sets, length)),
    # Flatten all gene vectors into a single vector to use as column
    gene = unlist(gene_sets, use.names = FALSE)
  )
  
  term2name = data.frame(
    term = names(gene_sets),
    name = sub('^[^:]+:\\s*', '', names(gene_sets))
  )
  
  # Filter gene universe to desired genes
  sig_mapped = sig_genes[sig_genes %in% term2gene$gene]
  bg_mapped = background[background %in% term2gene$gene]
  
  if (length(sig_mapped) < 3) {
    cat('Too few mapped genes for ORA\n')
    return(NULL)
  }
  
  cat(sprintf('  ORA: %d of %d sig genes mapped, background: %d\n',
              length(sig_mapped), length(sig_genes), length(bg_mapped)))
  
  clusterProfiler::enricher(
    gene = sig_mapped,
    universe = bg_mapped,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = alpha,
    pAdjustMethod = 'BH'
  )
}

# Summarise ORA results for logging
# Inputs: ORA enrichResult object, alpha, label
# Outputs: Summary to console, invisible significant pathways df
summarise_ora = function(ora_res, alpha = config$alpha, label) {
  if (is.null(ora_res)) {
    cat(sprintf('%s: No ORA results (too few mapped genes)\n', label))
    return(invisible(NULL))
  }
  
  res_df = as.data.frame(ora_res)
  sig = res_df %>% filter(p.adjust < alpha)
  
  cat(sprintf('%s: %d / %d terms significant (padj < %.2f)\n',
              label, nrow(sig), nrow(res_df), alpha))
  
  if (nrow(sig) > 0) {
    top = head(sig, 5)
    cat(sprintf('Top %d:\n', nrow(top)))
    for (i in seq_len(nrow(top))) {
      cat(sprintf('%s (padj=%.3g, %s/%s genes)\n',
                  top$Description[i], top$p.adjust[i],
                  top$Count[i], top$BgRatio[i]))
    }
  }
  
  invisible(sig)
}
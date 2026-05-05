# R/annotation_functions.R
# Functions for building annotated gene tables and merging with DESeq2 results

source('R/packages.R')

# Parse feature table into an annotation data frame
# Input: Path to feature_table.txt.gz
# Output: Annotated data frame of genes
parse_feature_table = function(ft_path) {
  ft = read_tsv(ft_path, show_col_types = FALSE, progress = FALSE)
  
  # Standardise column names
  names(ft) = gsub('^#\\s*', '', names(ft))
  names(ft) = gsub(' ', '_', names(ft))
  
  # Create gene table with descriptions & symbols
  genes = ft %>%
    filter(feature == 'gene', !is.na(GeneID)) %>%
    transmute(ncbi_gene_id = as.character(GeneID),
              symbol = symbol,
              description = name) %>%
    distinct(ncbi_gene_id, .keep_all = TRUE)
  
  # CDS (coding sequence) table with protein accessions (1 per gene)
  cds = ft %>%
    filter(feature == 'CDS', !is.na(product_accession), !is.na(GeneID)) %>%
    transmute(ncbi_gene_id = as.character(GeneID),
              protein_accession = product_accession) %>%
    distinct(ncbi_gene_id, .keep_all = TRUE)
  
  # Join dataframes
  annotated_df = genes %>%
    left_join(cds, by = 'ncbi_gene_id') %>%
    mutate(gene_id = if_else(!is.na(symbol),
                             symbol, paste0('LOC', ncbi_gene_id))) %>%
    relocate(gene_id)
  
  annotated_df
}

# Parse gene2go file & filter for P. vannamei
# Input: Path to gene2go.gz, taxid
# Output: Data frame with 1 row per gene matched with GO terms
parse_gene2go = function(g2g_path, taxid) {
  cat('Reading gene2go (filtering for taxid', taxid, ')...\n')
  
  # Read g2g in chunks of 500,000 to prevent zip-bombing your PC
  g2g = read_tsv_chunked(
    g2g_path,
    callback = DataFrameCallback$new(function(chunk, pos) {
      names(chunk) = gsub('^#\\s*', '', names(chunk))
      names(chunk) = trimws(names(chunk))
      chunk[chunk$tax_id == taxid, ]
    }),
    chunk_size = 500000,
    show_col_types = FALSE,
    progress = FALSE,
    na = c('-', ''))
  
  if (nrow(g2g) == 0) {
    cat('No GO annotations found for taxid', taxid)
    return(tibble(gene_id = character(), go_id = character(),
                  go_categories = character(), go_terms = character()))
  }
  
  g2g = g2g %>%
    mutate(ncbi_gene_id = as.character(GeneID)) %>%
    group_by(ncbi_gene_id) %>%
    summarise(go_id = paste(GO_ID, collapse = '; '),
              go_categories = paste(Category, collapse = '; '),
              go_terms = paste(GO_term, collapse = '; '),
              .groups = 'drop') %>%
    select(ncbi_gene_id, go_id, go_categories, go_terms)
  
  g2g
}

# Combine feature table & gene2go table into a single annotation lookup table
# Inputs: Paths to feature table and gene2go table, taxid
# Output: Data frame with both annotation sets
build_annotation_table = function(ft_path, g2g_path, taxid) {
  ft = parse_feature_table(ft_path)
  go = parse_gene2go(g2g_path, taxid)
  anno = left_join(ft, go, by = 'ncbi_gene_id')
  
  anno = anno %>%
    mutate(annotated = !is.na(description) &
             !grepl('^uncharacterized LOC|^hypothetical protein|^LOW QUALITY PROTEIN',
                    description, ignore.case = TRUE))
  
  cat(sprintf('Annotation table: %d genes total, %d with functional annotation (%.1f%%)\n',
              nrow(anno), sum(anno$annotated), 100 * mean(anno$annotated)))
  
  anno
}

# Summarise annotation coverage
# Inputs: Annotation table, character vector of gene IDs
# Output: Printed summary
summarise_annotation_coverage = function(annotation_df, gene_ids) {
  in_data = filter(annotation_df, gene_id %in% gene_ids)
  missing = setdiff(gene_ids, annotation_df$gene_id)
  has_go = sum(!is.na(in_data$go_id))
  
  cat(sprintf('Annotation coverage for %d genes:\n', length(gene_ids)))
  cat(sprintf('Matched to NCBI: %d (%.1f%%)\n',
              nrow(in_data), 100 * nrow(in_data) / length(gene_ids)))
  cat(sprintf('Functional name: %d (%.1f%%)\n', sum(in_data$annotated),
              100 * sum(in_data$annotated) / length(gene_ids)))
  cat(sprintf('Uncharacterised: %d\n', sum(!in_data$annotated)))
  cat(sprintf('Not in reference: %d\n', length(missing)))
  cat(sprintf('With GO terms: %d (%.1f%%)',
              has_go, 100 * has_go / length(gene_ids)))
  
  invisible(in_data)
}

# Join annotation information to any results data frame with a 'gene' column
# Inputs: Results df, annotation table
# Outputs: Results df with label column
annotate_results = function(results_df, annotation_df) {
  results_df %>%
    left_join(
      annotation_df %>% select(gene_id, symbol, description, annotated),
      by = c('gene' = 'gene_id')
    ) %>%
    mutate(label = case_when(
      annotated ~ description,
      !is.na(symbol) ~ symbol,
      TRUE ~ gene
    ))
}

# Build a gene_id -> display label named vector from an annotation table
# Inputs: Annotation table from build_annotation_table
# Output: Named character vector (names = gene_id, values = display label)
build_gene_labels = function(annotation_df) {
  labels = annotation_df %>%
    mutate(label = case_when(
      annotated ~ description,
      !is.na(symbol) & symbol != gene_id ~ symbol,
      TRUE ~ gene_id
    ))
  
  setNames(labels$label, labels$gene_id)
}

# Resolve gene IDs to display labels using a gene_labels vector
# Falls back to the raw gene ID when gene_labels is NULL or the gene is absent
# Inputs: Character vector of gene IDs, gene_labels named vector (or NULL)
# Output: Character vector of display labels
resolve_labels = function(gene_ids, gene_labels = NULL) {
  if (is.null(gene_labels)) return(gene_ids)
  ifelse(gene_ids %in% names(gene_labels), gene_labels[gene_ids], gene_ids)
}
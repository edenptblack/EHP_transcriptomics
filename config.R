#config.R

config = list(
  # Organism
  TAXID = 6689,
  ASSEMBLY_TAG = 'GCF_042767895.1_ASM4276789v1',
  G2G_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz',
  ref_level = 'NaCl',
  
  # Aesthetics
  # Project colour palette
  group_colours = c(
    NaCl   = "green",
    dsPTP2 = "blue",
    dsEHP  = "red",
    EHP    = "orange"
  ),
  
  # Group labels
  group_labels = c(
    NaCl   = "NaCl (control)",
    dsPTP2 = "dsRNA only",
    dsEHP  = "dsRNA + EHP",
    EHP    = "EHP only"
  ),
  
  # Project up/down regulation palette
  direction_colours = c(
    up = 'red',
    down = 'blue',
    NS = 'grey'
  ),
  
  # Group order for charts
  group_order = c('NaCl', 'dsPTP2', 'EHP', 'dsEHP'),
  
  # Plot settings
  plot_dpi = 300,
  plot_width = 8,
  plot_height = 5,
  
  # Filtering
  min_count = 10,
  min_samples = 3,
  
  # DESeq2
  alpha = 0.05,
  pca_n_top = 500,
  
  # Data paths
  raw_dir = 'data/raw',
  ref_dir = 'data/reference',
  proc_dir = 'data/processed',
  
  # Pairwise contrasts
  pairwise_contrasts = list(
    EHP_vs_NaCl     = c('group', 'EHP',    'NaCl'),
    dsEHP_vs_EHP    = c('group', 'dsEHP',  'EHP'),
    dsPTP2_vs_NaCl  = c('group', 'dsPTP2', 'NaCl'),
    dsEHP_vs_NaCl   = c('group', 'dsEHP',  'NaCl'),
    dsEHP_vs_dsPTP2 = c('group', 'dsEHP',  'dsPTP2')
  ),
  
  # GSEA contrasts vector
  gsea_contrasts = c(
    'interaction',
    'EHP_vs_NaCl',
    'dsPTP2_vs_NaCl',
    'dsEHP_vs_EHP'
  )
)

# GSEA contrast list (derived from vector)
config$gsea_stat_cols = setNames(
  paste0('stat_', config$gsea_contrasts),
  config$gsea_contrasts
)
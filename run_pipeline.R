# Run full analysis pipeline

source("scripts/00_data_preparation.R")
source("scripts/01_qc_exploratory_analysis.R")
source("scripts/02_build_subsets.R")
source("scripts/03_differential_expression.R")
source("scripts/04_collate_results.R")
source("scripts/05_pathway_enrichment.R")
source("scripts/06_candidate_pathways.R")
source("scripts/07_mechanism_interrogation.R")
# R/metadata.R

source('config.R')

# Metadata preparation for project-wide use
# Inputs: Vector of sample names, reference group for DESeq2
# Output: Data frame of metadata
build_metadata = function(sample_names, ref_level = config$ref_level) {
  groups = sub('_[0-9]+$', '', sample_names)
  
  metadata = data.frame(
    sample = sample_names,
    group = factor(groups),
    vaccinated = factor(ifelse(groups %in% c("dsEHP", "dsPTP2"), "yes", "no")),
    challenged = factor(ifelse(groups %in% c("dsEHP", "EHP"), "yes", "no")),
    row.names  = sample_names
  )
  
  metadata$group = relevel(metadata$group, ref = ref_level)
  metadata
}
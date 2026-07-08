# EHP transcriptomics - dsRNA protection in *Penaeus vannamei*

Bulk RNA-seq analysis investigating whether dsRNA (*EhPTP2*) protection against the microsporidian parasite *Ecytonucleospora hepatopenaei* (EHP) is driven by immune priming in the pacific white shrimp *P. vannamei*. 

## Description

A reproducible pipeline for a 2x2 factorial design pilot study (vaccine x challenge, n = 12) running from a featureCounts matrix to biological interpretation of differentially expressed genes and pathway analysis. Steps ran from biological QC and outlier handling -> differential expression analysis -> pathway enrichment analysis -> mechanism interpretation.

## Findings

The analysis identified 1,463 differentially expressed genes (DEGs) in challenged animals and 259 DEGs in the multifactor interaction term capturing vaccine-modified challenge responses. Identified immune genes made up only 3.9% of interaction term DEGs. The dominant transcriptomic signature was a systematic reversal of infection-induced changes in the interaction term, and pathway analysis revealed strong reversals towards baseline expression in endoplasmic reticulum (ER) stress, protein chaperone, cytoskeletal, and adhesion pathways. No hallmarks of invertebrate immune priming, such as epigenetic reprogramming, antimicrobial peptide (AMP) production, or pro-glycolytic metabolic shifts were detected, although the sample with near-complete clearance showed a prophenoloxidase (proPO) cascade signature. These findings are more consistent with a model of the treated shrimp undergoing less cellular stress due to a lower parasite burden than a model of an enhanced immune response and therefore support the hypothesis that dsEhPTP2 protection against EHP is mediated by RNAi.

## Stack & methods

R - renv - DESeq2 - ashr - tidyverse - clusterProfiler - fgsea - vsn - forcats - ggplot2 - RColorBrewer - GGally - ggrepel - ggVennDiagram - patchwork - magick - reshape2 - pheatmap

## Repository layout

    config.R          # all paths and analysis parameters
    run_pipeline.R    # runs the full analysis
    R/                # reusable functions
    scripts/          # numbered pipeline stages
    thesis_figures/   # figure-generation code not included in the main pipeline
    renv.lock         # locked package versions
    data/raw/         # (git-ignored) place input data here
    logs/, results/   # (generated) readable summaries and full outputs

## Reproducibility

1.  **Packages:** Install renv and run renv::restore() to install the exact locked versions from renv.lock. Some packages seem to interfere with this process; if necessary, all necessary packages are listed in packages.R

2.  **Data:** Create data/raw/ and add EHP_featurecounts.csv: first column Geneid (gene names from reference genome GCF_042767895.1), then 12 sample columns labelled dsEHP_1/2/3, dsPTP2_1/2/3, EHP_1/2/3, NaCl_1/2/3. The dataset is not yet published; if you wish to use the same data used in my thesis, please enquire with Dr. Jamie Bojko via J.bojko@tees.ac.uk

3.  **Run:** Exectute run_pipeline.R. The first run downloads and caches the *P. vannamei* reference genome and gene2go table, increasing run time; later runs reuse the cache. Key results are stored in logs/, full outputs are stored in results/

## Data availability

Raw sequencing data is currently withheld to preserve confidentiality.
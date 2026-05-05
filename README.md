This repository contains all code required to reproduce the results of my master's thesis, 'Transcriptomic evidence for reduced parasite burden rather than immune priming in Ecytonucleospora hepatopenaei-infected Penaeus vannamei treated with dsRNA', with all other elements removed to maintain data confidentiality.

Package preparation: All packages required can be found in 'packages.R' located in the root folder. To install these packages, install the renv package and use the command renv::restore, or install manually

Data preparation: From the root folder, create the filepath /data/raw. In this folder, place the raw featureCounts data from the study. The first column should be labelled 'Geneid', and contain strings of gene names from the P. vannamei reference genome GCF_042767895.1. The subsequent 12 columns should contain the sample data, and should be labelled dsEHP_1/2/3, dsPTP2_1/2/3, EHP_1/2/3, and NaCl_1/2/3. The file should be saved as EHP_featurecounts.csv

Once these steps are prepared, run the analysis using 'run_pipeline.R' in the root folder. Note that the first run will take longer due to the need to download the reference genome and the gene2go table; any subsequent runs will use the cached tables. Key results are saved in a readable format in the logs folder; detailed outputs are saved in the results folder.

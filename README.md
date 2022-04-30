# TCGA Glycogene Metaanalysis
Scripts and Pipelines used to analyze the TCGA for patterns in Glycogene Expression across Cancer Types

## Data Gathering Using TCGABiolinks:

The "TCGAbl_gather_data.R" gathers TCGA data if passed a valid TCGA code.  If in an HPC environment all datasets can be downloaded simultaneously using independent processes using a slurm script formatted in a manner similarly to "slurm_TCGAblRun" script.  The TCGAbl_gather_data.R script utilizes functions in the "TCGAbl_functions.R" file.

## Directories:

### glycoBio_data:

This directory contains the GlycoOnto.owl file, as well as R scripts used to create pathway and function lists used for the enrichment analyses.  Running the "make_pathLists.sh" script will regenerate the "pathwayList.rda" and "functionList.rda" lists.  The script runs the following in order:

1. glycoOnto_Graph.R will do the following to get the ontology into R:
  - Converts ontology into a JSON file using the "ROBOT" package
  - JSON ontology is read into R, and glycogene and pathway relationships are stored into a nested list structure.  This nested list structure is then saved for use in step 2.

2. create_pathLists.R will read the nested list structure to create path lists, where each pathway has a list of glycogenes which fit into each pathway classification.

The output from this pipeline will be two files: "pathwayList.rda" and "functionList.rda", which are used for enrichments after differential expression analyses.

### analyses:

This directory contains analysis scripts for each section of the manuscript:

- t-SNE embeddings
- Consensus Clustering
- Differential expression analysiss
- Enrichment analyses
- Breast Cancer Glyco50 classifier code
- Glyco50 Cancer Progression Maps

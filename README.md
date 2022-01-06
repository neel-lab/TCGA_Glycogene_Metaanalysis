# TCGA Glycogene Metaanalysis
Scripts and Pipelines used to analyze the TCGA for patterns in Glycogene Expression across Cancer Types

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
- Differential expression analysiss
- Enrichment analyses
- Breast Cancer Glyco30 classifier code

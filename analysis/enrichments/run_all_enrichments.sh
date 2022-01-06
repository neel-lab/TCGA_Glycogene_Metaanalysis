#!/bin/sh

#Go to root directory of repo:
cd ../../

#Script to run enrichments:
script=analysis/enrichments/TCGA_glycogeneEnrich_Fisher.R
files=$(ls ./analysis/DE_Analysis/*_DEData.rda)

for f in ${files}
do
	ctype=$(echo $f | rev | cut -d'/' -f1 | rev | cut -d'_' -f1)
	echo "Analysis of ${ctype}"
	Rscript ${script} $f
done

#Merge the results into a .tsv file:
mergeScript=analysis/enrichments/TCGA_mergeEnrich.R

Rscript ${mergeScript}

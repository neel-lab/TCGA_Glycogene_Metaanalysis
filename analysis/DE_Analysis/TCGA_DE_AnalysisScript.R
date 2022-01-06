library(dplyr);library(tidyr);library(tibble)
library(SummarizedExperiment)
library(edgeR)
library(limma);library(variancePartition)#dream package
library(parallel)

#Source Helper Functions:
source('./analysis/DE_Analysis/normUtils.R')
source('./analysis/DE_Analysis/DEModelRun.R')


#Function Input Variables:
load('./data_sets/reference/cancer_subtypes_mRNA.rda')

#Process inputs:
args=commandArgs(trailing=T)
inFile=args[1]
load(inFile)
#Parse relevant data:
RNASeq_Expt <- TCGAdata_list[['RNASeq']]
purity <- TCGAdata_list[['purity']]
#Clear large TCGA data objects:
rm(TCGAdata_list)
#Parse cancer-specific metadata:
ctype <- sub('(TCGA)_(\\D+)\\_.*','\\1-\\2',basename(inFile))
cancerSpecificSubtypes <- cancer_subtypes_mRNA[[ctype]]



# -------------RUN Differential Expression Analysis---------------

DEResults <- run_differential_expression(RNASeq_Expt,purity,cancerSpecificSubtypes)

if (!is.null(DEResults)){
	save_results(ctype,DEResults)
}

library(dplyr);library(tidyr);library(tibble)
library(SummarizedExperiment)
library(edgeR)
library(limma);library(variancePartition)#dream package
library(parallel)

dream_model <- function(e_dta,modelMatrix,fe,re){
	#Fits Gene linear models based on provided modelMatrix:
	# modelMatrix has random effects included in the design.
	
	#Create formula:
	fmla <- paste('~',paste(fe,collapse='+'),'+',paste(paste('(1|',re,')',sep=''),collapse='+'),sep='')
	#Filter out any rows that have any missing values for fixed or random
	# effects from the expression and model matrices:

	keepInds <- which(apply(modelMatrix[,c(fe,re)],1,function(x) all(!is.na(x))))
	e_dta <- e_dta[,keepInds];modelMatrix <- modelMatrix[keepInds,]
	
	param=(SnowParam(10,'SOCK',progressbar=T))
	vst <- voomWithDreamWeights(e_dta,fmla,modelMatrix,BPPARAM=param)
	m <- dream(vst,fmla,modelMatrix,BPPARAM=param)
	return(m)
}

save_results <- function(ctype,DEResults){
	#Create File Name:
	fname <- paste(ctype,'DEData.rda',sep='_')
	#Specify Location:
	fLoc <- file.path('/projects/academic/neel/tgroth/TCGA_Glycogene_Metaanalysis/analysis/DE_Analysis',fname)
	save(file=fLoc,DEResults)
}

run_differential_expression <- function(RNASeq_Expt,purity,cancerSpecificSubtypes){
	# Runs differential expression analyses across 
	# provided metadata for each cancer type:

	# Create model matrix:
	# If NULL returned, means there aren't enough 
	# Adjacent normal tissues to measure. 
	# Exit Program:
	modelMatrix <- create_model_matrix(RNASeq_Expt,purity,cancerSpecificSubtypes)
	if (is.null(modelMatrix)){
		return(NULL)
	}
	#Align RNASeq with modelMatrix:
	RNASeq_Expt<-RNASeq_Expt[,modelMatrix$barcode]

	# Create DGEList:
	e_dta <- DGEList(counts=assay(RNASeq_Expt))
	#Filtering of expression data:
	e_dta <- filter_and_norm(e_dta,modelMatrix)
	
	#Run differential expression analyses:

	# Create list of variables to test:
	# If there are no cancerSpecificSubtypes, filter out empties:
	#Fix cancerSpecificSubtypes:
	cancerSpecificSubtypes<-cancerSpecificSubtypes %>% gsub('[\\:|\\-]','_',.) %>% gsub('\\ ','.',.)
	fe_cols <- c('shortLetterCode',cancerSpecificSubtypes) %>% .[.!=""]

	#Filter the relevant random effects:
	re_cols <- colnames(modelMatrix)[colnames(modelMatrix) %in% c('age','gender','center')]
	DEResults <- lapply(X=fe_cols,FUN=function(fe){
		#Fit Model for FE variable:
		if ('CPE' %in% colnames(modelMatrix)){
			fe_vals <- c(fe,'CPE')
		} else {
			fe_vals <- fe
		}
		m <- dream_model(e_dta,modelMatrix,fe_vals,re_cols)
		#Extract effect and purity effects:

		# Main fixed effect comparisions (With respect to normal):
		mainNames=colnames(m$design)[!colnames(m$design) %in% c('(Intercept)','CPE')]
		#Get purity and main effects into two different 
		#Sometimes tumor purity may be missing:
	       	# If it's missing, make "de_CPE" NULL, otherwise get Purity effects:	
		if ('CPE' %in% fe_vals){
			de_CPE<-topTable(m,coef='CPE',n=Inf) %>% as.data.frame() %>% 
				mutate(geneName=rowData(RNASeq_Expt)[row.names(.),'external_gene_name'])
		} else {
			de_CPE<-NULL
		}
		#Parse all effects into a de_main list
		de_main <- lapply(mainNames,function(nme){
			tbl <- topTable(m,coef=nme,n=Inf) %>% as.data.frame() %>% 
			mutate(geneName=rowData(RNASeq_Expt)[row.names(.),'external_gene_name'])
			return(tbl)
		});names(de_main)<-mainNames
		outList <- list('DE_CPE'=de_CPE,
				'DE_Main'=de_main
				)
		return(outList)
	})
	names(DEResults) <- fe_cols
	return(DEResults)
}


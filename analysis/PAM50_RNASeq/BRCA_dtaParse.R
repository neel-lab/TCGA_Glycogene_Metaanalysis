library(dplyr);library(tidyr);library(tibble)
library(SummarizedExperiment)
library(edgeR)
library(limma);library(variancePartition)#dream package
library(parallel)


processFEData <- function(coldta,purity,cancerSpecificSubtypes){
	#Process fixed effects from RNASeq data:
	# -Converts non Tumor/Stroma columns to "Normal" for Normal Tissue

	colnames(coldta) <- sapply(colnames(coldta),function(x) gsub('[\\:|\\-]','_',x))
	#Extract fixed effect columns from ColData
	# for Differential Expression Analysis:
	fe_cols <- c('barcode','sample','shortLetterCode',cancerSpecificSubtypes)
	fe_cols <- sapply(fe_cols,function(x) gsub('[\\:|\\-]','_',x))
	#Fix column names if they are broken:
	fe_df <- coldta[,fe_cols] %>% as.data.frame()
	# Fix short letter codes (N: Normal , T: Tumor)
	fe_df <- fe_df %>% 
		mutate(
		       shortLetterCode=sapply(shortLetterCode,function(x){
				if (grepl('^T',x)){
					return('T')
				} else if (grepl('^N',x)){
					return('N')
				}
		       })
		) %>%
		mutate(shortLetterCode=as.factor(shortLetterCode))
	# Add purity data:
	if (nrow(purity)==0){
		message("No purity data, skipping purity assessment")
	} else {
		fe_df <- fe_df %>% left_join(purity,by=c('sample'='Sample.ID')) %>%
			select(-Cancer.type,-ESTIMATE,-ABSOLUTE,-LUMP,-IHC) %>% 
			mutate(CPE=as.numeric(sub('\\,','.',CPE)))
		#Fix purity:
		if (sum(is.na(fe_df$CPE))/nrow(fe_df) >=0.5){
			message('Not enough data for purity analysis, removing from model matrix')
			fe_df <- fe_df %>% select(-CPE)
		} else {
			fe_df <- fe_df %>% 
				mutate(CPE=ifelse(shortLetterCode=='N',0,CPE))
		}
	}
	# Fix subtype variables:
	fe_df <- fe_df %>% 
		mutate_at(vars(matches('paper_')),function(x){
			ifelse(.$shortLetterCode=='N','N',as.character(x))
		})

	return(fe_df)
}

processREData <- function(coldta){
	#Process random effects from RNASeq data:
	# 
	# List of Random Effects:
	# - Center ID
	# - Age
	# - Sex
	re_cols <- c('sample','age_at_index','gender')
	re_df <- coldta[,re_cols] %>% as.data.frame()
	# Final mutations on matrix:
	re_df <- re_df %>% mutate(
			center=sapply(sample,function(x) unlist(strsplit(split='\\-',x=x))[2]),
			age=case_when(

				age_at_index<25 ~ 1,
				age_at_index>25 & age_at_index<=50 ~ 2,
				age_at_index>50 & age_at_index<=75 ~ 3,
				age_at_index>75 ~ 4
			) %>% factor(.,levels=seq(1,4))
		) %>% select(-age_at_index)

	#Remove gender if all same sex:
	if (length(unique(re_df$gender[!is.na(re_df$gender)]))==1){
		#remove gender:
		re_df <- re_df %>% select(-gender)
	}
	return(re_df)
}

create_model_matrix <- function(RNASeq_Expt,purity,cancerSpecificSubtypes){
	colDta <- colData(RNASeq_Expt)
	#Gather Fixed and Random Effects:
	fe_dta <- processFEData(colDta,purity,cancerSpecificSubtypes)
	re_dta <- processREData(colDta)
	#Merge fixed and random effects together:
	total_dta <- fe_dta %>% left_join(re_dta,by='sample') %>% distinct()
	#Handling rank issues in RE:
	# Remove samples belonging to groups with low 
	# center and age representation:
	# Age:
	ageDist <- table(total_dta$age)
	removeAge <- names(ageDist)[which(ageDist==1)]
	total_dta <- total_dta[!(total_dta$age %in% removeAge),]
	# Center:
	centerDist <- table(total_dta$center)
	removeCenter <- names(centerDist)[which(centerDist==1)]
	total_dta <- total_dta[!(total_dta$center %in% removeCenter),]

	#Pre-run checks
	# Check if cancer type is eligible for DE analysis:
	#  Normals>=10
	if (sum(total_dta$shortLetterCode=='N')<10){
		message('Insufficient Normals, returning NULL')
		return(NULL)
	}
	# Make reference level for all fixed effects "N":
	#Fix names:
	cancerSpecificSubtypes<-cancerSpecificSubtypes %>% gsub('[\\:|\\-]','_',.) %>% gsub('\\ ','.',.)
	for (c in c('shortLetterCode',cancerSpecificSubtypes)){
		nonNa <- unique(total_dta[,c]) %>% .[!(is.na(.) | .=="NA")]
		total_dta[,c]=factor(total_dta[,c],levels=nonNa)
		total_dta[,c]=relevel(total_dta[,c],ref='N')
	}
	return(total_dta)
}


#Library Filtering & Normalization Methods:

gm_mean <- function(x, na.rm=TRUE){
	  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

median_sf <- function(e_dta){
	#Median Ratio Normalization:
	# Geometric mean of each gene:
	gene_GM <- apply(e_dta,1,gm_mean)
	sf <- apply(e_dta,2,function(x) median(x/gene_GM))
	return(sf)
}

filter_and_norm <- function(e_dta,modelMatrix){
	#Filter lowly expressed genes:
	# Use the Tumor and Normal classifications to filter genes
	#  If gene is low in both tumor and normal, remove it
	tn_design <- modelMatrix %>% pull(shortLetterCode)
	e_dta <- e_dta[filterByExpr(e_dta,group=tn_design),]
	# Compute size factors and assign in "e_dta":
	e_dta$samples$norm.factors <-median_sf(e_dta$counts)
	return(e_dta)
}

dream_vst <- function(e_dta,modelMatrix,fe,re){
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
	return(vst)
}

#Function Input Variables:
load('./data_sets/reference/cancer_subtypes_mRNA.rda')

#Process inputs:
load('./data_sets/TCGA_BRCA_dataList.rda')
#Parse relevant data:
RNASeq_Expt <- TCGAdata_list[['RNASeq']]
purity <- TCGAdata_list[['purity']]
#Clear large TCGA data objects:
rm(TCGAdata_list)
#Parse cancer-specific metadata:
ctype <- sub('(TCGA)_(\\D+)\\_.*','\\1-\\2',basename(inFile))
cancerSpecificSubtypes <- cancer_subtypes_mRNA[[ctype]]

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

cancerSpecificSubtypes<-cancerSpecificSubtypes %>% gsub('[\\:|\\-]','_',.) %>% gsub('\\ ','.',.)
fe_cols <- c('shortLetterCode',cancerSpecificSubtypes) %>% .[.!=""]

#Fit Model for FE variable:
if ('CPE' %in% colnames(modelMatrix)){
	fe_vals <- c(fe,'CPE')
} else {
	fe_vals <- fe
}
vst <- dream_vst(e_dta,modelMatrix,fe_vals,re_cols)

#Parse glycogene-specific info and pam50 info


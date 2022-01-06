library(dplyr);library(tidyr);library(ggplot2);library(Rtsne)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(variancePartition)#dream package

#--------- FUNCTIONS -------------
processFEData <- function(coldta,purity){
	#Process fixed effects from RNASeq data:
	# -Converts non Tumor/Stroma columns to "Normal" for Normal Tissue

	colnames(coldta) <- sapply(colnames(coldta),function(x) sub('[\\:|\\-]','_',x))
	#Extract fixed effect columns from ColData
	# for Differential Expression Analysis:
	fe_cols <- c('barcode','sample','shortLetterCode',cancerSpecificSubtypes)
	fe_cols <- sapply(fe_cols,function(x) sub('[\\:|\\-]','_',x))
	#Fix column names if they are broken:
	fe_df <- coldta[,fe_cols]
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
	if (length(unique(re_df$gender))==1){
		#remove gender:
		re_df <- re_df %>% select(-gender)
	}
	return(re_df)
}



gm_mean <- function(x, na.rm=TRUE){
	  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

median_sf <- function(e_dta){
	gene_GM <- apply(e_dta,1,gm_mean)
	sf <- apply(e_dta,2,function(x) median(x/gene_GM))
	return(sf)
}

get_sf <- function(e_dta,modelVec){
	#Get vector of main effect:
	# Handle different forms of modelVec:
	if (is.matrix(modelVec)){
		modelVec<-as.vector(modelVec[,colnames(modelVec)!='CPE'])
	}
	#Compute size factors:
	sf <- median_sf(e_dta$counts)
	return(sf)
}

normalize_dta <- function(e_dta,modelVec){
	#Compute size factors:
	e_dta$samples$norm.factors <- get_sf(e_dta,modelVec)
	#Variance stabilizing transformation:
	vst <- voom(e_dta,model.matrix(~.,data=modelVec))
	#vst <- voom(e_dta)
	return(vst)
}


run_tSNE <- function(RNASeq_Expt,metadta,glycogenes){
	#Parse data:
	e_dta <- DGEList(counts=assay(RNASeq_Expt))
	geneData<-rowData(RNASeq_Expt)
	glycogeneNames <- geneData %>% as.data.frame() %>% filter(external_gene_name %in% glycogenes) %>% pull(ensembl_gene_id)
	#Filter metadata:
	metadta<-metadta[apply(metadta,1,function(x) all(!is.na(x))),]
	e_dta<-e_dta[,metadta$barcode]
	metadta_filt<-metadta[,colnames(metadta)!='barcode']
	#Normalize data:
	e_dta <- normalize_dta(e_dta,metadta_filt)
	#Transpose:
	e_dta <- t(e_dta$E)
	#Filter:
	e_dta <- e_dta[,colnames(e_dta) %in% glycogeneNames]
	#Run tSNE:
	tsne <- Rtsne(X=e_dta,dims=2,initial_dims=150,perplexity=nrow(e_dta)/10,max_iter=2000,verbose=T,num_cores=8,check_duplicates=F)	
	#Make dataFrame:
	mainEffect=metadta_filt[,colnames(metadta_filt)!='CPE']
	tsne_df <- data.frame(tSNE1=tsne$Y[,1],tSNE2=tsne$Y[,2],label=mainEffect)
	return(tsne_df)
}



#Script Inputs:

#Function Input Variables:
load('./data_sets/reference/cancer_subtypes_mRNA.rda')
args=commandArgs(trailing=T)
inFile=args[1]
load(inFile)
RNASeq_Expt<-TCGAdata_list$RNASeq
coldata<-colData(RNASeq_Expt)
purity<-TCGAdata_list$purity
#Parse cancer-specific metadata:
ctype <- sub('(TCGA)_(\\D+)\\_.*','\\1-\\2',basename(inFile))
cancerSpecificSubtypes <- cancer_subtypes_mRNA[[ctype]]
gps <- c('shortLetterCode',cancerSpecificSubtypes) %>% .[.!=""] %>% 
	sub('\\ ','.',.) %>% sub('[\\:|\\-]','_',.)
rm(TCGAdata_list)

#For all the "cancerSpecificSubtypes" and "shortLetterCode", make tSNE data:
#Get metadata:
classDta <- processFEData(coldata,purity) 
reDta<-processREData(coldata)
totalClassDta<-merge(classDta,reDta,by='sample')
tSNE_dta <- lapply(gps,function(g){
	message(g)
	#Gather relevant metadata:
	if ('CPE' %in% colnames(classDta)){
		mdta <- totalClassDta[,c('barcode',g,'CPE')]
	} else {
		mdta <- totalClassDta[,c('barcode',g)]
	}
	#Check if homogeneous:
	# All tumor samples, mostly
	if (length(unique(mdta[,g]))==1){
		message(paste(g,'is homogeneous, returning NULL for this data'))
		return(NULL)
	}
	#Run through tSNE algorithm:
	res<-run_tSNE(RNASeq_Expt,mdta,glycogenes)
	return(res)
});names(tSNE_dta)<-gps

save(file=file.path('/projects/academic/neel/tgroth/TCGA_Glycogene_Metaanalysis/analysis/tSNE_figures',paste(ctype,'tSNE_dta.rda',sep='_')),tSNE_dta)

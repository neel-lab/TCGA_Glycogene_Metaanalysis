library(dplyr);library(tidyr);library(ggplot2);library(Rtsne)
library(SummarizedExperiment)
library(edgeR)
library(limma)


#--------- FUNCTIONS -------------

gm_mean <- function(x, na.rm=TRUE){
	  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

median_sf <- function(e_dta){
	gene_GM <- apply(e_dta,1,gm_mean)
	sf <- apply(e_dta,2,function(x) median(x/gene_GM))
	return(sf)
}

get_sf <- function(e_dta,modelVec){
	e_dta <- e_dta[filterByExpr(e_dta,group=modelVec),]
	#Compute size factors:
	sf <- median_sf(e_dta$counts)
	return(sf)
}

normalize_dta <- function(e_dta,modelVec){
	#Compute size factors:
	e_dta$samples$norm.factors <- get_sf(e_dta,modelVec)
	#Filter out samples with abnormally low or high sfs:
	e_dta_filt<-e_dta[,e_dta$samples$norm.factors<2 & e_dta$samples$norm.factors>0.1]
	#Variance stabilizing transformation:
	modelVec_filt<-modelVec[e_dta$samples$norm.factors<2 & e_dta$samples$norm.factors>0.1]
	mm<-model.matrix(~.,as.data.frame(modelVec_filt))
	vst <- voom(e_dta_filt,mm)
	return(vst)
}


agg_dta <- function(glycogenes){
	#Get glycogenes:
	load(file.path('data_sets','TCGA_READ_dataList.rda'))
	geneData<-rowData(TCGAdata_list$RNASeq)
	rm(TCGAdata_list)
	glycogeneNames <- geneData %>% as.data.frame() %>%
		filter(external_gene_name %in% glycogenes) %>% pull(ensembl_gene_id)
	
	#Parse all glycogene data from tumor TCGA tissues:
	totalDta <- do.call(cbind,lapply(list.files('./data_sets',pattern='*.rda'),function(x){
			load(file.path('data_sets',x))
			rnaseq <- TCGAdata_list$RNASeq
			tumorInds <- which(colData(rnaseq)$shortLetterCode!="NT")
			e_dta <- assay(rnaseq)[,tumorInds]
			ctype_vec <- rep(sub('\\_dataList\\.rda','',basename(x)),ncol(e_dta))
			e_dta<-rbind(e_dta,'ctype'=ctype_vec)
			return(e_dta)
		}))
	ctype_vec<-totalDta['ctype',]
	e_dta <-totalDta[row.names(totalDta)!='ctype',]
	geneNames<-row.names(e_dta)
	e_dta<-apply(e_dta,2,as.numeric);row.names(e_dta)<-geneNames
	e_dta <- DGEList(counts=e_dta)
	normDta <- normalize_dta(e_dta,ctype_vec)
	save(file='TCGA_AllNormData.rda',normDta)
	#Transpose:
	dta <- t(normDta$E)
	dta=dta[,colnames(dta) %in% glycogeneNames]
	return(list('dta'=dta,'ctype_vec'=ctype_vec))
}

#Running t-SNE on All Tumor Types:

outDta <- agg_dta(glycogenes)
allDta <- outDta[['dta']] %>% distinct()
ctype_vec <- outDta[['ctype_vec']]

all_tSNE <- Rtsne(X=allDta,dims=2,initial_dims=150,perplexity=75,max_iter=2000,verbose=T,num_threads=8,check_duplicates=F)

#Construct DF:
all_tSNE_df <- data.frame(tSNE1=all_tSNE$Y[,1],tSNE2=all_tSNE$Y[,2],ctype=ctype_vec)

#Save dataframe:
save(file='/user/tgroth/neel/tgroth/TCGA_Glycogene_Metaanalysis/analysis/tSNE_figures/all_tSNE_df.rda',all_tSNE_df)

library(dplyr)
library(parallel)
library(rlang)
library(TCGAbiolinks)
library(cli)
library(BiocGenerics)
library(GenomicRanges)
library(utf8)
library(SummarizedExperiment)

source('TCGAbl_functions.R')

args=commandArgs(trailing=T)
proj_id=args[1]

#Construct queries for TCGA data:
# Get RNA-seq, Methylation, Copy Number Variation, miRNA-seq, and Clinical metadata

#Get molecular Subtypes for this cancer type:
subtypes=PanCancerAtlas_subtypes() %>% filter(cancer.type==sub('TCGA-','',proj_id))

#Parsing Function List:
f_list<-c(get_tx,get_miRNA,get_mth,get_cnv,get_clinic)

#Get all data based on project ID and store in a list
TCGAdata_list<- mclapply(f_list,function(f) f(proj_id),mc.cores=8)
names(TCGAdata_list)<-c('RNASeq','miRNASeq','meth_seq','cnv_dta','clinical')
#Adding subtypes and purity data:
TCGAdata_list[['subtypes']]<-subtypes
TCGAdata_list[['purity']]<-Tumor.purity %>%
	dplyr::filter(Cancer.type==unlist(strsplit(split='-',x=proj_id))[2])

#Save list of downloaded data:
#Objects are called "TCGAdata_list"
fname_save=paste(sub('-','_',proj_id),'dataList.rda',sep='_')
save(file=file.path('data_sets',fname_save),TCGAdata_list)


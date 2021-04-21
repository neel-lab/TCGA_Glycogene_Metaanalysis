library(dplyr)
library(parallel)
library(rlang)
library(TCGAbiolinks)
library(cli)
library(BiocGenerics)
library(GenomicRanges)
library(utf8)
library(SummarizedExperiment)


args=commandArgs(trailing=T)
proj_id=args[1]

#Construct queries for TCGA data:
# Get RNA-seq, Methylation, Copy Number Variation, miRNA-seq, and Clinical metadata

# Additionally, get data regarding tumor purity

#Transcriptome Profiling Parse:
get_tx<-function(project){
	query<-GDCquery(project,
			data.category='Transcriptome Profiling',
			data.type='Gene Expression Quantification',
			workflow.type='HTSeq - Counts')
	GDCdownload(query)
	dta<-GDCprepare(query=query)
	return(dta)
}

get_miRNA<-function(project){
	query<-GDCquery(project,
			data.category='Transcriptome Profiling',
			data.type='miRNA Expression Quantification')
	GDCdownload(query)
	dta<-GDCprepare(query=query)
	return(dta)
}

get_mth<-function(project){
	query<-GDCquery(project,
			data.category='DNA Methylation',
			platform='Illumina Human Methylation 450')
	#GDCdownload(query)
	tryCatch(GDCdownload(query,files.per.chunk = 10),
		     error = function(e) GDCdownload(query,files.per.chunk=5))
	dta<-GDCprepare(query=query)
	return(dta)
}

get_cnv<-function(project){
	query<-GDCquery(project,
			data.category='Copy Number Variation',
			data.type='Copy Number Segment')
	GDCdownload(query)
	dta<-GDCprepare(query=query)
	return(dta)
}

#get_ATAC<-function(project){
#	#Go to ATAC-Seq Directory and gather data:
#	ctype_name=unlist(strsplit(split='TCGA-',x=project))[1]
#	fname=grep(project,file.path('ATACSeq_dta'),value=T)
#	if (is.null(fname)){
#		return()
#	else {
#		return(read.table(file=file.path('ATACSeq_dta',fname),sep='\t',header=T))
#	}
#}


get_clinic<-function(project){
	query<-GDCquery(project,
			data.category='Clinical',
			data.type='Clinical Supplement',
			data.format='BCR Biotab'
	)
	GDCdownload(query)
	dta<-GDCprepare(query=query)
	return(dta)
}
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

library(dplyr);library(SummarizedExperiment);library(tidyr)
library(tibble)

#Specify input file:
#args=commandArgs(trailing=T)
#inFile=args[1]
#load(inFile)


#Utility functions:

get_patient_id<-function(barcodeVec){
	return(sapply(barcodeVec,function(x) sub('(TCGA-\\w{2}-\\w{4}).*','\\1',x)))
}

#miRNA data processing:
miRNA_proc<-function(TCGAdata_list,inFile){
	#Fix column labels:
	miRNA_dta<-TCGAdata_list$miRNA %>%
		select(-contains('cross-mapped')) %>%
		pivot_longer(colnames(.)[colnames(.)!='miRNA_ID']) %>%
		mutate(
			dtaType=sapply(name,function(x) unlist(strsplit(split='_TCGA',x=x))[1]),
			sampleID=sapply(name,function(x) unlist(strsplit(split='(read_count_|reads_per_million_miRNA_mapped_)',x=x))[2])
		      ) %>% 
		select(sampleID,dtaType,miRNA_ID,value) %>% mutate(row=row_number()) %>% 
		pivot_wider(id_cols=c('sampleID','miRNA_ID'),names_from='dtaType',values_from='value',values_fn=sum)
	#Create read count and reads per million count tables [miRNA x patient]
	miRNA_reads<-miRNA_dta %>%
		select(sampleID,miRNA_ID,read_count) %>%
		mutate(row=row_number()) %>%
		pivot_wider(names_from='sampleID',values_from='read_count',id_cols='miRNA_ID') %>% 
		column_to_rownames('miRNA_ID')
	miRNA_rpm<-miRNA_dta %>%
		select(sampleID,miRNA_ID,reads_per_million_miRNA_mapped) %>%
		mutate(row=row_number()) %>%
		pivot_wider(names_from='sampleID',values_from='reads_per_million_miRNA_mapped',id_cols='miRNA_ID') %>%
		column_to_rownames('miRNA_ID')

	#Parse relevant clinical metadata for patients:
	clinic_name<-paste('clinical_patient_',tolower(unlist(strsplit(split='\\_',x=inFile))[2]),sep='')
	print(clinic_name)
	miRNA_colData<-TCGAdata_list$clinical_data[[1]][[clinic_name]]
	print(miRNA_colData)
	#Merge sample info with clinical information:
	miRNA_colData<-data.frame(sampleID=colnames(miRNA_reads)[grepl('^TCGA',colnames(miRNA_reads))]) %>%
		mutate(bcr_patient_barcode=get_patient_id(sampleID)) %>% 
		left_join(miRNA_colData,by='bcr_patient_barcode') %>% 
		column_to_rownames('sampleID')
	#Create SummarizedExperiment
	miRNA_expt<-SummarizedExperiment(assays=list(miRNA_reads,miRNA_rpm),colData=miRNA_colData)
	return(miRNA_expt)
}





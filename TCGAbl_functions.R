library(TCGAbiolinks)

#Transcriptome Profiling Parse:
get_tx<-function(project){
	query<-GDCquery(project,
			data.category='Transcriptome Profiling',
			data.type='Gene Expression Quantification',
			workflow.type='HTSeq - Counts')
	GDCdownload(query,directory='./GDCdata')
	dta<-GDCprepare(query=query)
	return(dta)
}

#MicroRNA Parse
get_miRNA<-function(project){
	query<-GDCquery(project,
			data.category='Transcriptome Profiling',
			data.type='miRNA Expression Quantification')
	GDCdownload(query,directory='./GDCdata')
	dta<-GDCprepare(query=query)
	return(dta)
}

#Methylation Parse:
get_mth<-function(project){
	query<-GDCquery(project,
			data.category='DNA Methylation',
			platform='Illumina Human Methylation 450')
	#GDCdownload(query,directory='GDCdata')
	tryCatch(GDCdownload(query,files.per.chunk = 10,directory='./GDCdata'),
		     error = function(e) GDCdownload(query,files.per.chunk=5,directory='./GDCdata'))
	dta<-GDCprepare(query=query)
	return(dta)
}

#Copy Number Variation Parse:
get_cnv<-function(project){
	query<-GDCquery(project,
			data.category='Copy Number Variation',
			data.type='Copy Number Segment')
	GDCdownload(query,directory='./GDCdata')
	dta<-GDCprepare(query=query)
	return(dta)
}

#Clinical data parse
get_clinic<-function(project){
	query<-GDCquery(project,
			data.category='Clinical',
			data.type='Clinical Supplement',
			data.format='BCR Biotab'
	)
	GDCdownload(query,directory='./GDCdata')
	dta<-GDCprepare(query=query)
	return(dta)
}

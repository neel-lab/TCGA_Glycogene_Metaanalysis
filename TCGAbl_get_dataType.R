source('TCGAbl_functions.R')

args=commandArgs(trailing=T)
proj_id=args[1]
dta_type=args[2]

f_list<-list(
'get_tx'=get_tx,
'get_miRNA'=get_miRNA,
'get_mth'=get_mth,
'get_cnv'=get_cnv,
'get_clinic'=get_clinic
)

if (!(dta_type) %in% names(f_list)){
	message('Invalid datatype, must select from the following:\n')
	message(paste(names(f_list),collapse='\n'))
	stop()
}

fileOutName<-paste(sub('-','_',proj_id),dta_type,'_result.rda',sep='_')

#Run getter:
reparsedData<-f_list[[dta_type]](proj_id)

#Save parsed result:
save(file=fileOutName,reparsedData)

library(dplyr);library(tidyr);library(tibble)
library(parallel)

#List of enrichment data:
enrichFiles<-list.files('./analysis/enrichments/',pattern='enrichments.rda')

#Merge everything:
allEnrichData<-lapply(enrichFiles,function(f){
	load(file.path('./analysis/enrichments/',f))
	return(enrichments)
})
names(allEnrichData)<-sapply(enrichFiles,function(f){
	ctype=unlist(strsplit(split='_',x=f))[1]
	ctype=sub('\\-','_',ctype)
	return(ctype)
})
	
#Aggregator:

aggDta<-function(enrichList,pathType,universe,level,direction){
	dta<-do.call(rbind,lapply(names(enrichList),function(ct){
		df<-enrichList[[ct]][[pathType]][[universe]][[level]][[direction]] %>%
	      mutate(ctype=ct,pathType=pathType,universe=universe,coefType=level,direction=direction)
		return(df)
		}))      
	return(dta)
}

#Aggregate data:
allGroups<-expand.grid(pathType=c('Pathways','Functions'),
		       universe=c('genome','glycogene'),
		       direction=c('up','down'),
		       level=c('tumorVnormal','purity')
		       )

allEnrichDF<-do.call(rbind,apply(allGroups,1,function(x){
		allDta<-do.call(rbind,lapply(names(allEnrichData),function(n){
			return(aggDta(allEnrichData,x['pathType'],x['universe'],x['level'],x['direction']))
}))
		return(allDta)
})) %>% distinct()

#Write out enrichments:
write.table(file='./analysis/enrichments/allEnrichData.tsv',
	    sep='\t',
	    quote=F,
	    row.names=F,
	    allEnrichDF)

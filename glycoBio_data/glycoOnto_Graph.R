library(rjson)
library(dplyr)

#Run "robot" command to convert glycoOnto.owl to glycoOnto.json
#system('robot convert --input glycoOnto/glycoOnto_individuals_reactions.rdf --output glycoOnto.json')
system('robot convert --input GlycoOntology/GlycoOnto.owl --output GlycoOnto.json')

#Read graph as ontology:
# Element 1, edges property is graph:
glycoPathwayObj <- fromJSON(file='./GlycoOnto.json')$graphs[[1]]$edges

glycoPathwayObj  <- lapply(glycoPathwayObj,function(elt){
	elt <- lapply(elt,function(e){
		return(sub('.*\\#(\\w.+)$','\\1',e))
		})	
	names(elt) <- c('sub','pred','obj')
	return(elt)
})

#Split pathway relationships and glycogene relationships:
#Parse for GlycoOnto-specific associations, not GO relationships:
pathwayRels<-glycoPathwayObj[sapply(glycoPathwayObj,function(x) x$pred=='is_a' & !grepl('\\/GO\\_.+$',x$obj))]

#Find all glycogenes under the "enzyme" class:
glycogenes<-glycoPathwayObj[sapply(glycoPathwayObj,function(x) grepl('enzyme',x$obj))]
#Gene name list:
glycogenes<-sapply(glycogenes,function(x) x$sub)
#Getting pathway associations:
glycogenes<-glycoPathwayObj[sapply(glycoPathwayObj,function(x) x$sub %in% glycogenes & x$pred=='type' & !grepl('GO\\_.+?$',x$obj) & x$obj!='enzyme')]

glycoPathwayObj <- list('pathways'=pathwayRels,'glycogenes'=glycogenes)

save(file='glycoPathwayObj.rda',glycoPathwayObj)

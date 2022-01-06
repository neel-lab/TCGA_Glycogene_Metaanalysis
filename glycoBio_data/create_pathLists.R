library(dplyr);library(tidyr)
load('./glycoPathwayObj.rda')

#Functions to traverse the ontology:
get_gg_rels<-function(geneName,glycogeneObj){
	funs <- glycogeneObj[sapply(glycogeneObj,function(x) x$sub==geneName)]
	return(sapply(funs,function(x) x$obj))
}

get_init_pathway <- function(geneName,ggPathObj){
	pths <- ggPathObj[sapply(ggPathObj,function(x) x$sub==geneName)]
	return(sapply(pths,function(x) x$obj))
}

traverse_graph <- function(term,pathwayRelObj,pathway_list){
	#Find parent term:	
	parentInd <- which(sapply(pathwayRelObj,function(x) x$sub==term))
	if (length(parentInd)==0){
		return(pathway_list)
	}
	newTerm <- pathwayRelObj[[parentInd]]$obj
	pathway_list <- c(pathway_list,newTerm)
	#Call again:
	traverse_graph(newTerm,pathwayRelObj,pathway_list)
}

find_related_pathways <- function(geneName,ggene_rels,pathway_rels){
	#Get pathways that glycogene is in:
	pths<-get_gg_rels(geneName,ggene_rels)
	#Find the rest of the related pathways:
	totalPaths <- c()
	for (p in pths){
		pathChain<-c(p,traverse_graph(p,pathway_rels,c()))
		totalPaths<-c(totalPaths,pathChain)
	}
	#Remove potential duplicates:
	totalPaths <- unique(totalPaths)
	return(totalPaths)
}


makePathLists<-function(glycoPathwayObj){
	#Parse out the different lists:
	pathwayRelList <- glycoPathwayObj$pathways
	glycogenes <- glycoPathwayObj$glycogenes
	#Get all Unique Glycogenes:
	glycogeneNames <- unique(sapply(glycogenes,function(x) x$sub))
	#Initialize pathway and function lists:
	pathwayList=list()
	functionList=list()
	for (g in glycogeneNames){
		# Find related Pathways and Functions:
		relClasses<-find_related_pathways(g,glycogenes,pathwayRelList)
		relPathways<-relClasses[which(sapply(relClasses,function(x) 'Pathway' %in% traverse_graph(x,pathwayRelList,c())))]
		relFunctions<-relClasses[which(sapply(relClasses,function(x) 'Glycogene' %in% traverse_graph(x,pathwayRelList,c())))]
		#relPathways <- find_related_pathways(g,pathwayRelList,get_init_pathway)
		#relFunctions <- find_related_pathways(g,pathwayRelList,get_gg_rels)
		# Add to lists if names not there:
		for (p in relPathways){
			if (!(p %in% names(pathwayList))){
				pathwayList[[p]]<-c(g)
			} else {
				pathwayList[[p]]<-c(pathwayList[[p]],g)
			}
		}
		for (f in relFunctions){
			if (!(f %in% names(functionList))){
				functionList[[f]]<-c(g)
			} else {
				functionList[[f]] <-c(functionList[[f]],g)
			}
		}
	}
	outList<-list('pathwayList'=pathwayList,'functionList'=functionList,'glycogenes'=glycogeneNames)
	return(outList)
}

outList<-makePathLists(glycoPathwayObj)
glycogenes<-outList$glycogenes;save(file='glycogenes.rda',glycogenes)
pathwayList<-outList$pathwayList;save(file='pathwayList.rda',pathwayList)
functionList<-outList$functionList;save(file='functionList.rda',functionList)

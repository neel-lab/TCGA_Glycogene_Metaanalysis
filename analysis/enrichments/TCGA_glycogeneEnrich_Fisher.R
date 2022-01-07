library(dplyr);library(tidyr);library(tibble)
library(parallel)

#--------Enrichment Functions ------------

get_tumorVnormal<-function(DEResults){
	#Parse the DE data List:
	return(DEResults$shortLetterCode$DE_Main[[1]])
}


get_tumorPurityEffects<-function(DEResults){
	#Parse the DE data List:
	return(DEResults$shortLetterCode$DE_purity)
}

get_glycogene_dedata<-function(dedata,glycogenes){
	return(dedata %>% filter(geneName %in% glycogenes))
}

numGenes<-function(dedata,logFC_thresh,adjP_thresh,pathList,pass=T,inPath=T,high=T){
	#Parses the number of genes that pass thresholds and belong to pathway list:
	#dplyr pipes:
	n_genes<-dedata %>% filter(if (inPath) geneName %in% pathList else !(geneName %in% pathList)) %>% 
	filter(case_when(
		 # Passes a high threshold:
		 pass & high ~ (logFC>logFC_thresh & adj.P.Val<=adjP_thresh),
		 # Passes a low threshold
		 pass & !high ~ (logFC<logFC_thresh & adj.P.Val<=adjP_thresh),
		 # NOT Passes a high threshold:
		 !pass & high ~ (logFC<=logFC_thresh | adj.P.Val>adjP_thresh),
		 # NOT Passes a low threshold:
		 !pass & !high ~ (logFC>=logFC_thresh | adj.P.Val>adjP_thresh)
	 )) %>% 
		 #Quantify
		 nrow()
	return(n_genes)
}


numGenes_bin<-function(dedata,logFC_ecdf,adjP_thresh,pathList,pass=T,inPath=T,high=T){
	#Binary version of numGenes:
	# Used for processing purity model coefficients
	#Parses the number of genes that pass thresholds and belong to pathway list:
	#dplyr pipes:
	ecdf_thresh<-ifelse(high,0.975,0.025)

	n_genes<-dedata %>% filter(if (inPath) geneName %in% pathList else !(geneName %in% pathList)) %>% 
	filter(case_when(
		 # Passes a high threshold:
		 pass & high ~ (logFC_ecdf(logFC)>ecdf_thresh & adj.P.Val<=adjP_thresh),
		 # Passes a low threshold
		 pass & !high ~ (logFC_ecdf(logFC)<ecdf_thresh & adj.P.Val<=adjP_thresh),
		 # NOT Passes a high threshold:
		 !pass & high ~ (logFC_ecdf(logFC)<=ecdf_thresh | adj.P.Val>adjP_thresh),
		 # NOT Passes a low threshold:
		 !pass & !high ~ (logFC_ecdf(logFC)>=ecdf_thresh | adj.P.Val>adjP_thresh)
	 )) %>% 
		 #Quantify
		 nrow()
	return(n_genes)
}


make_contingency<-function(dedata,logFC_thresh,adjP_thresh,pathList,high){
	#Matrix:
	cMat<-matrix(rep(0,4),nrow=2)
	#cMat[1,1]: Genes IN path AND pass thresholds
	cMat[1,1]<-numGenes(dedata,logFC_thresh,adjP_thresh,pathList,pass=T,inPath=T,high)
	#cMat[1,2]: Gene IN path AND do NOT pass thresholds:
	cMat[1,2]<-numGenes(dedata,logFC_thresh,adjP_thresh,pathList,pass=F,inPath=T,high)
	#cMat[2,1]: Genes NOT IN path AND pass thresholds
	cMat[2,1]<-numGenes(dedata,logFC_thresh,adjP_thresh,pathList,pass=T,inPath=F,high)
	#cMat[1,2]: Gene IN path AND do NOT pass thresholds:
	cMat[2,2]<-numGenes(dedata,logFC_thresh,adjP_thresh,pathList,pass=F,inPath=F,high)
	#Return the contingency matrix:
	return(cMat)
}

make_contingency_bin<-function(dedata,logFC_ecdf,adjP_thresh,pathList,high){
	#Binary version of contingency matrix
	#Matrix:
	cMat<-matrix(rep(0,4),nrow=2)
	#cMat[1,1]: Genes IN path AND pass thresholds
	cMat[1,1]<-numGenes_bin(dedata,logFC_ecdf,adjP_thresh,pathList,pass=T,inPath=T,high)
	#cMat[1,2]: Gene IN path AND do NOT pass thresholds:
	cMat[1,2]<-numGenes_bin(dedata,logFC_ecdf,adjP_thresh,pathList,pass=F,inPath=T,high)
	#cMat[2,1]: Genes NOT IN path AND pass thresholds
	cMat[2,1]<-numGenes_bin(dedata,logFC_ecdf,adjP_thresh,pathList,pass=T,inPath=F,high)
	#cMat[1,2]: Gene IN path AND do NOT pass thresholds:
	cMat[2,2]<-numGenes_bin(dedata,logFC_ecdf,adjP_thresh,pathList,pass=F,inPath=F,high)
	#Return the contingency matrix:
	return(cMat)
}


run_test<-function(cMat){
	#Tests for overrepresentation of geneList:
	return(fisher.test(cMat,alternative='greater')$p.value)
}

test_pathways<-function(dedata,logFC_thresh,adjP_thresh,pathways_list,high){
	#Tests every path in pathways_list:

	enrichmentTests<-mclapply(X=pathways_list,FUN=function(p){
			#1. Make Contingency Matrix:
			# high argument means is LFC greater than or equal to thresh.
			# True tests for upregulation, False tests for downregulation.
			cMat<-make_contingency(dedata,logFC_thresh,adjP_thresh,p,high)
			#2. Test
			pval<-run_test(cMat)
			return(pval)
	});names(enrichmentTests)<-names(pathways_list)
	enrichmentTests<-data.frame(pathways=names(pathways_list),pvals=unlist(enrichmentTests))
	#Correct p-values for false discovery:
	enrichmentTests$adjP<-p.adjust(enrichmentTests$pvals,method='BH')
	return(enrichmentTests)
}


test_pathways_bin<-function(dedata,logFC_ecdf,adjP_thresh,pathways_list,high){
	#Tests every path in pathways_list:

	enrichmentTests<-mclapply(X=pathways_list,FUN=function(p){
			#1. Make Contingency Matrix:
			# high argument means is LFC greater than or equal to thresh.
			# True tests for upregulation, False tests for downregulation.
			cMat<-make_contingency_bin(dedata,logFC_ecdf,adjP_thresh,p,high)
			#2. Test
			pval<-run_test(cMat)
			return(pval)
	});names(enrichmentTests)<-names(pathways_list)
	enrichmentTests<-data.frame(pathways=names(pathways_list),pvals=unlist(enrichmentTests))
	#Correct p-values for false discovery:
	enrichmentTests$adjP<-p.adjust(enrichmentTests$pvals,method='BH')
	return(enrichmentTests)
}

#Enricher Main:

do_enrich<-function(dedata,glycogenes,geneSet){
	#Define logFC and padj thresholds:
	logFC_lowthresh=log(1/2)
	logFC_highthresh=log(2)
	p.adj.thresh=0.05
	#Parse glycogene Data:
	glycogene_dedata<-dedata %>% dplyr::filter(geneName %in% glycogenes)
	# Gene Set Enrichments:
	# - Entire Genome:
	pathwayEnrich_wholeGenome_upreg<-test_pathways(dedata,logFC_highthresh,p.adj.thresh,geneSet,high=T)
	pathwayEnrich_wholeGenome_downreg<-test_pathways(dedata,logFC_lowthresh,p.adj.thresh,geneSet,high=F)
	# - Glycogene-Centered:
	pathwayEnrich_glycogeneCentered_upreg<-test_pathways(glycogene_dedata,logFC_highthresh,p.adj.thresh,geneSet,high=T)
	pathwayEnrich_glycogeneCentered_downreg<-test_pathways(glycogene_dedata,logFC_lowthresh,p.adj.thresh,geneSet,high=F)
	
	#Organize into list:
	dta_out<-list(
		      'genome'=list('up'=pathwayEnrich_wholeGenome_upreg,'down'=pathwayEnrich_wholeGenome_downreg),
		      'glycogene'=list('up'=pathwayEnrich_glycogeneCentered_upreg,'down'=pathwayEnrich_glycogeneCentered_downreg)
	      )
	return(dta_out)
}


do_enrich_purity<-function(dedata,glycogenes,geneSet){
	#Define logFC and padj thresholds:
	logFC_ecdf=ecdf(dedata$logFC)
	p.adj.thresh=0.05
	#Parse glycogene Data:
	glycogene_dedata<-dedata %>% dplyr::filter(geneName %in% glycogenes)
	# Gene Set Enrichments:
	# - Entire Genome:
	pathwayEnrich_wholeGenome_upreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,geneSet,high=T)
	pathwayEnrich_wholeGenome_downreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,geneSet,high=F)
	# - Glycogene-Centered:
	pathwayEnrich_glycogeneCentered_upreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,geneSet,high=T)
	pathwayEnrich_glycogeneCentered_downreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,geneSet,high=F)
	
	#Organize into list:
	dta_out<-list(
		      'genome'=list('up'=pathwayEnrich_wholeGenome_upreg,'down'=pathwayEnrich_wholeGenome_downreg),
		      'glycogene'=list('up'=pathwayEnrich_glycogeneCentered_upreg,'down'=pathwayEnrich_glycogeneCentered_downreg)
	      )
	return(dta_out)
}


# All Group Wrapper:

DE_data_enrichWrapper<-function(DEResults,glycogenes,pathwayList,functionList){
	#Loops through each set of differential expression analysis
	# performed on each cancer type and does enrichment.
	# Each cancer must have a "shortLetterCode" set, which is the 
	# tumor vs normal differential expression data
	# Other cancers may have more molecular subtypes/clinical metadata to
	# test, and will be stored.
	enrich_results<-lapply(DEResults,function(deset){
		#Main effects:
		main_data<-deset$DE_Main
		main_data_names<-names(main_data)
		main_enrichments<-lapply(main_data,function(m_deset){
				enrichSet_pathway<-do_enrich(m_deset,glycogenes,pathwayList)
				enrichSet_function<-do_enrich(m_deset,glycogenes,functionList)
				enrichSet=list('Pathways'=enrichSet_pathway,'Functions'=enrichSet_function)
				return(enrichSet)
		 });names(main_enrichments)<-main_data_names
		#Purity effects:
		purity_data<-deset$DE_purity
		purity_enrichments_pathway<-do_enrich_purity(purity_data,glycogenes,pathwayList)
		purity_enrichments_function<-do_enrich_purity(purity_data,glycogenes,functionList)
		purity_enrichments<-list('purity'=list('Pathways'=purity_enrichments_pathway,'Function'=purity_enrichments_function))
		#Concatenate enrichment lists:
		total_enrich<-c(main_enrichments,purity_enrichments)
		return(total_enrich)
      });names(enrich_results)<-names(DEResults)
      return(enrich_results)
}

# -------- Execution ----------

args=commandArgs(trailing=T)
inFile=args[1]
#Load Differential Expression Analysis:
load(inFile)
ctype <- sub('(TCGA\\-\\D+)\\_.*','\\1',basename(inFile))
enrichments<-DE_data_enrichWrapper(DEResults,glycogenes,pathwayList,functionList)
save(file=file.path('analysis/enrichments',paste(ctype,'enrichments.rda',sep='_')),enrichments)

#Gather DE Data:
#dedata<-get_tumorVnormal(DEResults)
#dedata_purity<-get_tumorPurityEffects(DEResults)
##Glycogene-specific DE data:
#glycogene_dedata<-get_glycogene_dedata(dedata,glycogenes)
#glycogene_dedata_purity<-get_glycogene_dedata(dedata_purity,glycogenes)
#
##Define logFC and padj thresholds:
#logFC_lowthresh=log(1/2)
#logFC_highthresh=log(2)
#logFC_ecdf=ecdf(dedata_purity$logFC)
#p.adj.thresh=0.05
#
## Glycosylation Pathway Enrichments:
## - Entire Genome:
#pathwayEnrich_wholeGenome_upreg<-test_pathways(dedata,logFC_highthresh,p.adj.thresh,pathwayList,high=T)
#pathwayEnrich_wholeGenome_downreg<-test_pathways(dedata,logFC_lowthresh,p.adj.thresh,pathwayList,high=F)
#pathwayEnrich_wholeGenome_purity_upreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,pathwayList,high=T)
#pathwayEnrich_wholeGenome_purity_downreg<-test_pathways_bin(dedata,logFC_ecdf,p.adj.thresh,pathwayList,high=F)
## - Glycogene-Centered:
#pathwayEnrich_glycogeneCentered_upreg<-test_pathways(glycogene_dedata,logFC_highthresh,p.adj.thresh,pathwayList,high=T)
#pathwayEnrich_glycogeneCentered_downreg<-test_pathways(glycogene_dedata,logFC_lowthresh,p.adj.thresh,pathwayList,high=F)
#pathwayEnrich_glycogeneCentered_purity_upreg<-test_pathways_bin(glycogene_dedata_purity,logFC_ecdf,p.adj.thresh,pathwayList,high=T)
#pathwayEnrich_glycogeneCentered_purity_downreg<-test_pathways_bin(glycogene_dedata_purity,logFC_ecdf,p.adj.thresh,pathwayList,high=F)
#
## Glycogene Function Enrichments:
## - Entire Genome:
#functionEnrich_wholeGenome_upreg<-test_pathways(dedata,logFC_highthresh,p.adj.thresh,functionList,high=T)
#functionEnrich_wholeGenome_downreg<-test_pathways(dedata,logFC_lowthresh,p.adj.thresh,functionList,high=F)
#functionEnrich_wholeGenome_purity_upreg<-test_pathways_bin(dedata_purity,logFC_ecdf,p.adj.thresh,functionList,high=T)
#functionEnrich_wholeGenome_purity_downreg<-test_pathways_bin(dedata_purity,logFC_ecdf,p.adj.thresh,functionList,high=F)
## - Glycogene-Centered:
#functionEnrich_glycogeneCentered_upreg<-test_pathways(glycogene_dedata,logFC_highthresh,p.adj.thresh,functionList,high=T)
#functionEnrich_glycogeneCentered_downreg<-test_pathways(glycogene_dedata,logFC_lowthresh,p.adj.thresh,functionList,high=F)
#functionEnrich_glycogeneCentered_purity_upreg<-test_pathways_bin(glycogene_dedata_purity,logFC_ecdf,p.adj.thresh,functionList,high=T)
#functionEnrich_glycogeneCentered_purity_downreg<-test_pathways_bin(glycogene_dedata_purity,logFC_ecdf,p.adj.thresh,functionList,high=F)
#
##Create Results List:
## Pathways:
## | --> genome --> tumorVnormal/purity --> up/down
## | --> glycogene --> tumorVnormal/purity --> up/down
## Functions:
## | --> genome --> tumorVnormal/purity --> up/down
## | --> glycogene --> tumorVnormal/purity --> up/down
#
#enrichments<-list(
#	'Pathways'=list('genome'=list(
#				'tumorVnormal'=list('up'=pathwayEnrich_wholeGenome_upreg,
#				     'down'=pathwayEnrich_wholeGenome_downreg
#				     ),
#				'purity'=list('up'=pathwayEnrich_wholeGenome_purity_upreg,
#				     'down'=pathwayEnrich_wholeGenome_purity_downreg
#				     )),
#	       		'glycogene'=list(
#				'tumorVnormal'=list('up'=pathwayEnrich_glycogeneCentered_upreg,
#					 'down'=pathwayEnrich_glycogeneCentered_downreg
#					 ),
#				'purity'=list('up'=pathwayEnrich_glycogeneCentered_purity_upreg,
#					 'down'=pathwayEnrich_glycogeneCentered_purity_downreg
#					 )
#			)),
#
#	'Functions'=list('genome'=list(
#				'tumorVnormal'=list('up'=functionEnrich_wholeGenome_upreg,
#				     'down'=functionEnrich_wholeGenome_downreg
#				     ),
#				'purity'=list('up'=functionEnrich_wholeGenome_purity_upreg,
#				     'down'=functionEnrich_wholeGenome_purity_downreg
#				     )),
#	       		'glycogene'=list(
#				'tumorVnormal'=list('up'=functionEnrich_glycogeneCentered_upreg,
#					 'down'=functionEnrich_glycogeneCentered_downreg
#					 ),
#				'purity'=list('up'=functionEnrich_glycogeneCentered_purity_upreg,
#					 'down'=functionEnrich_glycogeneCentered_purity_downreg
#					 ))
#			)
#)
#

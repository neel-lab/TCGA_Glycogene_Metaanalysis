library(dplyr);library(tidyr);library(tibble)
library(parallel)

#--------Enrichment Functions ------------

get_tumorVnormal<-function(DEResults){
	#Parse the DE data List:
	return(DEResults$shortLetterCode$DE_Main[[1]])
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
		 pass & high ~ (logFC>=logFC_thresh & adj.P.Val<=adjP_thresh),
		 # Passes a low threshold
		 pass & !high ~ (logFC<=logFC_thresh & adj.P.Val<=adjP_thresh),
		 # NOT Passes a high threshold:
		 !pass & high ~ (logFC<logFC_thresh | adj.P.Val>adjP_thresh),
		 # NOT Passes a low threshold:
		 !pass & !high ~ (logFC>logFC_thresh | adj.P.Val>adjP_thresh)
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

# -------- Execution ----------
#Load Differential Expression Analysis:
load('./analysis/DE_Analysis/TCGA-BRCA_DEData.rda')

#Define logFC and padj thresholds:
logFC_lowthresh=log(1/2)
logFC_highthresh=log(2)
p.adj.thresh=0.05

#PAM50 subtypes analyzed:
pam50_names<-names(DEResults$paper_BRCA_Subtype_PAM50$DE_Main)

pam50_enrichment_dta<-lapply(pam50_names,function(n){
	dedata<-DEResults$paper_BRCA_Subtype_PAM50$DE_Main[[n]]
	#Glycogene-specific DE data:
	glycogene_dedata<-get_glycogene_dedata(dedata,glycogenes)
	# Glycosylation Pathway Enrichments:
	# - Entire Genome:
	pathwayEnrich_wholeGenome_upreg<-test_pathways(dedata,logFC_highthresh,p.adj.thresh,pathwayList,high=T)
	pathwayEnrich_wholeGenome_downreg<-test_pathways(dedata,logFC_lowthresh,p.adj.thresh,pathwayList,high=F)
	# - Glycogene-Centered:
	pathwayEnrich_glycogeneCentered_upreg<-test_pathways(glycogene_dedata,logFC_highthresh,p.adj.thresh,pathwayList,high=T)
	pathwayEnrich_glycogeneCentered_downreg<-test_pathways(glycogene_dedata,logFC_lowthresh,p.adj.thresh,pathwayList,high=F)

	# Glycogene Function Enrichments:
	# - Entire Genome:
	functionEnrich_wholeGenome_upreg<-test_pathways(dedata,logFC_highthresh,p.adj.thresh,functionList,high=T)
	functionEnrich_wholeGenome_downreg<-test_pathways(dedata,logFC_lowthresh,p.adj.thresh,functionList,high=F)
	# - Glycogene-Centered:
	functionEnrich_glycogeneCentered_upreg<-test_pathways(glycogene_dedata,logFC_highthresh,p.adj.thresh,functionList,high=T)
	functionEnrich_glycogeneCentered_downreg<-test_pathways(glycogene_dedata,logFC_lowthresh,p.adj.thresh,functionList,high=F)

	enrichments<-list(
		'Pathways'=list('genome'=
					list('up'=pathwayEnrich_wholeGenome_upreg,
					     'down'=pathwayEnrich_wholeGenome_downreg
					     ),
				'glycogene'=list('up'=pathwayEnrich_glycogeneCentered_upreg,
						 'down'=pathwayEnrich_glycogeneCentered_downreg
						 )
				),

		'Functions'=list('genome'=
					list('up'=functionEnrich_wholeGenome_upreg,
					     'down'=functionEnrich_wholeGenome_downreg
					     ),
				'glycogene'=list('up'=functionEnrich_glycogeneCentered_upreg,
						 'down'=functionEnrich_glycogeneCentered_downreg
						 )
				)
	)
	return(enrichments)
});names(pam50_enrichment_dta)<-pam50_names
save(file='pam50_enrichment_dta.rda',pam50_enrichment_dta)

pam50_enrichment_df<-do.call(rbind,lapply(names(pam50_enrichment_dta),function(x){
		pathDta<-rbind(
	pam50_enrichment_dta[[x]][['Pathways']][['genome']][['up']] %>% 
		mutate(universe='genome',direction='up'),
	pam50_enrichment_dta[[x]][['Pathways']][['genome']][['down']] %>% 
		mutate(universe='genome',direction='down'),
	pam50_enrichment_dta[[x]][['Pathways']][['glycogene']][['up']] %>% 
		mutate(universe='glycogene',direction='up'),
	pam50_enrichment_dta[[x]][['Pathways']][['glycogene']][['down']] %>% 
		mutate(universe='glycogene',direction='down')
	) %>% mutate(pathClass='Pathway')
		funDta<-rbind(
	pam50_enrichment_dta[[x]][['Functions']][['genome']][['up']] %>% 
		mutate(universe='genome',direction='up'),
	pam50_enrichment_dta[[x]][['Functions']][['genome']][['down']] %>% 
		mutate(universe='genome',direction='down'),
	pam50_enrichment_dta[[x]][['Functions']][['glycogene']][['up']] %>% 
		mutate(universe='glycogene',direction='up'),
	pam50_enrichment_dta[[x]][['Functions']][['glycogene']][['down']] %>% 
		mutate(universe='glycogene',direction='down')
	) %>% mutate(pathClass='Function')
		totalDta<-rbind(pathDta,funDta) %>% mutate(grouping=x)
		return(totalDta)
})) %>% mutate(grouping=sub('paper_BRCA_Subtype_PAM50','',grouping))
save(file='pam50_enrichment_df.rda',pam50_enrichment_df)

library(dplyr);library(tidyr);library(tibble)
library(SummarizedExperiment)
library(edgeR)
library(limma);library(variancePartition)#dream package
library(parallel)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(RobustRankAggreg)

#List files:

dedata_files<-list.files('./analysis/DE_Analysis/',pattern='DEData.rda')

#Aggregate Tumor v Normal data:

allTumorVNormal<-lapply(dedata_files,function(f){
	load(file.path('./analysis/DE_Analysis',f))
	return(DEResults$shortLetterCode$DE_Main$shortLetterCodeT)
});names(allTumorVNormal)<-sapply(dedata_files,function(x) unlist(strsplit(split='_',x=x))[1])


allpurity<-lapply(dedata_files,function(f){
	load(file.path('./analysis/DE_Analysis',f))
	return(DEResults$shortLetterCode$DE_purity)
});names(allpurity)<-sapply(dedata_files,function(x) unlist(strsplit(split='_',x=x))[1])

glycogene_dedata_df<-do.call(rbind,lapply(names(allTumorVNormal),function(n){
	dta<-allTumorVNormal[[n]] %>%
		right_join(data.frame(geneName=glycogenes),by='geneName') %>%
		select(geneName,logFC,adj.P.Val) %>%
		mutate(logFC=ifelse(is.na(logFC),0,logFC),
		       adj.P.Val=ifelse(is.na(adj.P.Val),1,adj.P.Val),
		       ctype=n)
	return(dta)
})) %>% right_join(data.frame(geneName=glycogenes),by='geneName')


glycogene_de_puritydata_df<-do.call(rbind,lapply(names(allpurity),function(n){
	dta<-allpurity[[n]] %>%
		right_join(data.frame(geneName=glycogenes),by='geneName') %>%
		select(geneName,logFC,adj.P.Val) %>%
		mutate(logFC=ifelse(is.na(logFC),0,logFC),
		       adj.P.Val=ifelse(is.na(adj.P.Val),1,adj.P.Val),
		       ctype=n)
	return(dta)
})) %>% right_join(data.frame(geneName=glycogenes),by='geneName')

#Create group annotations:
glycogene_dedata_df<-glycogene_dedata_df %>% 
	mutate(pathway=case_when(
		geneName %in% pathwayList[['donor_synthesis']] ~ 'Donor Synthesis',
		geneName %in% c(pathwayList[['Glycosylation']],pathwayList[['Glycan_Degradation']]) ~ 'Glycosylation and Glycan Degradation',
		geneName %in% pathwayList[['Metabolic_GlcAT']] ~ 'Glucuronidation',
		geneName %in% c(pathwayList[['Sugar_Metabolism']],pathwayList[['Transporters']]) ~ 'Sugar Metabolism and Transporters',
		TRUE ~ 'Other'
))

#Create subgroup annotations:
glycogene_dedata_df<-glycogene_dedata_df %>% 
	mutate(subpathway=case_when(
		#Glycosylation & Degradation:
		geneName %in% pathwayList[['C-Man']] ~ "C-Mannose",
		geneName %in% pathwayList[['GAG']] ~ "GAG",
		geneName %in% pathwayList[['Lipid-linked']] ~ "Lipid-linked",
		geneName %in% pathwayList[['N-linked']] ~ "N-linked",
		geneName %in% pathwayList[['O-GalNAc']] ~ "O-GalNAc",
		geneName %in% pathwayList[['O-Fuc']] ~ "O-Fuc",
		geneName %in% pathwayList[['O-Glc']] ~ "O-Glc",
		geneName %in% pathwayList[['O-Man']] ~ "O-Man",
		geneName %in% pathwayList[['POGLUT1_Type']] ~ "POGLUT1",
		geneName %in% pathwayList[['LacdiNAc']] ~ "LacdiNAc",
		geneName %in% pathwayList[['Type_I_LacNAc']] ~ "Type I LacNAc",
		geneName %in% pathwayList[['Type_II_LacNAc']] ~ "Type II LacNAc",
		geneName %in% pathwayList[['Type_III_LacNAc']] ~ "Type III LacNAc",
		geneName %in% pathwayList[['Cytosol_Degradation']] ~ "Cytosolic Degradation",
		geneName %in% pathwayList[['Extracellular_Degradation']] ~ "Extracellular Degradation",
		geneName %in% pathwayList[['Lysosome_Degradation']] ~ "Lysosomal Degradation",
		#donor_synthesis:
		geneName %in% pathwayList[['Nucleotide_Sugar_Synthesis']] ~ "Nucleotide Sugar Synthesis",
		geneName %in% pathwayList[['Nucleotide_Sugar_Transport']] ~ "Nucleotide Sugar Transport",
		#Glucuronidation
		geneName %in% pathwayList[['Metabolic_GlcAT']] ~ 'Glucuronidation',
		#Sugar metabolism and transport:
		geneName %in% c(pathwayList[['Sugar_Metabolism']],pathwayList[['Transporters']]) ~ 'Sugar Metabolism and Transporters',
		TRUE ~ 'Other'
))


annot_dta<-glycogene_dedata_df %>% 
	filter(pathway!="Other" & subpathway!='Other') %>% 
	select(geneName,pathway,subpathway) %>% distinct() %>%
	arrange(pathway,subpathway) %>% 
	column_to_rownames('geneName')

#Filter desired annotation data here:
annot_dta<-annot_dta %>% 
	select(subpathway) %>% 
	dplyr::rename(pathway=subpathway)

glycogene_dedata_classifMat<-glycogene_dedata_df %>% 
	filter(pathway!="Other") %>% 
	mutate(logFC=ifelse(is.na(logFC),0,logFC)) %>%
	pivot_wider(id_cols='ctype',
		    names_from='geneName',values_from='logFC') %>% 
	column_to_rownames('ctype') %>% t()

#Make the heatmap:
breaks<-seq(-15,15,1/10)
colorVec<-colorRampPalette(c('blue','white','red'))(length(breaks))
#Label genes with 0 lfc with black to show lowly expressed
colorVec[(length(colorVec)/2)+1]<-'black'
p<-pheatmap(glycogene_dedata_classifMat[row.names(annot_dta),],color=colorVec,cluster_rows=F,annotation_row=annot_dta,labels_row='',breaks=breaks,cellheight=2,cellwidth=2)

ggsave(filename='Figure3_logFC_heatmap.eps',device='eps',plot=p,width=8,height=5.5,units='in')



# ---------- Ranked Bar Plots --------------

glycogene_dedata_df<-do.call(rbind,lapply(names(allTumorVNormal),function(n){
	dta<-allTumorVNormal[[n]] %>%
		right_join(data.frame(geneName=glycogenes),by='geneName') %>%
		select(geneName,logFC,adj.P.Val) %>%
		#mutate(logFC=ifelse(is.na(logFC),0,logFC),
		#       adj.P.Val=ifelse(is.na(adj.P.Val),1,adj.P.Val),
		#       ctype=n)
		mutate(ctype=n)
	return(dta)
})) %>% right_join(data.frame(geneName=glycogenes),by='geneName')
save(file='./analysis/DE_Analysis/glycogene_dedata_df.rda',glycogene_dedata_df)


get_path_barData<-function(glycogene_dedata_df,pathway,ct,pathwayList){
	#Gather data:
	dta<-glycogene_dedata_df %>% 
		filter(geneName %in% pathwayList[[pathway]] & ctype==ct) %>% 
		arrange(logFC) %>% 
		mutate(geneName=factor(geneName,levels=geneName))
	print(dta)
	# Make barplot:
	p<-ggplot(dta,aes(x=geneName,y=logFC,fill=logFC)) +
		geom_bar(stat='identity',color='black',size=1.1) + 
		coord_flip() + 
		scale_fill_gradientn(colors=c('blue','white','red')) + 
		labs(title=paste(ct,pathway,sep=': '),y='log fold change',x='Glycogene') + 
		theme_classic()
	return(p)
}



# -------- Robust Rank Aggregation --------

get_ranked_lists<-function(glycogene_dedata_df,pathway,pathwayList,direction='up'){
	lists<-glycogene_dedata_df %>% 
		filter(geneName %in% pathwayList[[pathway]]) %>% 
		group_by(ctype) %>% 
		group_map(~{
			.x %>% filter(!is.na(logFC) & !is.na(adj.P.Val)) %>% 
				filter(abs(logFC)>=log(2) & adj.P.Val<=0.05) %>% 
				arrange(
					case_when(
						  direction=='up' ~ -logFC,
						  direction=='down' ~ logFC
						  )
					) %>%  pull(geneName)
		    })
	return(lists)
}

run_rra<-function(glycogene_dedata_df,pathwayList){
	#Scope the pathways :
	pthList<-pathwayList
	#pthList<-pathwayList[which(sapply(pathwayList,function(x) length(x)>=3))]
	#Run RRA:
	upResults<-lapply(names(pthList),function(x){
		message(x)
		lsts<-get_ranked_lists(glycogene_dedata_df,x,pthList,direction='up')
		if (all(sapply(lsts,function(x) length(x)==0))){
			return(NULL)
		} else {
			return(aggregateRanks(lsts,N=length(pthList[[x]])))
		}
	    });names(upResults)<-names(pthList)
	downResults<-lapply(names(pthList),function(x){
		message(x)
		lsts<-get_ranked_lists(glycogene_dedata_df,x,pthList,direction='down')
		if (all(sapply(lsts,function(x) length(x)==0))){
			return(NULL)
		} else {
			return(aggregateRanks(lsts,N=length(pthList[[x]])))
		}
	    });names(downResults)<-names(pthList)
	resList<-list('upResults'=upResults,'downResults'=downResults)
	return(resList)
}

RRA_res<-run_rra(glycogene_dedata_df,pathwayList)
RRA_res$upResults[sapply(RRA_res$upResults,is.null)]<-NULL
RRA_res$downResults[sapply(RRA_res$downResults,is.null)]<-NULL

#Function:
RRA_res_fun<-run_rra(glycogene_dedata_df,functionList)
RRA_res_fun$upResults[sapply(RRA_res_fun$upResults,is.null)]<-NULL
RRA_res_fun$downResults[sapply(RRA_res_fun$downResults,is.null)]<-NULL


#DataFrame for pathways:
RRA_res_df<-do.call(rbind,lapply(names(RRA_res),function(n){
	dta<-do.call(rbind,lapply(names(RRA_res[[n]]),function(p){
		tbl<-RRA_res[[n]][[p]] %>% 
			mutate(pathway=p,direction=n)
		return(tbl)
	}))
	return(dta)
}))

#DataFrame for functions:
RRA_res_fun_df<-do.call(rbind,lapply(names(RRA_res_fun),function(n){
	dta<-do.call(rbind,lapply(names(RRA_res_fun[[n]]),function(p){
		tbl<-RRA_res_fun[[n]][[p]] %>% 
			mutate(fun=p,direction=n)
		return(tbl)
	}))
	return(dta)
}))

#RRA_pathwayData<-RRA_res_df %>% filter(Score<=0.001 & !(pathway %in% c('Pathway','Glycogene')))
#RRA_functionData<-RRA_res_fun_df %>% filter(Score<=0.001 & !(fun %in% c('Pathway','Glycogene')))

#Save all results into objects:
RRA_pathwayData<-RRA_res_df
RRA_functionData<-RRA_res_fun_df
message('Writing out table')
write.table(file='RRA_statistics.tsv',rbind(
				RRA_pathwayData %>% mutate(enrichType='Pathway') %>% rename(geneSet=pathway),
				RRA_functionData %>% mutate(enrichType='Function' %>% rename(geneSet=fun))
				),sep='\t',quote=F,row.names=F)

save(file='./RRA_pathwayData.rda',RRA_pathwayData)
save(file='./RRA_functionData.rda',RRA_functionData)

#Illustrate degree of consistent up/down regulation:
rra_plotFun<-function(glycogene_dedata_df,RRA_res_sig,pathwayList){
	plotList<-lapply(unique(RRA_res_sig$pathway),function(x){
		pthGenes<-pathwayList[[x]]
		message(x)
		hits<-RRA_res_sig[RRA_res_sig$pathway==x,] %>% 
			group_by(pathway,Name) %>% 
			filter(Score==min(Score)) %>% 
			#Remove ties:
			filter(n()<2)
		if (nrow(hits)==0){
			return(NULL)
		}
		highGenes<-hits[hits$direction=='upResults','Name']$Name
		lowGenes<-hits[hits$direction=='downResults','Name']$Name
		#Filter out genes in both high and low:
		orderList<-c(pthGenes[pthGenes %in% lowGenes],
			     setdiff(pthGenes,c(highGenes,lowGenes)),
		       pthGenes[pthGenes %in% highGenes]
		       )
		#Parse dedata:
		dta<-glycogene_dedata_df %>% 
			filter(geneName %in% pthGenes) %>% 
			filter(abs(logFC)>=log(2) & adj.P.Val<=0.05) %>% 
			mutate(geneName=factor(geneName,levels=orderList)) %>% 
			mutate(Enrichment=case_when(
					geneName %in% highGenes ~ 'High',
					geneName %in% lowGenes ~ 'Low',
					TRUE ~ 'N.S.'
					)) %>% 
			mutate(Significant=case_when(
					abs(logFC)>=log(2) & adj.P.Val<=0.05 ~ 'Significant',
				       TRUE ~ 'Not Significant'))	
		# Report how many logFC values 
		# are present for each gene:
		repDta<-dta %>% filter(!is.na(logFC)) %>% 
			group_by(geneName) %>% 
			summarize(tally=n()) %>% 
			mutate(labels=paste(geneName,paste('[N=',tally,']',sep=''),sep=':')) %>% 
			select(-tally) %>% 
			mutate(labels=factor(labels,levels=labels))
		dta<-dta %>% left_join(repDta,by='geneName')
		#Make a plot:
		# Plot constants:
		p<-ggplot(dta) +
		  geom_boxplot(aes(x=labels,y=logFC,fill=Enrichment),outlier.size=0) + 
		  geom_jitter(aes(x=labels,y=logFC),size=0.85,color='black') + 
		  theme_classic() + 
		  labs(title=x,x='Glycogene',y='logFC') + 
		  scale_fill_manual(values=c('High'='firebrick1','Low'='royalblue1','N.S.'='white')) + 
		  scale_color_manual(values=c('darkgrey','black')) + 
		  guides(fill=guide_legend(title='RRA Significant?'),colour=guide_legend(title='DE Significant?')) + 
		  theme(axis.text.x=element_text(hjust=1,angle=45))
		return(p)
	    });names(plotList)<-unique(RRA_res_sig$pathway)
	return(plotList)
}

rra_pathway_plotLists<-rra_plotFun(glycogene_dedata_df,RRA_pathwayData,pathwayList)
rra_function_plotLists<-rra_plotFun(glycogene_dedata_df,RRA_functionData,functionList)
#rra_plotLists[sapply(rra_plotLists,function(x) is.null(x))]<-NULL
#rra_pathway_plotLists[sapply(rra_pathway_plotLists,function(x) is.null(x))]<-NULL
#rra_function_plotLists[sapply(rra_function_plotLists,function(x) is.null(x))]<-NULL
save(file='rra_pathway_plotLists.rda',rra_pathway_plotLists)
save(file='rra_function_plotLists.rda',rra_function_plotLists)

# ---------- RRA organized gene panel: ---------

#Make Grob Objects:
#make_grob_set<-function(plotList,nrow,ncol){
#	layout_matrix<-matrix(seq(1,length(plotList)),nrow=nrow,ncol=ncol,byrow=T)
#	grobs<-lapply(plotList,function(x) ggplotGrob(x))
#	g<-grid.arrange(plotList,grobs=grobs,layout_matrix=layout_matrix)
#	return(g)
#}
#
## N-linked:
#nlinkedPlots<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('Branching','Dolichol_Pathway','Lysosomal_Targeting','Processing'))]
#nlinked_grid<-make_grob_set(nlinkedPlots,2,1)
#ggsave(filename='N-linked_glycogeneRRA.eps',device='eps',plot=nlinked_grid,width=8.5,height=6,units='in')
#
## Lipid-linked(none):
##llPlots<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('GalCer','GlcCer','GPI_Linked'))]
##ll_grid<-make_grob_set(llPlots,2,1)
##ggsave(filename='N-linked_glycogeneRRA.eps',device='eps',plot=nlinked_grid,width=8.5,height=6,units='in')
#
## O-Linked(none):
#OlinkedPlots<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('O-Fuc','Core1','Core2','Core3','Core4','O-Glc','O-Man'))]
#Olinked_grid<-make_grob_set(OlinkedPlots,2,2)
#
## GAG:
#GAGPlots<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('Hyaluronan','Keratan_Sulfate','CS&DS_Polymerization','CS&DS_Sulfation','HS_Polymerization','HS_Sulfation'))]
#GAG_grid<-make_grob_set(GAGPlots,3,1)
#ggsave(filename='GAG_glycogeneRRA.eps',device='eps',plot=GAG_grid,width=8.5,height=6,units='in')
#
## Capping:
#CappingPlot<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('Capping'))]
#Capping_grid<-make_grob_set(CappingPlot,3,1)
#
##Glycan Degradation:
#GlyDegPlot<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('Lysosome_Degradation','Extracellular_Degradation','Cytosol_Degradation'))]
#GlyDeg_grid<-make_grob_set(GlyDegPlot,2,1)
#ggsave(filename='GlyDeg_glycogeneRRA.eps',device='eps',plot=GlyDeg_grid,width=8.5,height=6,units='in')
#
#
##Elongation and Branching (none):
##EBPlot<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('LacdiNAc','Type_I_LacNAc','Type_II_LacNAc','Type_III_LacNAc'))]
##EB_grid<-make_grob_set(EBPlot,3,1)
#
#
##Sugar Metabolims
#SugarMetabolismPlot<-rra_pathway_plotLists[which(names(rra_pathway_plotLists) %in% c('Sugar_Metabolism'))]
#SugarMetabolism_grid<-make_grob_set(SugarMetabolismPlot,1,1)

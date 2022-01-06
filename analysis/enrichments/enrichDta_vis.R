library(dplyr);library(tidyr);library(ggplot2);library(tibble)


#Read in enrichment data:
enrichdta<-read.table('./analysis/enrichments/allEnrichData.tsv',sep='\t',header=T)

p<-enrichdta %>% group_by(ctype,pathways) %>% filter(any(adjP<=0.1)) %>% 
	mutate(enrichDir=case_when(adjP<=0.1 & direction=='up' ~ 'Upregulated',
				   adjP<=0.1 & direction=='down' ~ 'Downregulated',
				   TRUE ~ 'N.S.'
				   )) %>% 
	filter(enrichDir!='N.S.',pathways!='Pathway',pathways!='Glycogene') %>% 
	mutate(enrichDir=factor(enrichDir,levels=c('Upregulated','Downregulated'))) %>%
	ggplot(aes(x=ctype,y=pathways,fill=enrichDir)) +
	geom_tile(color='black',size=0.85) +
	theme_classic() +
	facet_wrap(~pathType + universe,scale='free') +
	guides(fill=guide_legend(title='Enrichment')) + 
	labs(x='Cancer Type',y='Glycogene Sets') + 
       	theme(axis.text.x=element_text(angle=45,hjust=1,size=5.75),
	axis.text.y=element_text(angle=15,size=5.5))	

ggsave(filename='enrichmentTestResults_all.eps',device='eps',plot=p,width=8,height=5.5,units='in')

p_genome<-enrichdta %>% group_by(ctype,pathways) %>% filter(any(adjP<=0.1) & universe=='genome') %>% 
	mutate(enrichDir=case_when(adjP<=0.1 & direction=='up' ~ 'Upregulated',
				   adjP<=0.1 & direction=='down' ~ 'Downregulated',
				   TRUE ~ 'N.S.'
				   )) %>% 
	filter(enrichDir!='N.S.',pathways!='Pathway',pathways!='Glycogene') %>% 
	mutate(enrichDir=factor(enrichDir,levels=c('Upregulated','Downregulated'))) %>%
		ggplot(aes(x=ctype,y=pathways,fill=enrichDir)) +
		geom_tile(color='black',size=0.85) +
		theme_classic() +
		facet_wrap(~pathType,scale='free') +
		guides(fill=guide_legend(title='Enrichment')) + 
		labs(x='Cancer Type',y='Glycogene Sets') + 
		theme(axis.text.x=element_text(angle=45,hjust=1,size=5.75),
		axis.text.y=element_text(angle=15,size=5.5))	

ggsave(filename='enrichmentTestResults.eps',device='eps',plot=p_genome,width=8,height=5.5,units='in')


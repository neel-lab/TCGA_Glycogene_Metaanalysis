library(dplyr);library(tidyr);library(ggplot2)
setwd('~/tgroth_2/Programs/TCGA_Analysis/BC_Progression/Manuscript_Analysis/')

#Glycogene pseudotime data:
load('./glycogene_pseudotime_normal_to_luminal.rda')
load('./glycogene_pseudotime_normal_to_basal.rda')
load('./pam50_pseudotime_normal_to_luminal.rda')
load('./pam50_pseudotime_normal_to_basal.rda')

#Finding common nodes:
commonGenes<-function(dta1,dta2){
  gene1=dta1 %>% pull(gene) %>% as.character()  %>% unique()
  gene2=dta2 %>% pull(gene) %>% as.character() %>% unique()
  return(length(intersect(gene1,gene2)))
}

comparePathSOMS<-function(set1,set2){
  set1_dta<-lapply(set1,function(x) x$data)
  set2_dta<-lapply(set2,function(x) x$data)
  pairGrid<-expand.grid(n1=seq(1,length(set1_dta)),n2=seq(1,length(set2_dta)))
  commonScores<-apply(pairGrid,1,function(x) commonGenes(set1_dta[[x['n1']]],set2_dta[[x['n2']]]))
  topMatch<-pairGrid[which(commonScores==max(commonScores)),]
  print(topMatch)
  pairedSOM<-list('s1_plot'=set1[[topMatch[1,'n1']]],'s2_plot'=set2[[topMatch[1,'n2']]])
  return(pairedSOM)
}

#Top common Nomral to basal SOM nodes between glycogene and PAM50 progression model:
commonBasal<-comparePathSOMS(glycogene_pseudotime_normal_to_basal,pam50_pseudotime_normal_to_basal)
layout_matrix<-matrix(seq(1,2),nrow=1,ncol = 2,byrow = TRUE)
commonBasalgrobs<-lapply(commonBasal,function(x) ggplotGrob(x))
commonBasalGrid<-arrangeGrob(commonBasal,grobs=commonBasalgrobs,layout_matrix=layout_matrix)
ggsave(filename='commonBasal_SOM_Nodes.eps',device='eps',plot=commonBasalGrid,width=6,height=3,units='in')
#Top common Nomral to luminal SOM nodes between glycogene and PAM50 progression model:
commonLuminal<-comparePathSOMS(glycogene_pseudotime_normal_to_luminal,pam50_pseudotime_normal_to_luminal)
layout_matrix<-matrix(seq(1,2),nrow=1,ncol = 2,byrow = TRUE)
commonLuminalgrobs<-lapply(commonLuminal,function(x) ggplotGrob(x))
commonLuminalGrid<-arrangeGrob(commonLuminal,grobs=commonLuminalgrobs,layout_matrix=layout_matrix)
ggsave(filename='commonLuminal_SOM_Nodes.eps',device='eps',plot=commonLuminalGrid,width=6,height=3,units='in')

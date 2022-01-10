#Breast Cancer Progression Model Analysis
#Find Expression trends along pseudotime
library(openxlsx)
library(kohonen)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)
setwd('~/tgroth_2/Programs/TCGA_Analysis/BC_Progression/Manuscript_Analysis/')
#There are 3 Paths in the xlsx files
BC_GG_pathdata<-list()
BC_GG_pseudotime<-list()
for (i in seq(1,3)){
  BC_GG_pathdata[[i]]<-read.xlsx(xlsxFile = 'BC_pseudotime_data.xlsx',
                            colNames = FALSE,
                            rowNames = TRUE,
                            sheet = i)
  names(BC_GG_pathdata)[i]<-paste('Path',i,sep='_')
  BC_GG_pseudotime[[i]]<-read.xlsx(xlsxFile = 'BC_pseudotime.xlsx',
                                 colNames = FALSE,
                                 sheet = i)
  names(BC_GG_pseudotime)[i]<-paste('Path',i,sep='_')
}

#Iterate through each path, find significantly correlated genes, then create a SOM
#to summarize the gene expression
sigcor_pathgenes<-list()
pathPlot_objs<-list()
pathPlot_rawPlots<-list()
for (p in 1:length(BC_GG_pathdata)){
  #Extract data
  data<-BC_GG_pathdata[[p]]
  pst<-t(BC_GG_pseudotime[[p]])
  gene_cordata<-c()
  for (g in 1:dim(data)[1]){
    genedata<-t(data[g,])
    #Corelate the pseudotime and gene information
    genecor<-cor(pst,genedata)
    gene_cordata<-rbind(gene_cordata,genecor)
  }
  row.names(gene_cordata)<-row.names(data)
  
  #Extract data that significantly correlates with pseudotime
  data_sig<-data[row.names(data) %in% names(gene_cordata[abs(gene_cordata)>=0.3,]),]
  #Normalize data across pseudotime:
  data_sig<-as.data.frame(t(apply(data_sig,1,function(x) (x-mean(x))/sd(x))))
  sigcor_pathgenes[[p]]<-row.names(data_sig)
  names(sigcor_pathgenes)[p]<-paste('Path',p,sep='_')
  #Create a square SOM to summarize trends in gene expression data
  genesom<-som(as.matrix(data_sig),somgrid(3,3,'rectangular',neighbourhood.fct = 'gaussian'))
  #Create plots for each SOM node, showing general trend in glycogene expression:
  #First, find max and min ranges for all plots:
  #Melting data:
  dta_pseudo<-as.data.frame(cbind(t(data_sig),t(BC_GG_pseudotime[[p]])))
  colnames(dta_pseudo)[dim(dta_pseudo)[2]]<-'pseudotime'
  dta_melt<-melt(dta_pseudo,id.vars='pseudotime')
  
  maxval<-dta_melt %>% group_by(pseudotime) %>% summarize(mean(value)) %>% ungroup(pseudotime) %>% max
  minval<-dta_melt %>% group_by(pseudotime) %>% summarize(mean(value)) %>% ungroup(pseudotime) %>% min
  
  plotList<-lapply(unique(genesom$unit.classif)[order(unique(genesom$unit.classif))],function(x){
    dta<-as.data.frame(t(rbind(as.matrix(data_sig[which(genesom$unit.classif==x),]),t(as.matrix(pst)))))
    colnames(dta)[dim(dta)[2]]<-'pseudotime'
    dta_melt<-melt(dta,id.vars = 'pseudotime')
    colnames(dta_melt)<-c('pseudotime','gene','value')
    # p<-ggplot(dta_ave,aes(x=pseudotime,y=zscore),color='black') + 
    #   geom_smooth(size=0.5,se=FALSE) + geom_hline(yintercept = 0,size=1.75,color='black',alpha=0.35) + 
    #   labs(x='',y='',title=paste('Node',x)) + 
    #   scale_x_continuous(limits = c(0,max(dta_ave$pseudotime))) + 
    #   scale_y_continuous(limits = c(minval,maxval)) + 
    #   theme_light() + 
    #   theme(axis.text=element_text(size=7),axis.title = element_text(size=9),
    #         title=element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    #         legend.text = element_text(size=8),legend.key.size = unit(0.01,'cm'),legend.position = 'bottom')
    
    
    p<-ggplot(dta_melt,aes(x=pseudotime,y=value,colour=gene)) + 
      geom_hline(yintercept = 0,size=0.85,color='black') + geom_smooth(se = FALSE,size=0.5) +
      # geom_smooth(aes(colour=gene),se = TRUE,size=0.5,alpha=0) + geom_hline(yintercept = 0,size=1.75,color='black') + 
      labs(x='',y='',title=paste('Node',x)) + guides(colour=guide_legend(title = '',nrow=6,ncol = 5)) +
      scale_x_continuous(limits = c(0,max(dta_melt$pseudotime))) +
      scale_y_continuous(limits = c(minval,maxval)) +
      theme_light() +
      theme(axis.text=element_text(size=7),axis.title = element_text(size=9),
            title=element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            legend.text = element_text(size=6),legend.key.size = unit(0.05,'cm'),
            legend.position = 'bottom',legend.margin = margin(t=-25,l=-20))
    return(p)
  })
  #Make the grid:
  layout_matrix<-matrix(seq(1,9),nrow=3,ncol = 3,byrow = TRUE)
  grobs<-lapply(plotList,function(x) ggplotGrob(x))
  plotGrid<-arrangeGrob(plotList,grobs=grobs,layout_matrix=layout_matrix)
  #Save the grid:
  pathPlot_objs[[paste('Path',p,sep = '_')]]<-plotGrid
  pathPlot_rawPlots[[paste('Path',p,sep = '_')]]<-plotList
  
  fname<-paste('Path',p,'GlycoGene_SOM.png',sep='_')
  png(filename = fname,width=720,height = 720)
  plot(genesom,main=paste('Path',p,'Glycogene SOM'))
  dev.off()
  print(paste('Path',p,'complete',sep=' '))
}
glycogene_pseudotime_normal_to_luminal<-pathPlot_rawPlots[['Path_1']];save(file='glycogene_pseudotime_normal_to_luminal.rda',glycogene_pseudotime_normal_to_luminal)
glycogene_pseudotime_normal_to_basal<-pathPlot_rawPlots[['Path_3']];save(file='glycogene_pseudotime_normal_to_basal.rda',glycogene_pseudotime_normal_to_basal)
ggsave(filename = 'NormalLuminal_SOMNodes.jpg',plot = pathPlot_objs[["Path_1"]],width = 8.5,height=8.5,units = 'in')
ggsave(filename = 'NormalBasal_SOMNodes.jpg',plot = pathPlot_objs[["Path_3"]],width = 8.5,height=8.5,units = 'in')
save(file='BC_sigcor_path_glycogenes.rda',sigcor_pathgenes)
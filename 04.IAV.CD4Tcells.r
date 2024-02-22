
setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

#IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/IAV_new_V1.rds")

Tcells<-subset(IAV,idents = c('2','3','4','7','8','9','12','13','14','16','17','18','19','20',
                              '21','22','23','24','25','27','28','30','31','32','34','35','36',
                              '37','38','39','40','45','47','48','49','50','55'))


#Tcells <- SCTransform(Tcells, vars.to.regress = "percent.mt", verbose = FALSE) 
CD4T<-subset(Tcells, subset = CD4 > 1) 
CD8T<-subset(Tcells, subset = CD8B > 1) 
##说明实际上有很多即不表达CD4也不表达CD8的。###生信技能树的博客已经说明了，即使是示范的例子也是如此。
#CD4T<-SCTransform(CD4T, vars.to.regress = "percent.mt", verbose = FALSE) 

#果然CD8比CD4多很多，与预期不符合，是不是很多文章都常见？查一下。

Tcell_genes = c(
  "PTPRC","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","IFNG",
  "TNF","EOMES","CCR7","SELL","CD27","CD28","CXCR3","IL2",
  "IL12A","IL12B","IL18","STAT4","STAT1","TBX21","CCR4",
  "IL4","IL5","GATA3","CCR3","CCR6","SPI1","IL17A","CCR10",
  "IL22","FOXP3","IL2RA","IL7R","CTLA4","TGFB1","IL10",
  "STAT5A","STAT5B","CXCR5","CD40LG","ICOS","IL21","BCL6",
  "IL23R","RORC","TRAV24","TRBV11-1"
)



CD4_Tcell_genes = c(
  "CCR7","SELL","CD27","CD28","IFNG","LEF1","LTB","TCF7","TNF","CXCR3","IL2",
  "IL12A","IL12B","IL18","STAT4","STAT3","STAT1",
  "TBX21","CCR4","IL4","IL5","GATA3","CCR3","CCR6",
  "SPI1","CCR10","FOXP3","IL2RA","IL7R","CTLA4","TGFB1","IL10",
  "STAT5A","STAT5B","CXCR5","CD40LG","ICOS","IL21",
  "BCL6","IL23R","RORC","TRAV24"
)

CD4_Tcell_plot = c("PTPRC","CD4",
  "CCR7","SELL","CD27","CD28","IFNG","CXCR3","IL2",
  "IL12A","IL12B","IL18","STAT4","STAT3","STAT1",
  "TBX21","CCR4","IL4","IL5","GATA3","CCR3","CCR6",
  "SPI1","CCR10","FOXP3","IL2RA","CTLA4","TGFB1","IL10",
  "STAT5A","STAT5B","CXCR5","CD40LG","ICOS","IL21",
  "BCL6","IL23R","RORC","TRAV24","LEF1","LTB","TCF7") ##IL13, IL3只有RNAassay, TNF基本都无表达


CD4T <- RunPCA(CD4T,features = CD4_Tcell_genes,verbose = FALSE)
CD4T <- RunUMAP(CD4T, dims = 1:14, verbose = FALSE)
CD4T <- FindNeighbors(CD4T, dims = 1:14, verbose = FALSE)
CD4T <- FindClusters(CD4T, resolution=0.4, verbose = FALSE)
DimPlot(CD4T, label = TRUE) + NoLegend() 
###在Dimplot图案上面，特点就是，用一些共有的高表达基因，容易把细胞聚集在一起。然而IL7R不合适，因为它是特异性的基因。
###所以理论上加上CD4，PTRPC，会好看一些。但是实践发现改善不大。


DotPlot(CD4T, features = CD4_Tcell_plot,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))



###################注释###############################
library("plyr")

CD4T@meta.data$self_cluster<-NULL
self_cluster = CD4T@meta.data$seurat_clusters
self_cluster<-revalue(self_cluster, c("0" = "Naive T", 
                                      "1" = "Th2",
                                      "2" = "Th1", 
                                      "3" = "Tfh",
                                      "4" = "Memory T", 
                                      "5" = "Naive T",
                                      "6" = "Memory T", 
                                      "7" = "Th17",
                                      "8" = "Treg1", 
                                      "9" = "Th1",
                                      "10" ="Treg2",
                                      "11" ="Naive T",
                                      "12" ="Unknown",
                                      "13" ="Th1",
                                      "14" ="Th1",
                                      "15" ="γδ T"
))


self_cluster<-factor(self_cluster,levels=c("Naive T","Th1","Th2","Th17","Tfh","Treg1","Treg2","γδ T",
                                           "Memory T","Unknown"))       

CD4T <- AddMetaData(object = CD4T,              
                          metadata = self_cluster,               
                          col.name = "self_cluster")    


DimPlot(CD4T, label = TRUE, group.by="self_cluster") 

DotPlot(CD4T, features = CD4_Tcell_plot,
        assay = 'SCT', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))



###################提取需要的cluster####################

Idents(CD4T) <- "self_cluster"                
CD4T_plot<-subset(CD4T,idents = c("Naive T","Th1","Th2","Th17","Tfh","Treg1",
                                  "Treg2","γδ T","Memory T"))                  

DimPlot(CD4T_plot, label = TRUE, group.by="self_cluster",raster=F) 


DotPlot(CD4T_plot, features = CD4_Tcell_plot,
        assay = 'SCT', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


###############统计cluster数目#################

table(CD4T@meta.data$self_cluster)

CD4T_SampleCell_number<-table(CD4T_plot@meta.data$Sample_ID)

IAV_SampleCell_number<-table(IAV@meta.data$Sample_ID)

cbind(IAV_SampleCell_number,CD4T_SampleCell_number)

  

###############各细胞占比#################
library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####太好了，果然是包冲突了导致dplyr没有起作用

CD4T_metadata<-CD4T_plot@meta.data
CD4T_metadata <- CD4T_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

CD4T_metadata <- CD4T_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

CD4T_metadata <- CD4T_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

CD4T_metadata_select=CD4T_metadata [ , c("Sample_ID",
                                             "Immune_Response_Level",
                                             "self_cluster",
                                             "Sample_Celltype_Sum",
                                             "Sample_Celltype_Pre",
                                             "Sample_Sum") ]

CD4T_metadata_select<-CD4T_metadata_select[!duplicated(CD4T_metadata_select),]


library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))

ggplot(CD4T_metadata_select,aes(x = Immune_Response_Level, 
                                  y = Sample_Celltype_Pre,
                                  fill = Immune_Response_Level)) +
  facet_wrap(~self_cluster,scales = "free",ncol=3,strip.position = "top")+
  # 小提琴图层
  geom_violin(position = position_dodge(0.9),alpha = 0.5,
              width = 1.2,trim = T,
              color = NA) +
  # 箱线图图层http://127.0.0.1:30467/graphics/plot_zoom_png?width=1124&height=834
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  # 主题调整
  #  theme_bw(base_size = 14) +
  #  theme(axis.text.x = element_text(angle = 0,color = 'black'),#,hjust = 1
  #        legend.position = 'top',
  #        aspect.ratio = 0.8) + ##图形长宽比例
  # 颜色设置
  scale_fill_manual(values = c('Low'='#398AB9','High'='red','Normal'='grey'),
                    name = '') +
  # 添加显著性标记
  stat_compare_means(comparisons = comparelist,vjust=0.8,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif"
                     
  ) +
  #ylim(0.5,3.5)
  ###添加标题
  labs(x="Differences in antibody titers and cell proportions among Inactivated vaccines")+
  theme(axis.title.x = element_text(color="black", size=12, face="bold")
  )



##############################差异基因####################################
##############################CD4 Naive_T######################################


library("limma")
CD4_NaiveT<-subset(CD4T_plot,idents = c('Naive T'))
Idents(CD4_NaiveT) <- "Immune_Response"
CD4_NaiveT_LowNormal_marker <- FindMarkers(CD4_NaiveT, ident.1 = "Low", 
                                        ident.2 = "Normal", verbose = FALSE)

CD4_NaiveT_NormalHigh_marker <- FindMarkers(CD4_NaiveT, ident.1 = "Normal", 
                                         ident.2 = "High", verbose = FALSE)

write.table(CD4_NaiveT_LowNormal_marker,
            file = "CD4_NaiveT_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(CD4_NaiveT_NormalHigh_marker,
            file = "CD4_NaiveT_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#######################Menory T cells##################################
CD4_MenoryT<-subset(CD4T_plot,idents = c('Memory T'))
Idents(CD4_MenoryT) <- "Immune_Response"
CD4_MenoryT_LowNormal_marker <- FindMarkers(CD4_MenoryT, ident.1 = "Low", 
                                         ident.2 = "Normal", verbose = FALSE)

CD4_MenoryT_NormalHigh_marker <- FindMarkers(CD4_MenoryT, ident.1 = "Normal", 
                                          ident.2 = "High", verbose = FALSE)

write.table(CD4_MenoryT_LowNormal_marker,
            file = "CD4_MenoryT_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CD4_MenoryT_NormalHigh_marker,
            file = "CD4_MenoryT_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Th1####################
Th1<-subset(CD4T_plot,idents = c('Th1'))
Idents(Th1) <- "Immune_Response"
Th1_LowNormal_marker <- FindMarkers(Th1, ident.1 = "Low", 
                                        ident.2 = "Normal", verbose = FALSE)

Th1_NormalHigh_marker <- FindMarkers(Th1, ident.1 = "Normal", 
                                         ident.2 = "High", verbose = FALSE)

write.table(Th1_LowNormal_marker,
            file = "Th1_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Th1_NormalHigh_marker,
            file = "Th1_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Th2####################
Th2<-subset(CD4T_plot,idents = c('Th2'))
Idents(Th2) <- "Immune_Response"
Th2_LowNormal_marker <- FindMarkers(Th2, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

Th2_NormalHigh_marker <- FindMarkers(Th2, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(Th2_LowNormal_marker,
            file = "Th2_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Th2_NormalHigh_marker,
            file = "Th2_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#####################Th17####################
Th17<-subset(CD4T_plot,idents = c('Th17'))
Idents(Th17) <- "Immune_Response"
Th17_LowNormal_marker <- FindMarkers(Th17, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

Th17_NormalHigh_marker <- FindMarkers(Th17, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(Th17_LowNormal_marker,
            file = "Th17_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Th17_NormalHigh_marker,
            file = "Th17_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Tfh####################
Tfh<-subset(CD4T_plot,idents = c('Tfh'))
Idents(Tfh) <- "Immune_Response"
Tfh_LowNormal_marker <- FindMarkers(Tfh, ident.1 = "Low", 
                                     ident.2 = "Normal", verbose = FALSE)

Tfh_NormalHigh_marker <- FindMarkers(Tfh, ident.1 = "Normal", 
                                      ident.2 = "High", verbose = FALSE)

write.table(Tfh_LowNormal_marker,
            file = "Tfh_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Tfh_NormalHigh_marker,
            file = "Tfh_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)



#####################CD4_Treg1####################
CD4_Treg1<-subset(CD4T_plot,idents = c('Treg1'))
Idents(CD4_Treg1) <- "Immune_Response"
CD4_Treg1_LowNormal_marker <- FindMarkers(CD4_Treg1, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

CD4_Treg1_NormalHigh_marker <- FindMarkers(CD4_Treg1, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(CD4_Treg1_LowNormal_marker,
            file = "CD4_Treg1_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CD4_Treg1_NormalHigh_marker,
            file = "CD4_Treg1_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################CD4_Treg2####################
CD4_Treg2<-subset(CD4T_plot,idents = c('Treg2'))
Idents(CD4_Treg2) <- "Immune_Response"
CD4_Treg2_LowNormal_marker <- FindMarkers(CD4_Treg2, ident.1 = "Low", 
                                          ident.2 = "Normal", verbose = FALSE)

CD4_Treg2_NormalHigh_marker <- FindMarkers(CD4_Treg2, ident.1 = "Normal", 
                                           ident.2 = "High", verbose = FALSE)

write.table(CD4_Treg2_LowNormal_marker,
            file = "CD4_Treg2_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CD4_Treg2_NormalHigh_marker,
            file = "CD4_Treg2_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


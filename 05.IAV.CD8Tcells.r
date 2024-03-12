
setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

CD8AT<-subset(Tcells, subset = CD8A > 1) 
CD8BT<-subset(Tcells, subset = CD8B > 1)

#######################合并CD8A和CD8B################

library(dplyr)

# Pull cell names
cells1 <- colnames(CD8AT)
cells2 <- colnames(CD8BT)

# Identify shared cells between both objects
shared_cells <- intersect(x = cells1, y = cells2)

# Create Obj1 minus shared cells
CD8AT_minus <- subset(CD8AT, cells = shared_cells, invert = TRUE)

# Merge objects
CD8T <- merge(CD8AT_minus, CD8BT)

#CD8T<-SCTransform(CD8T, vars.to.regress = "percent.mt", verbose = FALSE)

############################分cluster#################################

CD8_Tcell_genes = c(
  "CD3D","CD3E","CD3G","CD8A","CD8B","IFNG","TNF","EOMES","CCR7","SELL",
  "CD27","CD28","LEF1","LTB","TCF7","NCAM1","TNFRSF8","TBX21",
  "PRDM1","IRF4","GATA3","IL2RA","FOXP3","STAT4","STAT3","CD58","CXCR3",
  "PDCD1","HAVCR2","ZNF683",
  
  "IL2",
  "IL12A","IL12B",
  "IL18",
  "STAT1")

CD8_Tcell_plot = c(
  "CD3D","CD3E","CD3G","CD8A","CD8B","IFNG","TNF","EOMES","CCR7","SELL",
  "CD27","CD28","LEF1","LTB","TCF7","NCAM1","TNFRSF8","CD69","TBX21",
  "PRDM1","IRF4","GATA3","IL2RA","FOXP3","STAT4","STAT3","CD58","CXCR3",
  "PDCD1","HAVCR2","ZNF683",

  "IL2",
  "IL12A","IL12B",
  "IL18",
  "STAT1"
)

CD8T <- RunPCA(CD8T,features = CD8_Tcell_genes,verbose = FALSE)
CD8T <- RunUMAP(CD8T, dims = 1:10, verbose = FALSE)
CD8T <- FindNeighbors(CD8T, dims = 1:10, verbose = FALSE)
CD8T <- FindClusters(CD8T, resolution=0.4, verbose = FALSE)
DimPlot(CD8T, label = TRUE,raster = T) + NoLegend() 

DotPlot(CD8T, features = CD8_Tcell_plot,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))

######################注释###############################
library("plyr")

CD8T@meta.data$self_cluster<-NULL
self_cluster = CD8T@meta.data$seurat_clusters
self_cluster<-revalue(self_cluster, c("0" = "Tc1", 
                                      "1" = "Tc2",
                                      "2" = "CTL", 
                                      "3" = "Naive T",
                                      "4" = "Naive T", 
                                      "5" = "Naive T",
                                      "6" = "Naive T", 
                                      "7" = "Naive T",
                                      "8" = "Naive T", 
                                      "9" = "Exhausted T",
                                      "10" ="CTL",
                                      "11" ="Memory T",
                                      "12" ="Unknown",
                                      "13" ="Memory T",
                                      "14" ="Memory T",
                                      "15" ="CTL",
                                      "16" = "Treg"
))


self_cluster<-factor(self_cluster,levels=c("Naive T","Memory T",
                                           "CTL","Tc1","Tc2","Treg",
                                           "Exhausted T","Unknown"))   

CD8T <- AddMetaData(object = CD8T,              
                    metadata = self_cluster,               
                    col.name = "self_cluster")    


DimPlot(CD8T, label = TRUE, group.by="self_cluster") 

DotPlot(CD8T, features = CD8_Tcell_plot,
        assay = 'SCT', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


###################提取需要的cluster####################

Idents(CD8T) <- "self_cluster"                
CD8T_plot<-subset(CD8T,idents = c("Naive T","Memory T","CTL",
                                  "Tc1","Tc2","Treg","Exhausted T"))                  

DimPlot(CD8T_plot, label = TRUE, group.by="self_cluster",raster = T) 


DotPlot(CD8T_plot, features = CD8_Tcell_plot,
        assay = 'SCT', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


########################统计cluster数目#################

table(CD8T@meta.data$self_cluster)

CD8T_SampleCell_number<-table(CD8T_plot@meta.data$Sample_ID)

IAV_SampleCell_number<-table(IAV@meta.data$Sample_ID)

cbind(IAV_SampleCell_number,CD8T_SampleCell_number)



#####################各细胞占比#####################
library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####包冲突了导致dplyr没有起作用

CD8T_metadata<-CD8T_plot@meta.data
CD8T_metadata <- CD8T_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

CD8T_metadata <- CD8T_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

CD8T_metadata <- CD8T_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

CD8T_metadata_select=CD8T_metadata [ , c("Sample_ID",
                                         "Immune_Response_Level",
                                         "self_cluster",
                                         "Sample_Celltype_Sum",
                                         "Sample_Celltype_Pre",
                                         "Sample_Sum") ]

CD8T_metadata_select<-CD8T_metadata_select[!duplicated(CD8T_metadata_select),]
CD8T_metadata_select$Immune_Response_Level<-factor(CD8T_metadata_select$Immune_Response_Level,
                                                   levels=c("Low","Normal","High"))  

library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))

ggplot(CD8T_metadata_select,aes(x = Immune_Response_Level, 
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
##############################CD8 Naive_T######################################


library("limma")
CD8_NaiveT<-subset(CD8T_plot,idents = c('Naive T'))
Idents(CD8_NaiveT) <- "Immune_Response"
CD8_NaiveT_LowNormal_marker <- FindMarkers(CD8_NaiveT, ident.1 = "Low", 
                                           ident.2 = "Normal", verbose = FALSE)

CD8_NaiveT_NormalHigh_marker <- FindMarkers(CD8_NaiveT, ident.1 = "Normal", 
                                            ident.2 = "High", verbose = FALSE)

write.table(CD8_NaiveT_LowNormal_marker,
            file = "CD8_NaiveT_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(CD8_NaiveT_NormalHigh_marker,
            file = "CD8_NaiveT_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#######################Menory T cells##################################
CD8_MenoryT<-subset(CD8T_plot,idents = c('Memory T'))
Idents(CD8_MenoryT) <- "Immune_Response"
CD8_MenoryT_LowNormal_marker <- FindMarkers(CD8_MenoryT, ident.1 = "Low", 
                                            ident.2 = "Normal", verbose = FALSE)

CD8_MenoryT_NormalHigh_marker <- FindMarkers(CD8_MenoryT, ident.1 = "Normal", 
                                             ident.2 = "High", verbose = FALSE)

write.table(CD8_MenoryT_LowNormal_marker,
            file = "CD8_MenoryT_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CD8_MenoryT_NormalHigh_marker,
            file = "CD8_MenoryT_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Tc1####################
Tc1<-subset(CD8T_plot,idents = c('Tc1'))
Idents(Tc1) <- "Immune_Response"
Tc1_LowNormal_marker <- FindMarkers(Tc1, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

Tc1_NormalHigh_marker <- FindMarkers(Tc1, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(Tc1_LowNormal_marker,
            file = "Tc1_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Tc1_NormalHigh_marker,
            file = "Tc1_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Tc2####################
Tc2<-subset(CD8T_plot,idents = c('Tc2'))
Idents(Tc2) <- "Immune_Response"
Tc2_LowNormal_marker <- FindMarkers(Tc2, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

Tc2_NormalHigh_marker <- FindMarkers(Tc2, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(Tc2_LowNormal_marker,
            file = "Tc2_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Tc2_NormalHigh_marker,
            file = "Tc2_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#####################CTL####################
CTL<-subset(CD8T_plot,idents = c('CTL'))
Idents(CTL) <- "Immune_Response"
CTL_LowNormal_marker <- FindMarkers(CTL, ident.1 = "Low", 
                                     ident.2 = "Normal", verbose = FALSE)

CTL_NormalHigh_marker <- FindMarkers(CTL, ident.1 = "Normal", 
                                      ident.2 = "High", verbose = FALSE)

write.table(CTL_LowNormal_marker,
            file = "CTL_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CTL_NormalHigh_marker,
            file = "CTL_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Exhausted T####################
Exhausted_T<-subset(CD8T_plot,idents = c('Exhausted T'))
Idents(Exhausted_T) <- "Immune_Response"
Exhausted_T_LowNormal_marker <- FindMarkers(Exhausted_T, ident.1 = "Low", 
                                    ident.2 = "Normal", verbose = FALSE)

Exhausted_T_NormalHigh_marker <- FindMarkers(Exhausted_T, ident.1 = "Normal", 
                                     ident.2 = "High", verbose = FALSE)

write.table(Exhausted_T_LowNormal_marker,
            file = "Exhausted_T_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Exhausted_T_NormalHigh_marker,
            file = "Exhausted_T_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)



#######################CD8_Treg#######################
CD8_Treg<-subset(CD8T_plot,idents = c('Treg'))
Idents(CD8_Treg) <- "Immune_Response"
CD8_Treg_LowNormal_marker <- FindMarkers(CD8_Treg, ident.1 = "Low", 
                                          ident.2 = "Normal", verbose = FALSE)

CD8_Treg_NormalHigh_marker <- FindMarkers(CD8_Treg, ident.1 = "Normal", 
                                           ident.2 = "High", verbose = FALSE)

write.table(CD8_Treg_LowNormal_marker,
            file = "CD8_Treg_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(CD8_Treg_NormalHigh_marker,
            file = "CD8_Treg_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)




VlnPlot(CD8T, features = "CD8A",slot = "scale.data")

#https://github.com/satijalab/seurat/issues/3761


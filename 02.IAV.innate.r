

setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/IAV_new_V1.rds")


Innate<-subset(IAV,idents = c('0','1','10','26','29','37','41','42','43','44','51','57','58'))  ###37不应该加进去
#Innate_1<-subset(IAV, subset = CD3D < 1) 
#Innate_2<-subset(Innate_1, subset = CD3E < 1) 
#Innate_3<-subset(Innate_2, subset = CD3G < 1) 
#Innate_4<-subset(Innate_3, subset = CD19 < 1) 
#Innate_5<-subset(Innate_4, subset = CD4 < 1) 
#Innate_6<-subset(Innate_5, subset = CD8A < 1) 
#Innate_7<-subset(Innate_6, subset = CD8B < 1) 
#Innate_1<-NULL
#Innate_2<-NULL
#Innate_3<-NULL
#Innate_4<-NULL
#Innate_5<-NULL
#Innate_6<-NULL
#gc()

Innate_genes = c(
  "PTPRC",
  "CD68",
  "CSF1R",
  "NCAM1",
  "ITGA2",
  "CD14",
  "FCGR3A",
  "ITGAX",
  "PPBP",
  "PF4",
  "SIGLEC10",
  "BST2",
  "CD34",
  
  "ITGAM",
  "KLRF1",
  "FLT3",
  "CD1C",
  "CLEC10A",
  "HLA-DRA",
  "CD163",
  "MRC1"
  )
  
#"CD4",
#"CD8A",
#"CD4B",
#"CD19"


#######CD14单核 CD16单核 巨噬细胞 干细胞  血小板 NK细胞 pDC cDC######

#Innate_SCT<-Innate
Innate<-Innate_SCT
#Innate <- SCTransform(Innate, vars.to.regress = "percent.mt", verbose = FALSE) 
Innate <- RunPCA(Innate,features = Innate_genes, verbose = FALSE)
Innate <- RunUMAP(Innate, dims = 1:10, verbose = FALSE)
Innate <- FindNeighbors(Innate, dims = 1:10, verbose = FALSE)
Innate <- FindClusters(Innate, resolution=0.2, verbose = FALSE)
DimPlot(Innate, label = TRUE) + NoLegend()


DotPlot(Innate, features = Innate_genes,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))

###################注释###############################
library("plyr")

Innate@meta.data$self_cluster<-NULL
self_cluster = Innate@meta.data$seurat_clusters
self_cluster<-revalue(self_cluster, c("0" = "Natural Killer", 
                                      "1" = "Non-classical Monocytes",
                                      "2" = "Macrophage", 
                                      "3" = "Conventional DC1",
                                      "4" = "Non-classical Monocytes", 
                                      "5" = "Non-classical Monocytes",
                                      "6" = "Non-classical Monocytes", 
                                      "7" = "Conventional DC2",
                                      "8" = "Unknown", 
                                      "9" = "Unknown",
                                      "10" = "Non-classical Monocytes", 
                                      "11" = "Classical Monocytes",
                                      "12" = "Plasmacytoid DC", 
                                      "13" = "Unknown",
                                      "14" = "Platelets", 
                                      "15" = "Conventional DC1"))
  
self_cluster<-factor(self_cluster,levels=c("Natural Killer","Classical Monocytes",
                                           "Non-classical Monocytes",
                                           "Macrophage","Conventional DC1",
                                           "Conventional DC2",
                                           "Plasmacytoid DC","Platelets","Unknown"))       

Innate <- AddMetaData(object = Innate,              
                       metadata = self_cluster,               
                       col.name = "self_cluster")    

DimPlot(Innate, label = TRUE)                    
DimPlot(Innate, label = TRUE, group.by="self_cluster") 

###################提取需要的cluster####################

Idents(Innate) <- "self_cluster"                
Innate_plot<-subset(Innate,idents = c("Natural Killer","Classical Monocytes","Non-classical Monocytes",
                                      "Macrophage","Conventional DC1","Conventional DC2",
                                      "Plasmacytoid DC","Platelets"
                                      ))                  

DimPlot(Innate_plot, label = TRUE, group.by="self_cluster") 
             
DotPlot(Innate_plot, features = Innate_genes,
        assay = 'SCT', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


###############统计cluster数目#################

table(Innate_plot@meta.data$self_cluster)

Innate_SampleCell_number<-table(Innate_plot@meta.data$Sample_ID)

IAV_SampleCell_number<-table(IAV@meta.data$Sample_ID)

cbind(IAV_SampleCell_number,Innate_SampleCell_number)


###############各细胞占比#################
library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####太好了，果然是包冲突了导致dplyr没有起作用

Innate_metadata<-Innate_plot@meta.data
Innate_metadata <- Innate_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

Innate_metadata <- Innate_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

Innate_metadata <- Innate_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

Innate_metadata_select=Innate_metadata [ , c("Sample_ID",
                                             "Immune_Response_Level",
                                             "self_cluster",
                                             "Sample_Celltype_Sum",
                                             "Sample_Celltype_Pre",
                                             "Sample_Sum") ]
      

Innate_metadata_select<-Innate_metadata_select[!duplicated(Innate_metadata_select),]


library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))


ggplot(Innate_metadata_select,aes(x = Immune_Response_Level, 
                                  y = Sample_Celltype_Pre,
                                  fill = Immune_Response_Level)) +
  facet_wrap(~self_cluster,scales = "free",ncol=4,strip.position = "top")+
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
                     label = "p.signif") +
  #ylim(0.5,3.5)
  ###添加标题
  labs(x="Differences in antibody titers and cell proportions among Inactivated vaccines")+
  theme(axis.title.x = element_text(color="black", size=12, face="bold")
  )





#########################差异基因####################################

library("limma")
Innate_NK<-subset(Innate_plot,idents = c('Natural Killer'))
Idents(Innate_NK) <- "Immune_Response"
Innate_NK_LowNormal_marker <- FindMarkers(Innate_NK, ident.1 = "Low", 
                                          ident.2 = "Normal", verbose = FALSE)

Innate_NK_NormalHigh_marker <- FindMarkers(Innate_NK, ident.1 = "Normal", 
                                           ident.2 = "High", verbose = FALSE)

write.table(Innate_NK_LowNormal_marker,
            file = "Innate_NK_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(Innate_NK_NormalHigh_marker,
            file = "Innate_NK_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#######################cMono##################################
Innate_cMono<-subset(Innate_plot,idents = c('Classical Monocytes'))
Idents(Innate_cMono) <- "Immune_Response"
Innate_cMono_LowNormal_marker <- FindMarkers(Innate_cMono, ident.1 = "Low", 
                                          ident.2 = "Normal", verbose = FALSE)

Innate_cMono_NormalHigh_marker <- FindMarkers(Innate_cMono, ident.1 = "Normal", 
                                           ident.2 = "High", verbose = FALSE)

write.table(Innate_cMono_LowNormal_marker,
            file = "Innate_cMono_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Innate_cMono_NormalHigh_marker,
            file = "Innate_cMono_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################NcMono####################
Innate_NcMono<-subset(Innate_plot,idents = c('Non-classical Monocytes'))
Idents(Innate_NcMono) <- "Immune_Response"
Innate_NcMono_LowNormal_marker <- FindMarkers(Innate_NcMono, ident.1 = "Low", 
                                             ident.2 = "Normal", verbose = FALSE)

Innate_NcMono_NormalHigh_marker <- FindMarkers(Innate_NcMono, ident.1 = "Normal", 
                                              ident.2 = "High", verbose = FALSE)

write.table(Innate_NcMono_LowNormal_marker,
            file = "Innate_NcMono_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Innate_NcMono_NormalHigh_marker,
            file = "Innate_NcMono_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)



#####################Macrophage##############################

Innate_Macrophage<-subset(Innate_plot,idents = c('Macrophage'))
Idents(Innate_Macrophage) <- "Immune_Response"
Innate_Macrophage_LowNormal_marker <- FindMarkers(Innate_Macrophage, ident.1 = "Low", 
                                              ident.2 = "Normal", verbose = FALSE)

Innate_Macrophage_NormalHigh_marker <- FindMarkers(Innate_Macrophage, ident.1 = "Normal", 
                                               ident.2 = "High", verbose = FALSE)

write.table(Innate_Macrophage_LowNormal_marker,
            file = "Innate_Macrophage_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Innate_Macrophage_NormalHigh_marker,
            file = "Innate_Macrophage_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################cDC1##############################

Innate_cDC1<-subset(Innate_plot,idents = c('Conventional DC1'))
Idents(Innate_cDC1) <- "Immune_Response"
Innate_cDC1_LowNormal_marker <- FindMarkers(Innate_cDC1, ident.1 = "Low", 
                                                  ident.2 = "Normal", verbose = FALSE)

Innate_cDC1_NormalHigh_marker <- FindMarkers(Innate_cDC1, ident.1 = "Normal", 
                                                   ident.2 = "High", verbose = FALSE)

write.table(Innate_cDC1_LowNormal_marker,
            file = "Innate_cDC1_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Innate_cDC1_NormalHigh_marker,
            file = "Innate_cDC1_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################cDC2##############################

Innate_cDC2<-subset(Innate_plot,idents = c('Conventional DC2'))
Idents(Innate_cDC2) <- "Immune_Response"
Innate_cDC2_LowNormal_marker <- FindMarkers(Innate_cDC2, ident.1 = "Low", 
                                            ident.2 = "Normal", verbose = FALSE)

Innate_cDC2_NormalHigh_marker <- FindMarkers(Innate_cDC2, ident.1 = "Normal", 
                                             ident.2 = "High", verbose = FALSE)

write.table(Innate_cDC2_LowNormal_marker,
            file = "Innate_cDC2_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(Innate_cDC2_NormalHigh_marker,
            file = "Innate_cDC2_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

###################Plasmacytoid DC#######################

Innate_pDC<-subset(Innate_plot,idents = c('Plasmacytoid DC'))
Idents(Innate_pDC) <- "Immune_Response"
Innate_pDC_LowNormal_marker <- FindMarkers(Innate_pDC, ident.1 = "Low", 
                                            ident.2 = "Normal", verbose = FALSE)

Innate_pDC_NormalHigh_marker <- FindMarkers(Innate_pDC, ident.1 = "Normal", 
                                             ident.2 = "High", verbose = FALSE)

write.table(Innate_pDC_LowNormal_marker,
            file = "Innate_pDC_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


write.table(Innate_pDC_NormalHigh_marker,
            file = "Innate_pDC_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


########################Platelets##################################


Innate_Platelets<-subset(Innate_plot,idents = c('Platelets'))
Idents(Innate_Platelets) <- "Immune_Response"
Innate_Platelets_LowNormal_marker <- FindMarkers(Innate_Platelets, ident.1 = "Low", 
                                           ident.2 = "Normal", verbose = FALSE)

Innate_Platelets_NormalHigh_marker <- FindMarkers(Innate_Platelets, ident.1 = "Normal", 
                                            ident.2 = "High", verbose = FALSE)

write.table(Innate_Platelets_LowNormal_marker,
            file = "Innate_Platelets_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(Innate_Platelets_NormalHigh_marker,
            file = "Innate_Platelets_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)





################如何展示差异基因？下面这些方法都不好，最好还是用dotplot###############
VlnPlot(Innate_NK, features = "JUN",group.by = "Sample_ID", 
        split.by = "Immune_Response",slot = "data")

library(dplyr)
detach('package:plyr') ####太好了，果然是包冲突了导致dplyr没有起作用
Innate_NK_metadata<-Innate_NK@meta.data
Innate_NK_metadata <- Innate_NK_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))
Innate_NK_metadata <- Innate_NK_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Mean = mean(nCount_SCT))

Innate_NK_metadata_select=Innate_NK_metadata [ , c("Sample_ID",
                                             "Immune_Response",
                                             "Sample_Sum",
                                             "Sample_Mean") ]
Innate_NK_metadata_select<-Innate_NK_metadata_select[!duplicated(Innate_NK_metadata_select),]

#install.packages("remotes")
#remotes::install_github("xmc811/Scillus")

write.table(Innate_NK_metadata_select,
            file = "Innate_NK_metadata_select.txt",
            quote = F,sep = "\t",
            row.names = T)
######################以上方法都不好，下面这个好#############

NK_DGE<-c("IL32",
          "GZMH",
          "KLRC2",
          "GBP1",
          "HLA-DRB5",
          "CCL5",
          "CD52",
          "TXNIP",
          "IL7R",
          "ZFP36L2",
          "CHST2",
          "FOS",
          "JUN",
          "DUSP1",
          "RPS10-NUDT3",
          "CTSW",
          "GNLY",
          "TNFAIP3",
          "EIF1",
          "DDX3Y",
          "PTPRCAP",
          "HLA-B",
          "NFKBIA",
          "KIR2DL1",
          "CORO1B",
          "MT-CO1",
          "BTN3A2",
          "MYOM2"
)


DotPlot(Innate_NK, features = NK_DGE,
        assay = 'SCT', group.by = "Immune_Response") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


library(dittoSeq)

DotPlot(Innate_plot, features = "MT-CO1",cols = dittoColors(),
        assay = 'SCT', group.by = "self_cluster",split.by = "Immune_Response") + 
  coord_flip() + RotatedAxis()
 

#################解决了分面的问题，但是不同细胞的趋势被拉平了#################

p1<-dittoDotPlot(Innate_NK, c("MT-CO1","GBP1","CD52","JUN",
                            "FOS","DUSP1","RPS10-NUDT3","GNLY",
                            "HLA-B","NFKBIA"),
             group.by = "Immune_Response",split.by = "self_cluster")+####太好了，完美
  coord_flip()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


p2<-dittoDotPlot(Innate_cDC1, c("MT-CO1","GBP1","CD52","JUN",
                              "FOS","DUSP1","RPS10-NUDT3","GNLY",
                              "HLA-B","NFKBIA"),
                 group.by = "Immune_Response",split.by = "self_cluster")+####太好了，完美
  coord_flip()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


p1+p2




p1<-VlnPlot(Innate_pDC,cols = NULL,features="HLA-DQB1",
        pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+xlab("Plasmacytioid DC")+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 



p2<-VlnPlot(Innate_cDC1,cols = NULL,features="FCER1G",
            pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+xlab("Conventional DC1")+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 



p3<-VlnPlot(Innate_cDC2,cols = NULL,features="HLA-DQB1",
        pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+xlab("Conventional DC2")+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 




VlnPlot(Innate_Macrophage,cols = NULL,features="RPL10",
        pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +xlab("Macrophage")



VlnPlot(Innate_Macrophage,cols = NULL,features="RNF213",
        pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 



VlnPlot(Innate_Macrophage,cols = NULL,features="LTA4H",
        pt.size = 0,group.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,6)+
  NoLegend()+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 




p1+p2+p3
plot_grid(p1,p2,p3, ncol=3, nrow = 1)


VlnPlot(Innate_plot,cols = NULL,features="MT-CO1",
        pt.size = 0,group.by = "self_cluster",split.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,5)+
  NoLegend()+
  geom_boxplot(width = .05,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') 
 

#############################################################################################

# 加载数据对象
#pbmc <- readRDS("pbmc_2k_1.rds")
#identity <- readRDS("pbmc_2k_2.rds")

features <- c("MT-CO1")

# 提取数据框子集
Innate_plot_data <- Innate_plot[,features]

# 添加细胞name和细胞群身份
Innate_plot_data$Cell <- rownames(Innate_plot)


# melt数据变长
pbmc <- reshape2::melt(Innate_plot_data,
                       id.vars = c("Cell","Self_cluster"),
                       measure.vars = features,
                       variable.name = "Self_cluster",
                       value.name = "Expr")
# 查看数据前10行
head(pbmc, 10)


















#########################################

# fig 2D 

p1=DotPlot(Innate_cMono, features =  "TMSB10", group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Classical Monocytes")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+
          ylab("")+xlab("")
p1

p2=DotPlot(Innate_Macrophage, features =  c("RPL10","LTA4H","IFITM2","MT-CO1","CCL5","HLA-DRB1","IFITM3","RNF213")
           , group.by = "Immune_Response"
           ) + 
  coord_flip()+labs(title = "Macrophage")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab("")+xlab("")
  p2

p3=DotPlot(Innate_cDC1, features =  c("FCER1G","RPS10-NUDT3"), group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Conventional DC1")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+ 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+ylab("")+xlab("")
  p3

p4=DotPlot(Innate_cDC2, features =  c("HLA-DQB1","DUSP1"), group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Conventional DC2")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+ylab("")+xlab("")
  p4

p5=DotPlot(Innate_pDC, features = "HLA-DQB1", group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Plasmacytioid DC")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+ 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+ylab("")+xlab("")
  p5

p6=DotPlot(Innate_Platelets, features =  "B2M", group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Platelets")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+ylab("")+xlab("")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))
  p6

  
library(cowplot)
plot_grid(p1,p5,p6,p3,p4,p2, ncol=3, nrow = 2,rel_heights=c(1.5,2))


#######################Dotplot的对比还是明显很多#############################


DotPlot(Innate_cMono, features =  "TMSB10", group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Classical Monocytes")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab("")+xlab("")


VlnPlot(Innate_cMono, cols = NULL, features = "TMSB10", pt.size = 0, 
        group.by = "Immune_Response") +
  geom_boxplot(width=.05,col="black",fill="white",outlier.size =0)+ 
  NoLegend() 




VlnPlot(Innate_Macrophage, cols = NULL,pt.size = 0,
        features =  c("RPL10","LTA4H","IFITM2","MT-CO1","CCL5","HLA-DRB1","IFITM3","RNF213"), 
        group.by = "Immune_Response")+
          geom_boxplot(width=.05,col="black",fill="white",outlier.size =0)+ 
          NoLegend() 
        

###Dotplot的对比还是明显很多，小提琴图看不出什么差别，应该查Dotplot和vlnplot的差别。

#########https://github.com/satijalab/seurat/issues/2798##########
#########差别其实就是Dotplot做了Z检验转化，明显一些#############
###########有些基因比如RNF213的确还是可以从图中看出来的#############



##################早期测试代码,不再使用#####################################
#####################使用高变基因而不是自定义基因来聚类分析#######################


Innate <- RunPCA(Innate, verbose = FALSE)http://10.16.68.5:8787/graphics/plot_zoom_png?width=919&height=613
Innate <- RunUMAP(Innate, dims = 1:10, verbose = FALSE)
Innate <- FindNeighbors(Innate, dims = 1:10, verbose = FALSE)
Innate <- FindClusters(Innate, resolution=0.2, verbose = FALSE)
DimPlot(Innate, label = TRUE) + NoLegend()


DotPlot(Innate, features = Innate_genes,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


####################scytype用于辅助判断############################
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData =  Innate[["SCT"]]@scale.data, 
                      scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(Innate@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Innate@meta.data[Innate@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), 
                  scores = es.max.cl, ncells = sum(Innate@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


Innate@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Innate@meta.data$customclassif[Innate@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}



DimPlot(Innate, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  





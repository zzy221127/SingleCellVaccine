



setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/IAV_new_V1.rds")

Bcells<-subset(IAV,idents = c('5','6','11','15','33','46','52','53','54','56','59','60','61','62')) 
Bcells_new<-subset(IAV, subset = CD19 > 1) 
#Bcells_new <- SCTransform(Bcells_new, vars.to.regress = "percent.mt", verbose = FALSE) 

Bcell_genes = c(
  "CD19","CD22","BCL6",
  "CD27","CD38","SDC1",
  "HBB","SLAMF7","IL6",
  "MS4A1","CD40","CD80",
  "PDCD1LG2","CXCR3","CXCR4",
  "CXCR5","CXCR6","IGHM","IGHG1",
  "IGHG2","IGHG3","IGHG4","IGHA1",
  "IGHA2","HLA-DRA","TNFRSF13B",
  "PAX5","FAS","PTPRJ","IGHD",
  "CD1A","CD1B","CD1C","CD1D","CD5",
  "CR2","CD24","TLR4","IL10","TGFB1"
)

Bcell_genes_plot = c(
  "CD19","CD22","BCL6",
  "CD27","CD38","SDC1",
  "HBB","SLAMF7","IL6",
  "MS4A1","CD40","CD80",
  "PDCD1LG2","CXCR3","CXCR4",
  "CXCR5","CXCR6","IGHM","IGHG1",
  "IGHG2","IGHG3","IGHG4","IGHA1",
  "IGHA2","HLA-DRA","TNFRSF13B",
  "PAX5","FAS","PTPRJ","IGHD",
  "CD1A","CD1B","CD1C","CD1D","CD5",
  "CR2","CD24","TLR4","IL10",
  "TGFB1","IFNA1","IFNB1","IFNG"
)



################################


#Bcells <- SCTransform(Bcells, vars.to.regress = "percent.mt", verbose = FALSE) 

Bcells_new <- RunPCA(Bcells_new,features = Bcell_genes, verbose = FALSE)
Bcells_new <- RunUMAP(Bcells_new, dims = 1:15, verbose = FALSE)
Bcells_new <- FindNeighbors(Bcells_new, dims = 1:15, verbose = FALSE)
Bcells_new <- FindClusters(Bcells_new, resolution=0.2, verbose = FALSE)
DimPlot(Bcells_new, label = TRUE) + NoLegend()

DotPlot(Bcells_new, features = Bcell_genes,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


Bcells <- RunPCA(Bcells,features = Bcell_genes, verbose = FALSE)
Bcells <- RunUMAP(Bcells, dims = 1:15, verbose = FALSE)
Bcells <- FindNeighbors(Bcells, dims = 1:15, verbose = FALSE)
Bcells <- FindClusters(Bcells, resolution=0.2, verbose = FALSE)
DimPlot(Bcells, label = TRUE) + NoLegend()

DotPlot(Bcells, features = Bcell_genes_plot,
        assay = 'SCT', group.by = "seurat_clusters") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))



###################注释###############################
library("plyr")

Bcells_new@meta.data$self_cluster<-NULL
self_cluster = Bcells_new@meta.data$seurat_clusters
self_cluster<-revalue(self_cluster, c("0" = "Naïve B cells", 
                                      "1" = "Naïve B cells",
                                      "2" = "Menory B cells", 
                                      "3" = "Menory B cells",
                                      "4" = "Plasma cells1", 
                                      "5" = "Plasma cells2",
                                      "6" = "Unknown", 
                                      "7" = "Unknown",
                                      "8" = "Regulatory B cells", 
                                      "9" = "Plasma cells3"
                                  ))
                          

self_cluster<-factor(self_cluster,levels=c("Naïve B cells",
                                           "Menory B cells",
                                           "Plasma cells1",
                                           "Plasma cells2",
                                           "Plasma cells3","Regulatory B cells","Unknown"))       

Bcells_new <- AddMetaData(object = Bcells_new,              
                      metadata = self_cluster,               
                      col.name = "self_cluster")    


#DimPlot(Bcells_new, label = TRUE)                    
DimPlot(Bcells_new, label = TRUE, group.by="self_cluster") 


###################提取需要的cluster####################

Idents(Bcells_new) <- "self_cluster"                
Bcells_plot<-subset(Bcells_new,idents = c("Naïve B cells","Menory B cells","Plasma cells1",
                                      "Plasma cells2","Plasma cells3","Regulatory B cells"
))                  

DimPlot(Bcells_plot, label = TRUE, group.by="self_cluster") 


DotPlot(Bcells_plot, features = Bcell_genes_plot,
        assay = 'RNA', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))  ###########注意，这里的assay 从SCT变成了RNA。


IFNI_genes<-c("IRF1","STAT1","UBA52","UBC","XAF1","DUSP1",
              "FOS","HLA-DPA1","NFKBIA","ZFP36","GBP4",
              "JUNB","PARP14","CD83")


DotPlot(Bcells_plot, features = IFNI_genes,
        assay = 'RNA', group.by = "self_cluster") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))  ###########注意，这里的assay 从SCT变成了RNA。



DotPlot(Plasma2, features = "HLA-DQA2",
        assay = 'SCT', group.by = "Immune_Response") + 
  coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))  





VlnPlot(Plasma2,cols = NULL,features="HLA-DQA2",
        pt.size = 0,group.by = "Immune_Response")+
  ylim(0,5)+NoLegend()+geom_boxplot(width = .05,show.legend = F,
                                    position = position_dodge(0.9),
                                    color = 'grey20',alpha = 0.5,
                                    outlier.color = 'grey50') 



VlnPlot(Plasma1,cols = NULL,features="FOSB",
            pt.size = 0,group.by = "Immune_Response",split.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,5)+NoLegend()+geom_boxplot(width = .05,show.legend = F,
                                    position = position_dodge(0.9),
                                    color = 'grey20',alpha = 0.5,
                                    outlier.color = 'grey50') 






########################差异基因小提琴图##############################

# 加载数据对象


# 提取数据框子集
Bcells_plot_IFNIgenes <- Bcells_plot[,IFNI_genes]

# 添加细胞name和细胞群身份
pbmc$Cell <- rownames(pbmc)
pbmc$Idents <- identity

# melt数据变长
pbmc <- reshape2::melt(pbmc,
                       id.vars = c("Cell","Idents"),
                       measure.vars = features,
                       variable.name = "Feat",
                       value.name = "Expr")
# 查看数据前10行
head(pbmc, 10)



ggplot(Bcells_plot, aes(Feat, Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0),
                     position="right",
                     labels = function(x)
                       c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(6, 6, 0, 6, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("Stack violin-plot with annotation version 1") +
  ylab("Expression Level")







###############统计cluster数目#################

table(Bcells_plot@meta.data$self_cluster)

Bcells_SampleCell_number<-table(Bcells_plot@meta.data$Sample_ID)

IAV_SampleCell_number<-table(IAV@meta.data$Sample_ID)

cbind(IAV_SampleCell_number,Bcells_SampleCell_number)


###############各细胞占比#################
library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####包冲突了导致dplyr没有起作用

Bcells_metadata<-Bcells_plot@meta.data
Bcells_metadata <- Bcells_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

Bcells_metadata <- Bcells_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

Bcells_metadata <- Bcells_metadata %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

Bcells_metadata_select=Bcells_metadata [ , c("Sample_ID",
                                             "Immune_Response_Level",
                                             "self_cluster",
                                             "Sample_Celltype_Sum",
                                             "Sample_Celltype_Pre",
                                             "Sample_Sum") ]

Bcells_metadata_select<-Bcells_metadata_select[!duplicated(Bcells_metadata_select),]


library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))

ggplot(Bcells_metadata_select,aes(x = Immune_Response_Level, 
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
##############################Naïve_B######################################


library("limma")
Naïve_B<-subset(Bcells_plot,idents = c('Naïve B cells'))
Idents(Naïve_B) <- "Immune_Response"
Naïve_B_LowNormal_marker <- FindMarkers(Naïve_B, ident.1 = "Low", 
                                          ident.2 = "Normal", verbose = FALSE)

Naïve_B_NormalHigh_marker <- FindMarkers(Naïve_B, ident.1 = "Normal", 
                                           ident.2 = "High", verbose = FALSE)

write.table(Naïve_B_LowNormal_marker,
            file = "Naive_B_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(Naïve_B_NormalHigh_marker,
            file = "Naive_B_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

#######################Menory B cells##################################
Menory_B<-subset(Bcells_plot,idents = c('Menory B cells'))
Idents(Menory_B) <- "Immune_Response"
Menory_B_LowNormal_marker <- FindMarkers(Menory_B, ident.1 = "Low", 
                                             ident.2 = "Normal", verbose = FALSE)

Menory_B_NormalHigh_marker <- FindMarkers(Menory_B, ident.1 = "Normal", 
                                              ident.2 = "High", verbose = FALSE)

write.table(Menory_B_LowNormal_marker,
            file = "Menory_B_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Menory_B_NormalHigh_marker,
            file = "Menory_B_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Plasma cells1####################
Plasma1<-subset(Bcells_plot,idents = c('Plasma cells1'))
Idents(Plasma1) <- "Immune_Response"
Plasma1_LowNormal_marker <- FindMarkers(Plasma1, ident.1 = "Low", 
                                              ident.2 = "Normal", verbose = FALSE)

Plasma1_NormalHigh_marker <- FindMarkers(Plasma1, ident.1 = "Normal", 
                                               ident.2 = "High", verbose = FALSE)

write.table(Plasma1_LowNormal_marker,
            file = "Plasma1_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Plasma1_NormalHigh_marker,
            file = "Plasma1_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)



#################Plasma cells2#################################

Plasma2<-subset(Bcells_plot,idents = c('Plasma cells2'))
Idents(Plasma2) <- "Immune_Response"
Plasma2_LowNormal_marker <- FindMarkers(Plasma2, ident.1 = "Low", 
                                                  ident.2 = "Normal", verbose = FALSE)

Plasma2_NormalHigh_marker <- FindMarkers(Plasma2, ident.1 = "Normal", 
                                                   ident.2 = "High", verbose = FALSE)

write.table(Plasma2_LowNormal_marker,
            file = "Plasma2_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Plasma2_NormalHigh_marker,
            file = "Plasma2_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Plasma cells3##############################

Plasma3<-subset(Bcells_plot,idents = c('Plasma cells3'))
Idents(Plasma3) <- "Immune_Response"
Plasma3_LowNormal_marker <- FindMarkers(Plasma3, ident.1 = "Low", 
                                            ident.2 = "Normal", verbose = FALSE)

Plasma3_NormalHigh_marker <- FindMarkers(Plasma3, ident.1 = "Normal", 
                                             ident.2 = "High", verbose = FALSE)

write.table(Plasma3_LowNormal_marker,
            file = "Plasma3_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)
write.table(Plasma3_NormalHigh_marker,
            file = "Plasma3_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)


#####################Regulatory B cells##############################

Regulatory_B<-subset(Bcells_plot,idents = c('Regulatory B cells'))
Idents(Regulatory_B) <- "Immune_Response"
Regulatory_B_LowNormal_marker <- FindMarkers(Regulatory_B, ident.1 = "Low", 
                                            ident.2 = "Normal", verbose = FALSE)

Regulatory_B_NormalHigh_marker <- FindMarkers(Regulatory_B, ident.1 = "Normal", 
                                             ident.2 = "High", verbose = FALSE)

write.table(Regulatory_B_LowNormal_marker,
            file = "Regulatory_B_LowNormal_marker.txt",
            quote = F,sep = "\t",
            row.names = T)

write.table(Regulatory_B_NormalHigh_marker,
            file = "Regulatory_B_NormalHigh_marker.txt",
            quote = F,sep = "\t",
            row.names = T)



#saveRDS(Bcells_new,file="/data/01.singlecell_vaccine/07.newbin/10.RDSfile/Bcells_new.rds")

#############################B cells vlnplot################################

#https://github.com/satijalab/seurat/issues/3761

#VlnPlot(Bcells_new, features = "CD19")
#VlnPlot(Bcells_new, features = "CD19",slot = "scale.data")
#saveRDS(Bcells_new,file="/data/01.singlecell_vaccine/07.newbin/10.RDSfile/Bcells_new.rds")

###Bcells_new<-NormalizeData(Bcells_new)  ###normalized 
####其实是因为单一指标替代了，所以scale.data没有东西了，
#####然后SCT的结果小提琴图本来就分层，而Dotplot做了Z检验的转化所以看起来比较和谐
##Note that this single command replaces NormalizeData() , ScaleData() , and FindVariableFeatures() 
###https://github.com/satijalab/seurat/issues/2023 问题基本无解
###https://github.com/satijalab/seurat/issues/7699


#Bcells_new2 <- NormalizeData(Bcells_new)  ###试过了，不用normalize就是还是分层
#all.features <- rownames(Bcells_new2)
#Bcells_new2 <- ScaleData(Bcells_new2, features = all.features)
#VlnPlot(Bcells_new2, features = "CD19",slot = "scale.data")  ###额，也不好看,data,counts都不行

#VlnPlot(Bcells_new2, cols = NULL, features = "CD19", pt.size = 0, 
#                  group.by = "self_cluster") +
#  geom_boxplot(width=.05,col="black",fill="white",outlier.size =0)+ 
#  NoLegend() #############这样看起来还行，也仅仅是还行而已。



library(ggplot2)

#c("HLA-DQA2","HLA-DRB5","RPS10-NUDT3","MT-CO1",
#           "FOSB","DUSP1","FOS","TXNIP","RPL36A-HNRNPH2",
#           "JUND","HLA-DQB1","PPP1R15A","EIF1","JUN",
#           "BTNL9","PSAP","IGHA1")

p1<-VlnPlot(Naïve_B,cols = NULL,features="HLA-DQA2",
        pt.size = 0,group.by = "Immune_Response",split.by = "Immune_Response")+
  stat_compare_means(comparisons = comparelist,label = "p.signif")+
  ylim(0,5)+NoLegend()+geom_boxplot(width = .05,show.legend = F,
                                    position = position_dodge(0.9),
                                    color = 'grey20',alpha = 0.5,
                                    outlier.color = 'grey50') 


p1
library(patchwork)
p1 + p2







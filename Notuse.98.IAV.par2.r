
setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/IAV_new_V1.rds")
meta.data<-IAV@meta.data
write.table(meta.data,file = "IAV.meta.data.new.txt",quote = F,sep = "\t",row.names = F)



#########结果2部分#######################################
#########提取APC的细胞组###############################



HLA_genes<-c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
             "HLA-DRA","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
             "HLA-DRB5","HLA-DQA2","HLA-DMB","HLA-DMA","HLA-DOB","HLA-DOA"
)


CD14mono <- subset(IAV, idents = c('42'))
CD16mono  <- subset(IAV, idents = c('10','26','51'))
cDC <- subset(IAV, idents = c('57'))
pDC<- subset(IAV, idents = c('41'))
NaiveB <- subset(IAV, idents = c('5','6','11','15','33'))
MemoryB <- subset(IAV, idents = c('46','52','53','54','56','59','60','61','62'))


# fig 2A 

p1=DotPlot(CD14mono, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "CD14+ Mono")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
p1

p2=DotPlot(CD16mono, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "CD16+ Mono")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+

p2

p3=DotPlot(cDC, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Classical DC")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+

p3

p4=DotPlot(pDC, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Plasmacytoid DC")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
p4

p5=DotPlot(NaiveB, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Naïve B")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+

p5

p6=DotPlot(MemoryB, features =  HLA_genes, group.by = "Immune_Response") + 
  coord_flip()+labs(title = "Memory B")+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+

p6

library(cowplot)
Fig2A = plot_grid(p1,p2,p3,p4,p5,p6, ncol=3, nrow = 2)
Fig2A


########################差异分析####################################

levels(IAV)
Idents(IAV) <- "Immune_Response"
levels(IAV)

IAV_LowNormal_marker <- FindMarkers(IAV, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
IAV_LowHigh_marker <- FindMarkers(IAV, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
IAV_NormalHigh_marker <- FindMarkers(IAV, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)

write.table(IAV_LowNormal_marker,file = "IAV_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(IAV_LowHigh_marker,file = "IAV_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(IAV_NormalHigh_marker,file = "IAV_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)


######################################先天免疫细胞#######################################
DC<-subset(IAV,idents = c('41','57'))
Idents(DC) <- "Immune_Response"
DC_LowNormal_marker <- FindMarkers(DC, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
DC_LowHigh_marker <- FindMarkers(DC, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
DC_NormalHigh_marker <- FindMarkers(DC, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(DC_LowNormal_marker,file = "DC_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(DC_LowHigh_marker,file = "DC_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(DC_NormalHigh_marker,file = "DC_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)



cDC<-subset(IAV,idents = c('57'))
Idents(cDC) <- "Immune_Response"
cDC_LowNormal_marker <- FindMarkers(cDC, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
cDC_LowHigh_marker <- FindMarkers(cDC, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
cDC_NormalHigh_marker <- FindMarkers(cDC, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(cDC_LowNormal_marker,file = "cDC_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(cDC_LowHigh_marker,file = "cDC_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(cDC_NormalHigh_marker,file = "cDC_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)


pDC<-subset(IAV,idents = c('41'))
Idents(pDC) <- "Immune_Response"
pDC_LowNormal_marker <- FindMarkers(pDC, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
pDC_LowHigh_marker <- FindMarkers(pDC, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
pDC_NormalHigh_marker <- FindMarkers(pDC, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(pDC_LowNormal_marker,file = "pDC_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(pDC_LowHigh_marker,file = "pDC_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(pDC_NormalHigh_marker,file = "pDC_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)




#HLA-DQB1,DUSP1

#VlnPlot(DC, features = c("HLA-DQB1","DUSP1"))
#DC_check_genes<-c("HLA-DQB1","DUSP1")
vln.dat=FetchData(DC,c(HLA_genes,"Immune_Response","Sample_ID"))

vln.dat$Cell <- rownames(vln.dat)
#宽转长######################这段代码有问题，但是原因不明，之前的细胞数目也是计算错了。
library(dplyr)
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","Immune_Response","Sample_ID"), 
                               measure.vars = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
                                                "HLA-DRA","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
                                                "HLA-DRB5","HLA-DQA2","HLA-DMB","HLA-DMA","HLA-DOB","HLA-DOA"),
                               variable.name = "HLA_genes", 
                               value.name = "SCT_counts")  %>% 
  group_by(HLA_genes,Sample_ID) %>% #分组 
  mutate(Sample_SCT_counts=mean(SCT_counts)) #计算均值


vln.dat.melt.used<-unique(subset(vln.dat.melt,
                                 select = c("Sample_ID","Immune_Response","Sample_SCT_counts","HLA_genes")))


library("ggpubr")

vln.dat.melt.used$Immune_Response<-factor(vln.dat.melt.used$Immune_Response,levels=c("Low","Normal","High"))
comparelist<-list(c("Low", "Normal"), 
                  c("Low", "High"), c("Normal", "High"))


# plot
ggplot(vln.dat.melt.used,aes(x = Immune_Response, y = Sample_SCT_counts,fill = Immune_Response)) +
  facet_wrap(~HLA_genes,scales = "free",strip.position = "top")+
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
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0,color = 'black'),#,hjust = 1
        legend.position = 'top',
        aspect.ratio = 0.8) + ##图形长宽比例
  # 颜色设置
  scale_fill_manual(values = c('Low'='#398AB9','High'='red','Normal'='grey'),
                    name = '') +
  # 添加显著性标记
  stat_compare_means(comparisons = comparelist,vjust=0.8,method="t.test"
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                     #                symbols = c("***", "**", "*", "ns")),
                     #label = "p.signif"
                     ) +
  #ylim(0.5,3.5)
  ###添加标题
  labs(x="DC genes")+
  theme(axis.title.x = element_text(color="black", size=12, face="bold")
  )




######################找差异基因######################################


CD16_Mono<-subset(IAV, idents = c('10','26','51'))
Idents(CD16_Mono) <- "Immune_Response"
CD16_Mono_LowNormal_marker <- FindMarkers(CD16_Mono, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
CD16_Mono_LowHigh_marker <- FindMarkers(CD16_Mono, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
CD16_Mono_NormalHigh_marker <- FindMarkers(CD16_Mono, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(CD16_Mono_LowNormal_marker,file = "CD16_Mono_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(CD16_Mono_LowHigh_marker,file = "CD16_Mono_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(CD16_Mono_NormalHigh_marker,file = "CD16_Mono_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)


CD14_Mono<-subset(IAV, idents = c('42'))
Idents(CD14_Mono) <- "Immune_Response"
CD14_Mono_LowNormal_marker <- FindMarkers(CD14_Mono, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
CD14_Mono_LowHigh_marker <- FindMarkers(CD14_Mono, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
CD14_Mono_NormalHigh_marker <- FindMarkers(CD14_Mono, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(CD14_Mono_LowNormal_marker,file = "CD14_Mono_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(CD14_Mono_LowHigh_marker,file = "CD14_Mono_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(CD14_Mono_NormalHigh_marker,file = "CD14_Mono_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)



Macrophage<-subset(IAV, idents = c('29','43'))  ###跟在mono的后面分析####
Idents(Macrophage) <- "Immune_Response"
Macrophage_LowNormal_marker <- FindMarkers(Macrophage, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
Macrophage_LowHigh_marker <- FindMarkers(Macrophage, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
Macrophage_NormalHigh_marker <- FindMarkers(Macrophage, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(Macrophage_LowNormal_marker,file = "Macrophage_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(Macrophage_LowHigh_marker,file = "Macrophage_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(Macrophage_NormalHigh_marker,file = "Macrophage_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)


NK<-subset(IAV, idents = c('0','1','44'))
Idents(NK) <- "Immune_Response"
NK_LowNormal_marker <- FindMarkers(NK, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
NK_LowHigh_marker <- FindMarkers(NK, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
NK_NormalHigh_marker <- FindMarkers(NK, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(NK_LowNormal_marker,file = "NK_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(NK_LowHigh_marker,file = "NK_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(NK_NormalHigh_marker,file = "NK_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)


Platelets<-subset(IAV, idents = c('58'))  
Idents(Platelets) <- "Immune_Response"
Platelets_LowNormal_marker <- FindMarkers(Platelets, ident.1 = "Low", ident.2 = "Normal", verbose = FALSE)
Platelets_LowHigh_marker <- FindMarkers(Platelets, ident.1 = "Low", ident.2 = "High", verbose = FALSE)
Platelets_NormalHigh_marker <- FindMarkers(Platelets, ident.1 = "Normal", ident.2 = "High", verbose = FALSE)
write.table(Platelets_LowNormal_marker,file = "Platelets_LowNormal_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(Platelets_LowHigh_marker,file = "Platelets_LowHigh_marker.txt",quote = F,sep = "\t",row.names = T)
write.table(Platelets_NormalHigh_marker,file = "Platelets_NormalHigh_marker.txt",quote = F,sep = "\t",row.names = T)



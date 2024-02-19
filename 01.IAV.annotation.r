
setwd("/data/01.singlecell_vaccine")
setwd("/data/01.singlecell_vaccine/07.newbin")


library(Seurat)
library(ggplot2)

IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v5.rds")

#先统计输出，再画boxplot:不需要，信息都在metadata里面
#meta.data<-IAV@meta.data
#write.table(meta.data,file = "IAV.meta.data.txt",quote = F,sep = "\t",row.names = F)


############手动注释#####################################
###################使用自定义的marker###################################

IAV_new<-IAV

#IAV_new <- SCTransform(IAV_new, vars.to.regress = "percent.mt", verbose = FALSE)
IAV_new <- RunPCA(IAV_new,features = gene1, verbose = FALSE)
IAV_new <- RunUMAP(IAV_new, dims = 1:20, verbose = FALSE)##测试过了dim为1：25，效果不好。
IAV_new <- FindNeighbors(IAV_new, dims = 1:20, verbose = FALSE)
IAV_new <- FindClusters(IAV_new, verbose = FALSE)
DimPlot(IAV_new, label = TRUE) + NoLegend()

##############为了画doHeatmap加入的代码##################
######https://zhuanlan.zhihu.com/p/574125880
IAV_new <- NormalizeData(IAV_new)
IAV_new <- FindVariableFeatures(IAV_new)
IAV_new <- ScaleData(IAV_new,verbose = FALSE,features = rownames(IAV_new),
                     vars.to.regress = 'self_cluster')

DoHeatmap(IAV_new,features = gene1,group.by="self_cluster",
          cells = 1:5000) #####2023年10月28，FigS1使用。

##+scale_fill_gradientn(colors=c("gold2","white","coral1"))


##############################################################

gene1 = c("PTPRC",#CD45，广泛表达于免疫细胞中
          "ITGAM",#尤其是在单核细胞和巨噬细胞中，ITGAM的表达非常显著，是其特异性标记之一。
          "CD68",#主要表达于单核细胞和巨噬细胞，是巨噬细胞特异性标记物
          "CD163",#CD163是巨噬细胞的标志物，在树突状细胞中不常见。
          "MRC1",#CD206（或MRC1）是巨噬细胞的标志物，在树突状细胞中不常见。
          "CSF1R",#CSF1R是巨噬细胞的标志物，树突状细胞中不常见。
          "FLT3",#FLT3是树突状细胞的标志物，在巨噬细胞中不常见。
          "NCAM1",#CD56，在外周血中主要表达于自然杀伤细胞（NK细胞）
          "KLRF1",#参与激活自然杀伤细胞
          "KLRB1",#CD161，可用于鉴定多种免疫细胞，包括自然杀伤细胞和一部分T细胞亚群
          "CD14",#单核细胞和巨噬细胞
          "FCGR3A",#CD16，主要表达于自然杀伤细胞、巨噬细胞和一部分淋巴细胞
          "CD1C",#主要在树突状细胞（dendritic cells，DCs）中表达 #B-1b细胞和一部分非调节性B-2细胞
          "CLEC10A",#主要在树突状细胞表面表达
          "PPBP",#主要表达于血小板和中性粒细胞
          "PF4",#主要表达于血小板
          "ITGA2",#CD49b
          
          "CD19",#主要在B细胞表面表达，是B细胞特异性标记物
          "MS4A1",#编码CD20，主要在B细胞中表达，是B细胞的标记物之一
          "CD38",#在免疫细胞中广泛表达，可用于鉴定淋巴细胞、浆细胞等。
          "CD27",#CD27是一种共刺激分子，主要表达于记忆性B细胞和Plasma cell上。
          "CR2",#编码CD21，广泛表达于B细胞和某些T细胞亚群，是B细胞的标记物之一
          "IGHM",#在成熟B细胞中表达
          "IGHG1",#
          "IGHG2",#
          "IGHG3",#
          "IGHG4",#
          "IGHA1",#
          "IGHA2",#
          "SDC1",#CD138 区分plasma B
          "SLAMF7",#CD319广泛表达于多种免疫细胞，尤其在浆细胞（Plasma cell）和某些淋巴细胞亚群中表达较高。
          "TNFRSF17",#BCMA，浆细胞和记忆B细胞
          "CXCR5",#CXCR5在中间记忆B细胞中表达
          "HLA-DRA",#HLA-DRA，滤泡B细胞，记忆
          "TNFRSF13B",#TACI，CD267，滤泡B细胞，记忆
          "PAX5",#滤泡B细胞，记忆
          "FAS",#CD95，记忆B
          "PTPRJ",#CD148，记忆B         
          
          
          "CD3D",#
          "CD3E",#
          "CD3G",#
          "CD4",#
          "CD8A",#主要表达于细胞毒性T细胞（CTLs）。
          "CD8B",#主要表达于细胞毒性T细胞（CTLs）。
          "CCR7",#CCR7是一种受体，参与T细胞的迁移至淋巴结。Naive CD8+T细胞通常表达CCR7，而已激活或记忆的CD8+T细胞则表达较低水平或不表达CCR7。。
          "SELL",# 编码选择素L（CD62L），已激活或记忆的CD8+T细胞通常在表面上不表达CD62L。
          "FOXP3",#是调节性T细胞（Treg细胞）特异性转录因子，在Treg细胞中表达。
          "CTLA4",#是T细胞的共刺激分子，主要在激活的T细胞表面表达
          "IL2RA",#编码IL-2受体α链（CD25），主要在激活的T细胞和调节性T细胞中表达。
          "STAT3",#是信号转导和转录激活因子3，广泛在免疫细胞中表达。
          "STAT4",#是信号转导和转录激活因子4，主要在Th1型辅助T细胞中表达。
          "TBX21",#编码T-box转录因子21，主要在Th1型辅助T细胞中表达。
          "CD28",#是T细胞的共刺激分子，主要在激活的T细胞表面表达。
          "CD58",#编码淋巴细胞功能相关抗原3（LFA-3），在多种免疫细胞中表达。
          "BCL6B",#是B细胞淋巴瘤/白血病6家族成员B，在多种免疫细胞中表达。
          "DPP4",#编码二肽基肽酶4（CD26），广泛在多种细胞类型中表达。
          "ITGAE",#编码整合素αE亚单位（CD103），在树突状细胞和一部分淋巴细胞表面表达。
          "CXCR3",#编码CXC趋化因子受体3，参与细胞迁移和定向，主要在多种免疫细胞中表达。
          "ZNF683")#编码核转录因子ZNF683，在Th17型辅助T细胞中表达。



IAV_gene_plot = c("PTPRC",#CD45，广泛表达于免疫细胞中
          "ITGAM",#尤其是在单核细胞和巨噬细胞中，ITGAM的表达非常显著，是其特异性标记之一。
          "CD68",#主要表达于单核细胞和巨噬细胞，是巨噬细胞特异性标记物
          "CD163",#CD163是巨噬细胞的标志物，在树突状细胞中不常见。
          "MRC1",#CD206（或MRC1）是巨噬细胞的标志物，在树突状细胞中不常见。
          "CSF1R",#CSF1R是巨噬细胞的标志物，树突状细胞中不常见。
          "FLT3",#FLT3是树突状细胞的标志物，在巨噬细胞中不常见。
          "NCAM1",#CD56，在外周血中主要表达于自然杀伤细胞（NK细胞）
          "KLRF1",#参与激活自然杀伤细胞
          "KLRB1",#CD161，可用于鉴定多种免疫细胞，包括自然杀伤细胞和一部分T细胞亚群
          "CD14",#单核细胞和巨噬细胞
          "FCGR3A",#CD16，主要表达于自然杀伤细胞、巨噬细胞和一部分淋巴细胞
          "CLEC10A",#主要在树突状细胞表面表达
          "PPBP",#主要表达于血小板和中性粒细胞
          "PF4",#主要表达于血小板

          "CD19",#主要在B细胞表面表达，是B细胞特异性标记物
          "MS4A1",#编码CD20，主要在B细胞中表达，是B细胞的标记物之一
          "IGHM",#在成熟B细胞中表达
          "IGHG1",#
          "IGHG2",#
          "IGHG3",#
          "IGHG4",#
          "IGHA1",#
          "IGHA2",#
 
          "CD3D",#
          "CD3E",#
          "CD3G",#
          "CD4",#
          "CD8A",#主要表达于细胞毒性T细胞（CTLs）。
          "CD8B",#主要表达于细胞毒性T细胞（CTLs）。
          "CCR7",#CCR7是一种受体，参与T细胞的迁移至淋巴结。Naive CD8+T细胞通常表达CCR7，而已激活或记忆的CD8+T细胞则表达较低水平或不表达CCR7。。
          "SELL",# 编码选择素L（CD62L），已激活或记忆的CD8+T细胞通常在表面上不表达CD62L。
          "FOXP3",#是调节性T细胞（Treg细胞）特异性转录因子，在Treg细胞中表达。
          "CTLA4",#是T细胞的共刺激分子，主要在激活的T细胞表面表达
          "IL2RA",#编码IL-2受体α链（CD25），主要在激活的T细胞和调节性T细胞中表达。
          "ZNF683")#编码核转录因子ZNF683，在Th17型辅助T细胞中表达。


DotPlot(IAV, features = unique(gene1),
        assay = 'SCT', group.by = "seurat_clusters") + coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


DimPlot(IAV_new, label = TRUE) 

self_cluster = IAV_new@meta.data$seurat_clusters

library("plyr")

self_cluster<-revalue(self_cluster, c("0" = "Natural killer", 
                                      "1" = "Natural killer",
                                      "2" = "CD8+ Naïve T", 
                                      "3" = "CD4+ NKT",
                                      "4" = "CD8+ Naïve T", 
                        "5" = "Naïve B",
                        "6" = "Naïve B", 
                        "7" = "CD8+ CTL",
                        "8" = "CD4+ Treg", 
                        "9" = "CD4+ Naïve T",
                        "10" = "CD16+ Mono", 
                        "11" = "Naïve B",
                        "12" = "CD4+ Memory", 
                        "13" = "CD4+ Memory",
                        "14" = "CD4+ Memory", 
                        "15" = "Naïve B",
                        "16" = "CD8+ NKT", 
                        "17" = "CD8+ Naïve T",
                        "18" = "DN Tregs", 
                        "19" = "CD4+ Memory",
                        "20" = "CD4+ Naïve T", 
                        "21" = "DN Tregs",
                        "22" = "CD4+ Memory", 
                        "23" = "CD4+ Naïve T",
                        "24" = "DN Tregs", 
                        "25" = "CD4+ Th1",
                        "26" = "CD16+ Mono", 
                        "27" = "DN NKT",
                        "28" = "CD4+ Memory", 
                        "29" = "Macrophage",
                        "30" = "DN NKT",
                        "31" = "CD4+ Naïve T",
                        "32" = "Unknown", 
                        "33" = "Naïve B",
                        "34" = "CD8+ SP", 
                        "35" = "CD4+ Naïve T",
                        "36" = "CD4+ Naïve T", 
                        "37" = "CD4+ Memory",
                        "38" = "CD4+ Treg", 
                        "39" = "CD4+ Naïve T",
                        "40" = "CD4+ Naïve T", 
                        "41" = "Dendritic cells",
                        "42" = "CD14+ Mono", 
                        "43" = "Macrophage",
                        "44" = "Natural killer", 
                        "45" = "CD8+ Naïve T",
                        "46" = "Memory B", 
                        "47" = "CD4+ Treg",
                        "48" = "CD4+ Naïve T", 
                        "49" = "CD4+ Naïve T",
                        "50" = "CD4+ Naïve T", 
                        "51" = "CD16+ Mono",
                        "52" = "Memory B", 
                        "53" = "Memory B",
                        "54" = "Memory B", 
                        "55" = "CD4+ Naïve T",
                        "56" = "Memory B",
                        "57" = "Dendritic cells",
                        "58" = "Platelets", 
                        "59" = "Memory B",
                        "60" = "Memory B", 
                        "61" = "Memory B",
                        "62" = "Memory B"
             ))


self_cluster<-factor(self_cluster,levels=c("CD14+ Mono", "CD16+ Mono","CD4+ Memory",
                                         "CD4+ Naïve T","CD4+ NKT","CD4+ Th1","CD4+ Treg","CD8+ CTL",
                                         "CD8+ Naïve T","CD8+ NKT","CD8+ SP","Dendritic cells",
                                         "DN NKT","DN Tregs","Macrophage","Memory B",
                                         "Naïve B","Natural killer","Platelets","Unknown"))


IAV_new <- AddMetaData(object = IAV_new,              
                   metadata = self_cluster,               
                   col.name = "self_cluster")    



library(RColorBrewer)
DimPlot(IAV_new, label = TRUE,group.by = 'self_cluster') 
DotPlot(IAV_new, features = unique(gene1),
        assay = 'SCT', group.by = "self_cluster") + coord_flip() + RotatedAxis()+
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))


gene_for_Fig1D<-c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FCGR3A","NCAM1",
                  "CD14","PPBP","CD79A","MS4A1","SERPINF1","CD1C","IGHM")

FeaturePlot(IAV_new, features =gene_for_Fig1D,
            cols=c("grey",(brewer.pal(n=9,name="Reds"))))
  

FeaturePlot(IAV, features =gene_for_Fig1D,ncol = 3,
            cols=c("grey",(brewer.pal(n=9,name="Reds"))))


Sample_ID = IAV_new@meta.data$SampleID
Sample_ID<-revalue(Sample_ID, c(
  "20S6575980"="IAV-005",
  "20S6575981"="IAV-024",
  "20S6575983"="IAV-006",
  "20S6575984"="IAV-007",
  "20S6575985"="IAV-008",
  "20S6575987"="IAV-009",
  "20S6575989"="IAV-010",
  "20S6575992"="IAV-011",
  "20S6575993"="IAV-012",
  "20S6575994"="IAV-013",
  "PS191478902"="IAV-014",
  "PS191478904"="IAV-015",
  "PS191478905"="IAV-025",
  "PS191478906"="IAV-016",
  "PS191478907"="IAV-017",
  "PS191478908"="IAV-018",
  "PS191478909"="IAV-019",
  "PS191478910"="IAV-001",
  "PS191478911"="IAV-020",
  "PS191478913"="IAV-021",
  "PS191478914"="IAV-002",
  "PS191478915"="IAV-003",
  "PS191478916"="IAV-026",
  "PS191478917"="IAV-027",
  "PS191478919"="IAV-022",
  "PS191478924"="IAV-023",
  "PS191478925"="IAV-004"))

Sample_ID<-factor(Sample_ID,levels=c(
  "IAV-001",
  "IAV-002",
  "IAV-003",
  "IAV-004",
  "IAV-005",
  "IAV-006",
  "IAV-007",
  "IAV-008",
  "IAV-009",
  "IAV-010",
  "IAV-011",
  "IAV-012",
  "IAV-013",
  "IAV-014",
  "IAV-015",
  "IAV-016",
  "IAV-017",
  "IAV-018",
  "IAV-019",
  "IAV-020",
  "IAV-021",
  "IAV-022",
  "IAV-023",
  "IAV-024",
  "IAV-025",
  "IAV-026",
  "IAV-027"
  ))


IAV_new <- AddMetaData(object = IAV_new,              
                       metadata = Sample_ID,               
                       col.name = "Sample_ID")    



Immune_Response = IAV_new@meta.data$SampleID
Immune_Response<-revalue(Immune_Response, c(
  "20S6575980"="Normal",
  "20S6575981"="High",
  "20S6575983"="Normal",
  "20S6575984"="Normal",
  "20S6575985"="Normal",
  "20S6575987"="Normal",
  "20S6575989"="Normal",
  "20S6575992"="Normal",
  "20S6575993"="Normal",
  "20S6575994"="Normal",
  "PS191478902"="Normal",
  "PS191478904"="Normal",
  "PS191478905"="High",
  "PS191478906"="Normal",
  "PS191478907"="Normal",
  "PS191478908"="Normal",
  "PS191478909"="Normal",
  "PS191478910"="Low",
  "PS191478911"="Normal",
  "PS191478913"="Normal",
  "PS191478914"="Low",
  "PS191478915"="Low",
  "PS191478916"="High",
  "PS191478917"="High",
  "PS191478919"="Normal",
  "PS191478924"="Normal",
  "PS191478925"="Low"
  ))

Immune_Response<-factor(Immune_Response,levels=c("Low","Normal","High"))

IAV_new <- AddMetaData(object = IAV_new,              
                       metadata = Immune_Response,               
                       col.name = "Immune_Response_Level")    



###############FIG S1A#################
library(reshape2)
library(tidyverse)
#提取样本和细胞数据，并且进行长宽数据转换

plotC<-IAV@meta.data
plotC <- plotC %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))
plotC <- plotC %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))
plotC <- plotC %>% group_by(Sample_ID,self_cluster) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

plotC1=plotC [ , c("Sample_ID","Immune_Response_Level", "Sample_Celltype_Pre","self_cluster") ]
plotC1<-plotC1[!duplicated(plotC1),]



#绘制每个组织中细胞数目柱状图
ggplot(data = plotC1, aes(x = Sample_ID, y = Sample_Celltype_Pre, fill = self_cluster)) +
  geom_bar(stat = "identity", width=0.8,aes(group=self_cluster),position="stack")+
  #  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))

#绘制每个组织中细胞比例柱状图

ggplot(data = plotC1, aes(x = Sample_ID, y = Sample_Celltype_Pre, fill = self_cluster)) +
  geom_bar(stat = "identity", width=0.8,aes(group=self_cluster),position="fill")+
  #  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  #  scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))#让横轴上的标签倾斜45度


saveRDS(IAV_new, file = "/data/01.singlecell_vaccine/05.rds.files/IAV_new_V1.rds", compress = F)
#结果1部分RDS保存


#############细胞占比小提琴图############################


ggplot(plotC1,aes(x = Immune_Response_Level, y = Sample_Celltype_Pre,fill = Immune_Response_Level)) +
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
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0,color = 'black'),#,hjust = 1
        legend.position = 'top',
        aspect.ratio = 0.8) + ##图形长宽比例
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


#####################cell classifications#################


cell_classifications = IAV@meta.data$seurat_clusters

library("plyr")

cell_classifications<-revalue(cell_classifications, c("0" = "Innate immune cells", 
                                      "1" = "Innate immune cells",
                                      "2" = "T cells", 
                                      "3" = "T cells",
                                      "4" = "T cells", 
                                      "5" = "B cells",
                                      "6" = "B cells", 
                                      "7" = "T cells",
                                      "8" = "Innate immune cells", 
                                      "9" = "T cells",
                                      "10" = "Innate immune cells", 
                                      "11" = "B cells",
                                      "12" = "T cells", 
                                      "13" = "T cells",
                                      "14" = "T cells", 
                                      "15" = "B cells",
                                      "16" = "T cells", 
                                      "17" = "T cells",
                                      "18" = "T cells", 
                                      "19" = "T cells",
                                      "20" = "T cells", 
                                      "21" = "Innate immune cells",
                                      "22" = "T cells", 
                                      "23" = "T cells",
                                      "24" = "Innate immune cells", 
                                      "25" = "T cells",
                                      "26" = "Innate immune cells", 
                                      "27" = "T cells",
                                      "28" = "T cells", 
                                      "29" = "Innate immune cells",
                                      "30" = "T cells",
                                      "31" = "T cells",
                                      "32" = "T cells", 
                                      "33" = "B cells",
                                      "34" = "T cells", 
                                      "35" = "T cells",
                                      "36" = "T cells", 
                                      "37" = "T cells",
                                      "38" = "Innate immune cells", 
                                      "39" = "T cells",
                                      "40" = "T cells", 
                                      "41" = "Innate immune cells",
                                      "42" = "Innate immune cells", 
                                      "43" = "Innate immune cells",
                                      "44" = "Innate immune cells", 
                                      "45" = "T cells",
                                      "46" = "B cells", 
                                      "47" = "Innate immune cells",
                                      "48" = "T cells", 
                                      "49" = "T cells",
                                      "50" = "T cells", 
                                      "51" = "Innate immune cells",
                                      "52" = "B cells", 
                                      "53" = "B cells",
                                      "54" = "B cells", 
                                      "55" = "T cells",
                                      "56" = "B cells",
                                      "57" = "Innate immune cells",
                                      "58" = "Innate immune cells", 
                                      "59" = "B cells",
                                      "60" = "B cells", 
                                      "61" = "B cells",
                                      "62" = "B cells"
))


cell_classifications<-factor(cell_classifications,levels=c("Innate immune cells", "B cells","T cells"))


IAV <- AddMetaData(object = IAV,              
                       metadata = cell_classifications,               
                       col.name = "cell_classifications")    



library(RColorBrewer)
DimPlot(IAV, label = T,group.by = 'cell_classifications') 
DotPlot(IAV, features = unique(gene1),
        assay = 'SCT', group.by = "cell_classifications") +
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1))


##################################重新画细胞占比图##################################


library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####太好了，果然是包冲突了导致dplyr没有起作用

IAV_metadata<-IAV@meta.data
IAV_metadata <- IAV_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

IAV_metadata <- IAV_metadata %>% group_by(Sample_ID,cell_classifications) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

IAV_metadata <- IAV_metadata %>% group_by(Sample_ID,cell_classifications) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

IAV_metadata_select=IAV_metadata [ , c("Sample_ID",
                                             "Immune_Response_Level",
                                             "cell_classifications",
                                             "Sample_Celltype_Sum",
                                             "Sample_Celltype_Pre",
                                             "Sample_Sum") ]


IAV_metadata_select<-IAV_metadata_select[!duplicated(IAV_metadata_select),]


library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))

ggplot(IAV_metadata_select,aes(x = Immune_Response_Level, 
                                  y = Sample_Celltype_Pre,fill = Immune_Response_Level)) +
  facet_wrap(~cell_classifications,scales = "free",ncol=3,strip.position = "top")+
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
  stat_compare_means(comparisons = comparelist,vjust=0.8,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif") +
  #ylim(0.5,3.5)
  ###添加标题
  labs(x="Differences in antibody titers and cell proportions among Inactivated vaccines")+
  theme(axis.title.x = element_text(color="black", size=12, face="bold")
  )



###########################################################################################
###############################第三次重新注释###########################################

cell_class3 = IAV@meta.data$seurat_clusters
#detach('package:dplyr')
library("plyr")

cell_class3<-revalue(cell_class3, c(
  '0'='Natural killer',
  '1'='Natural killer',
  '2'='CD8+ Naïve T',
  '3'='Th1',
  '4'='CD8+ Naïve T',
  '5'='Naïve B',
  '6'='Memory B',
  '7'='CTL',
  '8'='DC',
  '9'='CD4+ Naïve T',
  '10'='Monocyte',
  '11'='Naïve B',
  '12'='CD4+ Naïve T',
  '13'='CD4+ Naïve T',
  '14'='Th2',
  '15'='Naïve B',
  '16'='CTL',
  '17'='CD8+ Naïve T',
  '18'='Tc1',
  '19'='CD4+ Memory',
  '20'='Th17',
  '21'='Monocyte',
  '22'='CD4+ Memory',
  '23'='CD4+ Memory',
  '24'='Monocyte',
  '25'='CD8+ Memory',
  '26'='Monocyte',
  '27'='CD8+ Memory',
  '28'='CD4+ Memory',
  '29'='Macrophage',
  '30'='Tc1',
  '31'='CD4+ Naïve T',
  '32'='Exhausted T',
  '33'='Naïve B',
  '34'='CD8+ Memory',
  '35'='CD4+ Naïve T',
  '36'='CD4+ Naïve T',
  '37'='CD4+ Memory',
  '38'='DC',
  '39'='CD4+ Naïve T',
  '40'='CD4+ Naïve T',
  '41'='DC',
  '42'='Monocyte',
  '43'='Macrophage',
  '44'='Tc2',
  '45'='CD8+ Naïve T',
  '46'='Plasma',
  '47'='DC',
  '48'='CD4+ Naïve T',
  '49'='CD4+ Naïve T',
  '50'='Treg',
  '51'='Monocyte',
  '52'='Plasma',
  '53'='Plasma',
  '54'='Plasma',
  '55'='CD4+ Naïve T',
  '56'='Plasma',
  '57'='DC',
  '58'='Platelets',
  '59'='Memory B',
  '60'='Memory B',
  '61'='Memory B',
  '62'='Regulatory B'
))


cell_class3<-factor(cell_class3,levels=c(
"Natural killer","Monocyte","Macrophage","DC","Platelets","Naïve B","Memory B","Plasma","Regulatory B",
"CD4+ Naïve T","CD4+ Memory",	"Tfh","Th1","Th2","Th17",
"CD8+ Naïve T",	"CD8+ Memory","CTL","Tc1","Tc2","Treg","Exhausted T","Unknown"))


IAV <- AddMetaData(object = IAV,              
                   metadata = cell_class3,               
                   col.name = "cell_class3")    


library(RColorBrewer)
DimPlot(IAV, label = T,group.by = 'cell_class3') 
DotPlot(IAV, features = unique(gene1),
        assay = 'SCT', group.by = "cell_class3") +
  scale_colour_gradientn(colours=rev(brewer.pal(n=10,name="RdBu")))+
  coord_flip() + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1))


##################################重新画细胞占比图##################################


library(reshape2)
library(tidyverse)
library(dplyr)
detach('package:plyr') ####太好了，果然是包冲突了导致dplyr没有起作用

IAV_metadata<-IAV@meta.data
IAV_metadata <- IAV_metadata %>% group_by(Sample_ID)  %>% mutate(Sample_Sum = sum(nCount_SCT))

IAV_metadata <- IAV_metadata %>% group_by(Sample_ID,cell_class3) %>% 
  mutate(Sample_Celltype_Sum = sum(nCount_SCT))

IAV_metadata <- IAV_metadata %>% group_by(Sample_ID,cell_class3) %>% 
  mutate(Sample_Celltype_Pre = Sample_Celltype_Sum/Sample_Sum)

IAV_metadata_select=IAV_metadata [ , c("Sample_ID",
                                       "Immune_Response_Level",
                                       "cell_class3",
                                       "Sample_Celltype_Sum",
                                       "Sample_Celltype_Pre",
                                       "Sample_Sum") ]


IAV_metadata_select<-IAV_metadata_select[!duplicated(IAV_metadata_select),]


library(ggpubr)
comparelist<- list(c("Low", "Normal"), c("Normal", "High"), c("Low","High"))

ggplot(IAV_metadata_select,aes(x = Immune_Response_Level, 
                               y = Sample_Celltype_Pre,fill = Immune_Response_Level)) +
  facet_wrap(~cell_class3,scales = "free",ncol=4,strip.position = "top")+
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
    #    aspect.ratio = 0.8
        ) + ##图形长宽比例
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










  
  
  

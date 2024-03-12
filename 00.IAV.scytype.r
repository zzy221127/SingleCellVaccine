

setwd("/data/01.singlecell_vaccine/07.newbin")
library(tidyverse)
library(Seurat)

#IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge.rds")
#control <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/control_v1.rds")
#mRNA<- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/mRNAP1_all.combined1.rds")
#IAV3 <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/vaccine.rds")

IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v2.rds")

################0611###########合并metadata,已完成，不需要重复##############################

# From https://github.com/satijalab/seurat/issues/3260

# Change orig.ident column to factor so that it can be joined later
IAV$orig.ident <- as.factor(IAV$orig.ident)
IAV$Library<-IAV$library
IAV$library<-NULL

# Pull existing meta data where samples are specified by orig.ident and remove everything but orig.ident
IAV_meta <- IAV@meta.data %>% 
  rownames_to_column("Barcodes")

SampleInfo<-read.table("/data/01.singlecell_vaccine/08.pub.analysis/01.IAV/lib-sample-mathc.txt",header = T,sep = "\t")

full_new_meta <- full_join(x = IAV_meta, y = SampleInfo,by="Library" )  %>% 
  column_to_rownames("Barcodes") %>% 
  select(-orig.ident)

IAV <- AddMetaData(object = IAV, metadata = full_new_meta)
IAV$Barcodes<-NULL
IAV[["percent.mt"]] <- PercentageFeatureSet(object = IAV, pattern = "^MT-")
saveRDS(IAV, file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v2.rds", compress = F)

###############确定注释哪些细胞，如何注释###################################


#setp1:过滤
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(IAV, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(IAV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

VlnPlot(object = IAV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IAV <- subset(IAV, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
VlnPlot(object = IAV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


IAV <- SCTransform(IAV, vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(IAV, file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v3.rds", compress = F)

IAV <- RunPCA(IAV, verbose = FALSE)
IAV <- RunUMAP(IAV, dims = 1:30, verbose = FALSE)

IAV <- FindNeighbors(IAV, dims = 1:30, verbose = FALSE)
IAV <- FindClusters(IAV, verbose = FALSE)
DimPlot(IAV, label = TRUE) + NoLegend()

saveRDS(IAV, file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v4.rds", compress = F)



####################注释######################################################
IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v4.rds")

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
es.max = sctype_score(scRNAseqData =  IAV[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(IAV@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(IAV@meta.data[IAV@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), 
                  scores = es.max.cl, ncells = sum(IAV@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


IAV@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  IAV@meta.data$customclassif[IAV@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}



DimPlot(IAV, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
DimPlot(IAV, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',split.by="Antibody.Titer.LeveL")  

################################
###得出每个cluster的marker：

top10 <- head(VariableFeatures(IAV), 10)
plot1 <- VariableFeaturePlot(IAV)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

library(future)
plan()
plan("multiprocess", workers = 80)
plan()
levels(IAV)
IAV.markers <- FindAllMarkers(IAV, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)

write.table(IAV.markers,
            file="IAV.marker.txt",quote=F,sep="\t")

write.table(cL_resutls,file = "IAV.cluster.txt",quote = F,sep = "\t")


###################区分B细胞####################

VlnPlot(IAV, features = c("CD27"), slot = "counts", log = TRUE)  ##SCT 方法储存在 scale.data slot 另外不定义slot的时候默认应该用的是SCT的值

###230817,最新IAV的rds文件
saveRDS(IAV, file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v5.rds")
IAV2<-readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/81libs_rmdoublet_merge_v5.rds")

###########################不分LOW median High，此段不用####################################
######################LOW##################################

obj.list <- SplitObject(IAV, split.by = "Antibody.Titer.LeveL")

IAV_low<-obj.list$LOW
IAV_median<-obj.list$Median
IAV_high<-obj.list$High


es.max.low = sctype_score(scRNAseqData =  IAV_low[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls_low = do.call("rbind", lapply(unique(IAV_low@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(IAV_low@meta.data[IAV_low@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), 
                  scores = es.max.cl, ncells = sum(IAV_low@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores_low = cL_resutls_low %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores_low$type[as.numeric(as.character(sctype_scores_low$scores)) < sctype_scores_low$ncells/4] = "Unknown"
print(sctype_scores_low[,1:3])


IAV_low@meta.data$customclassif = ""
for(j in unique(sctype_scores_low$cluster)){
  cl_type = sctype_scores_low[sctype_scores_low$cluster==j,]; 
  IAV_low@meta.data$customclassif[IAV_low@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


DimPlot(IAV_low, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  


######################median##################################################

es.max.median = sctype_score(scRNAseqData =  IAV_median[["SCT"]]@scale.data, scaled = TRUE, 
                             gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls_median = do.call("rbind", lapply(unique(IAV_median@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(IAV_median@meta.data[IAV_median@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), 
                  scores = es.max.cl, ncells = sum(IAV_median@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores_median = cL_resutls_median %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores_median$type[as.numeric(as.character(sctype_scores_median$scores)) < sctype_scores_median$ncells/4] = "Unknown"
print(sctype_scores_median[,1:3])


IAV_median@meta.data$customclassif = ""
for(j in unique(sctype_scores_median$cluster)){
  cl_type = sctype_scores_median[sctype_scores_median$cluster==j,]; 
  IAV_median@meta.data$customclassif[IAV_median@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(IAV_median, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

################################high#########################################################

es.max.high= sctype_score(scRNAseqData =  IAV_high[["SCT"]]@scale.data, scaled = TRUE, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls_high = do.call("rbind", lapply(unique(IAV_high@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(IAV_high@meta.data[IAV_high@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), 
                  scores = es.max.cl, ncells = sum(IAV_high@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores_high = cL_resutls_high %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores_high$type[as.numeric(as.character(sctype_scores_high$scores)) < sctype_scores_high$ncells/4] = "Unknown"
print(sctype_scores_high[,1:3])


IAV_high@meta.data$customclassif = ""
for(j in unique(sctype_scores_high$cluster)){
  cl_type = sctype_scores_high[sctype_scores_high$cluster==j,]; 
  IAV_high@meta.data$customclassif[IAV_high@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(IAV_high, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

###########################不分LOW median High，此段不用####################################
####################测试scCATCH###########################################
#install.packages("scCATCH")
library("scCATCH")

obj_low <- createscCATCH(data = IAV_low[['SCT']]@scale.data, cluster = as.character(Idents(IAV_low)))



obj_low <- findmarkergene(object = obj_low,species = "Human",marker = cellmatch,
                      tissue = "Peripheral blood")

obj_low <- findcelltype(obj_low)


obj_low_1 <- findmarkergene(object = obj_low,species = "Human",marker = cellmatch,
                          tissue = "Peripheral blood",use_method = "1",comp_cluster = 1)

obj_low_1 <- findcelltype(obj_low_1)
bb<-obj_low_1@celltype

obj_low_2 <- findmarkergene(object = obj_low,species = "Human",marker = cellmatch,
                            tissue = "Peripheral blood",use_method = "2")

obj_low_2 <- findcelltype(obj_low_2)

###################总而言之，效果不好###############################

obj.sample <- SplitObject(IAV, split.by = "SampleID")
IAV.20S6575985<-obj.sample$'20S6575985'

SCT_20S6575985<-GetAssayData(object = IAV.20S6575985, assay = "SCT", slot = "scale.data")
write.table(SCT_20S6575985,file="SCT_20S6575985",sep="\t",quote=F)




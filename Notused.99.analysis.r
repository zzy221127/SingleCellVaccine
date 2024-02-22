##/etc/rstudio/rserver.conf


library("dplyr")
library("Seurat")

control <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/control_v1.rds")
IAV <- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/vaccine.rds")
mRNA<- readRDS(file = "/data/01.singlecell_vaccine/05.rds.files/mRNAP1_all.combined1.rds")



####常规marker

#################################Series Markers###################################
################################循环外面定义combined matrix########################

mRNA.sample.id<-unique(mRNA@meta.data$orig.ident)
control.sample.id<-unique(control@meta.data$sample)
IAV.sample.id<-unique(IAV@meta.data$donor)

##CD4(IFN-γ,TNF,IL-2,IL-4,IL-5),CD8(IFN-γ,TNF,IL-2)
###th2:(IL-4,IL-5,IL-13),th1:(IFN-γ)


###########到底是用表达量还是用表达细胞，附录里证明，表达量和表达细胞的关联性是很强的。都可以用。#########
#######################记得每次loop前先清空dataframe##########################

set1.markers.combined<--data.frame()
set1.markers<-c('CD4+ T cells IFNγ','CD4+ T cells TNF',
                'CD4+ T cells IL2','CD4+ T cells IL4','CD4+ T cells IL5',
                'CD8+ T cells IFNγ ','CD8+ T cells TNF','CD8+ T cells IL2')



CD4.AIM.markers.combined1<-data.frame()
CD4.AIM.markers.combined2<-data.frame()
CD4.AIM.markers.combined3<-data.frame()
CD4.AIM.markers<-c('CD4+ T cells','Non-naive CD4+ T cell',
                'CD4+ Tcm cell','CD4+ Tem cell','CD4+ Temra cell',
                'Memory CD4+ T cell','Th1 cell')


CD8.AIM.markers.combined<-data.frame()
CD8.AIM.markers<-c('CD8+ T cells','Non-naive CD8+ T cell',
                   'CD8+ Tcm cell','CD8+ Tem cell','CD8+ Temra cell',
                   'Memory CD8+ T cell')



###########################################################################################


########################开始循环################################################

#for (i in c(1:length(mRNA.sample.id)) ) {
  
#  ID<-mRNA.sample.id[i]
#  name<-subset(x=mRNA,subset = orig.ident ==ID)   ####按样品名subset
 
################################################################################  
   
  for (i in c(1:length(control.sample.id)) ) {
    
    ID<-control.sample.id[i]
    name<-subset(x=control,subset = sample ==ID)   ####按样品名subset  
  
################################################################################

# for (i in c(1:length(IAV.sample.id)) ) {

#    ID<-IAV.sample.id[i]
#    name<-subset(x=IAV,subset = donor ==ID)   ####按样品名subset  

    
  ###################set1##########################  
  CD19.ID.c<-grep('^CD19$', rownames(name@assays$RNA), value = FALSE)  
  CD14.ID.c<-grep('^CD14$', rownames(name@assays$RNA), value = FALSE)  
  CD3D.ID.c<-grep('^CD3D$', rownames(name@assays$RNA), value = FALSE)  
  CD4.ID.c<-grep('^CD4$', rownames(name@assays$RNA), value = FALSE)  
  CD8A.ID.c<-grep('^CD8A$', rownames(name@assays$RNA), value = FALSE)
  IFNG.ID.c<-grep('^IFNG$', rownames(name@assays$RNA), value = FALSE)  
  TNF.ID.c<-grep('^TNF$', rownames(name@assays$RNA), value = FALSE)  
  IL2.ID.c<-grep('^IL2$', rownames(name@assays$RNA), value = FALSE)  
  IL4.ID.c<-grep('^IL4$', rownames(name@assays$RNA), value = FALSE) 
  IL5.ID.c<-grep('^IL5$', rownames(name@assays$RNA), value = FALSE) 
  IL13.ID.c<-grep('^IL13$', rownames(name@assays$RNA), value = FALSE) 
  CCR7.ID.c<-grep('^CCR7$', rownames(name@assays$RNA), value = FALSE)  
  
    
  ####################set2########################
  PTPRC.ID.c<-grep('^PTPRC$', rownames(name@assays$RNA), value = FALSE)  ###CD45RA
  CCR7.ID.c<-grep('^CCR7$', rownames(name@assays$RNA), value = FALSE) 
  CD27.ID.c<-grep('^CD27$', rownames(name@assays$RNA), value = FALSE)
  
  CXCR5.ID.c<-grep('^CXCR5$', rownames(name@assays$RNA), value = FALSE)  
  CXCR3.ID.c<-grep('^CXCR3$', rownames(name@assays$RNA), value = FALSE)  
  CCR6.ID.c<-grep('^CCR6$', rownames(name@assays$RNA), value = FALSE)
  
  TNFRSF4.ID.c<-grep('^TNFRSF4$', rownames(name@assays$RNA), value = FALSE) ###CD134,OX40
  TNFRSF9.ID.c<-grep('^TNFRSF9$', rownames(name@assays$RNA), value = FALSE) ###CD137
  CD69.ID.c<-grep('^CD69$', rownames(name@assays$RNA), value = FALSE)
 # OX40.ID.c<-grep('^OX40$', rownames(name@assays$RNA), value = FALSE) ###CD134
 
   ####################set3########################  
  #TBX21.ID.c<-grep('^TBX21$', rownames(name@assays$RNA), value = FALSE)  
  
  ###################set4##########################
  #LAMP1.ID.c<-grep('^LAMP1$', rownames(name@assays$RNA), value = FALSE) 
  
  ###################set5##########################
 

  
############ count how many cells express the gene (all non-zero expressing cells)
  
  markerall.c<-length(colnames(name@assays$RNA@counts)) ##all
  
#####CD4+Tcell#################其实可以把CD4和CD8的分别subset出来，代码就不会这么冗余。就是把name改掉(有些问题，晚点)####
  
  marker05.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0
  )) ##CD19- & CD14- CD3D+ CD4+ CD8A-


  
############CD4+ AIM############################
##########CD4+ T cell#########  
  aim.cd4.t.c1<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                            name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.cd4.t.c2<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.cd4.t.c3<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  
  
######Non-naive CD4+ T cell############ 
  
  aim.nonnaive.cd4t.c1.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                              name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c1.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c1.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c1<-aim.nonnaive.cd4t.c1.a+aim.nonnaive.cd4t.c1.b+aim.nonnaive.cd4t.c1.c
  
  
   
  aim.nonnaive.cd4t.c2.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
 
  aim.nonnaive.cd4t.c2.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  
  aim.nonnaive.cd4t.c2.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c2<-aim.nonnaive.cd4t.c2.a+aim.nonnaive.cd4t.c2.b+aim.nonnaive.cd4t.c2.c
  
  
   
  aim.nonnaive.cd4t.c3.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                              name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  

  aim.nonnaive.cd4t.c3.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                      name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c3.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                      name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.nonnaive.cd4t.c3<-aim.nonnaive.cd4t.c3.a+aim.nonnaive.cd4t.c3.b+aim.nonnaive.cd4t.c3.c
  
 
  
######CD4+ Tcm cell############  
  

  aim.Tcm.cd4t.c1<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                             name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                             name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.Tcm.cd4t.c2<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                             name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                             name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.Tcm.cd4t.c3<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                             name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                             name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  
  
 # aim.Tcm.t.c<- aim.Tcm.t.c1+ aim.Tcm.t.c2+ aim.Tcm.t.c3
  
  
  
######CD4+ Tem cell############  
  
  
  aim.Tem.cd4t.c1<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                               name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.Tem.cd4t.c2<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                               name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.Tem.cd4t.c3<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                               name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  
  
  #aim.Tem.t.c<- aim.Tem.t.c1+ aim.Tem.t.c2+ aim.Tem.t.c3
  
  
######CD4+ Term cell############  
  
  
  aim.Term.cd4t.c1<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] != 0 &
                                name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.Term.cd4t.c2<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] != 0 &
                                name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.Term.cd4t.c3<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                name@assays$RNA@counts[CD8A.ID.c,] == 0 & name@assays$RNA@counts[PTPRC.ID.c,] != 0 &
                                name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  #aim.Term.t.c<- aim.Term.t.c1+ aim.Term.t.c2+ aim.Term.t.c3
  
  
######Memory CD4+ T cell############ 
  
  aim.memory.cd4t.c1.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd4t.c1.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd4t.c1.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd4t.c1<-aim.memory.cd4t.c1.a+aim.memory.cd4t.c1.b+aim.memory.cd4t.c1.c
  
  
  
  aim.memory.cd4t.c2.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.memory.cd4t.c2.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.memory.cd4t.c2.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.memory.cd4t.c2<-aim.memory.cd4t.c2.a+aim.memory.cd4t.c2.b+aim.memory.cd4t.c2.c
  
  
  aim.memory.cd4t.c3.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  
  aim.memory.cd4t.c3.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                       name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.memory.cd4t.c3.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                       name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                                       name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                                       name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                                       name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0)) 
  
  aim.memory.cd4t.c3<-aim.memory.cd4t.c3.a+aim.memory.cd4t.c3.b+aim.memory.cd4t.c3.c


######Th1 cell############ 
  
  aim.th1.c1.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c1.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  
  aim.th1.c1.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c1<-aim.th1.c1.a+aim.th1.c1.b+aim.th1.c1.c
  
  
  aim.th1.c2.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c2.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c2.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                               name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  
  aim.th1.c2<-aim.th1.c2.a+aim.th1.c2.b+aim.th1.c2.c
  
  
  aim.th1.c3.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  
  aim.th1.c3.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c3.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                               name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] != 0 &
                               name@assays$RNA@counts[CD8A.ID.c,] == 0 & 
                               name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                               name@assays$RNA@counts[TNFRSF9.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF4.ID.c,] != 0 &
                               name@assays$RNA@counts[CXCR5.ID.c,] == 0 & name@assays$RNA@counts[CXCR3.ID.c,] != 0 &
                               name@assays$RNA@counts[CCR6.ID.c,] == 0) )
  
  aim.th1.c3<-aim.th1.c3.a+aim.th1.c3.b+aim.th1.c3.c
  
  
  
#########CD4+ (IL2 IL4 IL5  IFNG+ TNF+)###################################
  
  marker06.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 &
                             name@assays$RNA@counts[IFNG.ID.c,] != 0
  )) ##CD19- CD14- CD3D+ CD4+ CD8A- IFNG+
  
  marker07.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 &
                             name@assays$RNA@counts[TNF.ID.c,] != 0
  )) ##CD19- & CD14- CD3D+ CD4+ CD8A- TNF+
  
  marker08.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 &
                             name@assays$RNA@counts[IL2.ID.c,] != 0
  )) ##CD19- & CD14- CD3D+ CD4+ CD8A- IL2+


  marker09.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 &
                             name@assays$RNA@counts[IL4.ID.c,] != 0
  )) ##CD19- & CD14- CD3D+ CD4+ CD8A- IL4+
  
  marker10.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] != 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] == 0 &
                             name@assays$RNA@counts[IL5.ID.c,] != 0
  )) ##CD19- & CD14- CD3D+ CD4+ CD8A- IL5+


########################CD8#############################
  
  marker25.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] == 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] != 0
  )) ##CD19- CD14- CD3D+ CD4- CD8A+
  
  
############CD8+ AIM############################
##########CD8+ T cell#########  
  aim.cd8.t.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                            name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                            name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                            name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  
  
######Non-naive CD8+ T cell############ 
  
  aim.nonnaive.cd8t.c.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd8t.c.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd8t.c.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.nonnaive.cd8t.c<-aim.nonnaive.cd8t.c.a+aim.nonnaive.cd8t.c.b+aim.nonnaive.cd8t.c.c
  
  
######CD8+ Tcm cell############  
  
  aim.Tcm.cd8t.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                 name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                 name@assays$RNA@counts[CD8A.ID.c,] != 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                                 name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                 name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0))   
  
  
  ######CD8+ Tem cell############  
  
  aim.Tem.cd8t.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                 name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                 name@assays$RNA@counts[CD8A.ID.c,] != 0 & name@assays$RNA@counts[PTPRC.ID.c,] == 0 &
                                 name@assays$RNA@counts[CCR7.ID.c,] == 0 &
                                 name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  ######CD8+ Term cell############  
  
  
  aim.Term.cd8t.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                  name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                  name@assays$RNA@counts[CD8A.ID.c,] != 0 & name@assays$RNA@counts[PTPRC.ID.c,] != 0 &
                                  name@assays$RNA@counts[CCR7.ID.c,] != 0 &
                                  name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  
  
  ######Memory CD8+ T cell############ 
  
  aim.memory.cd8t.c.a<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd8t.c.b<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] != 0 & name@assays$RNA@counts[CD27.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd8t.c.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & name@assays$RNA@counts[CD14.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD3D.ID.c,] != 0 & name@assays$RNA@counts[CD4.ID.c,] == 0 &
                                      name@assays$RNA@counts[CD8A.ID.c,] != 0 & 
                                      name@assays$RNA@counts[PTPRC.ID.c,] == 0 & name@assays$RNA@counts[CD27.ID.c,] != 0 &
                                      name@assays$RNA@counts[CD69.ID.c,] != 0 & name@assays$RNA@counts[TNFRSF9.ID.c,] != 0)) 
  
  aim.memory.cd8t.c<-aim.memory.cd8t.c.a+aim.memory.cd8t.c.b+aim.memory.cd8t.c.c
  
  
    
  ###########CD8 (IFN-γ/ TNF-α/ IL-2/ CD107α+)
  
  marker26.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] == 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] != 0 &
                             name@assays$RNA@counts[IFNG.ID.c,] != 0
  )) ##CD19- CD14- CD3D+ CD4- CD8A+ IFNG+
  
  
  marker27.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] == 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] != 0 &
                             name@assays$RNA@counts[TNF.ID.c,] != 0
  )) ##CD19- CD14- CD3D+ CD4- CD8A+ TNF+
  
  marker28.c<-length(which(name@assays$RNA@counts[CD19.ID.c,] == 0  & 
                             name@assays$RNA@counts[CD14.ID.c,] == 0 &
                             name@assays$RNA@counts[CD3D.ID.c,] != 0 & 
                             name@assays$RNA@counts[CD4.ID.c,] == 0 &
                             name@assays$RNA@counts[CD8A.ID.c,] != 0 &
                             name@assays$RNA@counts[IL2.ID.c,] != 0
  )) ##CD19- CD14- CD3D+ CD4- CD8A+ IL2+
  
  
#######################
#######Sample matrix combined start##############
###set1###
  
  ##marker05.c  is CD4+ T cell 
  ### marker25.c is CD8+ T cell
  
##'CD4+ T cells IFNγ ','CD4+ T cells TNF','CD4+ T cells IL2','CD4+ T cells IL4','CD4+ T cells IL5','CD4+ th1/th2'
##'CD8+ T cells IFNγ ','CD8+ T cells TNF','CD8+ T cells IL2', 
  
  
  set1.pre<-c(marker06.c/marker05.c,marker07.c/marker05.c,marker08.c/marker05.c,marker09.c/marker05.c,
              marker10.c/marker05.c,marker26.c/marker25.c,marker27.c/marker25.c,marker28.c/marker25.c)
  
  sample.name1 <- data.frame(Marker=set1.markers,precentage=set1.pre,SampleID=rep(ID,8))   
  set1.markers.combined<-rbind(sample.name1,set1.markers.combined)

  
  
  aim.cd4.set1<-c(aim.cd4.t.c1/marker05.c, aim.nonnaive.cd4t.c1/marker05.c,aim.Tcm.cd4t.c1/marker05.c, 
                  aim.Tem.cd4t.c1/marker05.c,aim.Term.cd4t.c1/marker05.c,aim.memory.cd4t.c1/marker05.c,
                  aim.th1.c1/marker05.c)
  
  sample.name2 <- data.frame(Marker=CD4.AIM.markers,precentage=aim.cd4.set1,SampleID=rep(ID,7),
                             Set=rep("Set1",7))
  CD4.AIM.markers.combined1<-rbind(sample.name2,CD4.AIM.markers.combined1)


    
  aim.cd4.set2<-c(aim.cd4.t.c2/marker05.c, aim.nonnaive.cd4t.c2/marker05.c,aim.Tcm.cd4t.c2/marker05.c, 
                  aim.Tem.cd4t.c2/marker05.c,aim.Term.cd4t.c2/marker05.c,aim.memory.cd4t.c2/marker05.c,
                  aim.th1.c2/marker05.c)
  
  sample.name3 <- data.frame(Marker=CD4.AIM.markers,precentage=aim.cd4.set2,SampleID=rep(ID,7),
                             Set=rep("Set2",7))
  CD4.AIM.markers.combined2<-rbind(sample.name3,CD4.AIM.markers.combined2)
  
  

  aim.cd4.set3<-c(aim.cd4.t.c3/marker05.c, aim.nonnaive.cd4t.c3/marker05.c,aim.Tcm.cd4t.c3/marker05.c, 
                  aim.Tem.cd4t.c3/marker05.c,aim.Term.cd4t.c3/marker05.c,aim.memory.cd4t.c3/marker05.c,
                  aim.th1.c3/marker05.c)
  sample.name4 <- data.frame(Marker=CD4.AIM.markers,precentage=aim.cd4.set3,SampleID=rep(ID,7),
                             Set=rep("Set3",7))
  CD4.AIM.markers.combined3<-rbind(sample.name4,CD4.AIM.markers.combined3)
  
  
  
  
  aim.cd8<-c(aim.cd8.t.c/marker25.c, aim.nonnaive.cd8t.c/marker25.c,aim.Tcm.cd8t.c/marker25.c, 
                  aim.Tem.cd8t.c/marker25.c,aim.Term.cd8t.c/marker25.c,aim.memory.cd8t.c/marker25.c)
  sample.name5 <- data.frame(Marker=CD8.AIM.markers,precentage=aim.cd8,SampleID=rep(ID,6))   
  CD8.AIM.markers.combined<-rbind(sample.name5,CD8.AIM.markers.combined)
  
  
  
  }
################################Loop Over#######################################  
################################################################################




################mRNA matrix##########################
mRNA.set1.markers.combined<-set1.markers.combined
mRNA.CD4.AIM.markers.combined<-rbind(CD4.AIM.markers.combined1,CD4.AIM.markers.combined2,CD4.AIM.markers.combined3) ##明天跑这个
mRNA.CD8.AIM.markers.combined<-CD8.AIM.markers.combined


tmp_c<-strsplit(mRNA.set1.markers.combined$SampleID,"_")
mRNA.set1.markers.combined$Sample<-sapply(tmp_c, function(x) x[[1]])
mRNA.set1.markers.combined$Day<-sapply(tmp_c, function(x) x[[2]])
tmp_b<-strsplit(mRNA.set1.markers.combined$Marker," ")
mRNA.set1.markers.combined$Celltype<-paste(sapply(tmp_b, function(x) x[[1]]),
                                           sapply(tmp_b, function(x) x[[2]]),
                                           sapply(tmp_b, function(x) x[[3]]),sep=" ")
mRNA.set1.markers.combined$Pre<-mRNA.set1.markers.combined$precentage*100



tmp_c<-strsplit(mRNA.CD4.AIM.markers.combined$SampleID,"_")
mRNA.CD4.AIM.markers.combined$Sample<-sapply(tmp_c, function(x) x[[1]])
mRNA.CD4.AIM.markers.combined$Day<-sapply(tmp_c, function(x) x[[2]])
mRNA.CD4.AIM.markers.combined$Celltype<-mRNA.CD4.AIM.markers.combined$Marker
mRNA.CD4.AIM.markers.combined$Pre<-mRNA.CD4.AIM.markers.combined$precentage*100


tmp_c<-strsplit(mRNA.CD8.AIM.markers.combined$SampleID,"_")
mRNA.CD8.AIM.markers.combined$Sample<-sapply(tmp_c, function(x) x[[1]])
mRNA.CD8.AIM.markers.combined$Day<-sapply(tmp_c, function(x) x[[2]])
mRNA.CD8.AIM.markers.combined$Celltype<-mRNA.CD8.AIM.markers.combined$Marker
mRNA.CD8.AIM.markers.combined$Pre<-mRNA.CD8.AIM.markers.combined$precentage*100





###############control matrix###########
control.set1.markers.combined<-set1.markers.combined
control.CD4.AIM.markers.combined<-rbind(CD4.AIM.markers.combined1,CD4.AIM.markers.combined2,CD4.AIM.markers.combined3)
control.CD8.AIM.markers.combined<-CD8.AIM.markers.combined


control.set1.markers.combined$Sample<-control.set1.markers.combined$SampleID
control.set1.markers.combined$Day<-'Unvaccinated'
control.set1.markers.combined$Pre<-control.set1.markers.combined$precentage*100
tmp_c<-strsplit(control.set1.markers.combined$Marker," ")
control.set1.markers.combined$Celltype<-paste(sapply(tmp_c, function(x) x[[1]]),
                                           sapply(tmp_c, function(x) x[[2]]),
                                           sapply(tmp_c, function(x) x[[3]]),sep=" ")


control.CD4.AIM.markers.combined$Sample<-control.CD4.AIM.markers.combined$SampleID
control.CD4.AIM.markers.combined$Day<-'Unvaccinated'
control.CD4.AIM.markers.combined$Pre<-control.CD4.AIM.markers.combined$precentage*100
control.CD4.AIM.markers.combined$Celltype<-control.CD4.AIM.markers.combined$Marker

control.CD8.AIM.markers.combined$Sample<-control.CD8.AIM.markers.combined$SampleID
control.CD8.AIM.markers.combined$Day<-'Unvaccinated'
control.CD8.AIM.markers.combined$Pre<-control.CD8.AIM.markers.combined$precentage*100
control.CD8.AIM.markers.combined$Celltype<-control.CD8.AIM.markers.combined$Marker



#################IAV matrix################
IAV.set1.markers.combined<-set1.markers.combined
IAV.CD4.AIM.markers.combined<-rbind(CD4.AIM.markers.combined1,CD4.AIM.markers.combined2,CD4.AIM.markers.combined3)
IAV.CD8.AIM.markers.combined<-CD8.AIM.markers.combined


IAV.set1.markers.combined$Sample<-IAV.set1.markers.combined$SampleID
IAV.set1.markers.combined$Day<-'IAV vaccine'
IAV.set1.markers.combined$Pre<-IAV.set1.markers.combined$precentage*100

tmp_c<-strsplit(IAV.set1.markers.combined$Marker," ")
IAV.set1.markers.combined$Celltype<-paste(sapply(tmp_c, function(x) x[[1]]),
                                              sapply(tmp_c, function(x) x[[2]]),
                                              sapply(tmp_c, function(x) x[[3]]),sep=" ")


IAV.CD4.AIM.markers.combined$Sample<-IAV.CD4.AIM.markers.combined$SampleID
IAV.CD4.AIM.markers.combined$Day<-'IAV vaccine'
IAV.CD4.AIM.markers.combined$Pre<-IAV.CD4.AIM.markers.combined$precentage*100
IAV.CD4.AIM.markers.combined$Celltype<-IAV.CD4.AIM.markers.combined$Marker


IAV.CD8.AIM.markers.combined$Sample<-IAV.CD8.AIM.markers.combined$SampleID
IAV.CD8.AIM.markers.combined$Day<-'IAV vaccine'
IAV.CD8.AIM.markers.combined$Pre<-IAV.CD8.AIM.markers.combined$precentage*100
IAV.CD8.AIM.markers.combined$Celltype<-IAV.CD8.AIM.markers.combined$Marker


#####################all matrix###############################
library(dplyr)

all.set1.markers.combined<-rbind(control.set1.markers.combined,IAV.set1.markers.combined,mRNA.set1.markers.combined)
all.set1.markers.combined<-all.set1.markers.combined %>% distinct()

all.CD4.AIM.markers.combined<-rbind(control.CD4.AIM.markers.combined,
                                    IAV.CD4.AIM.markers.combined,mRNA.CD4.AIM.markers.combined)



all.CD8.AIM.markers.combined<-rbind(control.CD8.AIM.markers.combined,
                                    IAV.CD8.AIM.markers.combined,mRNA.CD8.AIM.markers.combined)


all.CD8.AIM.markers.combined$Set<-"Set4.cd8.AIM"

all.AIM.cd4.cd8<-rbind(all.CD8.AIM.markers.combined,all.CD4.AIM.markers.combined)

#################ALL plot ############



all.set1.markers.combined$Day1<-factor(all.set1.markers.combined$Day,
                                   levels=c("Unvaccinated","IAV vaccine",
                                            "Day0","Day1","Day2","Day7","Day21","Day22","Day28","Day42"))


all.CD4.AIM.markers.combined$Day1<-factor(all.CD4.AIM.markers.combined$Day,
                                       levels=c("Unvaccinated","IAV vaccine",
                                                "Day0","Day1","Day2","Day7","Day21","Day22","Day28","Day42"))


all.CD8.AIM.markers.combined$Day1<-factor(all.CD8.AIM.markers.combined$Day,
                                          levels=c("Unvaccinated","IAV vaccine",
                                                   "Day0","Day1","Day2","Day7","Day21","Day22","Day28","Day42"))


all.AIM.cd4.cd8$Day1<-factor(all.AIM.cd4.cd8$Day,levels=c("Unvaccinated","IAV vaccine",
                                                   "Day0","Day1","Day2","Day7","Day21","Day22","Day28","Day42"))



library(ggplot2)


ggplot(data = all.set1.markers.combined) + geom_boxplot(aes(x = Day1, y = Pre,color = factor(Marker)))+
  geom_jitter(aes(x = Day1, y = Pre, color = factor(Marker)), position = position_jitterdodge())+
  facet_wrap(Celltype~Marker,scales = "free",ncol=2)+
  theme(axis.text.x = element_text(angle = 35,vjust = 0.5,hjust = 0.5))+ xlab(NULL)+ylab(NULL)


ggplot(data = all.CD4.AIM.markers.combined) + geom_boxplot(aes(x = Day1, y = Pre,color = factor(Marker)))+
  geom_jitter(aes(x = Day1, y = Pre, color = factor(Marker)), position = position_jitterdodge())+
  facet_wrap(Set~Marker,scales = "free",nrow=3)+
  theme(axis.text.x = element_text(angle = 35,vjust = 0.5,hjust = 0.5))+ xlab(NULL)+ylab(NULL)



ggplot(data = all.CD8.AIM.markers.combined) + geom_boxplot(aes(x = Day1, y = Pre,color = factor(Marker)))+
  geom_jitter(aes(x = Day1, y = Pre, color = factor(Marker)), position = position_jitterdodge())+
  facet_wrap(.~Marker,scales = "free",nrow=3)+
  theme(axis.text.x = element_text(angle = 35,vjust = 0.5,hjust = 0.5))+ xlab(NULL)+ylab(NULL)



ggplot(data = all.AIM.cd4.cd8) + geom_boxplot(aes(x = Day1, y = Pre,color = factor(Marker)))+
  geom_jitter(aes(x = Day1, y = Pre, color = factor(Marker)), position = position_jitterdodge())+
  facet_wrap(Set~Marker,scales = "free",nrow=4)+
  theme(axis.text.x = element_text(angle = 35,vjust = 0.5,hjust = 0.5))+ xlab(NULL)+ylab(NULL)


#################################AIM告一段落###########################











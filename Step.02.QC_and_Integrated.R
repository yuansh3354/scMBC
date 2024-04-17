setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]
########################## Step.02 QC & Integrated All Samples ########################## 
if(F){
  # sce.merge.obj = JoinLayers(sce.merge.obj)
  # sce = JoinLayers(sce)
  
  sce = sce.merge.obj
  sce[["pMT"]] = PercentageFeatureSet(object = sce, pattern = "^MT-")
  # sce[["pRP"]] = PercentageFeatureSet(sce, pattern = "^RP[SL][[:digit:]]")
  summary(sce[["pMT"]])
  # summary(sce[["pRP"]])
  summary(sce[["nFeature_RNA"]])
  summary(sce[["nCount_RNA"]])
  nFeature_lower = 300
  nFeature_upper = 10000
  nCount_lower = 1000
  nCount_upper = 100000
  pMT_lower = 0
  pMT_upper = 25
  pRP_lower = 0
  
  sce = subset(sce, subset =  nFeature_RNA > nFeature_lower & 
                 nCount_RNA < nCount_upper & nCount_RNA > nCount_lower &
                 pMT < pMT_upper)
  meta = read.csv('MetaInfo.csv')
  sce@meta.data$CellID = rownames(sce@meta.data)
  combined_data <- merge(sce@meta.data, meta, by = "SampleID", all.x = T)
  rownames(combined_data) = combined_data$CellID
  sce@meta.data = combined_data
  
  table(sce$SampleID)
  table(sce$Type)
  sce.merge.obj = sce
}
if(F){
  # Run 1 
  if(F){
    sce = sce.merge.obj
    sce = NormalizeData(sce,verbose = F)
    sce = FindVariableFeatures(sce,verbose = F)
    sce = ScaleData(sce,verbose = F)
    sce = Seurat::RunPCA(sce,verbose = F)
    # sce = FindNeighbors(sce,  dims = 1:30,verbose = F)
    # sce = FindClusters(sce,verbose = F)
    # sce = RunUMAP(sce, dims = 1:30,verbose = F)
    options(future.globals.maxSize = 50 * 1024^3)
    sce <- IntegrateLayers(
      object = sce, method = RPCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.rpca",
      verbose = FALSE
    )
    sce <- FindNeighbors(sce, reduction = "integrated.rpca", dims = 1:50)
    sce <- FindClusters(sce, cluster.name = "rpca.cluster")
  }
  
  ids = names(which(table(sce$rpca.cluster)<500))
  ids = rownames(sce@meta.data[which(sce$rpca.cluster %in% ids),])
  ids = setdiff(Cells(sce),ids)
  sce.merge.obj = subset(sce, cells = ids)
  
  # Run 2 
  if(F){
    sce = sce.merge.obj
    sce = NormalizeData(sce,verbose = F)
    sce = FindVariableFeatures(sce,verbose = F)
    sce = ScaleData(sce,verbose = F)
    sce = Seurat::RunPCA(sce,verbose = F)
    # sce = FindNeighbors(sce,  dims = 1:30,verbose = F)
    # sce = FindClusters(sce,verbose = F)
    # sce = RunUMAP(sce, dims = 1:30,verbose = F)
    options(future.globals.maxSize = 50 * 1024^3)
    sce <- IntegrateLayers(
      object = sce, method = RPCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.rpca",
      verbose = FALSE
    )
    sce <- FindNeighbors(sce, reduction = "integrated.rpca", dims = 1:50)
    sce <- FindClusters(sce, cluster.name = "rpca.cluster")
  }
  sce <- RunTSNE(sce, reduction = "integrated.rpca", dims = 1:50, reduction.name = "tsne.rpca")
}
sce = JoinLayers(sce)
saveRDS(sce,'result/RDS/02.integrated.rpca.RDS')
########################## 
# INFO[SAVE FILENAME]: result/RDS/02.integrated.rpca.RDS
# INFO[END TIME]: 2024-04-06 11:07:31


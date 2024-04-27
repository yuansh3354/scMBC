########################## Step.02 QC & Integrated All Samples ########################## 
if(T){
  sce[["pMT"]] = PercentageFeatureSet(object = sce, pattern = "^MT-")
  summary(sce[["pMT"]])
  summary(sce[["nFeature_RNA"]])
  summary(sce[["nCount_RNA"]])
  nFeature_lower = 300
  nFeature_upper = 10000
  nCount_lower = 1000
  nCount_upper = 100000
  pMT_lower = 0
  pMT_upper = 25
  
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
saveRDS(sce,'result/RDS/02.integrated.rpca.RDS')


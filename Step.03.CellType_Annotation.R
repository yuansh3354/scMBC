setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]
########################## Step.03 Manually Annotation ########################## 
########################## Step.03.1 For All Cells ########################## 
sce = readRDS('result/RDS/02.integrated.rpca.RDS')

# Check Markers
CellMarkers <- list(
  `Mast cell` = c('KIT','ENPP3','CPA3'),
  `B cell` = c('MS4A1', 'CD79A','CD79B'), # 
  `Plasma cell` = c('MZB1', 'IGKC', 'JCHAIN'), # 
  `Neutrophils` = c('CSF3R','SORL1','S100A8'),
  `Fibroblast` = c('COL1A1','COL1A2','DCN'), # 
  `Endothelium` = c('PTPRB', 'PECAM1','RAMP2'), # 
  `Macrophage` = c('CD68','CD163','CD14'),
  `T cell` = c('CD3D', 'CD3E','CD3G'),
  `Epithelium` = c('EPCAM','KRT8','KRT18')# 
)
DotPlot(sce, features = CellMarkers,cluster.idents=T) + theme_bw() + RotatedAxis()

# Manually Annotate Cell Types
sce$celltype = NA
sce$celltype = ifelse(Idents(sce) %in% c(23),'Mast cell',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(21,5),'B cell',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(26,18,22),'Plasma cell',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(24),'Neutrophils',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(6,19,11),'Fibroblast',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(7,25),'Endothelium',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(9,20),'Macrophage',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(4,2,15,13),'T cell',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(14,10,16,8,0,1,17,12,3),'Epithelium',sce$celltype)


DimPlot(sce, group.by = 'celltype',label = T,label.size = 5,raster = F) + 
  scale_color_manual(values = cell.colors) + NoLegend()
ggsave('~/Desktop/Dim.pdf',height = 8,width = 10)

if(F){
  plot = DimPlot(sce[,which(sce$celltype %in% c('Epithelium') )], group.by = 'celltype',label = T,label.size = 5,pt.size = 1.3) + 
    scale_color_manual(values = colors) + NoLegend() + xlim(-25,25) + ylim(-25,25)
  cells = CellSelector(plot)
  
  ids = setdiff(Cells(sce),cells)
  sce = subset(sce, cells=ids)
  
  sce$celltype = ifelse(sce$CellID %in% cells,'B cell',sce$celltype)
  sce$celltype = ifelse(sce$CellID %in% cells,'Macrophage',sce$celltype)
  sce$celltype = ifelse(sce$CellID %in% cells,'T cell',sce$celltype)
  sce$celltype = ifelse(sce$CellID %in% cells,'Endothelium',sce$celltype)
  sce$celltype = ifelse(sce$CellID %in% cells,'Epithelium',sce$celltype)
  sce$celltype = ifelse(sce$CellID %in% cells,'Fibroblast',sce$celltype)
  FeaturePlot(sce,features = c('EPCAM','CD24','KRT8','KRT18'),min.cutoff = 0.5,reduction = 'tsne.rpca')
  DimPlot(sce,group.by=c('celltype'),reduction = 'tsne.rpca',label = T,label.size = 5) + scale_color_manual(values = colors)
  saveRDS(sce,'result/RDS/03.Seurat.annotation.RDS')
}

########################## Step.03.2 For CAFs Subtype ########################## 
# Get Cells
sce  = readRDS('result/RDS/03.Seurat.annotation.RDS')
Idents(sce) = sce$celltype
fibroblast = subset(x = sce, subset = celltype == "Fibroblast")
fibroblast = CreateSeuratObject(counts = fibroblast@assays$RNA$counts)
fibroblast = AddMetaData(fibroblast, metadata = sce@meta.data)

# Combine samples with low cell counts 
fibroblast$Groups = ifelse(fibroblast$PatientID %in% names(table(fibroblast$PatientID))[which(table(fibroblast$PatientID)<100)],
                           'G0',fibroblast$PatientID)
fibroblast <- split(fibroblast, f = fibroblast$Groups)

# Run Seurat Protocols
if(T){
  fibroblast <- NormalizeData(fibroblast,verbose = F)
  fibroblast <- FindVariableFeatures(fibroblast,verbose = F)
  fibroblast <- ScaleData(fibroblast,verbose = F)
  fibroblast = Seurat::RunPCA(fibroblast,verbose = F)
  options(future.globals.maxSize = 50 * 1024^3)
  # 整合
  fibroblast <- IntegrateLayers(
    fibroblast, k.weight=30,
    method = RPCAIntegration,
    orig.reduction = "pca", 
    new.reduction = "integrated.rpca",
    verbose = FALSE
  )
  fibroblast = JoinLayers(fibroblast)
  fibroblast <- FindNeighbors(fibroblast, 
                              reduction = "integrated.rpca",
                              dims = 1:30)
  fibroblast <- FindClusters(fibroblast,resolution=1, cluster.name = "rpca.cluster")
  fibroblast <- RunTSNE(fibroblast, reduction = "integrated.rpca", 
                        dims = 1:15, 
                        reduction.name = "tsne.rpca")
}

# Check Markers
DimPlot(fibroblast,raster = F,label = T,reduction = 'tsne.rpca',pt.size = 1) + NoLegend()
CellMarkers = list(
  'EMT-like CAFs' = c('KRT8','KRT18','KRT19'),
  'Antigene CAFs' = c('CD74','HLA-DRA','HLA-DRB1'),
  'Inflammatory CAFs' = c('CXCL12','CXCL14','IGF1','IGFBP6'),
  'Vascular CAFs' = c('MCAM','MYH11','ACTA2')
)
DotPlot(fibroblast, features = CellMarkers,cluster.idents = T) + RotatedAxis()

# Define Cell Subtypes
fibroblast$Subtype = NA
fibroblast$Subtype = ifelse(Idents(fibroblast) %in% c(10,19,6,17,2,7,8,13,22),
                            'Vascular CAFs',fibroblast$Subtype)
fibroblast$Subtype = ifelse(Idents(fibroblast) %in% c(9,18,0,4,3,15,1,11,5),
                            'Inflammatory CAFs',fibroblast$Subtype)
fibroblast$Subtype = ifelse(Idents(fibroblast) %in% c(20,21,16),
                            'Antigene CAFs',fibroblast$Subtype)
fibroblast$Subtype = ifelse(Idents(fibroblast) %in% c(14,12),
                            'EMT-like CAFs',fibroblast$Subtype)

DimPlot(fibroblast,raster = F,label = T,group.by = 'Subtype',
        reduction = 'tsne.rpca',pt.size = 1)
saveRDS(fibroblast,'result/RDS/03.1.fibroblast.RDS')

########################## Step.03.3 For T Cells Subtype ########################## 
# Get Cells
sce  = readRDS('result/RDS/03.Seurat.annotation.RDS')
Idents(sce) = sce$celltype
T.cell= subset(x = sce, subset = celltype == "T cell")
T.cell = CreateSeuratObject(counts = T.cell@assays$RNA$counts,
                            min.cells = 100)
T.cell = AddMetaData(T.cell, metadata = sce@meta.data)
T.cell <- split(T.cell, f = T.cell$PatientID)

# Run Seurat Protocols
if(T){
  T.cell <- NormalizeData(T.cell,verbose = F)
  T.cell <- FindVariableFeatures(T.cell,verbose = F)
  T.cell <- ScaleData(T.cell,verbose = F)
  T.cell = Seurat::RunPCA(T.cell,verbose = F)
  options(future.globals.maxSize = 50 * 1024^3)
  T.cell <- IntegrateLayers(
    object = T.cell, method = RPCAIntegration,k.weight = 30,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE,
  )
  T.cell = JoinLayers(T.cell)
  T.cell <- FindNeighbors(T.cell, 
                          reduction = "integrated.rpca",
                          dims = 1:50)
  T.cell <- FindClusters(T.cell,resolution=1, cluster.name = "rpca.cluster")
  T.cell <- RunTSNE(T.cell, reduction = "integrated.rpca", 
                    dims = 1:50, 
                    reduction.name = "tsne.rpca")
}

# Check Markers
CellMarkers = list(
  'Treg' = c('IL2RA','FOXP3','IKZF2','FOXP3'),
  'Proliferating' =c('MKI67','TOP2A','STMN1'),
  'NK/NKT' = c('GNLY','NKG7','CCL5','TRDC','FCER1G','KLRF1'),
  'Naive' = c('TCF7','LEF1','CCR7','SELL','IL7R'),
  'Inhibitory' = c('IKZF2','TIGIT','CTLA4','LGA3','PDCD1','HAVCR2'),
  'Exhausted' = c('TOX','HAVCR2','PDCD1'),
  'Cytotoxic' = c('S100A10','GPX4','IL2','GZMA','PRF1','GZMB','GZMK','IFNG')
)
if(T){
  options(future.globals.maxSize = 50 * 1024^3)
  sce$Subtype = NA
  sce = AddModuleScore(sce,features = CellMarkers)
  colnames(sce@meta.data)
  df = DotPlot(sce, features = paste0('Cluster',seq(1,length(names(CellMarkers)))),col.min=0.2,
               cols = c("white", "blue")) + theme_cleveland() + RotatedAxis()
  df = df$data
  df$score = df$avg.exp.scaled * df$pct.exp
  
  df = spread(df[,c('id','features.plot','score')], key = "features.plot",
              value = "score")
  rownames(df) = df$id
  df = df[,-1]
  df = df - apply(df,1,max)
  for(i in 1:dim(df)[1]){
    ids = df[i,,drop=F]
    cluster.id = as.numeric(gsub('Cluster','',names(ids[which(ids>=0)])))
    Subtype = names(CellMarkers)[cluster.id]
    sce$Subtype = ifelse(sce$seurat_clusters == as.character(rownames(ids)),
                         Subtype, sce$Subtype)
  }
  sce@meta.data = sce@meta.data[,-c(grep('Cluster',colnames(sce@meta.data)))]
}

CellMarkers = list(
  'CD8' = c('CD8A','CD8B'),
  'CD4' = c('CD4')
)
if(T){
  options(future.globals.maxSize = 50 * 1024^3)
  sce = AddModuleScore(sce,features = CellMarkers)
  colnames(sce@meta.data)
  df = DotPlot(sce, features = paste0('Cluster',seq(1,length(names(CellMarkers)))),col.min=0.2,
               cols = c("white", "blue")) + theme_cleveland() + RotatedAxis()
  df = df$data
  df$score = df$avg.exp.scaled * df$pct.exp
  
  df = spread(df[,c('id','features.plot','score')], key = "features.plot",
              value = "score")
  rownames(df) = df$id
  df = df[,-1]
  df = df - apply(df,1,max)
  for(i in 1:dim(df)[1]){
    ids = df[i,,drop=F]
    cluster.id = as.numeric(gsub('Cluster','',names(ids[which(ids>=0)])))
    Subtype = names(CellMarkers)[cluster.id]
    sce$Subtype = ifelse(sce$seurat_clusters == as.character(rownames(ids)),
                         paste0(Subtype,'_',sce$Subtype),sce$Subtype )
  }
  sce@meta.data = sce@meta.data[,-c(grep('Cluster',colnames(sce@meta.data)))]
}

# Check Result
DimPlot(sce,raster = F,label = T,pt.size = 0.75,
        group.by ='Subtype',reduction = 'tsne.rpca') + 
  theme(legend.position = 'bottom')  + NoLegend()+
  theme(
    axis.title = element_blank(),  #轴标题
    axis.text = element_blank(), # 文本
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('') + scale_color_simpsons()

CellMarkers = list(
  'CD4/8' = c('CD4','CD8A','CD8B'),
  'Proliferating' =c('MKI67','TOP2A','STMN1'),
  'NK/NKT' = c('GNLY','NKG7','CCL5','TRDC','FCER1G','KLRF1'),
  'Naive' = c('TCF7','LEF1','CCR7','SELL','IL7R'),
  'Exhausted' = c('TOX','HAVCR2','PDCD1'),
  'Cytotoxic' = c('S100A10','GPX4','IL2','GZMA','PRF1','GZMB','GZMK','IFNG'),
  'Treg' = c('IL2RA','FOXP3','IKZF2','FOXP3'),
  'Inhibitory' = c('IKZF2','TIGIT','CTLA4','PDCD1','HAVCR2')
)
x = jjDotPlot(object = sce,
              gene = unlist(CellMarkers),
              xtree = F,
              ytree = F,
              dot.col = c("#A6CEE3",'white',"#c10534"),
              id = 'Subtype',
              legend.position='bottom') + ggtitle('Marker Genes')
x
topptx(x,'plts/DotPlot.pptx')
saveRDS(sce,'result/RDS/03.1.Tcell.RDS')
########################## Step.03.4 Macrophage ########################## 
# Get Cells
sce  = readRDS('result/RDS/03.Seurat.annotation.RDS')
Idents(sce) = sce$celltype
Macrophage = subset(x = sce, subset = celltype == "Macrophage")
Macrophage = CreateSeuratObject(counts = Macrophage@assays$RNA$counts,
                                min.cells = 100)
Macrophage = AddMetaData(Macrophage, metadata = sce@meta.data)
Macrophage <- split(Macrophage, f = Macrophage$PatientID)

# Run Seurat Protocols
if(T){
  Macrophage <- NormalizeData(Macrophage,verbose = F)
  Macrophage <- FindVariableFeatures(Macrophage,verbose = F)
  Macrophage <- ScaleData(Macrophage,verbose = F)
  Macrophage = Seurat::RunPCA(Macrophage,verbose = F)
  options(future.globals.maxSize = 50 * 1024^3)
  Macrophage <- IntegrateLayers(
    object = Macrophage, method = RPCAIntegration,k.weight = 30,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE,
  )
  Macrophage = JoinLayers(Macrophage)
  Macrophage <- FindNeighbors(Macrophage,
                              reduction = "integrated.rpca",
                              dims = 1:30)
  Macrophage <- FindClusters(Macrophage,resolution=1, cluster.name = "rpca.cluster")
  Macrophage <- RunTSNE(Macrophage, reduction = "integrated.rpca",
                        dims = 1:15,
                        reduction.name = "tsne.rpca")
}
DimPlot(Macrophage,reduction = 'tsne.rpca',pt.size = 2)

# Check Markers
CellMarkers = c('SPP1','CXCL9','S100A9','CCL4','STMN1','MERTK')
DotPlot(Macrophage,features = CellMarkers,cluster.idents = T) + theme_cleveland()

FeaturePlot(Macrophage,features = CellMarkers)

Macrophage$Subtype = NA
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(0,1,11,12,16,17,6,2),'Mph SPP1+',Macrophage$Subtype)
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(5,14),'Mph S100A9+',Macrophage$Subtype)
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(3,10,4,9),'Mph CCL4+',Macrophage$Subtype)
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(13,20,15),'Mph STNM1+',Macrophage$Subtype)
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(19,18),'Mph CXCL9+',Macrophage$Subtype)
Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(8,7),'Mph MERTK+',Macrophage$Subtype)
# Macrophage$Subtype = ifelse(Idents(Macrophage) %in% c(4,2,9,6),'Mph FCER1G',Macrophage$Subtype)
# DotPlot(Macrophage[,which(is.na(Macrophage$Subtype))],features = CellMarkers,cluster.idents = T) + theme_cleveland()

DimPlot(Macrophage,raster = F,label = T,pt.size = 0.75,
        group.by ='Subtype',reduction = 'tsne.rpca') + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = colors) + NoLegend()+
  theme(
    axis.title = element_blank(),  #轴标题
    axis.text = element_blank(), # 文本
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
saveRDS(Macrophage,'result/RDS/03.1.Macrophage.RDS')

########################## Step.03.5 Epithelium ########################## 
# Get Cells
sce  = readRDS('result/RDS/03.Seurat.annotation.RDS')
Idents(sce) = sce$celltype
Epithelium = subset(x = sce, subset = celltype == "Epithelium")
Epithelium = CreateSeuratObject(counts = Epithelium@assays$RNA$counts,
                                min.cells = 100)
Epithelium = AddMetaData(Epithelium, metadata = sce@meta.data)
Epithelium <- split(Epithelium, f = Epithelium$SampleID)
# Run Seurat Protocols
if(T){
  Epithelium <- NormalizeData(Epithelium,verbose = F)
  Epithelium <- FindVariableFeatures(Epithelium,verbose = F)
  Epithelium <- ScaleData(Epithelium,verbose = F)
  Epithelium = Seurat::RunPCA(Epithelium,verbose = F)
  options(future.globals.maxSize = 50 * 1024^3)
  Epithelium <- IntegrateLayers(
    object = Epithelium, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE,
  )
  Epithelium = JoinLayers(Epithelium)
  Epithelium <- FindNeighbors(Epithelium,
                              reduction = "integrated.rpca",
                              dims = 1:50)
  Epithelium <- FindClusters(Epithelium,cluster.name = "rpca.cluster")
  Epithelium <- RunUMAP(Epithelium, reduction = "integrated.rpca",
                        dims = 1:50,
                        reduction.name = "tsne.rpca")
}

# Check Markers
CellMarkers = list(
  'Prolification-like' =  c(unlist(strsplit('SAT1,CD24,TNFRSF12A',','))),
  'Neuroendocrine-like' = c(unlist(strsplit('GPHN,GFRA1,YAP1',','))),
  "Metastasis-like" = c(unlist(strsplit('RNF181,MDK,BST2',','))),
  'Invasion-like' = c('MKI67','TOP2A','NDC80'),
  'Immune-like' = c(unlist(strsplit('CD74,PTPRC,PARP14',','))),
  'EMT-like' = c(unlist(strsplit('VIM,FN1,MMP2',',')))
)
DotPlot(Epithelium, features = CellMarkers,col.min=0.2,
        cols = c("white", "blue"),cluster.idents = T) + theme_cleveland() + RotatedAxis()

Epithelium$Subtype = 'u'
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(21,4,14,8,3,0),'Neuroendocrine-like',Epithelium$Subtype)
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(17,5,22),'Immune-like',Epithelium$Subtype)
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(10,19,16,20),'EMT-like',Epithelium$Subtype)
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(7,9,1,2),'Metastasis-like',Epithelium$Subtype)
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(11,6,15,18),'Prolification-like',Epithelium$Subtype)
Epithelium$Subtype = ifelse(Idents(Epithelium) %in% c(12,13),'Invasion-like',Epithelium$Subtype)
DotPlot(Epithelium, features = CellMarkers,group.by = 'Subtype') + theme_cleveland() + RotatedAxis()

DimPlot(Epithelium,group.by = 'Subtype',raster = F,label = T,label.size = 5)

jjDotPlot(object = Epithelium,
          gene = unlist(CellMarkers),
          xtree = F,
          ytree = F,
          dot.col = c("#A6CEE3",'white',"#c10534"),
          id = 'Subtype',
          rescale = T,
          rescale.min = -2,
          rescale.max = 2,
          point.geom = F,
          tile.geom = T)

saveRDS(Epithelium,'result/RDS/03.1.Epithelium.RDS')
########################## Step.03.6 Combine All ########################## 
sce  = readRDS('result/RDS/03.Seurat.annotation.RDS')
Tcell = readRDS('result/RDS/03.1.Tcell.RDS')
fibroblast = readRDS('result/RDS/03.1.fibroblast.RDS')
Macrophage = readRDS('result/RDS/03.1.Macrophage.RDS')
Epithelium = readRDS('result/RDS/03.1.Epithelium.RDS')

sce = Merge_MetaData(sce,Tcell,by = 'Subtype') %>% 
  Merge_MetaData(fibroblast,by = 'Subtype') %>% 
  Merge_MetaData(Macrophage,by = 'Subtype') %>% 
  Merge_MetaData(Epithelium,by = 'Subtype')

sce$Subtype = factor(sce$Subtype,levels = c(
  'Endothelium','Plasma cell','B cell','Mast cell',
  'EMT-like','Neuroendocrine-like','Metastasis-like','Prolification-like','Immune-like','Invasion-like',
  'Inhibitory T','NKT','Exhausted T','Naive T','Cytotoxic T',
  'Inflammatory CAFs','Vascular CAFs','Antigene CAFs','EMT-like CAFs',
  'Mph MERTK+','Mph SPP1+','Mph CCL4+','Mph STNM1+','Mph S100A9+',
  'Mph CXCL9+','Neutrophils'
))

# Check Result
df = as.data.frame.array(table(sce$Subtype,sce$celltype))
pheatmap::pheatmap(df,display_numbers=T,number_format = "%f",scale = 'row')
saveRDS(sce,'result/RDS/03.Seurat.Sub.annotation.RDS')

########################## 
# INFO[SAVE FILENAME]: result/RDS/03.Seurat.annotation.RDS
# INFO[SAVE FILENAME]: result/RDS/03.1.fibroblast.RDS
# INFO[SAVE FILENAME]: result/RDS/03.1.Tcell.RDS
# INFO[SAVE FILENAME]: result/RDS/03.1.Macrophage.RDS
# INFO[SAVE FILENAME]: result/RDS/03.1.Epithelium.RDS
# INFO[SAVE FILENAME]: result/RDS/03.Seurat.Sub.annotation.RDS

# # Step.00 ----------------------------------------------------------------------
rm(list = ls())
gc(reset = T)
source('script/Functions/Functions.R')
# unload packages
pkgs[which(load.packages == F)]
rpid = Sys.getpid()
options(future.globals.maxSize = 1000 * 1024^3)
################################################################################
sce.total = qread('Result/sce.total.qs',nthreads = 512)
sce.list = qread('Result/sce.list.total.qs',nthreads = 512)
gs = qread('Result/final_use_gs.qs',nthreads = 512)
hd.list = qread('Result/cellpose_obj.hd.sample.list.qs',nthreads = 512)
# ============================================================================

library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(parallel)
library(qs)
library(Banksy)

# ============================================================================
{
  obj.adata <- zellkonverter::readH5AD('Result/Xenium_with_annotations.h5ad', use_hdf5 = TRUE)
  obj <- CreateSeuratObject(obj.adata@assays@data$X)
  obj <- AddMetaData(obj, metadata = as.data.frame(obj.adata@colData))
  
  spatial_coords <- as.data.frame(reducedDim(obj.adata, "spatial"))
  colnames(spatial_coords) <- c("spatial_1", "spatial_2")
  obj[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(spatial_coords),
    key = "spatial_",
    assay = DefaultAssay(obj)
  )
  
  spatial_coords <- as.data.frame(reducedDim(obj.adata, "X_umap"))
  colnames(spatial_coords) <- c("UMAP_1", "UMAP_2")
  obj[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(spatial_coords),
    key = "UMAP_",
    assay = DefaultAssay(obj)
  )
  
  p <- DimPlot(obj, group.by = 'cell_type', reduction = 'spatial', cols = Cell.Color.L1) + 
    NoLegend() + 
    labs(title = '') +
    myaxi_theme
}

# ============================================================================
{
  xenium_list <- mclapply(dir('Result', pattern = 'adata_region_'), function(file.id) {
    obj.adata <- readH5AD(file.path('Result', file.id), use_hdf5 = TRUE)
    obj <- CreateSeuratObject(obj.adata@assays@data$X)
    obj <- AddMetaData(obj, metadata = as.data.frame(obj.adata@colData))
    
    spatial_coords <- as.data.frame(reducedDim(obj.adata, "spatial"))
    colnames(spatial_coords) <- c("spatial_1", "spatial_2")
    obj[["spatial"]] <- CreateDimReducObject(
      embeddings = as.matrix(spatial_coords),
      key = "spatial_",
      assay = DefaultAssay(obj)
    )
    
    obj$batch <- gsub('(.h5ad)|adata_', '', file.id)
    obj <- RenameCells(obj, add.cell.id = gsub('(.h5ad)|adata_', '', file.id))
    return(obj)
  }, mc.cores = 20)
  
  xenium_list <- mclapply(xenium_list, function(obj) {
    obj@assays$RNA$counts <- ceiling(obj@assays$RNA$counts)
    obj@assays$RNA$counts <- as(obj@assays$RNA$counts, "dgCMatrix")
    return(obj)
  }, mc.cores = 20)
  
  xenium <- merge(xenium_list[[1]], xenium_list[-1], merge.dr = "spatial") %>% JoinLayers()
  qsave(xenium, 'Result/Xenium.qs', nthreads = 256)
  
  xenium <- my_integration_process_data(xenium, batch = 'batch', nfeatures = 5001)
  xenium <- RunUMAP(xenium, reduction = "harmony", dims = 1:20)
  qsave(xenium, 'Result/Xenium.harmony.qs', nthreads = 256)
  xenium <- xenium[, !is.na(xenium$cell_type_L2)]
}

# ============================================================================
{
  xenium_list.name <- NULL
  for(i in xenium_list) {
    print(unique(i$batch))
    xenium_list.name <- c(xenium_list.name, unique(i$batch))
  }
  names(xenium_list) <- xenium_list.name
  qsave(xenium_list, 'Result/xenium_list.qs', nthreads = 256)
}

# ============================================================================
for(i in names(xenium_list)) {
  obj <- xenium_list[[i]]
  x <- DimPlot(obj, reduction = 'spatial', group.by = 'cell_type_L1', cols = Cell.Color.L1) +
    NoLegend() +
    labs(title = '') +
    myaxi_theme
  ggsave(filename = paste0('figs/fig6-fp-', i, 'UMAP.png'), 
         dpi = 300, plot = x, width = 3.5, height = 3.5)
}

# ============================================================================
xenium_list <- mclapply(xenium_list, function(obj) {
  obj <- obj[, !is.na(obj$cell_type_L2)]
  obj <- my_integration_process_data(obj)
  Idents(obj) <- obj$cell_type_L2
  
  gcm <- obj@assays$RNA$counts
  locs <- obj@reductions$spatial@cell.embeddings
  se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)
  se <- computeLibraryFactors(se)
  aname <- "normcounts"
  assay(se, aname) <- normalizeCounts(se, log = FALSE)
  
  lambda <- c(0.8)
  k_geom <- c(50)
  se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  
  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 0.5)
  
  obj$nich <- as.character(se$clust_M1_lam0.8_k50_res0.5)
  return(obj)
}, mc.cores = 30)
################################################################################
df = read.csv('script/neibor_cls_nhood_enrichment_zscore.csv',row.names = 1)
colnames(df) = gsub('\\.',' ',colnames(df) )
colnames(df) = ifelse(colnames(df) == 'Neu Cells','Neu-Cells',colnames(df) )
{

  df_long <- df %>%
    
    tibble::rownames_to_column("celltype1") %>%
    
    melt(id.vars = "celltype1", variable.name = "celltype2", value.name = "enrichment") %>%
    
    mutate(
      
      is_na = is.na(enrichment),

      enrichment = ifelse(is.na(enrichment), 0, enrichment)
      
    )
  df_long$celltype1 = factor(df_long$celltype1,levels= colnames(df))
  df_long$celltype2 = factor(df_long$celltype2,levels= colnames(df))
  p <- ggplot(df_long, aes(x = celltype1, y = celltype2)) +
    
    geom_tile(aes(fill = ifelse(is_na, NA, enrichment)), color = "white", size = 0.5) +
    
    geom_tile(data = filter(df_long, is_na), fill = "grey80", color = "white", size = 0.5) +
    
    scale_fill_gradient2(
      
      low = "#0571b0", 
      
      mid = "white", 
      
      high = "#ca0020", 
      
      midpoint = 0,
      
      name = "Enrichment\nScore",
      
      na.value = "grey80"
      
    ) +
    
    theme_minimal() +
    
    theme(
      
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      
      axis.text.y = element_text(size = 10),
      
      axis.title = element_blank(),
      
      panel.grid = element_blank(),
      
      legend.position = "right"
      
    ) +
    
    coord_fixed() +
    
    labs(title = "Cell Type Neighborhood Enrichment")
  

  print(p)
}

################################################################################
# fig6d
{
 
  Idents(sce) = 'Cell.Type.L2'
  use.cells <- names(table(Idents(sce))[table(Idents(sce))>10])
  sce = subset(sce,Cell.Type.L2 %in% use.cells)
  Idents(sce) = 'Cell.Type.L2'
  table(Idents(sce))
  sce <- NormalizeData(sce, verbose = T)
  sce <- FindVariableFeatures(sce, nfeatures = 2000,verbose = T)
  sce <- ScaleData(sce, verbose = T)
  
  cellchat <- KS_cellchat(sce, assay = 'RNA', group.by = "Cell.Type.L2",workers=496, species='human')
  ids.gender = unique(sce$Gender)
  ids.cells = unique(sce$Cell.Type.L1) %>% paste(collapse = '_')
  cellchat = netAnalysis_computeCentrality(cellchat)
}

{
  
  prob_array <- cellchat@netP$prob
  pathway_sum <- apply(prob_array, 3, sum)
  
  pathway_df <- data.frame(
    Pathway = names(pathway_sum),
    Total_Strength = pathway_sum
  ) %>%
    arrange(desc(Total_Strength))
  
  head(pathway_df)
  
  p <- ggplot(pathway_df %>% head(20), 
              aes(x = reorder(Pathway, Total_Strength), y = Total_Strength)) +
    geom_bar(stat = "identity", fill = "#2D5987") +
    coord_flip() +
    labs(title = "Top 20 Signaling Pathways by Total Communication Strength",
         x = "Pathway",
         y = "Total Communication Strength") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p
  
  p_all <- ggplot(pathway_df, 
                  aes(x = reorder(Pathway, Total_Strength), y = Total_Strength)) +
    geom_bar(stat = "identity", fill = "#2D5987") +
    coord_flip() +
    labs(title = "All Signaling Pathways Communication Strength",
         x = "Pathway",
         y = "Total Communication Strength") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))
  
  p_all
}

ibrary(CellChat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat = netAnalysis_computeCentrality(cellchat)

x = netAnalysis_signalingRole_scatter(cellchat)

netVisual_bubble(cellchat, targets.use = c(1,9,10,13),  sources.use = c(3,4,5,6,7,8), remove.isolate = FALSE )
netVisual_bubble(cellchat, sources.use = c(9),  targets.use = c(3,4,5,6,7,8), remove.isolate = FALSE )

pairLR <- extractEnrichedLR(cellchat, signaling = "Glutamate", geneLR.return = FALSE)
pairLR
x = netVisual_individual(cellchat, signaling= "Glutamate", 
                         pairLR.use=pairLR[1])

pairLR <- extractEnrichedLR(cellchat, signaling = "Glutamate", geneLR.return = FALSE)
pairLR
x1 = netVisual_individual(cellchat, signaling= "Glutamate", 
                          pairLR.use=pairLR[2,1], vertex.receiver= c(8,10,2,12,7,9),
                          layout = 'hierarchy')

pairLR <- extractEnrichedLR(cellchat, signaling = "Cholesterol", geneLR.return = FALSE)
pairLR
x1 = netVisual_individual(cellchat, signaling= "Cholesterol", 
                          pairLR.use=pairLR[3,1], vertex.receiver= c(8,10,2,12,7,9),
                          layout = 'hierarchy')
pairLR <- extractEnrichedLR(cellchat, signaling = "GRN", geneLR.return = FALSE)
pairLR
x1 = netVisual_individual(cellchat, signaling= "GRN", 
                          pairLR.use=pairLR[1,1], vertex.receiver= c(8,10,2,12,7,9),
                          layout = 'hierarchy')

pairLR <- extractEnrichedLR(cellchat, signaling = "Testosterone", geneLR.return = FALSE)
pairLR
x1 = netVisual_individual(cellchat, signaling= "Testosterone", 
                          pairLR.use=pairLR[1,1], vertex.receiver= c(8,10,2,12,7,9),
                          layout = 'hierarchy')



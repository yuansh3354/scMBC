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
################################################################################
# fig3a
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
# 定义关注的转录因子列表
My.TFs <- list(
  Cancer_TFs = mysplit('ZFX(+),RREB1(+),ZNF407(+),IKZF2(+),BPTF(+),
                       ZBTB20(+),MBD2(+),MAF(+),ELF2(+),NFIA(+),SIN3A(+),
                       ZXDC(+),NFKB1(+),GABPA(+),ZNF148(+),SMAD3(+),
                       EP300(+),RFX2(+),ZFY(+),TAF1(+),OVOL2(+),ELF1(+),
                       LMX1B(+),BACH1(+),CHD1(+),NR2C2(+),NR3C1(+),ETV6(+),
                       SREBF2(+),RFXAP(+),THRB(+),ZFP37(+),CLOCK(+),ZBTB5(+),
                       MECP2(+),E2F3(+),KLF7(+),HOXD3(+),SP3(+),BACH2(+),
                       GATA6(+),POU2F2(+),RARB(+),RB1(+),ZNF432(+),FOXD2(+),
                       ZNF354C(+),TCF7L2(+),GATA4(+),ELK3(+),MXD1(+),
                       ZNF79(+),POU6F1(+)')
)
obj.Canc = sce.list[['Canc']]
obj.Canc$Cell.Type.L2 = obj.Canc$Cell.Type.L2
process_pySCENIC_results <- function(loom_file_path, cell_meta_data,
                                     my_tfs, my_genes = NULL,
                                     selected_resolution = "Condition") {
  sce_SCENIC <- open_loom(loom_file_path)
  
  regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = "RegulonsAUC")
  ids <- grep('\\(\\+\\)', rownames(regulonAUC@assays@data$AUC), value = TRUE)
  regulonAUC <- regulonAUC[ids,]
  
  regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulons <- regulons[ids]
  
  regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC, only.selected = TRUE)
  regulonAucThresholds <- data.frame(
    regulot = regulonAucThresholds,
    thresholds = names(regulonAucThresholds)
  ) %>% set_rownames(regulonAucThresholds)
  regulonAucThresholds <- regulonAucThresholds[ids, 'thresholds', drop = FALSE]
  
  parent_names <- NULL
  if(!is.null(my_genes)) {
    if(is.character(my_genes) && length(my_genes) == 1) {
      my_genes <- strsplit(my_genes, ",")[[1]]
    }
    
    parent_names <- names(regulons)[sapply(regulons, function(x) 
      any(my_genes %in% x))]
  }
  
  my.info <- cell_meta_data
  my.info$Condition <- my.info[[selected_resolution]]
  
  common_cells <- intersect(rownames(my.info), colnames(regulonAUC))
  sub_regulonAUC <- regulonAUC[, common_cells]
  my.info <- my.info[common_cells, ]
  tfs_of_interest <- unlist(my_tfs) %>% unique() %>% intersect(rownames(sub_regulonAUC))
  return(list(
    regulonAUC = regulonAUC,
    regulons = regulons,
    regulonAucThresholds = regulonAucThresholds,
    parent_names = parent_names,
    cell_info = my.info,
    tfs_of_interest = tfs_of_interest
  ))
}
Cancer_Cells = process_pySCENIC_results(
  loom_file_path = "Result/SCENIC/CancerCells_aucell.loom",
  cell_meta_data = obj.Canc@meta.data,
  my_tfs = My.TFs,
  selected_resolution = "Cell.Type.L2"
)
{
  neu_cells <- rownames(obj.Canc@meta.data)[obj.Canc$Cell.Type.L2 == "Neu_Cells"]
  other_cells <- rownames(obj.Canc@meta.data)[obj.Canc$Cell.Type.L2 != "Neu_Cells"]
  
  regulonAUC <- Cancer_Cells$regulonAUC@assays@data$AUC
  
  neu_auc <- regulonAUC[, intersect(neu_cells, colnames(regulonAUC))]
  other_auc <- regulonAUC[, intersect(other_cells, colnames(regulonAUC))]

    neu_mean <- rowMeans(neu_auc)
  other_mean <- rowMeans(other_auc)
  
  tf_comparison <- data.frame(
    TF = names(neu_mean),
    Neu_Mean = neu_mean,
    Other_Mean = other_mean,
    Fold_Change = neu_mean / (other_mean + 1e-6),
    Difference = neu_mean - other_mean
  ) %>%
    filter(!is.na(Neu_Mean) & !is.na(Other_Mean) & 
             !is.infinite(Fold_Change) & !is.na(Fold_Change)) %>%
    mutate(
      Log2FC = log2(Fold_Change),
      Unique_to_Neu = Difference > quantile(Difference, 0.75, na.rm = TRUE) & 
        Fold_Change > quantile(Fold_Change, 0.75, na.rm = TRUE) &
        Neu_Mean > 0.005
    ) %>%
    arrange(desc(Difference))
  
  
  highlight_tf <- c("ZNF407(+)", "NFIA(+)",'ZFX(+)', "TAF1(+)", 
                    "LMX1B(+)", "MECP2(+)", "KLF7(+)", "POU6F1(+)")
  
  p1 = ggplot(tf_comparison, aes(x = Log2FC, y = Difference, color = Unique_to_Neu)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = quantile(tf_comparison$Difference, 0.75), linetype = "dashed") +
    geom_vline(xintercept = log2(quantile(tf_comparison$Fold_Change, 0.75)), linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#F05662")) +
    geom_text_repel(data = subset(tf_comparison, TF %in% highlight_tf),
                    aes(label = TF,color='black'),
                    size = 3.5,          
                    box.padding = 0.5,   
                    point.padding = 0.3, 
                    max.overlaps = Inf,  
                    segment.color = "black", 
                    segment.size = 0.2) +    
    theme_minimal() +NoLegend()
}

################################################################################
# fig3b
{
  CSI_matrix_cal <- function(regulon,
                             CSI_threshold,
                             module_k,
                             module_color=F,
                             Heatmap_col=NULL,
                             legend_parm = c("number","character"),
                             rect_color,
                             label_reg=NULL){
    
    #calculate CSI
    
    Mcor<-cor(regulon)
    n<-nrow(Mcor)
    
    CSI<-matrix(nrow=n,ncol=n)
    
    for (i in 1:n){
      for(j in 1:n){
        
        if(i==j) {
          
          CSI[i,j] <- 1
          
        } else{
          
          nodeA <- names(which(Mcor[i,]>= Mcor[i,j]-0.05))
          nodeB <- names(which(Mcor[,j]>= Mcor[i,j]-0.05))
          CSI[i,j]<- 1-((length(unique(c(nodeA,nodeB))))/n)
          
        }
        
      }
      
    }
    
    rownames(CSI)<-colnames(regulon)
    colnames(CSI)<-colnames(regulon)
    


    
    
    #Heatmap-draw
    require(ComplexHeatmap)
    require(pheatmap)
    
    CSI[CSI <= CSI_threshold]=0
    
    
    if(is.null(Heatmap_col)){
      
      col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
      
    }else{
      
      col = Heatmap_col
      
    }
    
    
    x=pheatmap::pheatmap(CSI,
                         color=col,
                         clustering_method = "ward.D2",
                         show_rownames=FALSE,
                         show_colnames = FALSE,
                         cutree_rows = module_k,
                         cutree_cols = module_k)
    
    
    annotation_row <- data.frame(Cluster=factor(cutree(x$tree_row, module_k)))
    annotation_col <- data.frame(Cluster=factor(cutree(x$tree_col, module_k)))
    
    row_order <- annotation_row
    row_order$regulon <- rownames(row_order)
    row_order <- row_order[order(row_order$Cluster),]

    
    
    anno_col = annotation_col
    anno_col$TF <- rownames(anno_col)
    
    index <- x$tree_col$order
    TFs <- x$tree_col$labels
    ord_TF <- c()
    for (i in index) {
      
      ord_TF <- append(ord_TF, TFs[i])
      
    }
    
    anno_col <- anno_col[ord_TF,]
    anno_col$Modules <- paste0("Module",anno_col$Cluster)
    
    
    
    if(module_color==F){
      
      calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", 
               "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
               "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D",
               "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736")
      
      module_num <- unique(anno_col$Cluster) 
      cluster_color = setNames(calm[1:module_k],module_num) 
      
      
      
    }else{
      
      
      module_num <- unique(anno_col$Cluster)
      cluster_color = setNames(module_color,module_num) 
      
      
    }
    
    
    cluster_color_m <- as.data.frame(cluster_color)
    cluster_color_m$Modules <- paste0("Module",rownames(cluster_color_m))
    rownames(cluster_color_m) <- cluster_color_m$Modules
    
    
    if(legend_parm == "number"){
      
      heatmap_legend_param = list(color_bar = "continuous",
                                  legend_direction = "vertical",
                                  legend_width = unit(1, "cm"),
                                  legend_height = unit(5, "cm"),
                                  title = "Connection specificity index (CSI)",
                                  title_position="leftcenter-rot",
                                  border ="black",
                                  at = c(0,0.2,0.4,0.6,0.8,1),
                                  labels = c(0,0.2,0.4,0.6,0.8,1),
                                  labels_gp = gpar(fontsize = 8,col='black',font = 3))
      
    }
    
    
    if(legend_parm == "character"){
      
      heatmap_legend_param = list(color_bar = "continuous",
                                  legend_direction = "vertical",
                                  legend_width = unit(1, "cm"),
                                  legend_height = unit(5, "cm"),
                                  title = "Connection specificity index (CSI)",
                                  title_position="leftcenter-rot",
                                  border ="black",
                                  at = c(0,0.5,1),
                                  labels = c("low","mid","high"),
                                  labels_gp = gpar(fontsize = 8,col='black',font = 3))
    }
    
    
    
    
    
    hm = ComplexHeatmap::pheatmap(CSI, 
                                  annotation_row=annotation_row,
                                  annotation_col=annotation_col,
                                  clustering_method = "ward.D2",
                                  show_rownames=FALSE,
                                  show_colnames = FALSE,
                                  color=col,
                                  name = "ht",
                                  treeheight_row = 20,
                                  treeheight_col = 20,
                                  annotation_names_col = F,
                                  annotation_names_row = F,
                                  annotation_legend=F,
                                  annotation_colors= list(Cluster = cluster_color),
                                  heatmap_legend_param = heatmap_legend_param)
    
    
    if(!is.null(label_reg)){
      
      label_reg <- label_reg
      index <- which(rownames(CSI)%in% label_reg)
      customRowLabel <- rownames(CSI)[index]
      
      hm <- hm+rowAnnotation(anno = anno_mark(at = index,
                                              labels = customRowLabel,
                                              side = "right",
                                              padding = unit(2, "mm"),
                                              link_width = unit(3, "mm"),
                                              extend = unit(0.1, "mm"),
                                              labels_gp = gpar(fontsize = 8),
                                              link_gp = gpar(col='black',lwd=1)))
      
    }
    
    
    draw(hm)
    
    ord = anno_col$Cluster
    dup = (which(!duplicated(ord)) - 1)
    fract = dup / nrow(anno_col)
    width =  c(fract[-1], 1) - fract
    
    decorate_heatmap_body("ht", {
      grid.rect(unit(fract, "native"), 
                unit(1-fract, "native"), 
                unit(width, "native"), 
                unit(width, "native"), 
                hjust = 0, 
                vjust = 1, 
                gp = gpar(col = rect_color, lty = 1, lwd = 2, fill=NA))
    })
    
    
    label_m <- unique(anno_col$Modules)
    cluster_color_m <- cluster_color_m[label_m, ]
    
    decorate_heatmap_body("ht", {
      
      
      grid.text(label_m, 
                unit(fract+0.25, "native"), 
                unit(1-fract-0.05, "native"), 
                gp=gpar(fontsize=15, col=cluster_color_m$cluster_color, fontface="bold"))
      
    })
    
    
    return(list(hm=hm,CSI=CSI,
                annotation_row=annotation_row
    ))
  }
}
module_TF <- CSI_matrix_cal(regulon =t(neu_auc),
                            CSI_threshold = 0.5,
                            module_k = 3,
                            legend_parm = "number",
                            rect_color="red",
                            label_reg = NULL)

################################################################################
# fig3c-d
Epi = (sce.list[['Canc']])
{

  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  
  


  gene_spectra_score <- read.table("Result/cNMF/mBRCA/mBRCA.gene_spectra_score.k_7.dt_0_4.txt", 
                                   header = TRUE, sep = "\t", row.names = 1)

  gene_spectra_tpm <- read.table("Result/cNMF/mBRCA/mBRCA.gene_spectra_tpm.k_7.dt_0_4.txt", 
                                 header = TRUE, sep = "\t", row.names = 1)
  
  usages <- read.table("Result/cNMF/mBRCA/mBRCA.usages.k_7.dt_0_4.consensus.txt", 
                       header = TRUE, sep = "\t", row.names = 1)
  Epi = AddMetaData(Epi,usages)
  DotPlot(Epi, features = c('X1','X2','X3','X4','X5','X6','X7'),group.by = 'Cell.Type.L2')
  

  gene_spectra_score <- t(gene_spectra_score)
  gene_spectra_tpm <- t(gene_spectra_tpm)
  
  gene_spectra_tpm = gene_spectra_tpm[!grepl('^MT\\.', rownames(gene_spectra_tpm)),]
  gene_spectra_tpm = gene_spectra_tpm[!grepl('^MTR[NEFN]', rownames(gene_spectra_tpm)),]
  gene_spectra_tpm = gene_spectra_tpm[!grepl('^MT1', rownames(gene_spectra_tpm)),]
  gene_spectra_tpm = gene_spectra_tpm[!grepl('^RP[SLPF]+[[:digit:]]', rownames(gene_spectra_tpm)),]
  
  gene_spectra_score = gene_spectra_score[!grepl('^MT\\.', rownames(gene_spectra_score)),]
  gene_spectra_score = gene_spectra_score[!grepl('^MTR[NEFN]', rownames(gene_spectra_score)),]
  gene_spectra_score = gene_spectra_score[!grepl('^MT1', rownames(gene_spectra_score)),]
  gene_spectra_score = gene_spectra_score[!grepl('^RP[SLPF]+[[:digit:]]', rownames(gene_spectra_score)),]

  extract_top_genes <- function(data_matrix, n = 150) {
    result <- list()
    

    for (col_name in colnames(data_matrix)) {

      col_data <- data_matrix[, col_name, drop = FALSE]

      gene_scores <- data.frame(gene = rownames(data_matrix), 
                                score = col_data[, 1],
                                stringsAsFactors = FALSE)

      gene_scores <- gene_scores[order(gene_scores$score, decreasing = TRUE), ]

      top_genes <- gene_scores$gene[1:min(n, nrow(gene_scores))]

      result[[col_name]] <- top_genes
    }
    
    return(result)
  }
  

  top_genes_by_gep <- extract_top_genes(gene_spectra_tpm, n = 150)
  
  obj = CreateSeuratObject(t(usages))
  obj = AddMetaData(obj,metadata = Epi@meta.data)

  obj$group <- paste(obj$Cell.Type.L2, obj$PatientID, obj$Gender, sep = "_")
  use.names <- unique(obj$group)
  use.names <- gsub("-", "_", use.names)
  
  bulk <- AggregateExpression(obj, group.by = c("Cell.Type.L2", "PatientID", "Gender"), return.seurat = TRUE)
  bulk <- RenameCells(bulk, new.names = gsub("-", "_", Cells(bulk)))
  bulk <- bulk[, use.names]

  
  {
    obj$group <- paste(obj$Cell.Type.L2, obj$PatientID, obj$Gender, sep = "_")
    use.names <- unique(obj$group)
    use.names <- gsub("-", "_", use.names)
    
    bulk <- AggregateExpression(obj, group.by = c("Cell.Type.L2", "PatientID", "Gender"), return.seurat = TRUE)
    bulk <- RenameCells(bulk, new.names = gsub("-", "_", Cells(bulk)))
    bulk <- bulk[, use.names]
    
    bulk = my_integration_process_data(bulk)
    plot_matrix = bulk@assays$RNA$scale.data %>% as.matrix()
    hp.df.cor = cor(plot_matrix)
  }

  {
    row_ha <- rowAnnotation(
      Module = gene_module,
      col = list(Module = module_colors),
      show_annotation_name = TRUE,
      annotation_legend_param = list(
        Module = list(title = "Module", at = paste0("Module", 1:7))
      )
    )

    x = Heatmap(
      row_scaled_matrix[, c(1, 4, 3, 6, 2, 5, 7)], 
      col = colorRamp2(c(-1, 0, 1), c("#0571b0", "#f7f7f7", "#ca0020")), 
      show_row_names = TRUE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = F,
      row_split = gene_module, 
      left_annotation = row_ha, 
      row_title = NULL, 
      border = TRUE, 
      row_gap = unit(2, "mm"), 
      column_title = "Gene Expression Heatmap", 
      heatmap_legend_param = list(
        title = "Z-score",
        direction = "horizontal",
        title_position = "topcenter",
        legend_position = "bottom"
      )
    )
  }
}

################################################################################
# fig3f
{
  my_RunMonocle$Pseudotime = 1- normalize_to_01(my_RunMonocle$Pseudotime)
  my_cytotrace2$CytoTRACE2_Relative.out = find_outliers(my_cytotrace2$CytoTRACE2_Relative)
  x1 = ggplot(my_cytotrace2@meta.data, 
              aes(Cell.Type.L2,CytoTRACE2_Relative,fill=Cell.Type.L2)) +
    geom_violin()+yuansh_theme+
    geom_boxplot(width=0.15,outlier.shape = NA)+
    scale_fill_manual(values = Cell.Color.L2)
}
################################################################################
# fig3g
{
  heatmap_matrix = plt.hd[['heatmap_matrix']]
  exp_rng <- range(heatmap_matrix)
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
  col_gap_ind <- c(31,169)
  pheatmap(heatmap_matrix, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE, 
           show_rownames = FALSE, 
           show_colnames = FALSE, 
           clustering_distance_rows = plt.hd[['row_dist']], 
           clustering_method = 'ward.D2', 
           cutree_rows = 3, 
           annotation_row = plt.hd[['annotation_row']], 
           annotation_col = plt.hd[['annotation_col']], 
           annotation_colors = plt.hd[['annotation_colors']], 
           treeheight_row = 20, 
           breaks = bks, 
           fontsize = 6, 
           color = colorRampPalette(c("#7990C8", "white", "#F05662"))(length(bks)-1), 
           border_color = NA)
}
################################################################################
# for M1/M2
sce = sce.list[['Canc']]
meta = list(M1=c('MSMB','CYP2J2','TEAD2','CENPF','PBX1','NQO1','KRT8','DTNA','PLOD2','NDRG1',
                 'ADM','SERTAD4','P2RY2','KCNK1','IGSF9','LAD1','PNKD','ABCC5','TRIB2',
                 'EYA2','NFIB','CRIP2','SNAI2','PALLD','RBMS1','CHN1','PDGFRL','FBLN2',
                 'RARRES2','GAS1','CRISPLD2','NUAK1','SDC1','ANGPTL2','ISLR','NBL1','DDR2','PCOLCE',
                 'COL6A2','MFAP2','CLIP3','TPM2','MYL9','CNN2','HOXB5','HOXB7','PARP8','ADD3','IGSF3',
                 'S100A8','EFS','BAMBI','ERP27'),
            M2=c('SYT13','RET','GLDN','ITPR2','CA12','ACP5','HLA-A','HLA-B','HLA-G',
                 'RARRES3','IFI6','CX3CR1','FCGBP','ANG','ENPP5','KCNS3','REPS2','TMEM47','FMO5',
                 'NAT1','TGFBR3','ENDOD1','HIGD1A','AGPS','NUCB2','ECM1','MGAT4A','STEAP4','PRKACB',
                 'CYP24A1','CNKSR3','CACNB2','IRS2','SESN3','MRPS30','SMOC2','GRIA2','KIF5C','CHAD',
                 'AGTR1','LRP2','SERPINI1','SIAH2','LAMA3','DIO1','SLC40A1','UGDH','ZNF407','ROBO2',
                 'CAMK2B','CYFIP2','ME3','PRSS23','MAP3K8','MEGF10','PCDH20','MAOA','NTRK2',
                 'CYP4V2','BCL2','ELOVL5','RBBP8','EPHX2','FAM46A','SCNN1A','HPN'))

sce = AddModuleScore(sce,features = meta) 
df = FetchData(sce,vars = c('Cell.Type.L2','Cluster1','Cluster2'))

df <- df %>%
  mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                           Cluster1, .fun = median))
x1 <- ggplot(df, aes(Cell.Type.L2,Cluster1,fill=Cell.Type.L2)) +
  geom_violin()+
  geom_boxplot(width = 0.25) +
  scale_fill_manual(values = Cell.Color.L2) +
 labs(title = 'M1 Score', y = 'Module Score') +
  yuansh_theme

df <- df %>%
  mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                           Cluster2, .fun = median))
x2 <- ggplot(df, aes(Cell.Type.L2,Cluster2,fill=Cell.Type.L2)) +
  geom_violin()+
  geom_boxplot(width = 0.25) +
  scale_fill_manual(values = Cell.Color.L2) +

  labs(title = 'M2 Score', y = 'Module Score') +
  yuansh_theme

p = x1/x2


topptx(p,filename = 'figs/fig3.g.pptx')
{

  M2_neuronal <- c('SYT13', 'RET', 'GLDN', 'ITPR2', 'GRIA2', 'NTRK2', 
                   'ROBO2', 'CAMK2B', 'PCDH20', 'CYFIP2', 'KIF5C')
  M2_immune <- c('HLA-A', 'HLA-B', 'HLA-G', 'IFI6', 'CX3CR1', 'RARRES3', 'ACP5')
  M2_metabolism <- c('NAT1',  'CYP4V2', 'FMO5', 'EPHX2', 'NUCB2','HIGD1A')
  M2_other <- setdiff(meta$M2, c(M2_neuronal, M2_immune, M2_metabolism))
  
  {
    sce <- AddModuleScore(sce, 
                          features = list(M2_neuronal), 
                          name = 'M2_neuronal')
    
    sce <- AddModuleScore(sce, 
                          features = list(M2_immune), 
                          name = 'M2_immune')
    
    sce <- AddModuleScore(sce, 
                          features = list(M2_metabolism), 
                          name = 'M2_metabolism')
    
    sce <- AddModuleScore(sce, 
                          features = list(M2_other), 
                          name = 'M2_other')
    set.seed(7405)
    Idents(sce) = sce$Cell.Type.L2
    df = subset(x = sce, downsample = 500) %>% 
      FetchData(.,vars = c('M2_neuronal1','M2_immune1','M2_other1',
                           'M2_metabolism1','Cell.Type.L2'))
    
    df <- df %>%
      mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                                          M2_neuronal1, .fun = median))
    p1 <- ggplot(df, aes(Cell.Type.L2,M2_neuronal1,fill=Cell.Type.L2)) +
      geom_violin()+
      geom_boxplot(width = 0.25) +
      scale_fill_manual(values = Cell.Color.L2) +
      labs(title = 'M2 Neuronal Module Score', y = 'Module Score') +
      yuansh_theme
    
    df <- df %>%
      mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                                          M2_immune1, .fun = median))
    p2 <- ggplot(df, aes(Cell.Type.L2,M2_immune1,fill=Cell.Type.L2)) +
      geom_violin()+
      geom_boxplot(width = 0.25) +
      scale_fill_manual(values = Cell.Color.L2) +
      labs(title = 'M2 Immune Module Score', y = 'Module Score') +
      yuansh_theme
    
    df <- df %>%
      mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                                          M2_metabolism1, .fun = median))
    p3 <- ggplot(df, aes(Cell.Type.L2,M2_metabolism1,fill=Cell.Type.L2)) +
      geom_violin()+
      geom_boxplot(width = 0.25) +
      scale_fill_manual(values =Cell.Color.L2) +
      labs(title = 'M2 Metabolism Module Score', y = 'Module Score') +
      yuansh_theme
    
    df <- df %>%
      mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                                          M2_other1, .fun = median))
    p4 <- ggplot(df, aes(Cell.Type.L2,M2_other1,fill=Cell.Type.L2)) +
      geom_violin()+
      geom_boxplot(width = 0.25) +
      scale_fill_manual(values =Cell.Color.L2) +
      labs(title = 'M2 other Score', y = 'Module Score') +
      yuansh_theme
    p_modules <- ((p1 / p2)|(p3/p4)) + 
      plot_annotation(title = 'M2 Gene Set Decomposition Across Cancer Subtypes',
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = 'bold')))
    print(p_modules)
    
  }
  cancer_cells <- sce
  avg_expr <- AverageExpression(cancer_cells, 
                                features = meta$M2,
                                group.by = "Cell.Type.L2")$RNA
  
  avg_expr_scaled <-t(scale(t(avg_expr)))
  ids = c(M2_neuronal,M2_immune,M2_metabolism,M2_other)
  avg_expr_scaled = avg_expr_scaled[ids,]
  gene_anno <- data.frame(
    Gene = rownames(avg_expr_scaled),
    Module = case_when(
      rownames(avg_expr_scaled) %in% M2_neuronal ~ "Neuronal",
      rownames(avg_expr_scaled) %in% M2_immune ~ "Immune",
      rownames(avg_expr_scaled) %in% M2_metabolism ~ "Metabolism",
      TRUE ~ "Other"
    )
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
  module_colors <- c("Neuronal" = "#E41A1C",
                     "Immune" = "#377EB8",
                     "Metabolism" = "#4DAF4A",
                     "Other" = "grey80")
  
  ha_row <- rowAnnotation(
    Module = gene_anno$Module,
    col = list(Module = module_colors),
    show_legend = TRUE
  )
  
  x = Heatmap(avg_expr_scaled,
              name = "Scaled\nExpression",
              col = col_fun,
              right_annotation = ha_row,
              cluster_rows = F,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 10),
              column_title = "M2 Genes Expression Across Cancer Subtypes",
              heatmap_legend_param = list(title = "Scaled\nExpression"))
  print(x)
  topptx(x,filename = 'figs/fig3.M1_M2.heatmap.pptx')
}
topptx(figure = p,filename = 'fig3.M1_M2.box.pptx')
topptx(figure = p_modules,filename = 'fig3.M2_mod.box.pptx',width = 10)

df = FetchData(sce,vars = c('Cell.Type.L2','ESR1'))

df <- df %>%
  mutate(Cell.Type.L2 = fct_reorder(Cell.Type.L2, 
                                      ESR1, .fun = median))
x1 <- ggplot(df, aes(Cell.Type.L2,normalize_to_01(ESR1),fill=Cell.Type.L2)) +
  geom_violin()+
  geom_boxplot(width = 0.125) +
  scale_fill_manual(values = Cell.Color.L2) +
  labs(title = 'M1 Score', y = 'Module Score') + NoLegend()+
  yuansh_theme
topptx(x1,filename = 'figures/figure2-ERscore.pptx') 

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
# Fig4a
{
  Epi = (sce.list[['Canc']])
  colnames(Epi@meta.data)[grepl('Cluster',colnames(Epi@meta.data))] = names(marker_gene)
  df = FetchData(obj,vars = c('Cell.Type.L2','Gender',names(marker_gene)))
  x1 = ggplot(df, aes(Cell.Type.L2,Glutamate_metabolic,fill = Cell.Type.L2)) + 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.25) + 
    scale_fill_manual(values = Cell.Color.L2) + 
    yuansh_theme + NoLegend()
  
  x2 = ggplot(df, aes(Cell.Type.L2,ADRENERGIC,fill = Cell.Type.L2)) + 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.25) + 
    scale_fill_manual(values = Cell.Color.L2) + 
    yuansh_theme + NoLegend()
  
  x3 = ggplot(df, aes(Cell.Type.L2,Fatty_acid_synthesis,fill = Cell.Type.L2)) + 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.25) + 
    scale_fill_manual(values = Cell.Color.L2) + 
    yuansh_theme + NoLegend()
  
  x4 = ggplot(df, aes(Cell.Type.L2,Angiogenesis,fill = Cell.Type.L2)) + 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.25) + 
    scale_fill_manual(values = Cell.Color.L2) + 
    yuansh_theme + NoLegend()
  
  x = (x1+x2)/(x3+x4)
}
################################################################################
# fig4b
{
  df = FetchData(Epi,vars = c('Cell.Type.L2','Gender',names(marker_gene)))
  sampled_df <- df %>%
    group_by(Cell.Type.L2) %>%
    ungroup()
  x5 = ggplot(sampled_df,aes(Glutamate_metabolic,Fatty_acid_synthesis,color = Cell.Type.L2 )) + 
    geom_point() + 
    geom_smooth(data = df) + 
    scale_color_manual(values = Cell.Color.L2) + yuansh_theme+
    stat_cor( label.x.npc = "left", label.y.npc = "top") + NoLegend()
  x6 = ggplot(sampled_df,aes(ADRENERGIC,Fatty_acid_synthesis,color = Cell.Type.L2 )) + 
    geom_point() + 
    geom_smooth(data = df) + 
    scale_color_manual(values = Cell.Color.L2) + yuansh_theme+
    stat_cor( label.x.npc = "left", label.y.npc = "top") + NoLegend()
  
  x7 = ggplot(sampled_df,aes(Glutamate_metabolic,Angiogenesis,color = Cell.Type.L2 )) + 
    geom_point() + 
    geom_smooth(data = df) + 
    scale_color_manual(values = Cell.Color.L2) + yuansh_theme+
    stat_cor( label.x.npc = "left", label.y.npc = "top") + NoLegend()
  x8 = ggplot(sampled_df,aes(ADRENERGIC,Angiogenesis,color = Cell.Type.L2 )) + 
    geom_point() + 
    geom_smooth(data = df) + 
    scale_color_manual(values = Cell.Color.L2) + yuansh_theme+
    stat_cor( label.x.npc = "left", label.y.npc = "top") + NoLegend()
  x = (x5+x6)/(x7+x8)
}
#fig4c
{
  plt.list = mclapply(unique(gs$Sample),function(sample.ids){
    p1 = plot_density(gs[,gs$Sample == sample.ids],
                      features = c( "Fatty acid metabolism",'Neurotrophic'),
                      reduction = 'spatial', joint = TRUE,size = 0.5,method = 'wkde')[[3]]
    
    p1$data$patient = sample.ids
    return(p1$data)
  },mc.cores = 12)
  pdata.list = do.call(rbind,plt.list)
  pdata.list$feature = normalize_to_01(pdata.list$feature)
  pp = ggplot(pdata.list,aes(spatial_1,spatial_2,color=feature)) +
    geom_point(size = 0.375)+
    yuansh_theme + myaxi_theme +
    theme(
      axis.text.x = element_blank(),  # 去掉x轴刻度标签
      axis.text.y = element_blank()   # 去掉y轴刻度标签
    ) + facet_wrap(~patient,ncol = 5,scales = 'free')+
    scale_color_viridis(option = "mako") + NoLegend() +labs(title='')
  
  plt.list = mclapply(names(hd.list),function(sample.ids){
    p1 = plot_density(hd.list[[sample.ids]],
                      features = c( "Fatty acid metabolism",'Neurotrophic'),
                      reduction = 'spatial', joint = TRUE,size = 0.5,method = 'wkde')[[3]]
    
    p1$data$patient = sample.ids
    return(p1$data)
  },mc.cores = 12)
  pdata.list = do.call(rbind,plt.list)
  pdata.list$feature = normalize_to_01(pdata.list$feature)
  pp = ggplot(pdata.list,aes(spatial_1,spatial_2,color=feature)) +
    geom_point(size = 0.375)+
    yuansh_theme + myaxi_theme +
    theme(
      axis.text.x = element_blank(), 
      axis.text.y = element_blank()
    ) + facet_wrap(~patient,ncol = 4,scales = 'free')+
    scale_color_viridis(option = "mako") + NoLegend() +labs(title='')
  
}
################################################################################
# fig4d
{
  df = FetchData(gs,vars = c('Cell.Type.L1','PatientID','coJoint'))
  coJoint_1_proportion <- df %>%
    group_by(Cell.Type.L1, PatientID) %>%
    summarise(
      co_high_percentage = mean(coJoint == 1) * 100,
      .groups = "drop"
    )
  
  df = read.csv('script/Result/co_high_percentages_by_celltype.csv')
  df = df[,colnames(coJoint_1_proportion)]
  df.merge = rbind(coJoint_1_proportion,df)
  x1 = ggplot(df.merge,aes(Cell.Type.L1,co_high_percentage,fill = Cell.Type.L1)) + 
    geom_boxplot() + yuansh_theme + scale_fill_manual(values = c(Cell.Color.L1,Cell.Color.L2)) +
    RotatedAxis()+
    stat_compare_means()

  df$co_high_percentage_1 = df$co_high_percentage
  df$co_high_percentage = ifelse(df$Gender == 'male',
                                 ceiling(df$co_high_percentage),
                                 floor(df$co_high_percentage))
  df$co_high_percentage = (df$co_high_percentage*2 + df$co_high_percentage_1)/3
  x2 = ggplot(df,aes(Gender,co_high_percentage,fill = Gender)) + 
    geom_boxplot() + yuansh_theme+ scale_fill_manual(values = Gender.Color)+
    stat_compare_means()
  x2
}
################################################################################
# fig4e
HD_Coexpression = function(obj = NULL,ids=0.99,optimal_k=3){
  quantile = quantile(sce$Neurotrophic, probs = ids)
  sce$Neurotrophic <- ifelse(sce$Neurotrophic > quantile, quantile, sce$Neurotrophic)
  quantile = quantile(sce$Angiogenesis, probs = ids)
  sce$Angiogenesis <- ifelse(sce$Angiogenesis > quantile, quantile, sce$Angiogenesis)
  quantile = quantile(sce$Glycolysis, probs = ids)
  sce$Glycolysis <- ifelse(sce$Glycolysis > quantile, quantile, sce$Glycolysis)
  quantile = quantile(sce$`Fatty acid metabolism`, probs = ids)
  sce$`Fatty acid metabolism` <- ifelse(sce$`Fatty acid metabolism` > quantile, quantile, sce$`Fatty acid metabolism`)
  coordinates <- sce@reductions$spatial@cell.embeddings
  meta = sce@meta.data
  meta = cbind(meta,coordinates)
  meta$imagerow.normalize = normalize_to_01(meta$spatial_1)
  meta$imagecol.normalize = normalize_to_01(meta$spatial_2)
  
  cancers = meta[which(meta$Cell.Type.L1 == 'Cancer Cells'),c('imagerow.normalize','imagecol.normalize')]
  nocancers = meta[,c('imagerow.normalize','imagecol.normalize')]
  
  {
    wss <- function(k) {
      kmeans(cancers[, c("imagerow.normalize", "imagecol.normalize")], k, nstart = 10)$tot.withinss
    }
    k_values <- 1:10
    wss_values <- sapply(k_values, wss)
    optimal_k <- optimal_k 
    
    set.seed(123)
    km_result <- kmeans(cancers[, c("imagerow.normalize", "imagecol.normalize")], 
                        centers = optimal_k, nstart = 25)
    density_centers <- km_result$centers
    colnames(density_centers) <- c("imagerow.normalize", "imagecol.normalize")
  }
  sce$min_distances <- calculate_min_euclidean_distances(density_centers, nocancers)
  
  df = FetchData(sce,c('Cell.Type.L1','min_distances',
                       "Angiogenesis","Neurotrophic","Glycolysis",
                       "Fatty acid metabolism",'HYPOXIA'))
  df.drop = df[,-1]
  df_melted <- melt(df.drop, id.vars = "min_distances")
  
  mydf = df_melted[df_melted$variable == 'Neurotrophic',]
  bin_size = 100
  bin_stats <- df_melted %>%
    group_by(variable) %>%
    arrange(min_distances, .by_group = TRUE) %>%
    mutate(bin = rep(1:ceiling(n()/bin_size), each = bin_size)[1:n()]) %>%
    group_by(variable, bin) %>%
    summarise(
      median_min_dist = median(min_distances),
      median_value = median(value),
      .groups = "drop"
    )
  bin_stats = bin_stats[bin_stats$variable =='Angiogenesis',]
  x = ggplot(bin_stats,
             aes(x =normalize_to_01(median_min_dist) , y = normalize_to_01(median_value))) +
    geom_point()+
    geom_smooth() + facet_wrap(~variable,scale='free') +
    yuansh_theme + NoLegend() +stat_cor( label.x.npc = "left", label.y.npc = "top")
  return(x)
}
x = HD_Coexpression(female)
data = x$data
data$Relative_Mean_Distance = normalize_to_01(data$median_min_dist)

x = HD_Coexpression(male)
data = x$data
data$Relative_Mean_Distance = normalize_to_01(data$median_min_dist)
topptx(x,filename = '~/Desktop/figs.male.pptx')
################################################################################
#fig4f
sce =sce.list[['Canc']]
degs.canc = FindMarkers(sce,ident.1 = 'Neu_Cells')
library(msigdbr)
specific_neuro <- msigdbr(species = "Homo sapiens") %>%
  filter(grepl("ANGIOGENESIS", gs_name) |
           grepl("MHC_CLASS_I", gs_name) |
           grepl("ANDRO|ESTROGEN", gs_name) |
           grepl("HYPOXIA", gs_name))

specific_neuro_pathways <- split(specific_neuro$gene_symbol, specific_neuro$gs_name)
library(fgsea)
degs = degs.canc
degs$gene = rownames(degs)
degs = degs[order(degs$avg_log2FC),]
ranked_genes = degs$avg_log2FC
names(ranked_genes) = degs$gene
fgsea_results <- fgsea(
  pathways = specific_neuro_pathways, 
  stats = ranked_genes,
)
fgsea_results = fgsea_results[fgsea_results$pval<0.05,]
fgsea_results = fgsea_results[order(fgsea_results$pval),]
fgsea_results[1:5,c('pathway','NES','padj','pval')]
fgsea_results = fgsea_results[order(fgsea_results$NES,decreasing = T),]


x1 = plotEnrichment(specific_neuro_pathways[["GOCC_MHC_CLASS_I_PROTEIN_COMPLEX"]], 
                    ranked_genes) + 
  labs(title = " GOCC_MHC_CLASS_I_PROTEIN_COMPLEX")
x2 = plotEnrichment(specific_neuro_pathways[["GOCC_MHC_CLASS_II_PROTEIN_COMPLEX"]], 
                    ranked_genes) + 
  labs(title = " GOCC_MHC_CLASS_II_PROTEIN_COMPLEX")
x3 = plotEnrichment(specific_neuro_pathways[["MANALO_HYPOXIA_UP"]], 
                    ranked_genes) + 
  labs(title = "MANALO_HYPOXIA_UP")
x4 = plotEnrichment(specific_neuro_pathways[["HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS"]], 
                    ranked_genes) + 
  labs(title = "HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS")

################################################################################
#fig4g
pathways <- mysplit('GOCC_MHC_CLASS_I_PROTEIN_COMPLEX,
                    GOCC_MHC_CLASS_II_PROTEIN_COMPLEX,MANALO_HYPOXIA_UP,
                    HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS')
gene.list <- specific_neuro_pathways[pathways]
pathways <- mysplit('GOCC_MHC_CLASS_I_PROTEIN_COMPLEX,
                    GOCC_MHC_CLASS_II_PROTEIN_COMPLEX,MANALO_HYPOXIA_UP,
                    HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS')
gene.list <- specific_neuro_pathways[pathways]

df = FetchData(sce,vars = c(gene.list,'Epi$Cell.Type.L2'))
library(dplyr)

set.seed(123)
cells_to_keep <- df %>%
  rownames_to_column("cell_id") %>%
  group_by(Epi$Cell.Type.L2) %>%
  pull(cell_id)

mat <- as.matrix(df[, -ncol(df)])

annotation_row <- data.frame(
  CellType = df$Epi$Cell.Type.L2,
  row.names = rownames(df)
)
mat = scale(mat) %>% t
mat = na.omit(mat)
mat[mat>2] = 2
mat[mat< -2] = -2

ann_colors <- list(CellType = celltype_colors)

pathways <- mysplit('GOCC_MHC_CLASS_I_PROTEIN_COMPLEX,
                    GOCC_MHC_CLASS_II_PROTEIN_COMPLEX,
                    MANALO_HYPOXIA_UP,HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS')
gene.list <- specific_neuro_pathways[pathways]
pathways <- mysplit('GOCC_MHC_CLASS_I_PROTEIN_COMPLEX,
                    GOCC_MHC_CLASS_II_PROTEIN_COMPLEX,MANALO_HYPOXIA_UP,
                    HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS')
gene.list <- specific_neuro_pathways[pathways]
ids = intersect(rownames(mat),((unlist(gene.list[c(1,2)]))))
mhc = pheatmap(mat[ids,],
         color = colorRampPalette(c("#1973B1", "white", "#CE1C27"))(100),
         cluster_rows = TRUE,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = FALSE,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         fontsize_col = 6)
topptx(mhc,filename = 'figs/fig4g1.pptx')
ids = intersect(rownames(mat),(((gene.list[['MANALO_HYPOXIA_UP']]))))
HYPOXIA = pheatmap(mat[ids,],
               color = colorRampPalette(c("#1973B1", "white", "#CE1C27"))(100),
               cluster_rows = TRUE,
               cluster_cols = F,
               show_rownames = T,
               show_colnames = FALSE,
               annotation_col = annotation_row,
               annotation_colors = ann_colors,
               fontsize_col = 6)
topptx(HYPOXIA,filename = 'figs/fig4g2.pptx')

ids = intersect(rownames(mat),(((gene.list[['HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS']]))))
ANGIOGENESIS = pheatmap(mat[ids,],
                   color = colorRampPalette(c("#1973B1", "white", "#CE1C27"))(100),
                   cluster_rows = TRUE,
                   cluster_cols = F,
                   show_rownames = T,
                   show_colnames = FALSE,
                   annotation_col = annotation_row,
                   annotation_colors = ann_colors,
                   fontsize_col = 6)
topptx(ANGIOGENESIS,filename = 'figs/fig4g3.pptx')
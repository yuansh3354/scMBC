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
# fig5b
{
  obj.Canc = sce.list[[Canc]]
  x1 = jjDotPlot(object = obj.Canc,
                 gene = ligand,
                 id = 'Cell.Type.L2',
                 xtree=F,ytree=F,
                 dot.col = Strength.cls,
                 point.geom =  F,
                 rescale = T,
                 legend.position='bottom',
                 tile.geom = T)
  ligand_df = x1$data[,c(1,4,2)] %>% set_colnames(c('cell_type','gene', 'expression'))
  
  x2 = jjDotPlot(object = obj.Myeloids,
                 gene = receptor,
                 id = 'Cell.Type.L2',
                 xtree=F,ytree=F,
                 dot.col = Strength.cls,
                 point.geom =  T,
                 rescale = T,
                 legend.position='bottom',
                 tile.geom = T)
  receptor_df = x2$data[,c(1,4,2)] %>% set_colnames(c('cell_type','gene', 'expression')) 
  
  result_reshape_alt <- dcast(x1$data, id ~ gene, value.var = "avg.exp", fill = 0)
  rownames(result_reshape_alt) = result_reshape_alt[,1]
  result_reshape_alt = result_reshape_alt[,-1]
  
  {
    interactions = CellChatDB.human$interaction[CellChatDB.human$interaction$ligand %in% 
                                                  ligand & CellChatDB.human$interaction$receptor %in% 
                                                  receptor,c('ligand','receptor')] %>% unique
    
    ligand_expr <- ligand_df %>%
      select(gene, cell_type, expression) %>%  
      distinct() %>% 
      pivot_wider(
        names_from = cell_type,
        values_from = expression,
        id_cols = gene
      ) %>%
      as.data.frame()
    rownames(ligand_expr) <- ligand_expr$gene
    ligand_expr = ligand_expr[,-1]
    receptor_expr <- receptor_df %>%
      select(gene, cell_type, expression) %>% 
      distinct() %>%  
      pivot_wider(
        names_from = cell_type,
        values_from = expression,
        id_cols = gene
      ) %>%
      as.data.frame()
    rownames(receptor_expr) <- receptor_expr$gene
    receptor_expr = receptor_expr[,-1]
    
    ligand_expr <- t(scale(t(ligand_expr)))
    receptor_expr <- t(scale(t(receptor_expr)))
    
    x = plot_ligand_receptor_interactions(
      ligand_data = ligand_expr,
      receptor_data = receptor_expr,
      interactions = interactions,
      widths = c(3, 3, 3),
      title = "CellChat Ligand-Receptor Interactions"
    )
  }
}

################################################################################
# fig5c
{
  gs$Cell.Type.L2
  sce = sce.total
  
  calculate_cell_proportions <- function(sce, cancer_subtype = "T Cells") {
    metadata <- sce@meta.data
    non_cancer_props <- metadata %>%
      group_by(PatientID, Cell.Type.L1) %>%
      summarise(n_cells = n(), .groups = 'drop') %>%
      group_by(PatientID) %>%
      mutate(proportion = n_cells / sum(n_cells)) %>%
      select('PatientID', 'Cell.Type.L1', 'proportion') %>% 
      pivot_wider(names_from = Cell.Type.L1, values_from = proportion, values_fill = 0) %>% 
      as.data.frame()
    cancer_props <- metadata[metadata$Cell.Type.L1 == 'Cancer Cells',] %>%
      group_by(PatientID, Cell.Type.L2) %>%
      summarise(n_cells = n(), .groups = 'drop') %>%
      group_by(PatientID) %>%
      mutate(proportion = n_cells / sum(n_cells)) %>%
      select('PatientID', 'Cell.Type.L2', 'proportion') %>% 
      pivot_wider(names_from = Cell.Type.L2, values_from = proportion, values_fill = 0) %>% 
      as.data.frame()
    
    rownames(non_cancer_props) <- non_cancer_props$PatientID
    rownames(cancer_props) <- cancer_props$PatientID
    common_ids <- intersect(non_cancer_props$PatientID, cancer_props$PatientID)
    
    non_cancer_props <- non_cancer_props[common_ids,]
    
    if(cancer_subtype %in% colnames(cancer_props)) {
      non_cancer_props[[cancer_subtype]] <- cancer_props[common_ids, cancer_subtype]
    }
    
    return(non_cancer_props)
  }
  df.sce = calculate_cell_proportions(sce)
  df.gs = calculate_cell_proportions(gs)
  df.sce = df.sce[,colnames(df.gs)]
  df.gs = df.gs[,colnames(df.gs)]
  mydf = rbind(df.sce,df.gs)
  
  x2 = ggplot(mydf, aes(Myeloids,Neu_Cells)) + 
    geom_point()+  geom_smooth(method = "lm", se = T) +  
    stat_cor(method = "pearson") + yuansh_theme+ 
    scale_color_manual(values = Gender.Color)
  x2
  
  x1 = ggplot(mydf, aes(`T Cells`,Neu_Cells)) + 
    geom_point()+  geom_smooth(method = "lm", se = T) +  
    stat_cor(method = "pearson") + yuansh_theme+ 
    scale_color_manual(values = Gender.Color)
  x1
}

{
  df.sce = calculate_cell_proportions(sce)
  df.gs = calculate_cell_proportions(gs)
  col.ids = intersect(colnames(df.gs),colnames(df.sce))
  df.sce = df.sce[,col.ids]
  df.gs = df.gs[,col.ids]
  mydf = rbind(df.sce,df.gs)
  
  x1 = ggplot(mydf, aes(Mph_APOE_CD163,Neu_Cells)) + 
    geom_point()+  geom_smooth(method = "lm", se = T) + 
    stat_cor(method = "pearson") + yuansh_theme+ 
    scale_color_manual(values = Gender.Color)
  x1
  x = (x2/x1)
}

{
  df.sce = calculate_cell_proportions(sce)
  mydf = df.sce
  x2 = ggplot(mydf, aes(Treg,Neu_Cells)) + 
    geom_point()+  geom_smooth(method = "lm", se = T) + 
    stat_cor(method = "pearson") + yuansh_theme+ 
    scale_color_manual(values = Gender.Color)
  x2
}

################################################################################
# fig5 f
{
  hd = hd.list[['C']]
  hd$spatial_cluster = as.character(hd$spatial_cluster)
  x = jjDotPlot(
    object = hd,
    gene = 'Cluster1',
    xtree = F,
    ytree = F, lwd = 0.2, bar.width = 3,
    dot.col = c("#0571b0", "#f7f7f7", "#ca0020"),
    id = "spatial_cluster",
    rescale = T, legend.position = "bottom",
    # point.geom = F,
    tile.geom = T
  )
  df.ids = x$data
  df.ids = df.ids[order(df.ids$avg.exp,decreasing = T),]
  df.ids$myid = paste0('R',1:9)
  df.ids = split(df.ids$myid,df.ids$id)
  hd$myid = NA
  for(ids in names(df.ids)){
    hd$myid = ifelse(as.character(hd$spatial_cluster) == ids,
                     df.ids[[ids]],hd$myid)
  }
  hd$Region = ifelse(hd$myid %in% mysplit('R2,R7'),'Region1','Region2')
  
  x = jjDotPlot(
    object = hd,
    gene = 'Cluster1',
    xtree = F,
    ytree = F, lwd = 0.2, bar.width = 3,
    dot.col = c("#0571b0", "#f7f7f7", "#ca0020"),
    id = "myid",
    rescale = T, legend.position = "bottom",
    tile.geom = T
  )
}

################################################################################
# fig5 g
{
  df = visualizeCellTypeDistribution(hd@meta.data,'myid','Cell.Type.L2')
  df = df$data
  df$Percentage =  df$Percentage * 100
  
  x = ggplot(df,aes(myid,Percentage,fill=myid,group =Cell.Type.L2 )) + 
    facet_wrap(~Cell.Type.L2,ncol = 5,scale='free') + 
    geom_bar(stat = 'identity')+
    yuansh_theme +
    scale_fill_npg()
  x
  
  df = hd@meta.data
  df = df[df$Cell.Type.L1 !="Cancer Cells",]
  df = visualizeCellTypeDistribution(df,'myid','Cell.Type.L2')
  df = df$data
  df$Percentage =  df$Percentage * 100
  x = ggplot(df,aes(myid,Percentage,fill=myid,group =Cell.Type.L2 )) + 
    facet_wrap(~Cell.Type.L2,ncol = 5,scale='free') + 
    geom_bar(stat = 'identity')+
    yuansh_theme +
    scale_fill_npg()
  x
  topptx(x,filename = 'figs/fi4.bar-g-2.pptx',width = 9)
  write_to_excel_sheet(
    file_path = "source datafig5.xlsx",
    data = x$data,
    sheet_name = "fig5g-mph",
    rowNames = TRUE
  )
  
}
################################################################################
# fig5i-j

{
  fgsea_results[1:5,c('pathway','NES','padj','pval')]
  fgsea_results = fgsea_results[order(fgsea_results$NES,decreasing = T),]
  df = as.data.frame(fgsea_results)
  df$pathway = gsub('HALLMARK_','',df$pathway)
  immune_pathways <- c(
    "INFLAMMATORY_RESPONSE",
    "IL2_STAT5_SIGNALING",
    "INTERFERON_ALPHA_RESPONSE",
    "COMPLEMENT",
    "INTERFERON_GAMMA_RESPONSE",
    "TNFA_SIGNALING_VIA_NFKB",
    "ALLOGRAFT_REJECTION",
    "IL6_JAK_STAT3_SIGNALING"
  )
  x = ggplot(df[df$pathway %in% immune_pathways,],
             aes(NES, reorder(pathway, NES), fill = -log10(padj))) +
    geom_bar(stat = 'identity') + 
    scale_fill_gradientn(
      colors = c("grey", "grey", "#0571b0", "#ca0020"),
      values = scales::rescale(c(0, 1.3, 1.3, 3.7)),
      breaks = c(0, 1.3, 2, 3),
      limits = c(0, 3.7)
    ) +
    yuansh_theme
  topptx(x,filename = 'figs/barplot.gsea.pptx')
  
}

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
# fig1a
DimPlot(sce.total,group.by = 'Cell.Type.L2') + NoLegend() + 
  scale_color_manual(values = c(Cell.Color.L2,Cell.Color.L1)) +
  NoLegend() +
  labs(title='') +
  myaxi_theme
################################################################################
# fig1b
patients = unique(gs$PatientID)
mclapply(patients,function(patient){
  x = DimPlot(gs[,gs$PatientID ==patient],
              reduction = 'spatial',
              group.by = 'Cell.Type.L1') + 
    yuansh_theme + myaxi_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + NoLegend() +
    scale_color_manual(values = c(Cell.Color.L1))
  ggsave(paste0("figs/fig1.GS.",patient,".UMAP.pdf"),plot = x,width = 5.52,height = 3.24)
})
mclapply(patients,function(patient){
  obj = gs[,gs$PatientID ==patient]
  obj$spatial_l2 = obj$Cell.Type.L2
  x = DimPlot(
    obj, raster = TRUE, pt.size = 1.8,
    reduction = "spatial",
    group.by = "spatial_l2"
  ) +
    yuansh_theme +
    myaxi_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    NoLegend()
  ggsave(paste0("figs/fig1.GS.sp.",patient,".UMAP.pdf"),plot = x,width = 5.52,height = 3.24)
})
################################################################################
# fig1c
mclapply(hd.list,
         function(obj){
           patient = unique(obj$patient_id)
           x = DimPlot(obj,raster = T,pt.size = 1.8,
                       reduction = 'spatial',
                       group.by = 'Cell.Type.L1') +
             yuansh_theme + 
             myaxi_theme +
             theme(
               axis.text.x = element_blank(),  
               axis.text.y = element_blank(),
               plot.title = element_blank(),
               plot.margin = margin(0, 0, 0, 0)
             ) + 
             NoLegend()
           ggsave(paste0("figs/fig1.HD.",patient,".UMAP.pdf"),plot = x,
                  width = 3.5,height = 2.8)
         },mc.cores = 128)

mclapply(hd.list, function(obj) {
  patient = unique(obj$patient_id)
  x = DimPlot(
    obj, raster = TRUE, pt.size = 1.8,
    reduction = "spatial",
    group.by = "spatial_l2"
  ) +
    yuansh_theme +
    myaxi_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    NoLegend() +
    scale_color_manual(values = c("MSTCs" = "#2D5987", "Others" = "grey90"))

}, mc.cores = 128)
mclapply(hd.list, function(obj) {
  patient = unique(obj$patient_id)
  x = DimPlot(
    obj, raster = TRUE, pt.size = 1.8,
    reduction = "spatial",
    group.by = "spatial_l2"
  ) +
    yuansh_theme +
    myaxi_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    NoLegend() +
    scale_color_manual(values = c("MSTCs" = "#2D5987",Mph_APOE='#7F2246', "Others" = "grey90"))

}, mc.cores = 128)
################################################################################
# fig1d
{
  df = visualizeCellTypeDistribution(sce.total@meta.data,'PatientID','Cell.Type.L1') + 
    scale_fill_manual(values = Cell.Color.L1)
  df = df$data
  male = unique(sce.total@meta.data[sce.total$Gender == 'male','PatientID'])
  df$Gender = ifelse(df$PatientID %in% male,'male','fmale')
  df = df[order(df$Gender),]
  ids = df[,c('PatientID','Gender')] %>% unique
  df$PatientID = factor(df$PatientID,levels= (ids$PatientID))
  x = ggplot(df,aes(PatientID,Percentage,fill=Cell.Type.L1))+
    geom_bar(stat = "identity", position = "stack") +
    yuansh_theme +
    scale_fill_manual(values = Cell.Color.L1) + RotatedAxis()  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  x
  
}
################################################################################
# fig1e
df = visualizeCellTypeDistribution(gs@meta.data,'PatientID','Cell.Type.L1') +
  scale_fill_manual(values = c(Cell.Color.L1))

df = visualizeCellTypeDistribution(hd.meta.total,'patient_id','Cell.Type.L1') +
  RotatedAxis() + scale_fill_manual(values = c(Cell.Color.L1,Cell.Color.L2))
################################################################################
# fig1f
{
  df = visualizeCellTypeDistribution(sce.total@meta.data[
    sce.total$Cell.Type.L1 == 'Cancer Cells',],'PatientID','Cell.Type.L2') +
    scale_fill_manual(values = Cell.Color.L1)
  df = df$data
  male = unique(sce.total@meta.data[sce.total$Gender == 'male','PatientID'])
  df$Gender = ifelse(df$PatientID %in% male,'male','fmale')
  df = df[order(df$Gender),]
  df$PatientID = factor(df$PatientID,levels= unique(df$PatientID))
  x = ggplot(df, aes(Cell.Type.L2, Percentage, fill = Gender)) +
    geom_violin(alpha = 0.7, position = position_dodge(width = 0.6)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.6)) +
    geom_jitter(size = 0.5, alpha = 1, color = "black",
                position = position_jitterdodge(dodge.width = 0.2)) +
    yuansh_theme +
    scale_fill_manual(values = Gender.Color) + 
    RotatedAxis() +
    stat_compare_means(method = 't.test') + 
    NoLegend()
  x
}
# fig1g
{
  df = sce.total@meta.data
  df = df[df$Cell.Type.L1 != "Cancer Cells" ,]
  df = visualizeCellTypeDistribution(df,'PatientID','Cell.Type.L1')$data
  male = unique(sce.total@meta.data[sce.total$Gender == 'male','PatientID'])
  df$Gender = ifelse(df$PatientID %in% male,'male','fmale')
  x = ggplot(df[df$Cell.Type.L1 == 'T Cells',], aes(Gender, Percentage, fill = Gender)) + 
    geom_violin(alpha = 0.7, position = position_dodge(width = 0.6)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.6)) +
    geom_jitter(size = 0.5, alpha = 1, color = "black",
                position = position_jitterdodge(dodge.width = 0.2)) +
    yuansh_theme +
    facet_wrap(~Cell.Type.L1, scale = 'free') +
    theme(legend.position = "none") +  
    stat_compare_means(method = 't.test') + 
    scale_fill_manual(values = Gender.Color)
}
{
  obj = sce.list[['Tcells']]
  df = obj@meta.data
  df = visualizeCellTypeDistribution(df,'PatientID','Cell.Type.L2')
  df = df$data
  male = unique(obj@meta.data[obj$Gender == 'male','PatientID'])
  df$Gender = ifelse(df$PatientID %in% male,'male','fmale')
  
  df$Percentage = df$Percentage *100
  x = ggplot(df, aes(Cell.Type.L2, Percentage, fill = Gender)) + 
    geom_violin(alpha = 0.7, position = position_dodge(width = 0.6)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.6)) +
    geom_jitter(size = 0.5, alpha = 1, color = "black",
                position = position_jitterdodge(dodge.width = 0.2)) +
    yuansh_theme +
    theme(legend.position = "none") +  
    stat_compare_means(method = 't.test') + 
    scale_fill_manual(values = Gender.Color) +
    RotatedAxis()
}
{
  obj = sce.list[['Myeloids']]
  df = obj@meta.data
  df = visualizeCellTypeDistribution(df,'PatientID','Cell.Type.L2')
  df = df$data
  male = unique(obj@meta.data[obj$Gender == 'male','PatientID'])
  df$Gender = ifelse(df$PatientID %in% male,'male','fmale')
  
  df$Percentage = df$Percentage *100
  x = ggplot(df, aes(Cell.Type.L2, Percentage, fill = Gender)) + 
    geom_violin(alpha = 0.7, position = position_dodge(width = 0.6)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.6)) +
    geom_jitter(size = 1, position = position_jitterdodge(dodge.width = 0.2)) +
    yuansh_theme +
    theme(legend.position = "none") +  
    stat_compare_means(method = 't.test') + 
    scale_fill_manual(values = Gender.Color) +
    RotatedAxis()
  x
}

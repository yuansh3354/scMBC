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
# fig2a
DimPlot(sce.list[["Canc"]],group.by = "Cell.Type.L2") +
  scale_color_manual(values = Cell.Color.L2)
jjDotPlot(
  object = sce.list[["Canc"]],
  gene = mysplit('CD24,CD44,PROM1,GRIA2,NEGR1,NF1,GFRA1,LSAMP,
                   MKI67,TOP2A,CCNB1,HSPB11,HSPE1,HSPD1'),
  xtree = F,
  ytree = F, lwd = 0.2, bar.width = 3,
  dot.col = c("#0571b0", "#f7f7f7", "#ca0020"),
  id = "Cell.Type.L2",
  rescale = T, legend.position = "bottom",
  # point.geom = F,
  tile.geom = T
)
################################################################################
# fig2b
df = FetchData(sce.list[['Canc']],vars = c('Cluster1','PGR','ERBB2','ESR1','Gender'))
for(ids.gender in c('male','fmale')){
  for(ids.gene in c('PGR','ERBB2','ESR1')){
    df.ids = df[df$Gender == ids.gender, ]
    df.ids$var = df.ids[[ids.gene]]
    x1 = ggplot(df.ids, aes(normalize_to_01(Cluster1),normalize_to_01(var))) +
      geom_point(size = 1) + 
      geom_smooth() + 
      stat_cor() + 
      yuansh_theme  +  NoLegend() +
      labs(title='') +
      myaxi_theme
    x1
    ggsave(plot= x1,filename = paste0('figs/',ids.gender,'-',ids.gene,'.png'),dpi = 450,height = 3,width = 3)
  }
}
################################################################################
# fig2e
data = read.csv('mIHC.HALO.csv')
ggplot(data, aes(Image_Tag, `% of Neuron-like epithelial cells`,color=Image_Tag)) + 
  geom_boxplot() +   geom_jitter()+
  yuansh_theme + scale_fill_simpsons() +
  stat_compare_means(method = 't.test') + 
  scale_color_manual(values = c('male' = '#8577B4',
                                'fmale' = '#B9888A'))

################################################################################
# fig2f
load('tcga.rdata')
x = ggplot(tcga,aes(Gender,score,fill = Gender)) + 
  geom_boxplot() + 
  stat_compare_means(method = 't.test') + 
  theme_bw()
x

load('gtex.rdata')
x = ggplot(gtex,aes(Gender,score,fill = Gender)) + 
  geom_boxplot() + 
  stat_compare_means(method = 't.test') + 
  theme_bw()
x
################################################################################
# fig2g
sfit1 <- survfit(Surv(time, event)~group, data=tcga) 
res1 = ggsurvplot(
  sfit1,
  pval = T,
  conf.int = T,
  palette = c( "#BF3C27","#0271B5") 
)
################################################################################
# fig2h
sfit1 <- survfit(Surv(time, event)~group, data=inhouse) 
res1 = ggsurvplot(
  sfit1,
  pval = T,
  conf.int = T,
  palette = c( "#BF3C27","#0271B5")  
)

################################################################################
{
  my_gsea_hall = qread('gsea_hall.qs',nthreads = 512)
  my_gsea_hall = mclapply(my_gsea_hall,function(x){
    x = intersect(x,rownames(gs))
  })
  gs = AddModuleScore(gs,features = my_gsea_hall)
  gs@meta.data[grep('Cluster',colnames(gs@meta.data))] = names(my_gsea_hall)
  
  plt.list = mclapply(unique(gs$Sample),function(sample.ids){
    p1 = plot_density(gs[,gs$Sample == sample.ids],
                      features = c( "Fatty acid metabolism",'Neurotrophic'),
                      reduction = 'spatial', joint = TRUE,size = 0.5,method = 'wkde')[[3]]
    
    p1$data$patient = sample.ids
    return(p1$data)
  },mc.cores = 12)
  pdata.list = do.call(rbind,plt.list)
  p1 = ggplot(pdata.list,aes(spatial_1,spatial_2,
                             color = normalize_to_01(feature),alpha = normalize_to_01(feature))) +
    geom_point(size = 0.125) +
    yuansh_theme + myaxi_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_gradient2(
      low = "grey50",  # 深蓝
      mid = "#3579A2FF",  # 浅蓝
      high = "#DEF5E5FF", # 黄绿
      midpoint = 0.5      # 中间值位置
    ) + facet_wrap(~patient,ncol = 5,scales = 'free')+
    NoLegend()

}
{
  pdata = mclapply(names(hd.list),
                   function(sample.id){
                     obj = hd.list[[sample.id]]
                     my_gsea_hall = mclapply(my_gsea_hall,function(x){
                       x = intersect(x,rownames(obj))
                     })
                     obj = AddModuleScore(obj,features = my_gsea_hall)
                     obj@meta.data[grep('Cluster',colnames(obj@meta.data))] = names(my_gsea_hall)
                     
                     p2 = plot_density(obj,
                                       features = c("Neurotrophic", "Fatty acid metabolism"),
                                       reduction = 'spatial', joint = TRUE,
                                       size = 0.5,method = 'wkde')[[3]]+
                       yuansh_theme + myaxi_theme +
                       theme(
                         axis.text.x = element_blank(),  # 去掉x轴刻度标签
                         axis.text.y = element_blank()   # 去掉y轴刻度标签
                       ) + facet_wrap(~patient,scale='free')+
                       scale_color_viridis(option = "mako") + NoLegend()
                     p2$data$patient = sample.id
                     
                     df = p2$data
                     p1 = ggplot(df,aes(spatial_1,spatial_2,
                                        color = normalize_to_01(feature),
                                        alpha = normalize_to_01(feature))) +
                       geom_point(size = 0.125) +
                       yuansh_theme + myaxi_theme +
                       theme(
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank()
                       ) +
                       scale_color_gradient2(
                         low = "grey50",  # 深蓝
                         mid = "#3579A2FF",  # 浅蓝
                         high = "#DEF5E5FF", # 黄绿
                         midpoint = 0.5      # 中间值位置
                       ) +
                       NoLegend()
                     p1
                     return(p1$data)
                   },
                   mc.cores = 12)
  }
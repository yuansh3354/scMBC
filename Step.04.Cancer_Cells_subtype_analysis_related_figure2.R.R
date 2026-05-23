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
  gene = my.use.genes,
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
ggplot(data, aes(Image_Tag, `% of MSTCs`,color=Image_Tag)) + 
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
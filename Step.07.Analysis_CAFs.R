setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]

########################## Step.07 CAFs ##########################
caf = readRDS('result/RDS/03.1.fibroblast.RDS')
use.cells = caf@meta.data[which(caf$Type !='Lymph'),] %>% rownames() 
caf = subset(caf, cells = use.cells)

########################## Step.07.1 UMAP ##########################
sce = caf
DimPlot(sce, group.by = 'Subtype',raster = F,order = T) + NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + 
  ggtitle('') + 
  scale_color_manual(values = cell.cls)  
ggsave('plts/umap-caf.pdf',height = 4,width = 4)

df = sce@meta.data
count_df <- df %>%
  group_by(Subtype,Type) %>%
  summarise(Total = n())

count_df <- count_df %>%
  group_by(Subtype) %>%
  mutate(Percentage = Total / sum(Total))

x = ggplot(count_df, aes(x =Subtype , y = Percentage, fill =Type )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "", fill = "", title = "") +
  theme_minimal() + scale_fill_manual(values = type.cls)
x

count_df <- df %>%
  group_by(Subtype) %>%
  summarise(Total = n()) %>%
  mutate(Percentage = Total / sum(Total) * 100)

x = ggplot(count_df, aes(x = reorder(Subtype,Total),fill=Subtype)) +
  geom_bar(aes(y = Total), stat = "identity") +
  geom_text(aes(y = Total, label = Total), hjust =1) +
  scale_y_continuous(sec.axis = sec_axis(~./sum(count_df$Total) * 100, name = "Percentage")) +
  labs(y = "Count", x = "", title = "") +
  theme_minimal() + coord_flip() + scale_fill_manual(values = cell.cls)
x

CellMarkers = list(
  'Vascular CAFs' = c('MCAM','MYH11','ACTA2'), 
  'Inflammatory CAFs' = c('CXCL12','CXCL14','IGF1','IGFBP6'),
  'EMT-like CAFs' = c('KRT8','KRT18','KRT19'),
  'Antigene CAFs' = c('CD74','HLA-DRA','HLA-DRB1','HLA-B','HLA-E')
)

x = jjDotPlot(object = caf,
              gene = unlist(CellMarkers),
              xtree = F,
              ytree = F,
              dot.col = c("#A6CEE3",'white',"#c10534"),
              id = 'Subtype',
              rescale = T,
              legend.position='bottom',
              # point.geom = F,
              tile.geom = T) + coord_flip()
x

########################## Step.07.2 scrabble for Meta module ##########################
sce = caf
degs = FindAllMarkers(sce, only.pos = T, logfc.threshold = 1,min.pct = 0.5)
sce = AddModuleScore(sce,features = split(degs$gene,degs$cluster))
df = scrabble::hierarchy(sce@meta.data[,c("Cluster1","Cluster2","Cluster3","Cluster4")])
df$Subtype = sce$Subtype
df$Type = sce$Type
groups = rownames(df)
names(groups) = df$Subtype %>% as.character()
scrabble::plot_hierarchy(df,legend.pos = "top",
                         quadrant.names = unique(sce$Subtype),
                         groups=groups,
                         group.cols = cell.cls,
                         legend = T)

ggplot(df, aes(x = X, y = Y, color = Subtype)) +
  geom_point() +
  scale_color_manual(values = cell.cls) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") + labs(x='',y='')

########################## Step.07.2 ssGSVA ##########################
res.list = list()
library(splatter)
library(Seurat)
library(msigdbr)
library(singleseqgset)
library(heatmap3)

for(cell in unique(sce$Subtype)){
  use.cells = rownames(sce@meta.data[which(sce$Subtype == cell),])
  temp = subset(sce, cells = use.cells)
  Idents(temp) = temp$Type

  h.human <- msigdbr(species="Homo sapiens",category="H")
  h.names <- unique(h.human$gs_name)
  h.sets <- vector("list",length=length(h.names))
  names(h.sets) <- h.names
  
  for (i in names(h.sets)) {
    h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
  }
  
  logfc.data <- logFC(cluster.ids=temp@meta.data$Type,
                      expr.mat=temp@assays$RNA$data)
  gse.res <- wmw_gsea(expr.mat=temp@assays$RNA$data,
                      cluster.cells=logfc.data[[1]],
                      log.fc.cluster=logfc.data[[2]],
                      gene.sets=h.sets)

  saveRDS(gse.res,'result/RDS/gse.epi.rds')
  res.stats <- gse.res[["GSEA_statistics"]]
  res.pvals <- gse.res[["GSEA_p_values"]]
  res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
  res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
  res.pvals[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets by p values
  rownames(res.stats) = gsub('HALLMARK_','',rownames(res.stats))
  rownames(res.stats) = gsub('_',' ',rownames(res.stats))
  
  res = list(res.stats = res.stats,res.pvals = res.pvals)
  res.list[[cell]] = res
}


gsva.list = list()
for (variable in names(res.list)) {
  res.stats = res.list[[variable]][[1]]
  top_local <- res.stats %>%
    arrange(desc(Local)) %>%
    slice_head(n = 5)
  
  top_near <- res.stats %>%
    arrange(desc(Near)) %>%
    slice_head(n = 5)
  
  
  combined_df <- bind_rows(top_local, top_near) %>%
    distinct() 
  colnames(combined_df) = paste0(colnames(combined_df),'(',variable,')')
  gsva.list[[variable]] = combined_df
}
gsva.list <- lapply(gsva.list, function(df) {
  df <- tibble::rownames_to_column(df, var = "RowName")
  return(df)
})
merged_df <- Reduce(function(x, y) {
  full_join(x, y, by = "RowName")
}, gsva.list)

rownames(merged_df) <- merged_df$RowName
merged_df$RowName <- NULL

print(merged_df)

x = pheatmap::pheatmap(merged_df,
                       scale = 'row',
                       cluster_rows = F,
                       cluster_cols = F,
                       na_col = 'grey',
                       cellwidth = 8,cellheight = 8,
                       color =colorRampPalette(c("#A6CEE3", "white", "firebrick3"))(64) )
x
topptx(x,filename = paste0('plts/CAFs-ssGSVA.pptx'),width = 12,height = 12)


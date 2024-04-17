setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]

########################## Step.06 Macrophage ##########################
mph = readRDS('result/RDS/03.1.Macrophage.RDS')

########################## Step.06.1 UMAP ##########################
sce = mph
DimPlot(sce, group.by = 'Subtype',raster = F,order = T) + NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + 
  ggtitle('') + 
  scale_color_manual(values = mph.colors)  
ggsave('plts/umap-mph.pdf',height = 4,width = 4)

df = sce@meta.data
count_df <- df %>%
  group_by(Type, Subtype) %>%
  summarise(Total = n())

count_df <- count_df %>%
  group_by(Type) %>%
  mutate(Percentage = Total / sum(Total))

x = ggplot(count_df, aes(x = Type, y = Percentage, fill = Subtype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "", fill = "", title = "") +
  theme_minimal() + scale_fill_manual(values = mph.colors)
x
topptx(x,'plts/mph-type.pptx')
source("script/distribution_Roe.R")

metainfo = FetchData(sce, vars = c('CellID','Type','Subtype'))
distribution_Roe(
  meta_data = metainfo,
  celltype_column = "Subtype",
  celltype_level = names(mph.colors) %>% rev(),
  condition_column = "Type",
  add_label = "sign",
  celltype_color = mph.colors,relative_width = 0.7,
  tile_color = NA
)





ggsave2('plts/roe.pdf',height = 6,width = 4)


count_df <- df %>%
  group_by(Subtype) %>%
  summarise(Total = n()) %>%
  mutate(Percentage = Total / sum(Total) * 100)

x = ggplot(count_df, aes(x = reorder(Subtype,Total),fill=Subtype)) +
  geom_bar(aes(y = Total), stat = "identity") +
  geom_text(aes(y = Total, label = Total), hjust =1) +
  scale_y_continuous(sec.axis = sec_axis(~./sum(count_df$Total) * 100, name = "Percentage")) +
  labs(y = "Count", x = "", title = "") +
  theme_minimal() + coord_flip() + scale_fill_manual(values = mph.colors)
x

degs = FindAllMarkers(sce, only.pos = T, logfc.threshold = 1,min.pct = 0.75)
CellMarkers = c('STMN1','TUBB','HMGB1',
                'SPP1','APOC1','APOE',
                'S100A9','LYZ','NAMPT',
                'MERTK','FRMD4B','LRMDA',
                'CXCL9','BTG1','ACTG1',
                'CCL4','CST3')

x = jjDotPlot(object = sce,
              gene = CellMarkers,
              xtree = F,
              ytree = F,
              dot.col = strength.cls,
              id = 'Subtype',
              rescale = T,
              legend.position='bottom',
              # point.geom = F,
              tile.geom = T)
########################## Step.06.2 Mph ssGSEA ##########################
sce = mph
if(F){
  deg.list = list()
  for(cell in unique(sce$Subtype)){
    use.cells = rownames(sce@meta.data[which(sce$Subtype == cell),])
    temp = subset(sce, cells = use.cells)
    Idents(temp) = temp$Type
    degs1 = FindMarkers(temp,ident.1 = 'Local',ident.2 = 'Near',min.pct = 0.5, logfc.threshold = 1,only.pos = T) %>% 
      tibble::rownames_to_column(var = "Gene") %>% 
      mutate(Group = "LocalvsNear")
    degs2 = FindMarkers(temp,ident.1 = 'Lymph',ident.2 = 'Near',min.pct = 0.5, logfc.threshold = 1,only.pos = T) %>% 
      tibble::rownames_to_column(var = "Gene") %>% 
      mutate(Group = "LymphvsNear")
    
    degs = rbind(degs1,degs2)
    degs = degs[!grepl('^LINC',degs$Gene),]
    degs = degs[!grepl('^MT-',degs$Gene),]
    degs = degs[!grepl('-',degs$Gene),]
    degs = degs[!grepl('\\.',degs$Gene),]
    degs = degs[!grepl('^RP[SL][[:digit:]]',degs$Gene),]
    deg.list[[cell]] = degs
  }
  library(ComplexHeatmap)
  genes =lapply(deg.list, function(deg){
    return(deg$Gene)
  })
}

library(splatter)
library(Seurat)
library(msigdbr)
library(singleseqgset)
library(heatmap3)

res.list = list()
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
  
  top_lymph <- res.stats %>%
    arrange(desc(Lymph)) %>%
    slice_head(n = 5)
  combined_df <- bind_rows(top_local, top_near, top_lymph) %>%
    distinct() 
  colnames(combined_df) = paste0(colnames(combined_df),'(',variable,')')
  gsva.list[[variable]] = combined_df
}

gsva.list <- lapply(gsva.list, function(df) {
  df <- tibble::rownames_to_column(df, var = "RowName")
  return(df)
})


x = pheatmap::pheatmap(merged_df,
                       scale = 'row',
                       cluster_rows = F,
                       cluster_cols = F,
                       na_col = 'grey',
                       cellwidth = 8,cellheight = 8,
                       color =colorRampPalette(strength.cls)(64) )



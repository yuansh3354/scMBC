########################## Step.04 T Cell ##########################
Tcell = readRDS('result/RDS/03.1.Tcell.RDS')
########################## Step.04.1 UMAP ##########################
sce = Tcell
DimPlot(sce, group.by = 'Subtype',raster = F,order = T) + NoLegend()+
  theme(
    axis.title = element_blank(),   
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + 
  ggtitle('') + 
  scale_color_manual(values = cell.cls) + stat_unchull(alpha = 0.5, color = "grey50", size = 0.5, lty = 2)

df = Tcell@meta.data
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
  theme_minimal() + scale_fill_manual(values = cell.cls)
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

metainfo = FetchData(sce, vars = c('CellID','Type','Subtype'))
distribution_Roe(
  meta_data = metainfo,
  celltype_column = "Subtype",
  celltype_level = names(cell.cls) %>% rev(),
  condition_column = "Type",
  add_label = "sign",
  celltype_color = cell.cls,relative_width = 0.7,
  tile_color = NA
)

########################## Step.04.2 DEGs ##########################
sce = Tcell
ids.type = unique(sce$Type)
sce = subset(sce, ident = ids.cells)
ids.genes = trimws(strsplit("PDCD1, CTLA4, LAG3, TIM3, TIGIT, IDO1 HLA-G, VISTA",',')[[1]])
sce@meta.data = sce@meta.data[,-c(grep('Cluster',colnames(sce@meta.data)))]
Idents(sce) = sce$Type
degs1 = FindMarkers(sce,ident.1 = 'Local',min.pct = 0.5, only.pos = T)
degs2 = FindMarkers(sce,ident.1 = 'Lymph',min.pct = 0.5, only.pos = T)
degs3 = FindMarkers(sce,ident.1 = 'Near',min.pct = 0.5, only.pos = T)

degs = rbind(degs1,degs2,degs3)
degs$gene = rownames(degs)
degs = degs[!grepl('LINC',degs$gene),]
degs = degs[!grepl('\\.',degs$gene),]
degs = degs[!grepl('^RP[SL][[:digit:]]',degs$gene),]
Subtype = "NK/NKT"
x = jjDotPlot(object = sce[,which(sce$Subtype == Subtype )],
              gene = degs$gene,
              xtree = F,
              ytree = F,
              dot.col = c("#A6CEE3",'white',"#c10534"),
              id = 'Type',
              rescale = T,legend.position='bottom',
              # point.geom = F,
              tile.geom = T)
genes3 = NULL
df = x$data
for(ids in c('Local','Near','Lymph')){
  temp = df[which(df$id == ids),]
  temp = df[which(df$gene %in% as.character(temp[which(temp$avg.exp.scaled==1),'gene'])),]
  temp =temp %>%
    group_by(gene) %>%
    summarise(total_avg_exp_scaled = sum(avg.exp.scaled, na.rm = TRUE)) %>%
    arrange(total_avg_exp_scaled) %>%
    slice_head(n = 10)
  genes3 =  unique(c(as.character(temp$gene),genes3))
}


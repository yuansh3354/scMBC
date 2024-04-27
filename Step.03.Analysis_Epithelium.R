########################## Step.03.1 UMAP & Markers ##########################
epi = readRDS('result/RDS/03.1.Epithelium.RDS')
sce = epi

# UMAP
DimPlot(sce, group.by = 'cancers',
        cols = epi.colors,
        raster = F,order = T) + NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
# ggsave('plts/umap-epi.pdf',height = 4,width = 4)

# Percentage Bar plot
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

# Counts Bar plot
df = sce@meta.data
count_df <- df %>%
  group_by(Subtype) %>%
  summarise(Total = n()) %>%
  mutate(Percentage = Total / sum(Total) * 100)
x = ggplot(count_df, aes(x = reorder(Subtype,Total),fill=Subtype)) +
  geom_bar(aes(y = Total), stat = "identity") +
  geom_text(aes(y = Total, label = Total), hjust =1) +
  scale_y_continuous(sec.axis = sec_axis(~./sum(count_df$Total) * 100, name = "Percentage")) +
  labs(y = "Count", x = "", title = "") +
  theme_minimal() + coord_flip() + scale_fill_manual(values = c(
    'EMT-like'='#F1787D',
    'Immune-like'='#F8D889',
    'Invasion-like'='#BD82BD',
    'Metastasis-like'='#EDA462',
    'Prolification-like'='#F6C4E6', 
    'Neuroendocrine-like'='#5EB7F1'))
x

# Feature plots 
x = jjDotPlot(object = epi,
              gene = unlist(CellMarkers),
              xtree = F,
              ytree = F,lwd=0.2,bar.width	=3,
              dot.col = c("#A6CEE3",'white',"#c10534"),
              id = 'cancers',
              rescale = T,legend.position='bottom',
              # point.geom = F,
              tile.geom = T)
x

########################## Step.03.2 dyno ##########################
library(dyno)
library(tidyverse)
dataset <- wrap_expression(
  expression = t(sce@assays$RNA$data),
  counts = t(sce@assays$RNA$counts)
)

guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected

# Set Prolifications as start cells, cause of no specific markers
dataset <- dataset %>% add_prior_information(start_id =
                                               rownames(sce@meta.data[which(sce$Subtype == 'Prolification-like'),])
)

# select paga_tree 
model <- infer_trajectory(dataset, 'paga_tree')
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
dataset <- add_dimred(
  dataset,
  sce@reductions$tsne.rpca@cell.embeddings
)
saveRDS(model,file = 'result/RDS/epi-dyno-model.RDS')
saveRDS(dataset,file = 'result/RDS/epi-dyno-dataset.RDS')

plot_dendro(model,  size_cells = 1,
            grouping =sce$Subtype) + 
  scale_color_manual(values =epi.colors ) + NoLegend()

########################## Step.03.3 Cell Distribution ##########################
sce = epi

meta.data = sce@meta.data
celltype_counts <- meta.data %>%
  group_by(Type, cancers) %>%
  summarise(count = n())
total_counts <- celltype_counts %>%
  group_by(Type) %>%
  summarise(total = sum(count))
celltype_counts <- celltype_counts %>%
  left_join(total_counts, by = "Type") %>%
  mutate(percentage = count / total * 100)
celltype_counts = celltype_counts[,c('Type','Subtype','percentage')]
celltype_counts$Subtype = as.character(celltype_counts$Subtype)
library(ggradar)
df = reshape2::dcast(celltype_counts, Type ~ Subtype, value.var = "percentage")
x1=ggradar(df,values.radar = c("0%", "25%", "50%"),
           grid.min = 0, 
           grid.mid = 25, 
           grid.max = 50,
           group.colours = type.cls,
           group.point.size = 1.5,
           group.line.width = 0.75, #
           background.circle.transparency = 0, 
           legend.position = 'bottom', 
           legend.text.size = 12, 
           fill = T,
)
x1

# Re/o
metainfo = FetchData(sce, vars = c('CellID','Type','Subtype'))
distribution_Roe(
  meta_data = metainfo,
  celltype_column = "Subtype",
  celltype_level = names(epi.colors) %>% rev(),
  condition_column = "Type",
  add_label = "sign",
  celltype_color = epi.colors,
  relative_width = 0.7,
  tile_color = NA
)
ggsave2('plts/roe.pdf',height = 6,width = 4)

for (cell.type in unique(Idents(sce))) {
  temp = subset(sce,idents=cell.type)
  temp$ShannonEntropy = estimateShannonEntropy(temp)
  df = FetchData(temp,c('Type','ShannonEntropy'))
  df$Type = factor(df$Type,levels = c("Local","Near","Lymph"))
  x = ggplot(df, aes(x = Type, y = ShannonEntropy, fill = Type)) +
    geom_boxplot(alpha=0.75) +
    theme_bw() +
    scale_fill_manual(values = type.cls) +
    stat_compare_means(aes(group = Type),
                       comparisons = list(c("Local", "Near"),c("Local", "Lymph")), label = "p") +
    theme(legend.position = "none")
}

########################## Step.03.4 CellCycleScoring ##########################
sce = epi
hg.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

# Select Genes
G1.first <- unique(annoGene(hg.pairs$G1$first, "ENSEMBL", 'human')$SYMBOL)
G1.second <- unique(annoGene(hg.pairs$G1$second, "ENSEMBL", 'human')$SYMBOL)
S.first <- unique(annoGene(hg.pairs$S$first, "ENSEMBL", 'human')$SYMBOL)
S.second <- unique(annoGene(hg.pairs$S$second, "ENSEMBL", 'human')$SYMBOL)
G2M.first <- unique(annoGene(hg.pairs$G2M$first, "ENSEMBL", 'human')$SYMBOL)
G2M.second <- unique(annoGene(hg.pairs$G2M$second, "ENSEMBL", 'human')$SYMBOL)

s.genes = setdiff(S.first,G1.first) %>% setdiff(.,G2M.first)
s.genes <- c(cc.genes$s.genes,s.genes) %>% unique
g2m.genes = setdiff(G2M.first,G1.first) %>% setdiff(.,S.first)
g2m.genes <- c(cc.genes$g2m.genes,g2m.genes) %>% unique

sce <- CellCycleScoring(sce, 
                        s.features = s.genes,
                        g2m.features = g2m.genes, 
                        set.ident = TRUE)

meta.data = sce@meta.data
celltype_counts <- meta.data %>%
  group_by(Subtype, Phase) %>%
  summarise(count = n())

total_counts <- celltype_counts %>%
  group_by(Subtype) %>%
  summarise(total = sum(count))

# 计算每个Type分组下每个celltype的比例
celltype_counts <- celltype_counts %>%
  left_join(total_counts, by = "Subtype") %>%
  mutate(percentage = count / total * 100)

celltype_counts$Phase = ifelse(celltype_counts$Phase =='G1','G1','S/G2M')
celltype_counts$Phase = factor(celltype_counts$Phase, levels = c('S/G2M','G1'))
x = ggplot(celltype_counts, aes(x = "", y = percentage, fill = Phase, label = paste0(round(percentage), "%"))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~Subtype) +
  labs(title = "", fill = "Celltype", y = "Percentage") +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_npg()

########################## Step.03.5 GenoType ##########################
# SComatic
if(T){
  file.list = dir(pattern = '_cellranger710')
  library(stringr)
  for (variable in file.list) {
    file.path = paste0(variable,'/Step4_VariantCalling/snpEff.anntag.vcf')
    mysnpEFF <- read.table(file.path,comment.char = '#',sep = "\t",header = T)
    mysnpEFF$X =str_split(mysnpEFF$ANN,'ANN=',simplify = T)[,2]
    # myann = str_split(mysnpEFF$X,'\\|',simplify = T)
    df <- data.frame(matrix(nrow = length(data), ncol = 21))
    df[1,] <- c('CHROM','Start','End','REF','ALT','FILTER','Cell_types',
                "Info", "Type", "Impact", "Gene", "Gene_ID", 
                "Feature_type", "Feature_ID", "Biotype", "Rank",
                "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "AA_pos")
    write.table(df, paste0('geneMut/',variable,'_geneMut.tsv'), append = TRUE, sep = "\t", row.names = FALSE, col.names = F)
    mysnpEFF = mysnpEFF[which(mysnpEFF$ALT %in% c("T","A","C","G" )),]
    for (i in 1:length(mysnpEFF$X)) {
      #if(i==1 | (i %% 1000 ==0)){    print(paste0(length(mysnpEFF$X),': ',i))}
      ids.row = mysnpEFF[i,1:7]
      ids = unlist(strsplit(mysnpEFF$X[i], ","))
      df = NULL
      for (id in ids) {
        id = unlist(strsplit(id, ";"))[1]
        info_parts <- unlist(strsplit(id, "\\|"))
        if(any(grepl("\\(",info_parts)) & length(info_parts)>15){print(i)}
        
        info_parts = c(as.character(ids.row),info_parts[1:14])
        
        df = rbind(df,info_parts)
      }
      write.table(df, paste0('geneMut/',variable,'_geneMut.tsv'), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
  for(variable in dir('result/OUT_SComatic/geneMut',full.names =T)){
    mut.ids = read.table(variable,header = T)
    ids = unlist(str_split(variable,'/'))[4]
    ids = gsub('_geneMut.tsv','',ids)
    SampleID = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',ids)
    mut.ids$SampleID = SampleID
    print(SampleID)
    assign(SampleID,mut.ids)
  }
  SComatic = SComatic[which(SComatic$Biotype =='protein_coding'),]
  SComatic = SComatic[which(!(SComatic$Type %in% c('synonymous_variant','splice_region_variant&synonymous_variant'))),]
  saveRDS(SComatic,file = 'result/RDS/scomatic.rds')
}
SComatic = readRDS('result/RDS/scomatic.rds')
# Remove lymph cells
tissue = sce@meta.data[which(sce$Type !='Lymph'),'SampleID'] %>% unique
df = SComatic
df = df[which(df$SampleID %in% tissue),]
result <- df %>%
  group_by(SampleID, Cell_types) %>%
  summarise(Count = log(n(),2))
result
result$Cell_types = factor(result$Cell_types,levels = c("EMT-like","Immune-like","Prolification-like",
                                                        "Invasion-like","Metastasis-like","Neuroendocrine-like"))

x = ggplot(result, aes(x = Count, y =Cell_types , fill = Cell_types)) +
  geom_boxplot() + scale_fill_npg(alpha = 0.8)+
  labs(x = NULL, y = "Log2(Mut budern)", fill = "Sample Type") + 
  scale_fill_manual(values = colors) + 
  theme_bw() 
x

result <- df %>%
  group_by(Cell_types,Gene) %>%
  summarise(Count = log(n(),2))
result = result[which(result$Gene %in% unlist(CellMarkers)),]
genes = names(which(table(result$Gene)>4))
result = result[which(result$Gene %in% genes),]
result$Cell_types = factor(result$Cell_types,levels = c("EMT-like","Immune-like","Prolification-like",
                                                        "Invasion-like","Metastasis-like","Neuroendocrine-like"))
x = ggplot(result, aes(x = Cell_types, y = Count, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  labs(x = "Cell_types", y = "Count", title = "Count of Genes in Different Cell_types") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_aaas()
x
topptx(x,'plts/gene-Mut_burdern.pptx')

df = SComatic
# df = df[which(df$SampleID %in% tissue),]
result <- df %>%
  group_by(SampleID, Cell_types,Gene) %>%
  summarise(Count = log(n(),2))
result

result <- df %>%
  group_by(Cell_types,Gene) %>%
  summarise(Count = log(n(),2))
result
result = result[which(result$Gene %in% c('TP53','UTY','NRS',
                                          'CTNNB1','PTEN','ATM','CDH1',
                                          'GATA3','MAP3K1','HERC2','KRAS','ERBB2')),]
ggplot(result, aes(Gene, Count,fill=Cell_types)) + 
  geom_bar(stat = 'identity',width = 0.75,position = 'dodge') + 
  # geom_line(aes(group = Cell_types), color = 'blue') +  # 添加折线
  # geom_point(aes(group = Cell_types), color = 'blue', size = 2) +
  theme_bw() + labs(x='',y='') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),  
        legend.position = "none"  
  ) + 
  scale_fill_manual(values = epi.colors)

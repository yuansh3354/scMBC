setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]
########################## Step.04 Epithelium ##########################
epi = readRDS('result/RDS/03.1.Epithelium.RDS')

########################## Step.04.1 UMAP & Markers ##########################
sce = epi

# UMAP
DimPlot(sce, group.by = 'cancers',
        cols = epi.colors,
        raster = F,order = T) + NoLegend()+
  theme(
    axis.title = element_blank(),  #轴标题
    axis.text = element_blank(), # 文本
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
CellMarkers = list(
  'Prolification-like' =  c(unlist(strsplit('CD24,TUBB3,TCIM',','))),
  'Neuroendocrine-like' = c(unlist(strsplit('GPHN,GFRA1,NEAT1',','))),
  "Metastasis-like" = c(unlist(strsplit('RNF181,MDK,BST2',','))),
  'Invasion-like' = c('MKI67','TOP2A','NDC80'),
  'Immune-like' = c(unlist(strsplit('CD74,PTPRC,PARP14',','))),
  'EMT-like' = c(unlist(strsplit('VIM,FN1,MMP2',','))),
  c(unlist(strsplit(paste(grep('^ASCL',rownames(sce),value = T),collapse = ','),',')))
)
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

########################## Step.04.2 dyno ##########################
library(dyno)
library(tidyverse)

# Select Local samples and sampling 20000 cells for calculated pseudo-time
use.cells = sce@meta.data[which(sce$Type =='Local'),] %>% rownames() 
sce = subset(sce, cells = use.cells)
Idents(sce) = sce$Subtype

use.cells = sample(Cells(sce),20000)
sce = subset(sce, cells = use.cells)

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

########################## Step.04.3 Cell Distribution ##########################
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

# ShannonEntropy
estimateShannonEntropy <- function(seurat_object){
  seurat_object <- FindNeighbors(seurat_object, k.param=30, return.neighbor=T)
  neighbors <- seurat_object@neighbors$RNA.nn@nn.idx
  cell_names <- seurat_object@neighbors$RNA.nn@cell.names
  sample <- seurat_object$Type
  entropy <- apply(neighbors, 1, function(x){
    data <- table(sample[x])
    if(length(unique(sample[x]))==1){
      percent <- 1
    }else{
      percent <- prop.table(data)
    }
    cell_entropy <- -1*sum(percent*log2(percent))
  })
  return(entropy)
}

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

########################## Step.04.4 Markers ##########################
########################## Step.04.5 GSEA ##########################
sce = epi
degs.std = FindAllMarkers(sce)
lapply(seq_along(unique(degs.std$cluster)), function(x){
  tmp <- degs.std |> 
    dplyr::filter(cluster == unique(degs.std$cluster)[x]) |>
    dplyr::arrange(desc(avg_log2FC))
  
  # gene list
  glist <- tmp$avg_log2FC
  names(glist) <- tmp$gene
  
  # enrichment
  ego <- clusterProfiler::gseGO(geneList = glist,
                                ont = "BP",
                                OrgDb = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                pvalueCutoff = 0.05)
  
  # to data.frame
  df <- data.frame(ego) |> 
    dplyr::arrange(pvalue,desc(NES))
  # df = df[which(df$NES > 0),]
  
  
  return(ego)
}) -> GSEA.list
names(GSEA.list) = unique(degs.std$cluster)

# Get Module for each Subtypes
C1 = c('mesenchyme development','angiogenesis')
C2 = c('negative regulation of immune response','negative regulation of T cell activation')
C3 = c("chromosome segregation","nuclear chromosome segregation")
C4 = c("positive regulation of cell migration","positive regulation of cell motility")
C5 = c('regulation of circadian rhythm','synaptic vesicle exocytosis')
C6 = c("organelle fission","nuclear division")

lapply(seq_along(names(GSEA.list)), function(x){
  
  # to data.frame
  df <- data.frame(GSEA.list[[x]]) |> 
    dplyr::arrange(pvalue,desc(NES))
  # df = df[which(df$NES > 0),]
  
  # plot
  lapply(c(1:2), function(y){
    GseaVis::gseaNb(object = GSEA.list[[x]],
                    subPlot = 2,
                    geneSetID = df[which(df$Description %in% get(paste0('C',x))),]$ID[y],
                    addPval = T,
                    pvalX = 0.9,
                    pvalY = 0.8)
  }) -> plist
  
  # combine
  pplist <- cowplot::plot_grid(plotlist = plist,nrow = 1,align = 'hv')
  
  return(pplist)
}) -> gglist

x1 = gglist[[1]]
x2 = gglist[[2]]
x3 = gglist[[3]]
x4 = gglist[[4]]
x5 = gglist[[5]]
x6 = gglist[[6]]

topptx(x1,'plts/EMT-GSEA.pptx',width = 8,height = 4)
topptx(x2,'plts/Immune-GSEA.pptx',width = 8,height = 4)
topptx(x3,'plts/Invasion-GSEA.pptx',width = 8,height = 4)
topptx(x4,'plts/Metastasis-GSEA.pptx',width = 8,height = 4)
topptx(x5,'plts/Neuroendocrine-GSEA.pptx',width = 8,height = 4)
topptx(x6,'plts/Prolification-GSEA.pptx',width = 8,height = 4)

features = c(C1,C2,C3,C4,C5,C6)
GSEA.list = lapply(GSEA.list, as.data.frame)

df1 = GSEA.list$`EMT-like`[which(GSEA.list$`EMT-like`$Description%in%features),c('Description','NES')]
df2 = GSEA.list$`Immune-like`[which(GSEA.list$`Immune-like`$Description%in%features),c('Description','NES')]
df3 = GSEA.list$`Invasion-like`[which(GSEA.list$`Invasion-like`$Description%in%features),c('Description','NES')]
df4 = GSEA.list$`Metastasis-like`[which(GSEA.list$`Metastasis-like`$Description%in%features),c('Description','NES')]
df5 = GSEA.list$`Neuroendocrine-like`[which(GSEA.list$`Neuroendocrine-like`$Description%in%features),c('Description','NES')]
df6 = GSEA.list$`Prolification-like`[which(GSEA.list$`Prolification-like`$Description%in%features),c('Description','NES')]
colnames(df1) = c('Description','df1') 
colnames(df2) = c('Description','df2')
colnames(df3) = c('Description','df3')
colnames(df4) = c('Description','df4')
colnames(df5) = c('Description','df5')
colnames(df6) = c('Description','df6')
df = join_all(list(df1,df2,df3,df4,df5,df6),by = 'Description',type = 'full')
rownames(df) = df[,1]
df = df[,-1]
df = df[features,]
colnames(df) = names(GSEA.list)
# df[is.na(df)] = -5
x = pheatmap::pheatmap(df, na_col = 'grey',
                       cluster_rows = F,
                       cluster_cols = F,
                       number_format = '%.1f',number_color = 'black',
                       cellwidth = 18,cellheight = 18,display_numbers = T,
                       color =colorRampPalette(c("navy", "white", "firebrick3"))(64) )

########################## Step.04.6 CellCycleScoring ##########################
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

########################## Step.04.7 GenoType ##########################
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
  
  SComatic = rbind(`CA-10-1`,`CA-10-2`,`CA-11-1`,`CA-11-2`,`CA1225-1`,`CA1225-2`,
                   `LN-10`,`LN-3-0913`,`M301_5T`,`MBDJ_2T`,`MLBA_1P`,`MLBA_1T`,
                   `MLYJ_4T`,`MLYZ_3L`,`MLYZ_3T`,`mQD2-P`,`mQD2-T`,`N-3-0913`,
                   `P-10`,`T-3-0913`) %>% as.data.frame()
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

CellMarkers = list(
  'Prolification-like' =  c(unlist(strsplit('CD24,TUBB3,TCIM',','))),
  'Neuroendocrine-like' = c(unlist(strsplit('GPHN,GFRA1,YAP1',','))),
  "Metastasis-like" = c(unlist(strsplit('RNF181,MDK,BST2',','))),
  'Invasion-like' = c('MKI67','TOP2A','NDC80'),
  'Immune-like' = c(unlist(strsplit('CD74,PTPRC,PARP14',','))),
  'EMT-like' = c(unlist(strsplit('VIM,FN1,MMP2',',')))
)
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
# 使用 ggplot 绘制柱状图，并添加折线、点和红色标记
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





# Numbat
if(F){
  for(ids in dir('result/OUT_Numbat',full.names = T)){
    num_files <- length(list.files(ids))
    if(num_files > 50){
      print(str_split(ids,'/',simplify = T)[3])
      temp =  Numbat$new(out_dir = ids)
      temp$clone_post$cell = paste0(sub('_cellranger710','',str_split(ids,'/',simplify = T)[3]),'_',
                                    temp$clone_post$cell)
      temp$joint_post$cell = paste0(sub('_cellranger710','',str_split(ids,'/',simplify = T)[3]),'_',
                                    temp$joint_post$cell)
      assign(str_split(ids,'/',simplify = T)[3],temp)
      rm(temp)
    }
  }
  
  Join_post = data.frame(cell = character(), CHROM = character(), cnv_state = character(),Score = numeric(),stringsAsFactors = FALSE)
  segs = NULL
  for (ids.n in ls(pattern = '_cellranger710$')) {
    print(ids.n)
    ids = get(ids.n)
    ids.feature = c("cell","CHROM","cnv_state",paste0('Z_',unique(ids$joint_post$cnv_state)))
    joint.pos = ids$joint_post[,..ids.feature]
    
    result = NULL
    for (stat in unique(joint.pos$cnv_state)) {
      print(stat)
      temp = subset(joint.pos, cnv_state == stat, select = c("cell", "CHROM", "cnv_state", paste0('Z_',stat)))
      colnames(temp) = c("cell","CHROM","cnv_state", "Score")
      result = rbind(result,temp)
    }
    
    Join_post = rbind(Join_post,result)
    seg = ids$segs_consensus
    
    seg = ids$bulk_clones %>%
      mutate(cnv_state = cnv_state_post) %>%
      distinct(CHROM, seg, cnv_state, seg_start, seg_end) %>%
      filter(cnv_state != 'neu') %>%
      mutate(length = seg_end - seg_start)
    seg$SampleID = ids.n
    segs = rbind(segs,seg)
  }
  
  ggplot(Join_post,aes(CHROM, Score)) +
    geom_jitter() + 
    facet_wrap(~cnv_state)
  
  
  
  segs = segs %>%
    mutate(group = SampleID) 
  ids = names(table(segs$SampleID))[table(segs$SampleID)>15]
  segs = segs[which(segs$SampleID %in% ids),]
  x=segs%>%
    cnv_heatmap() +
    theme(
      axis.title.y = element_text(size = 8, hjust = -5),
      plot.margin = margin(0, 0, 0, 0.5, unit = 'mm')
    )
  topptx(x,'plts/Numbat.pptx',width = 18)
  
}


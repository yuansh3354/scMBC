### ---------------
### Create: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
### Date:  2025-02-18
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Github: https://github.com/yuansh3354/
### Official Account: DeepBioinformatics
### Address:
###         1. Fujian Medical University. No. 1 Xue Yuan Road, University Town, 350122 FuZhou Fujian, China.
###         2. National Center for Nanoscience and Technology (NCNST). No.11 ZhongGuanCun BeiYiTiao, 100190 Beijing, China.
###         3. Department of Urology, National Cancer Center/National Clinical Research Center for Cancer/Cancer Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, 100021 China
### ---------------

# Step.00 ----------------------------------------------------------------------
# 配置函数
rm(list = ls())
gc(reset = T)
source('script/Functions/Functions.R')
# unload packages
pkgs[which(load.packages == F)]
library(parallel)
library(qs)
rpid = Sys.getpid()
options(scipen = 100)
options(future.globals.maxSize = 2000 * 1024^3)
#-------------------------------------------------------------------------------
My_cellchat <- function(input_obj,
                        assay= NULL,
                        group.by = NULL,
                        workers,
                        species=c('human','mouse'),
                        CellChatDB.use=NULL,
                        PPIuse=F,
                        type="triMean",
                        min.cells = 10
){
  

  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if(species=='human'){
    
    CellChatDB <- CellChatDB.human
    ppi = PPI.human
  }
  
  if(species =="mouse"){
    
    CellChatDB <- CellChatDB.mouse
    ppi = PPI.mouse
  }
  
  
  if(is.null(CellChatDB.use)){
    
    cellchat.obj@DB <- CellChatDB
    
  }else{
    
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj) 
  future::plan("multisession", workers = workers) 
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if(PPIuse==F){
    
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
    
  }else{
    
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use=F, type = type)
  }
  
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)


  cellchat.obj <- aggregateNet(cellchat.obj)
  
  return(cellchat.obj)
  
}

if(T){
  sce$myCellChat = ifelse(sce$Cell.Type.L1 %in% c('Cancer Cells','Macrophage'),sce$Cell.Type.L2,sce$Cell.Type.L1)

  Idents(ids.sce) = 'myCellChat'
  use.cells <- names(table(Idents(ids.sce))[table(Idents(ids.sce))>10])
  ids.sce = subset(ids.sce,myCellChat %in% use.cells)
  Idents(ids.sce) = 'myCellChat'  
  cellchat <- My_cellchat(ids.sce, assay = 'RNA', group.by = "myCellChat",workers=512, species='human')
}


### ---------------
### Create: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
### Date:  2024-04-06
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Github: https://github.com/yuansh3354/
### Official Account: DeepBioinformatics
### Address:
###         1. Fujian Medical University. No. 1 Xue Yuan Road, University Town, 350122 FuZhou Fujian, China.
###         2. National Center for Nanoscience and Technology (NCNST). No.11 ZhongGuanCun BeiYiTiao, 100190 Beijing, China.
###         3. Department of Urology, National Cancer Center/National Clinical Research Center for Cancer/Cancer Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, 100021 China
### ---------------

setwd("/Volumes/BANQ/")
rm(list = ls())
gc(reset = T)
source('script/00.Functions.R')
# unload packages
pkgs[which(load.packages == F)]
########################## Step.01 Load data ########################## 
# Load all CellRanger result
if(T){
  files = dir('data/scRNA/',pattern = '_cellranger710')
  for(file in files){
    file.path = paste('data/scRNA',file,'filtered_feature_bc_matrix',sep = '/')
    print(file.path)
    temp = CreateSeuratObject(Read10X(file.path),min.cells = 200,min.features = 300)
    temp$SampleID = gsub('_cellranger710','',file)
    temp$CellID = rownames(temp@meta.data)
    # temp = eliminate_genes(temp)
    print(file)
    print(temp)
    temp = RenameCells(object =  temp, add.cell.id = sub('_cellranger710','',file))
    assign(gsub('_cellranger710','',file),temp)
  }
  rm(temp)
  
  for(ids in ls(pattern = "_cellranger710$")){print(ids)}
  ids = ls(pattern = "_cellranger710$")
  sce.merge.obj = merge(get(ids[1]),
                        y=c(get(ids[2]),get(ids[3]),get(ids[4]),get(ids[5]),get(ids[6])
                           get(ids[7]),get(ids[8]),get(ids[9]),get(ids[10]),get(ids[11])
                           get(ids[12]),get(ids[13]),get(ids[14]),get(ids[15]),get(ids[16])
                           get(ids[17]),get(ids[18]),get(ids[19]),get(ids[20]))
  )
  saveRDS(sce.merge.obj,'result/RDS/01.Merge_RawData.RDS')
}

# INFO[SAVE FILENAME]: result/RDS/01.Merge_RawData.RDS
# INFO[END TIME]: 2024-04-06 11:04:53

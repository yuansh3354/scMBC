########################## Step.01 Load data ########################## 
# Load all CellRanger result
if(T){
  files = dir('data/scRNA/',pattern = '_cellranger')
  for(file in files){
    file.path = paste('data/scRNA',file,'filtered_feature_bc_matrix',sep = '/')
    print(file.path)
    temp = CreateSeuratObject(Read10X(file.path),min.cells = 200,min.features = 300)
    temp$SampleID = gsub('_cellranger','',file)
    temp$CellID = rownames(temp@meta.data)
    # temp = eliminate_genes(temp)
    print(file)
    print(temp)
    temp = RenameCells(object =  temp, add.cell.id = sub('_cellranger','',file))
    assign(gsub('_cellranger','',file),temp)
  }
  rm(temp)
  
  for(ids in ls(pattern = "_cellranger$")){print(ids)}
  ids = ls(pattern = "_cellranger$")
  merge_sce_by_ids <- function(ids, env = .GlobalEnv, project = "merged") {
    missing_ids <- ids[!sapply(ids, exists, envir = env)]
    obj.list <- mget(ids, envir = env)
    if (length(obj.list) == 1) {
      return(obj.list[[1]])
    }
    merged.obj <- merge(
      x = obj.list[[1]],
      y = obj.list[-1],
      add.cell.ids = ids,
      project = project
    )
    return(merged.obj)
  }
  sce.merge.obj <- merge_sce_by_ids(ids)
  saveRDS(sce.merge.obj,'result/RDS/01.Merge_RawData.RDS')
}

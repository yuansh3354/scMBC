### ---------------

### Create: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)

### Date:  2024-02-27

### Email: yuansh3354@163.com

### Github: https://github.com/yuansh3354/

### Official Account: DeepBioinformatics

### Address:

###         1. Fujian Medical University. No. 1 Xue Yuan Road, University Town, 350122 FuZhou Fujian, China.

###         2. National Center for Nanoscience and Technology (NCNST). No.11 ZhongGuanCun BeiYiTiao, 100190 Beijing, China.

###         3. Department of Urology, National Cancer Center/National Clinical Research Center for Cancer/Cancer Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, 100021 China

### ---------------
rm(list = ls())
gc(reset = T)
source('script/Functions/Functions.R')
pkgs[which(load.packages == F)]
rpid = Sys.getpid()
options(future.globals.maxSize = 1000 * 1024^3)
#-------------------------------------------------------------------------------
sce = readRDS('Result/sce_Cancer_only.RDS')
{
  obj_monocle = My_RunMonocle(sce, assay = "RNA", slot = "counts")
  obj_monocle = orderCells(obj_monocle)
  
  cytotrace_result <- cytotrace2(sce, 
                                   is_seurat = TRUE, 
                                   slot_type = "counts", 
                                   species = 'human',
                                   ncores = 64)
                                   
  obj_monocle$CytoTRACE2_Relative = cytotrace_result$CytoTRACE2_Relative
  obj_monocle$CytoTRACE2_Potency = cytotrace_result$CytoTRACE2_Potency
  obj_monocle$CytoTRACE2_Score = cytotrace_result$CytoTRACE2_Score
  obj_monocle$preKNN_CytoTRACE2_Score = cytotrace_result$preKNN_CytoTRACE2_Score
  obj_monocle$preKNN_CytoTRACE2_Potency = cytotrace_result$preKNN_CytoTRACE2_Potency  
  
  plot_cell_trajectory(obj_monocle, 
                       color_by="State", 
                       cell_size = 1,
                       cell_link_size = 1)
  obj_monocle = orderCells(obj_monocle, root_state = 3, num_paths = NULL, reverse = NULL)
  BEAM_res <- BEAM(obj_monocle, branch_point = ids,
                     cores = 496, 
                     progenitor_method = 'duplicate')
}

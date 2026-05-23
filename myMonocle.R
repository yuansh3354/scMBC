rm(list = ls())
gc(reset = T)
source('script/Functions/Functions.R')
pkgs[which(load.packages == F)]
rpid = Sys.getpid()
options(future.globals.maxSize = 1000 * 1024^3)
#-------------------------------------------------------------------------------
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
  obj_monocle = orderCells(obj_monocle)
  BEAM_res <- BEAM(obj_monocle, branch_point = ids,
                     cores = 496, 
                     progenitor_method = 'duplicate')
}

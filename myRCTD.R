### ---------------
### Create: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
### Date:  2024-02-27
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Github: https://github.com/yuansh3354/
### Official Account: DeepBioinformatics
### Address:
###         1. Fujian Medical University. No. 1 Xue Yuan Road, University Town, 350122 FuZhou Fujian, China.
###         2. National Center for Nanoscience and Technology (NCNST). No.11 ZhongGuanCun BeiYiTiao, 100190 Beijing, China.
###         3. Department of Urology, National Cancer Center/National Clinical Research Center for Cancer/Cancer Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, 100021 China
### ---------------

# ============================================================================
# Environment Setup
# ============================================================================
rm(list = ls())
gc(reset = T, verbose = F)
source('script/Functions/Functions.R')

pkgs[which(load.packages == F)]
rpid = Sys.getpid()
options(future.globals.maxSize = 1000 * 1024^3) 

ids.file = 'results/04.RCTD/'
dir.exists(ids.file) || dir.create(ids.file)

script_start_time <- Sys.time()
cat("Start:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# ============================================================================
# Load Libraries and Parse Arguments
# ============================================================================
library(argparse)
library(spacexr)

parser = ArgumentParser(description = "RCTD spatial deconvolution analysis")

parser$add_argument(
  "--input_dir",
  type = "character",
  help = "Working directory for input files [default: %(default)s]"
)

parser$add_argument(
  "--h5ad_files",
  type = "character",
  help = "Input h5ad file names [default: %(default)s]"
)

parser$add_argument(
  "--ref_file",
  type = "character",
  help = "Reference Seurat object file path [default: %(default)s]"
)

parser$add_argument(
  "--cell_type_column",
  type = "character",
  default = "Cell.Type.L2",
  help = "Column name for cell type annotation [default: %(default)s]"
)

parser$add_argument(
  "--min_cells",
  type = "integer",
  default = 50,
  help = "Minimum number of cells per cell type [default: %(default)s]"
)

args = parser$parse_args()

input_dir = args$input_dir
h5ad_files = mysplit(args$h5ad_files)
ref_file = args$ref_file
cell_type_column = args$cell_type_column
min_cells = args$min_cells

# ============================================================================
# Load and Filter Reference Data
# ============================================================================
ref = qread(ref_file)

use.cells = table(ref[[cell_type_column]])[table(ref[[cell_type_column]]) > min_cells] %>% names
ref = ref[, ref[[cell_type_column]] %in% use.cells]

# ============================================================================
# RCTD Analysis Loop
# ============================================================================
for(h5ad_file in h5ad_files){
  print(h5ad_file)
  
  object = schard::h5ad2seurat_spatial(file.path(input_dir, h5ad_file))
  object@images[[unique(object$sample)]]@scale.factors$spot = 20
  
  counts <- ref@assays$RNA$counts
  cluster <- as.factor(ref[[cell_type_column, drop = TRUE]])
  nUMI <- ref$total_counts
  cluster <- droplevels(cluster)
  
  reference <- Reference(counts, cluster, nUMI)
  
  counts_hd <- object[["Spatial"]]$counts
  object_cells_hd <- colnames(object[["Spatial"]])
  coords <- GetTissueCoordinates(object)[object_cells_hd, 1:2]
  
  query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
  
  RCTD <- create.RCTD(query, reference, max_cores = 100)
  myRCTD.d <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  qsave(myRCTD.d, 
        file = paste0(ids.file, h5ad_file, '_RCTD.', cell_type_column, '.qs'), 
        nthreads = 256)
}

# ============================================================================
# Script Completion
# ============================================================================
script_end_time <- Sys.time()
total_duration <- difftime(script_end_time, script_start_time, units = "mins")
cat("Total script runtime:", round(total_duration, 2), "minutes\n")

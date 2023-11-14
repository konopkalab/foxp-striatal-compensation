##-------------------------------------------------------
## SEURAT ANALYSIS | MODULES & LIBRARIES
##-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

##-------------------------------------------------------
## LOAD LIBRARIES
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(harmony)
library(scCustomize)
set.seed(10)


##-------------------------------------------------------
## FIND ALL MARKERS
rm(list = ls())
load("SEURAT_NA_P56_INTERGRATION_HARMONY_SPN.RData")
# seuObjSel

Idents(seuObjSel) <- "RNA_snn_res.1.2"

str.markers <- FindAllMarkers(seuObjSel, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

write.table(str.markers, "SEURAT_NA_P56_INTERGRATION_HARMONY_SPN_CLU_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
save(str.markers, file = "SEURAT_NA_P56_INTERGRATION_HARMONY_SPN_CLU_MARKERS.RData")



#-------------------------------------------------------
# END
#-------------------------------------------------------





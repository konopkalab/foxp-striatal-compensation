##-------------------------------------------------------
## LOAD MODULES AND LIBRARIES
#-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl



rm(list = ls())
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Matrix.utils))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(scCustomize))


##-------------------------------------------------------
## DEG | PSEUDOBULK
##-------------------------------------------------------
load("./../../SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt$CellType)
#  A_PROG B_NEURONAL_PROG          C_dSPN          D_iSPN          E_eSPN 
#    5181           13062           52124           43691            8140

table(seuObjFilt$Genotype)
#   CTL   P1CKO P1P2CKO   P2CKO 
# 33889   28990   28844   30475

table(seuObjFilt$GenoAgeSample)
#     CTL_P09_NA13     CTL_P09_NA14     CTL_P09_NA29   P1CKO_P09_NA32 
#             8647            12579            12663            10972 
#   P1CKO_P09_NA33   P1CKO_P09_NA42 P1P2CKO_P09_NA15 P1P2CKO_P09_NA20 
#             7006            11012             9760             7846 
# P1P2CKO_P09_NA24   P2CKO_P09_NA16   P2CKO_P09_NA22   P2CKO_P09_NA43 
#            11238             4212             7690            18573

Idents(seuObjFilt) <- "CellType"

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

seuObjFiltSel <- subset(seuObjFilt, ident = c("C_dSPN"))

table(seuObjFiltSel$CellTypeCluster)
#  C_dSPN_0  C_dSPN_1 C_dSPN_12  C_dSPN_5  C_dSPN_8 
#     16111     15045      3188      9748      8032

striosomeGenes <- c("Kcnip1", "Oprm1", "Sepw1", "Isl1", "Pdyn", "Lypd1", "Nat1", "Sorsc1")
matrixGenes <- c("Epha4", "Ebf1", "Rasgrp2", "Mef2c", "Brinp3", "Sgcz")
otherGenes <- c("Drd", "Tac1", "Drd2", "Penk", "Casz1", "Foxp1", "Foxp2")
mygenes <- c(striosomeGenes, matrixGenes, otherGenes)

plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjFiltSel, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", label = TRUE, raster = TRUE, label.size = 2)
ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER.pdf", plot = plotCluUMAP1, width = 8, height = 6, units = "in", dpi = 300)

plotDot2 <- DotPlot_scCustom(seurat_object = seuObjFiltSel, features = mygenes, group.by = "CellTypeCluster", flip_axes = T, remove_axis_titles = FALSE, x_lab_rotate = TRUE, col.min = 0, col.max = 3, dot.min = 0, dot.scale = 10)
ggsave(filename = "NA_SEURAT_DOTPLOT_STRIASOME_MATRIX_A.pdf", plot = plotDot2, width = 6, height = 9, units = "in", dpi = 300)

plotVln3 <- Stacked_VlnPlot(seurat_object = seuObjFiltSel, features = mygenes, group.by = "CellTypeCluster", x_lab_rotate = TRUE)
ggsave(filename = "NA_SEURAT_VLNPLOT_STRIASOME_MATRIX_B.pdf", plot = plotVln3, width = 6, height = 12, units = "in", dpi = 300)

# Matrix <- C_dSPN_0, C_dSPN_1, C_dSPN_5
# Striosome <- C_dSPN_8, C_dSPN_12 

seuObjFiltSel$SubType <- gsub("^C_dSPN_0$|^C_dSPN_1$|^C_dSPN_5$", "Matrix", seuObjFiltSel$CellTypeCluster)
seuObjFiltSel$SubType <- gsub("^C_dSPN_8$|^C_dSPN_12$", "Striosome", seuObjFiltSel$SubType)
table(seuObjFiltSel$SubType)
#  Matrix Striosome 
#   40904     11220

Idents(seuObjFiltSel) <- "Genotype"

seuObjFiltSel_CTL <- subset(seuObjFiltSel, ident = c("CTL"))

table(seuObjFiltSel_CTL$SubType)
#  Matrix Striosome 
#   11992      3518

str.selected <- seuObjFiltSel_CTL
Idents(str.selected) <- "SubType"

## Prepare data for MAST ##
mat <- str.selected@assays$RNA@data
meta <- str.selected@meta.data

## Keep genes with 10% expression in at least one species
cells1 <- WhichCells(str.selected, idents = 'Striosome')
cells2 <- WhichCells(str.selected, idents = 'Matrix')
pass1 = rowSums(mat[,cells1])/length(cells1) > 0.1
pass2 = rowSums(mat[,cells2])/length(cells2) > 0.1
mat = mat[pass1 | pass2,]

## Create MAST object
sca <- MAST::FromMatrix(exprsArray = as.matrix(x = mat),
              cData = meta,
              fData = data.frame(rownames(mat)))

## Scale number of detected genes (cngeneson. Same as nFeature_RNA only after filtering)
cdr2 <- colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)

## Calculate fold change
avg_logfc <- log(rowMeans(expm1(mat[,cells1])) + 1) - log(rowMeans(expm1(mat[,cells2])) + 1)

## MAST package with covariates and glm
mastfix = MAST::zlm(~SubType + cngeneson + Batch + Sex, sca, ebayes = F, method = 'glm')

summaryCond <- summary(object = mastfix, doLRT = 'SubTypeStriosome')
summaryDt <- summaryCond$datatable
p_val <- summaryDt[summaryDt$component == "H", 4]
genes.return <- summaryDt[summaryDt$component == "H", 1]
to.return <- data.frame(p_val, row.names = genes.return$primerid)
colnames(to.return)[1] = 'p_value'
fix_res = to.return
fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
fix_res$avg_logfc = avg_logfc
print(sum(abs(fix_res$avg_logfc) >= 0.1375 & fix_res$adj_p_value <= 0.05))

mast_allcov_glm = fix_res


## average gene expression
avg.exp <- as.data.frame(log1p(AverageExpression(str.selected, verbose = FALSE)$RNA))
avg.exp$Gene <- row.names(avg.exp)

degavg <- merge(avg.exp, fix_res, by = "row.names", all.y = TRUE)
row.names(degavg) <- degavg$Row.names
degavg$Row.names <- NULL

## write degavg table
write.table(degavg, "SEURAT_NA_DATA_CTL_numeratorSTRIOSOME_denominatorMATRIX_MASTGLM.txt", row.names = T, col.names = T, quote = F, sep = "\t")

## save dge table and avg gene expression
save(mast_allcov_glm, avg.exp, degavg, file = "SEURAT_NA_DATA_CTL_numeratorSTRIOSOME_denominatorMATRIX_MASTGLM.RData")



##-------------------------------------------------------
##-------------------------------------------------------
## END
##-------------------------------------------------------
##-------------------------------------------------------



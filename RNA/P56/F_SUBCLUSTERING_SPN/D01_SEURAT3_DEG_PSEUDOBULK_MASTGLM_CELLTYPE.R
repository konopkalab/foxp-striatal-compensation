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
# source("~/onlybiohpc/pr3/organized/utility_functions.R")


##-------------------------------------------------------
## DEG | PSEUDOBULK
##-------------------------------------------------------
load("./../../SEURAT_NA_P56_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt$CellType)
# A_dSPN B_iSPN C_eSPN 
#  14465  11181   6001

table(seuObjFilt$Genotype)
#     CTL   P1CKO P1P2CKO   P2CKO 
#   11143    6365    8579    5560

table(seuObjFilt$GenoAgeSample)
#     P56_CTL_NA03     P56_CTL_NA05     P56_CTL_NA07   P56_P1CKO_NA01 
#             3615             3941             3587             2621 
#   P56_P1CKO_NA06   P56_P1CKO_NA10 P56_P1P2CKO_NA02 P56_P1P2CKO_NA04 
#             2342             1402              803             3360 
# P56_P1P2CKO_NA45   P56_P2CKO_NA11   P56_P2CKO_NA25   P56_P2CKO_NA41 
#             4416              184             1475              185 
#   P56_P2CKO_NA44 
#             3716

Idents(seuObjFilt) <- "CellType"

table(seuObjFilt@active.ident)
# A_dSPN B_iSPN C_eSPN 
#  14465  11181   6001

seuObjFiltSel <- subset(seuObjFilt, ident = c("A_dSPN"))

table(seuObjFiltSel$CellTypeCluster)
#  A_dSPN_1 A_dSPN_10 A_dSPN_11  A_dSPN_2  A_dSPN_6  A_dSPN_8 
#      4481      1489      1397      3115      2209      1774

striosomeGenes <- c("Kcnip1", "Oprm1", "Sepw1", "Isl1", "Pdyn", "Lypd1", "Nat1", "Sorsc1")
matrixGenes <- c("Epha4", "Ebf1", "Rasgrp2", "Mef2c", "Brinp3", "Sgcz")
otherGenes <- c("Drd", "Tac1", "Drd2", "Penk", "Casz1", "Foxp1", "Foxp2")
mygenes <- c(striosomeGenes, matrixGenes, otherGenes)

plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjFiltSel, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", label = TRUE, raster = TRUE, label.size = 3)
ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER.pdf", plot = plotCluUMAP1, width = 8, height = 6, units = "in", dpi = 300)

plotDot2 <- DotPlot_scCustom(seurat_object = seuObjFiltSel, features = mygenes, group.by = "CellTypeCluster", flip_axes = T, remove_axis_titles = FALSE, x_lab_rotate = TRUE, col.min = 0, col.max = 3, dot.min = 0, dot.scale = 10)
ggsave(filename = "NA_SEURAT_DOTPLOT_STRIASOME_MATRIX_A.pdf", plot = plotDot2, width = 6, height = 9, units = "in", dpi = 300)

plotVln3 <- Stacked_VlnPlot(seurat_object = seuObjFiltSel, features = mygenes, group.by = "CellTypeCluster", x_lab_rotate = TRUE)
ggsave(filename = "NA_SEURAT_VLNPLOT_STRIASOME_MATRIX_B.pdf", plot = plotVln3, width = 6, height = 12, units = "in", dpi = 300)

# Matrix <- A_dSPN_1, A_dSPN_6, A_dSPN_8
# Striosome <- A_dSPN_2, A_dSPN_10, A_dSPN_11

seuObjFiltSel$SubType <- gsub("^A_dSPN_1$|^A_dSPN_6$|^A_dSPN_8$", "Matrix", seuObjFiltSel$CellTypeCluster)
seuObjFiltSel$SubType <- gsub("^A_dSPN_2$|^A_dSPN_10$|^A_dSPN_11$", "Striosome", seuObjFiltSel$SubType)
table(seuObjFiltSel$SubType)
#    Matrix Striosome 
#      8464      6001

Idents(seuObjFiltSel) <- "Genotype"

seuObjFiltSel_CTL <- subset(seuObjFiltSel, ident = c("CTL"))

table(seuObjFiltSel_CTL$SubType)
#    Matrix Striosome 
#      3344      1982

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
# mastfix = MAST::zlm(~SubType + cngeneson + LibPrep + SeqBatch + Sex, sca, ebayes = F, method = 'glm')
mastfix = MAST::zlm(~SubType + cngeneson + LibPrep + SeqBatch, sca, ebayes = F, method = 'glm') ## Sex as a co-variate was removed as all CTL samples are F.

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



# for (myclu in unique(sort(seuObjFiltSel@active.ident)))
# 	{
#   ## select a cell-type    
#   # myclu <- "E_eSPN"
# 	cluSelected <- myclu
# 	print(cluSelected)

#   ## subset for a cell-type
# 	str.selected <- subset(x = seuObjFilt, idents = c(cluSelected))
#   Idents(str.selected) <- "Genotype"

#   ## Prepare data for MAST ##
#   mat <- str.selected@assays$RNA@data
#   meta <- str.selected@meta.data

#   ## Keep genes with 10% expression in at least one species
#   cells1 <- WhichCells(str.selected, idents = 'CTL')
#   cells2 <- WhichCells(str.selected, idents = 'P1CKO')
#   pass1 = rowSums(mat[,cells1])/length(cells1) > 0.1
#   pass2 = rowSums(mat[,cells2])/length(cells2) > 0.1
#   mat = mat[pass1 | pass2,]

#   ## Create MAST object
#   sca <- MAST::FromMatrix(exprsArray = as.matrix(x = mat),
#                 cData = meta,
#                 fData = data.frame(rownames(mat)))

#   ## Scale number of detected genes (cngeneson. Same as nFeature_RNA only after filtering)
#   cdr2 <- colSums(assay(sca)>0)
#   colData(sca)$cngeneson <- scale(cdr2)

#   ## Calculate fold change
#   avg_logfc <- log(rowMeans(expm1(mat[,cells1])) + 1) - log(rowMeans(expm1(mat[,cells2])) + 1)

#   # ## MAST package all covariates #
#   # mastfix = MAST::zlm(~Genotype + cngeneson + Batch + Sex, sca)

#   # summaryCond <- summary(object = mastfix, doLRT = 'GenotypeP1CKO')
#   # summaryDt <- summaryCond$datatable
#   # p_val <- summaryDt[summaryDt$component == "H", 4]
#   # genes.return <- summaryDt[summaryDt$component == "H", 1]
#   # to.return <- data.frame(p_val, row.names = genes.return$primerid)
#   # colnames(to.return)[1] = 'p_value'
#   # fix_res = to.return
#   # fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
#   # fix_res$avg_logfc = avg_logfc
#   # sum(abs(fix_res$avg_logfc) >= 0.1375 & fix_res$adj_p_value <= 0.05)

#   # mast_allcov_bayes = fix_res


#   ## MAST package with covariates and glm
#   mastfix = MAST::zlm(~Genotype + cngeneson + Batch + Sex, sca, ebayes = F, method = 'glm')

#   summaryCond <- summary(object = mastfix, doLRT = 'GenotypeP1CKO')
#   summaryDt <- summaryCond$datatable
#   p_val <- summaryDt[summaryDt$component == "H", 4]
#   genes.return <- summaryDt[summaryDt$component == "H", 1]
#   to.return <- data.frame(p_val, row.names = genes.return$primerid)
#   colnames(to.return)[1] = 'p_value'
#   fix_res = to.return
#   fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
#   fix_res$avg_logfc = avg_logfc
#   print(sum(abs(fix_res$avg_logfc) >= 0.1375 & fix_res$adj_p_value <= 0.05))

#   mast_allcov_glm = fix_res


#   # ## MAST package with covariates and glmer and mixed model
#   # mastfix_glmer = MAST::zlm(~Genotype + (1|Sample) + Batch + Sex + cngeneson, sca, ebayes = F, method = 'glmer')
  
#   # summaryCond <- summary(object = mastfix_glmer, doLRT = 'GenotypeP1CKO')
#   # summaryDt <- summaryCond$datatable
#   # p_val <- summaryDt[summaryDt$component == "H", 4]
#   # genes.return <- summaryDt[summaryDt$component == "H", 1]
#   # to.return <- data.frame(p_val, row.names = genes.return$primerid)
#   # colnames(to.return)[1] = 'p_value'
#   # fix_res_glmer = to.return
#   # fix_res_glmer$adj_p_value = p.adjust(fix_res_glmer$p_value, method = 'BH')
#   # fix_res_glmer$avg_logfc = avg_logfc
#   # # sum(abs(fix_res_glmer$avg_logfc) >= 0.125 & fix_res_glmer$adj_p_value <= 0.05)
  
#   ## write dge table
#   write.table(mast_allcov_glm, paste("SEURAT_NA_DATA_", cluSelected, "_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

# 	## average gene expression
# 	str.selected$CellType <- paste(Idents(str.selected), str.selected$Genotype, sep = "_")
# 	Idents(str.selected) <- "CellType"

# 	avg.exp <- as.data.frame(log1p(AverageExpression(str.selected, verbose = FALSE)$RNA))
# 	avg.exp$Gene <- row.names(avg.exp)

#   # avgExp <- avg.exp[row.names(avg.exp) %in% row.names(fix_res_glmer),]
#   degavg <- merge(avg.exp, fix_res, by = "row.names", all.y = TRUE)
#   row.names(degavg) <- degavg$Row.names
#   degavg$Row.names <- NULL

#   ## write average gene expression table
# 	# write.table(avg.exp, paste("SEURAT_NA_DATA_", cluSelected, "_AVG_EXP.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

#   ## save dge table and avg gene expression
#   save(mast_allcov_glm, avg.exp, degavg, file = paste("SEURAT_NA_DATA_P1CKO_X_CTL_", cluSelected, "_PSEUDOBULK.RData", sep = ""))
# 	}


##-------------------------------------------------------
##-------------------------------------------------------
## END
##-------------------------------------------------------
##-------------------------------------------------------



# 
# 
# ## MAST package all covariates # Batch_cDNA
# mastfix = MAST::zlm(~Condition + Batch_cDNA, sca)
# 
# summaryCond <- summary(object = mastfix, doLRT = 'ConditionCS')
# summaryDt <- summaryCond$datatable
# p_val <- summaryDt[summaryDt$component == "H", 4]
# genes.return <- summaryDt[summaryDt$component == "H", 1]
# to.return <- data.frame(p_val, row.names = genes.return$primerid)
# colnames(to.return)[1] = 'p_value'
# fix_res = to.return
# fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
# fix_res$avg_logfc = avg_logfc
# sum(abs(fix_res$avg_logfc) > 0.25 & fix_res$adj_p_value < 0.05)
# 
# mast_allcov_bayes = fix_res
# 
# 
# 
# # MAST package with no covariates
# mastfix = MAST::zlm(~Species, sca)
# 
# summaryCond <- summary(object = mastfix, doLRT = 'Specieschimp')
# summaryDt <- summaryCond$datatable
# p_val <- summaryDt[summaryDt$component == "H", 4]
# genes.return <- summaryDt[summaryDt$component == "H", 1]
# to.return <- data.frame(p_val, row.names = genes.return$primerid)
# colnames(to.return)[1] = 'p_value'
# fix_res = to.return
# fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
# fix_res$avg_logfc = avg_logfc
# sum(abs(fix_res$avg_logfc) > 0.25 & fix_res$adj_p_value < 0.05)
# 
# mast_nocov_bayes = fix_res
# 
# 
# # MAST package without covariates but glm
# mastfix = MAST::zlm(~Species, sca, ebayes = F, method = 'glm')
# 
# summaryCond <- summary(object = mastfix, doLRT = 'Specieschimp')
# summaryDt <- summaryCond$datatable
# p_val <- summaryDt[summaryDt$component == "H", 4]
# genes.return <- summaryDt[summaryDt$component == "H", 1]
# to.return <- data.frame(p_val, row.names = genes.return$primerid)
# colnames(to.return)[1] = 'p_value'
# fix_res = to.return
# fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
# fix_res$avg_logfc = avg_logfc
# sum(abs(fix_res$avg_logfc) > 0.25 & fix_res$adj_p_value < 0.05)
# 
# mast_nocov_glm = fix_res
# 
# # MAST package with covariates and glm
# mastfix = MAST::zlm(~Species + cngeneson + log(human_age) + lib_batch + sex, sca, ebayes = F, method = 'glm')
# 
# summaryCond <- summary(object = mastfix, doLRT = 'Specieschimp')
# summaryDt <- summaryCond$datatable
# p_val <- summaryDt[summaryDt$component == "H", 4]
# genes.return <- summaryDt[summaryDt$component == "H", 1]
# to.return <- data.frame(p_val, row.names = genes.return$primerid)
# colnames(to.return)[1] = 'p_value'
# fix_res = to.return
# fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
# fix_res$avg_logfc = avg_logfc
# sum(abs(fix_res$avg_logfc) > 0.25 & fix_res$adj_p_value < 0.05)
# 
# mast_allcov_glm = fix_res
# 
# # MAST package with covariates and glmer and mixed model
# mastfix = MAST::zlm(~Species + (1|Sample) + cngeneson + log(human_age) + lib_batch + sex, sca, ebayes = F, method = 'glmer')
# 
# summaryCond <- summary(object = mastfix, doLRT = 'Specieschimp')
# summaryDt <- summaryCond$datatable
# p_val <- summaryDt[summaryDt$component == "H", 4]
# genes.return <- summaryDt[summaryDt$component == "H", 1]
# to.return <- data.frame(p_val, row.names = genes.return$primerid)
# colnames(to.return)[1] = 'p_value'
# fix_res = to.return
# fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
# fix_res$avg_logfc = avg_logfc
# sum(abs(fix_res$avg_logfc) > 0.25 & fix_res$adj_p_value < 0.05)
# 
# mast_allcov_glm_mixed = fix_res
# 
# 








# strSCE <- as.SingleCellExperiment(str.integrated)

# strGroups <- colData(strSCE)[, c("ConditionSample", "Condition")]
# strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
# strAggMat <- strPseudoSCE@assays@data$sum
# colnames(strAggMat) <- colData(strPseudoSCE)[['ConditionSample']]


# csExp = which(rowSums(strAggMat[, grepl('CS', colnames(strAggMat))] <= 10) == 0) %>% names
# sdExp = which(rowSums(strAggMat[, grepl('SD', colnames(strAggMat))] <= 10) == 0) %>% names
# expgns = Reduce(union, list(csExp, sdExp))

# DefaultAssay(str.integrated) <- "RNA"
# selSCE <- as.SingleCellExperiment(str.integrated)

# selGroups <- colData(selSCE)[, c("ConditionSample", "Condition", "Sample", "Batch")]
# selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
# selGroups <- colData(selPseudoSCE)[, c("ConditionSample", "Condition", "Sample", "Batch")]
# selGroups$ConditionSample <- factor(selGroups$ConditionSample)
# selGroups$Condition <- factor(selGroups$Condition)
# selGroups$Sample <- factor(selGroups$Sample)
# selGroups$Batch <- factor(selGroups$Batch)

# selAggMat <- selPseudoSCE@assays@data$sum
# colnames(selAggMat) = colData(selPseudoSCE)[['ConditionSample']]
# selAggMat <- selAggMat[expgns,]

# selDGEL <- DGEList(counts = selAggMat)
# selDGEL <- calcNormFactors(selDGEL)
# selDesign <- model.matrix(~Condition, data = selGroups)
# selDGEL <- estimateDisp(selDGEL, selDesign)

# selFit <- glmFit(selDGEL,selDesign)
# selLrt <- glmLRT(selFit)
# selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame
# DEGs <- subset(selRes, FDR <= 0.05)


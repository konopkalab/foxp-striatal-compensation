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
# source("~/onlybiohpc/pr3/organized/utility_functions.R")


##-------------------------------------------------------
## DEG | PSEUDOBULK
##-------------------------------------------------------
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/SEURAT_NA_P56_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

print(table(seuObjFilt$CellType))

print(table(seuObjFilt$Genotype))

print(table(seuObjFilt$GenoAgeSample))

Idents(seuObjFilt) <- "CellType"

print(table(seuObjFilt@active.ident))


for (myclu in unique(sort(seuObjFilt@active.ident)))
	{
  ## select a cell-type    
  # myclu <- "E_eSPN"
	cluSelected <- myclu
	print(cluSelected)

  ## subset for a cell-type
	str.selected <- subset(x = seuObjFilt, idents = c(cluSelected))
  Idents(str.selected) <- "Genotype"

  ## Prepare data for MAST ##
  mat <- str.selected@assays$RNA@data
  meta <- str.selected@meta.data

  ## Keep genes with 10% expression in at least one species
  cells1 <- WhichCells(str.selected, idents = 'CTL')
  cells2 <- WhichCells(str.selected, idents = 'P2CKO')
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

  # ## MAST package all covariates #
  # mastfix = MAST::zlm(~Genotype + cngeneson + Batch + Sex, sca)

  # summaryCond <- summary(object = mastfix, doLRT = 'GenotypeP1CKO')
  # summaryDt <- summaryCond$datatable
  # p_val <- summaryDt[summaryDt$component == "H", 4]
  # genes.return <- summaryDt[summaryDt$component == "H", 1]
  # to.return <- data.frame(p_val, row.names = genes.return$primerid)
  # colnames(to.return)[1] = 'p_value'
  # fix_res = to.return
  # fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
  # fix_res$avg_logfc = avg_logfc
  # sum(abs(fix_res$avg_logfc) >= 0.1375 & fix_res$adj_p_value <= 0.05)

  # mast_allcov_bayes = fix_res


  ## MAST package with covariates and glm
  mastfix = MAST::zlm(~Genotype + cngeneson + LibPrep + SeqBatch + Sex, sca, ebayes = F, method = 'glm')

  summaryCond <- summary(object = mastfix, doLRT = 'GenotypeP2CKO')
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
  
  ## write dge table
  write.table(mast_allcov_glm, paste("SEURAT_NA_DATA_", cluSelected, "_P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

	## average gene expression
	str.selected$CellTypeClusterGeno <- paste(str.selected$CellTypeCluster, Idents(str.selected), sep = "_")
	Idents(str.selected) <- "CellTypeClusterGeno"

	avg.exp <- as.data.frame(log1p(AverageExpression(str.selected, verbose = FALSE)$RNA))
	avg.exp$Gene <- row.names(avg.exp)

  # avgExp <- avg.exp[row.names(avg.exp) %in% row.names(fix_res_glmer),]
  degavg <- merge(avg.exp, fix_res, by = "row.names", all.y = TRUE)
  row.names(degavg) <- degavg$Row.names
  degavg$Row.names <- NULL

  ## write average gene expression table
	# write.table(avg.exp, paste("SEURAT_NA_DATA_", cluSelected, "_AVG_EXP.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

  ## save dge table and avg gene expression
  save(mast_allcov_glm, avg.exp, degavg, file = paste("SEURAT_NA_DATA_P2CKO_X_CTL_", cluSelected, "_PSEUDOBULK.RData", sep = ""))
	}


##-------------------------------------------------------
##-------------------------------------------------------
## END
##-------------------------------------------------------
##-------------------------------------------------------


##-------------------------------------------------------
## LOAD MODULES AND LIBRARIES
#-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(plyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(scran)
library(scater)
library(SingleCellExperiment)
library(edgeR)

##-------------------------------------------------------
## DEG PLOTS | PSEUDOBULK
##-------------------------------------------------------

filterDEGs <- function(degtable, pvalcutoff, l2fccutoff) 
	{
	deg.cs.sd.temp <- read.table(degtable, sep = "\t", header = TRUE)
	deg.cs.sd <- deg.cs.sd.temp
	deg.cs.sd.sig.up <- deg.cs.sd[deg.cs.sd$adj_p_value <= pvalcutoff & deg.cs.sd$avg_logfc >= l2fccutoff,] 
	deg.cs.sd.sig.dn <- deg.cs.sd[deg.cs.sd$adj_p_value <= pvalcutoff & deg.cs.sd$avg_logfc <= -l2fccutoff,]
	return(list("UP" = row.names(deg.cs.sd.sig.up), "DOWN" = row.names(deg.cs.sd.sig.dn)))
	}

pvalcutoff <- 0.05
l2fccutoff <- 0.25 #0.1375 # 10% difference

comp <- "P1P2CKO_x_CTL"
degfiles <- list.files(path = ".", pattern = paste("_", comp, "_", sep = ""))
degtables <- degfiles[grepl(".txt$", degfiles)]

deglist <- list()
degcountlist <- list()
for(i in 1:length(degtables))
 {
 print(degtables[[i]])
 prefx <- gsub(paste("_", comp, "_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = ""), "", gsub("SEURAT_NA_DATA_", "", degtables[[i]]))
 updnlist <- filterDEGs(degtables[[i]], pvalcutoff, l2fccutoff)
 deglist <- c(deglist, updnlist)
 upcount <- length(updnlist$UP)
 names(upcount) <- paste(prefx, "UP", sep = "_")
 dncount <- length(updnlist$DOWN)
 names(dncount) <- paste(prefx, "DN", sep = "_")
 degcountlist <- c(degcountlist, upcount, dncount)
 }

degup <- as.data.frame(t(as.data.frame(degcountlist[grepl("_UP", names(degcountlist))])))
degup$celltype <- gsub("_UP", "", row.names(degup))
degdn <- as.data.frame(t(as.data.frame(degcountlist[grepl("_DN", names(degcountlist))])))
degdn$celltype <- gsub("_DN", "", row.names(degdn))

degcounts <- merge(degup, degdn, by = "celltype", all.x = TRUE, incomparables = 0)
colnames(degcounts) <- c("CELLTYPE", "UP", "DOWN")
degcounts[is.na(degcounts)] <- 0
write.table(degcounts, paste("NA_DEG_", comp, "_L2FC_0.1375.tsv", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

p1 <- ggplot(degcounts) + 
	  geom_hline(yintercept = 0) + 
	  geom_segment(aes(x = CELLTYPE, xend = CELLTYPE, y = 0, yend = UP)) +
	  geom_point(data = degcounts, aes(x = CELLTYPE, y = UP, color = CELLTYPE)) + 
	  geom_segment(aes(x = CELLTYPE, xend = CELLTYPE, y = 0, yend = -DOWN)) +
	  geom_point(data = degcounts, aes(x = CELLTYPE, y = -DOWN, color = CELLTYPE)) + 	  
	  theme_bw() + 
	  ylim(-1000,1000) +  #ylim(-150,150) + 
	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
	  theme(legend.position = "none") + 
	  xlab("CellTypes") + 
	  ylab("# DEGs") + 
	  NULL
ggsave(filename = paste("NA_DEG_", comp, "_L2FC_0.1375.pdf", sep = ""), plot = p1, width = 4, height = 4, units = "in", dpi = 300, useDingbats = FALSE)

# save(deglist, file = paste("NA_DEG_", comp, "_L2FC_0.1375.RData", sep = ""))











# ### VOLCANO PLOT
# VOLplot <- ggplot(data = rg.deg.wt.ko, aes(x = avg_log2FC*(-1), y = -log10(p_val_adj))) +
# 			geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
# 			geom_point(data=rg.deg.wt.ko.sig.up, colour="green3", size=1.2, shape=16, alpha=0.8) +
# 			geom_point(data=rg.deg.wt.ko.sig.dn, colour="royalblue", size=1.2, shape=16, alpha=0.8) +
# 			annotate(geom = "text", x = 1.2, y = 300, label = nrow(rg.deg.wt.ko.sig.up)) +
# 			annotate(geom = "text", x = -1.2, y = 300, label = nrow(rg.deg.wt.ko.sig.dn)) +
# 			theme_bw() +
# 			geom_vline(xintercept = log2(1.20), colour="red", linetype=2) + geom_vline(xintercept = -log2(1.20), colour="red", linetype=2) +
# 			geom_hline(yintercept = -log10(0.05), colour="red", linetype=2) +
# 			# xlim(-(abs(max(resNew2$log2FoldChange) + 0.5)), abs(max(resNew2$log2FoldChange) + 0.5)) +
# 			xlim(-1.7,1.7) + 
# 			ylim(0, 320) +
# 			labs(x="Avg. log2FC", y="- log10 adj. P-value")
# ggsave(filename = "EP_RG_DEG_Volcano.pdf", plot = VOLplot, width = 4, height = 4, units = "in", dpi = 300, useDingbats = FALSE)





# load("SEURAT_NA_DATA_A_PROG_PSEUDOBULK.RData")


# load("./../SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt

# table(seuObjFilt$CellType)
# #  A_PROG B_NEURONAL_PROG          C_dSPN          D_iSPN          E_eSPN 
# #    5181           13062           52124           43691            4245

# table(seuObjFilt$Genotype)
# #     CTL   P1CKO P1P2CKO   P2CKO 
# #   33331   28132   28422   28418

# table(seuObjFilt$GenoAgeSample)
# #     CTL_P09_NA13     CTL_P09_NA14     CTL_P09_NA29   P1CKO_P09_NA32 
# #             8542            12440            12349            10408 
# #   P1CKO_P09_NA33   P1CKO_P09_NA42 P1P2CKO_P09_NA15 P1P2CKO_P09_NA20 
# #             6941            10783             9628             7763 
# # P1P2CKO_P09_NA24   P2CKO_P09_NA16   P2CKO_P09_NA22   P2CKO_P09_NA43 
# #            11031             4109             7358            16951

# Idents(seuObjFilt) <- "CellType"

# table(seuObjFilt@active.ident)
# # D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
# # 43691           52124            4245            5181           13062


# for (myclu in unique(sort(seuObjFilt@active.ident)))
# 	{
#     ## select a cell-type    
# 	cluSelected <- myclu
# 	print(cluSelected)

# 	seuObjSel <- subset(x = seuObjFilt, idents = c(cluSelected))

#     ## create SCE object for all genes
#     strSCE <- as.SingleCellExperiment(seuObjSel)

#     strGroups <- colData(strSCE)[, c("GenoAgeSample", "Genotype")]
#     strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
#     strAggMat <- strPseudoSCE@assays@data$sum
#     colnames(strAggMat) <- colData(strPseudoSCE)[['GenoAgeSample']]

#     ## filter genes
#     ctlExp = which(rowSums(strAggMat[, grepl('CTL_', colnames(strAggMat))] <= 10) == 0) %>% names
#     p1cExp = which(rowSums(strAggMat[, grepl('P1CKO', colnames(strAggMat))] <= 10) == 0) %>% names
#     expgns = Reduce(union, list(ctlExp, p1cExp))

#     ## create SCE object for filtered genes
#     DefaultAssay(seuObjSel) <- "RNA"
#     selSCE <- as.SingleCellExperiment(seuObjSel)

#     selGroups <- colData(selSCE)[, c("GenoAgeSample", "Genotype", "Batch")]
#     selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
#     selGroups <- colData(selPseudoSCE)[, c("GenoAgeSample", "Genotype", "Batch")]
#     selGroups$GenoAgeSample <- factor(selGroups$GenoAgeSample)
#     selGroups$Genotype <- factor(selGroups$Genotype)
#     selGroups$Batch <- factor(selGroups$Batch)

#     selAggMat <- selPseudoSCE@assays@data$sum
#     colnames(selAggMat) = colData(selPseudoSCE)[['GenoAgeSample']]
#     selAggMat <- selAggMat[expgns,]

#     ## perform differential gene expression
#     selDGEL <- DGEList(counts = selAggMat)
#     selDGEL <- calcNormFactors(selDGEL)
#     selDesign <- model.matrix(~Batch + Genotype, data = selGroups)
#     selDGEL <- estimateDisp(selDGEL, selDesign)

#     selFit <- glmFit(selDGEL,selDesign)
#     selLrt <- glmLRT(selFit)
#     selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame
#     degTableSel <- subset(selRes, FDR <= 0.05)

#     ## write dge table
#     write.table(selRes, paste("SEURAT_NA_DATA_", cluSelected, "_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

# 	## average gene expression
# 	seuObjSel$CellType <- paste(Idents(seuObjSel), seuObjSel$Genotype, sep = "_")
# 	Idents(seuObjSel) <- "CellType"

# 	avg.exp <- as.data.frame(log1p(AverageExpression(seuObjSel, verbose = FALSE)$RNA))
# 	avg.exp$Gene <- row.names(avg.exp)

#     # avgExp <- avg.exp[row.names(avg.exp) %in% row.names(selRes),]
#     degavg <- merge(avg.exp, selRes, by = "row.names", all.y = TRUE)
#     row.names(degavg) <- degavg$Row.names
#     degavg$Row.names <- NULL

#     ## write average gene expression table
# 	write.table(avg.exp, paste("SEURAT_STR_DATA_", cluSelected, "_AVG_EXP.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

#     ## save dge table and avg gene expression
#     save(selRes, avg.exp, degavg, file = paste("SEURAT_STR_DATA_", cluSelected, "_PSEUDOBULK.RData", sep = ""))
# 	}


##-------------------------------------------------------
##-------------------------------------------------------
## END
##-------------------------------------------------------
##-------------------------------------------------------


# strSCE <- as.SingleCellExperiment(seuObjFilt)

# strGroups <- colData(strSCE)[, c("ConditionSample", "Condition")]
# strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
# strAggMat <- strPseudoSCE@assays@data$sum
# colnames(strAggMat) <- colData(strPseudoSCE)[['ConditionSample']]


# csExp = which(rowSums(strAggMat[, grepl('CS', colnames(strAggMat))] <= 10) == 0) %>% names
# sdExp = which(rowSums(strAggMat[, grepl('SD', colnames(strAggMat))] <= 10) == 0) %>% names
# expgns = Reduce(union, list(csExp, sdExp))

# DefaultAssay(seuObjFilt) <- "RNA"
# selSCE <- as.SingleCellExperiment(seuObjFilt)

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


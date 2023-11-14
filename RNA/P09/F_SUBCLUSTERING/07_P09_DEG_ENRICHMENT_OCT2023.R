##-------------------------------------------------------
## SEURAT ANALYSIS | MODULES & LIBRARIES
##-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

##-------------------------------------------------------
## SEURAT ANALYSIS | ANNOTATE CLUSTERS
##-------------------------------------------------------
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(WGCNA)
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(ggpubr)
library(nichenetr) ## has function to convert human genes to mouse


# ##-------------------------------------------------------
# ## Load Annotated Seurat Object to get background gene list
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt


# DefaultAssay(seuObjFilt) <- "RNA"

# Idents(seuObjFilt) <- "Genotype"
# table(seuObjFilt@active.ident)
# #     CTL   P1CKO   P2CKO P1P2CKO 
# #   33889   28990   30475   28844

# na.ctl <- subset(seuObjFilt, idents = "CTL", slot = "counts")

# Idents(na.ctl) <- "CellType"

# table(na.ctl@active.ident)
# #  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
# #   12247           15510            1683            1288            3161

# DefaultAssay(na.ctl) <- "RNA"

# na.counts <- GetAssayData(na.ctl, slot = "counts")

# ## Genes with no expression across all cells
# na.counts2 <- na.counts[apply(na.counts, 1, function(x) !all(x==0)),]
# na.counts3 <- na.counts[apply(na.counts, 1, function(x) all(x==0)),]

# background_genes_1 <- row.names(na.counts2)
# write.table(background_genes_1, "Background_Genes_for_GO_List1.txt", row.names = F, col.names = F, quote = F, sep = "\t")
# # 17462 genes


# ##-------------------------------------------------------
# ## Load DEG Tables
# rm(list = ls())
# p9dspn1 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P1CKO_CTL/SEURAT_NA_DATA_C_dSPN_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)
# p9dspn2 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P2CKO_CTL/SEURAT_NA_DATA_C_dSPN_P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)
# p9dspn3 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P1P2CKO_CTL/SEURAT_NA_DATA_C_dSPN_P1P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)

# p56dspn1 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)
# p56dspn2 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)
# p56dspn3 <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P1P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = "\t", header = TRUE)

# background_genes_2_temp <- c(row.names(p9dspn1), row.names(p9dspn2), row.names(p9dspn3), row.names(p56dspn1), row.names(p56dspn2), row.names(p56dspn3))

# background_genes_2 <- unique(sort(background_genes_2_temp)) # 6101 genes
# write.table(background_genes_2, "Background_Genes_dSPN_for_GO_List.txt", row.names = F, col.names = F, quote = F, sep = "\t")
# # 6101 genes




# rm(list = ls())
# ##-------------------------------------------------------
# ## PREPARE REFERENCE GENE LISTS
# ## Load reference geneset
# load("GeneSets_Disorders.RData")

# asd13_human <- as.vector(GeneSets$`ASD_1-3`$Gene)
# asd_human <- as.vector(GeneSets$ASD$Gene)
# fmrp_human <- as.vector(GeneSets$FMRP$Gene)
# id_human <- as.vector(GeneSets$ID$Gene)
# syn_human <- as.vector(GeneSets$SYN$Gene)

# asd13_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asd13_human))  # 467 genes
# asd_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asd_human))  # 824 genes
# fmrp_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = fmrp_human))  # 773 genes
# id_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = id_human))  # 1010 genes
# syn_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = syn_human))  # 1809 genes


# ## create a list for all reference gene lists
# myref <- list("01_ASD_1-3" = asd13_mouse, 
#               "02_ASD" = asd_mouse, 
#               "03_FMRP" = fmrp_mouse, 
#               "04_ID" = id_mouse, 
#               "05_SYN" = syn_mouse)

# save(myref, file = "NA_ENRICHMENT_REFERENCE_GENE_LISTS.RData")

load("NA_ENRICHMENT_REFERENCE_GENE_LISTS.RData")

myref$`05_SYN` <- NULL

##-------------------------------------------------------
## READ & FILTER DEG TABLES
filterDEGs <- function(degtable, pvalcutoff, l2fccutoff) 
	{
	deg.cs.sd.temp <- read.table(degtable, sep = "\t", header = TRUE)
	deg.cs.sd <- deg.cs.sd.temp
	deg.cs.sd.sig.up <- deg.cs.sd[deg.cs.sd$adj_p_value <= pvalcutoff & deg.cs.sd$avg_logfc >= l2fccutoff,] 
	deg.cs.sd.sig.dn <- deg.cs.sd[deg.cs.sd$adj_p_value <= pvalcutoff & deg.cs.sd$avg_logfc <= -l2fccutoff,]
	return(list("UP" = row.names(deg.cs.sd.sig.up), "DOWN" = row.names(deg.cs.sd.sig.dn)))
	}

pvalcutoff <- 0.05
l2fccutoff <- 0.25

p9_dspn1_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P1CKO_CTL/SEURAT_NA_DATA_C_dSPN_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
p9_dspn2_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P2CKO_CTL/SEURAT_NA_DATA_C_dSPN_P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
p9_dspn3_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SUBCLUSTERING_SPNS_wP_wNP/PSEUDOBULK_DEG_CELLTYPE/MASTGLM_BATCH_SEX/P1P2CKO_CTL/SEURAT_NA_DATA_C_dSPN_P1P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)

# p56_dspn1_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P1CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
# p56_dspn2_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
# p56_dspn3_list <- filterDEGs("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/14_INTEGRATION_wHarmony/SUBCLUSTERING_SPN/SUBCLUSTERING_SPN2/PSEUDOBULK_DGE_CELLTYPE/MASTGLM_ALL_COVARIATES/SEURAT_NA_DATA_A_dSPN_P1P2CKO_x_CTL_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)


nadspndeg <-list("01_P09_dSPN_P1CKO_X_CTL_UP" = p9_dspn1_list$DOWN, "02_P09_dSPN_P1CKO_X_CTL_DOWN" = p9_dspn1_list$UP,
                 "03_P09_dSPN_P2CKO_X_CTL_UP" = p9_dspn2_list$DOWN, "04_P09_dSPN_P2CKO_X_CTL_DOWN" = p9_dspn2_list$UP,
                 "05_P09_dSPN_P1P2CKO_X_CTL_UP" = p9_dspn3_list$DOWN, "06_P09_dSPN_P1P2CKO_X_CTL_DOWN" = p9_dspn3_list$UP)


##-------------------------------------------------------
## Perform Super Exact Test
library(SuperExactTest)
library(ggplot2)
library(ggnewscale)
library(rlist)
library(wesanderson)

setout <- list()
setoutAll <- list()
i <- 0

for(mod in names(nadspndeg))
 {
 i <- i + 1
 newGenes <- list.append(myref, mod = nadspndeg[[mod]])
 names(newGenes)[[length(newGenes)]] <- mod
 setres <- supertest(newGenes, n = 6101, degree = 2)  # 6101 is background number of genes for dSPNs across P9 and P56
 setresData <- as.data.frame(summary(setres)$Table)
 setoutAll[[i]] <- setresData
 setresDataSel <- setresData[grep(mod, setresData$Intersections),c(1, 3, 5, 6)]
 setout[[i]] <- setresDataSel
 print(paste(mod, i, sep = " "))
 }

names(setoutAll) <- names(nadspndeg)
names(setout) <- names(nadspndeg)

setoutData <- Reduce("rbind", setout)

save(setout, setoutData, setoutAll, file = "OCT2023_NA_P09_ENRICHMENT_DGE_CLASS_SET.RData")


## REORGANIZE DATA FOR PLOTS
setoutDataTemp <- as.data.frame(matrix(unlist(strsplit(setoutData$Intersections, " & ")), byrow = T, ncol = 2))
colnames(setoutDataTemp) <- c("REF_GENES", "CELLTYPE_DEG")
setoutDataTemp$Intersections <- setoutData$Intersections

setoutData2 <- merge(setoutData, setoutDataTemp, by = "Intersections")
setoutData2$NegLog10 <- -log10(setoutData2$P.value + 10^-10)

write.table(setoutData2, "OCT2023_NA_P09_ENRICHMENT_DGE_CLASS_SET_FOR_PLOT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

setoutData3 <- setoutData2
# setoutData3$NewPval <- setoutData3$P.value
# setoutData3$NewPval[setoutData3$NewPval > 0.05] <- 1
# setoutData3$NegLog10 <- -log10(setoutData3$NewPval + 10^-10)



##-------------------------------------------------------
## HEATMAP
pAll <- ggplot(data = setoutData3) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "01_ASD_1-3",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        # scale_fill_viridis_c(option = "turbo") + 
        # scale_fill_gradientn(colours = c("red3", "gold", "green4")) +
        scale_fill_gradientn(colours = c("grey95", "green4", "gold", "red3")) +
        theme(panel.background = element_blank(), # bg of the panel
              panel.grid.major = element_blank(), # get rid of major grid
              panel.grid.minor = element_blank(), # get rid of minor grid
              axis.title=element_blank(),
              panel.grid=element_blank(),
              axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
              axis.ticks=element_blank(),
              axis.text.y=element_text(size=10)) +
        NULL
pAll <- pAll + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "02_ASD",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "03_FMRP",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "04_ID",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "05_SYN",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        NULL
ggsave(filename = "OCT2023_NA_P09_ENRICHMENT_DGE_CLASS_SET_FOR_PLOT.pdf", plot = pAll, width = 5, height = 6, units = "in", dpi = 300)


# pdf("NA_ENRICHMENT_DGE_CLASS_SET_FOR_BUBBLE_PLOT.pdf", width = 6, height = 6)
# ggscatter(setoutData3, 
#           x = "CELLTYPE_DEG",
#           y = "REF_GENES",
#           size = "FE",
#           color = "NegLog10",
#           alpha = 1,
#           xlab = "",
#           ylab = "",) +
#           theme_minimal() +
#           rotate_x_text(angle = 45)+
#           coord_flip()+
#           scale_size(range = c(2, 12))+ 
#           gradient_color(c("grey95", "green4", "gold", "red3"))
# dev.off()


pdf("OCT2023_NA_P09_ENRICHMENT_DGE_CLASS_SET_FOR_BUBBLE_PLOT.pdf", width = 8, height = 6)
ggscatter(setoutData3, 
          x = "CELLTYPE_DEG",
          y = "REF_GENES",
          size = "FE",
          color = "NegLog10",
          alpha = 1,
          xlab = "",
          ylab = "",) +
          theme_minimal() +
          rotate_x_text(angle = 45)+
          coord_flip()+
          scale_size(range = c(2, 12))+ 
          # gradient_color(c("grey95", "green4", "gold", "red3"))
          gradient_color(c("grey95", "blue2"))
dev.off()



##-------------------------------------------------------
## END
##-------------------------------------------------------

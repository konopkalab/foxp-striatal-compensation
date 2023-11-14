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


## https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html


##-------------------------------------------------------
## PLOTS
load("SEURAT_NA_P56_INTERGRATION_HARMONY.RData")
# seuObjFilt

# Perhaps 8, 11, 47, and 51 are progenitors
# SPNs: 1, 3, 14, 17, 18, 20, 22, 28, 32, 35, 41
# CTX (5, 6, 7, 9, 10, 25, 26, 31, 38, 46, 52, 57, 58

seuObjFilt$CellType <- gsub("^47$", "A_PROG", seuObjFilt$RNA_snn_res.1.2)
seuObjFilt$CellType <- gsub("^8$|^11$|^51$", "B_NEURONAL_PROG", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^1$|^28$|^35$", "C_dSPN", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^3$|^20$|^22$|^41$|^32$", "D_iSPN", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^14$|^17$|^18$", "E_eSPN", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^0$|^45$|^34$|^39$|^53$|^54$", "F_ASTRO", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^2$|^12$|^36$|^49$", "G_OLIG", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^15$|^50$|^55$", "H_MICRO", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^13$|^24$|^18$", "I_ENDO_EPEN", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^4$|^16$|^19$|^21$|^23$|^29$|^30$|^37$|^43$|^44$", "J_INTER", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^5$|^6$|^7$|^9$|^10$|^25$|^26$|^31$|^38$|^46$|^52$|^57$|^58$", "K_CORTICAL", seuObjFilt$CellType)
seuObjFilt$CellType <- gsub("^27$|^33$|^40$|^42$|^48$|^56$", "L_UNDET", seuObjFilt$CellType)


save(seuObjFilt, file = "SEURAT_NA_P56_INTERGRATION_HARMONY_ANNOTATED.RData")




Idents(seuObjFilt) <- "CellType"

seuObjSel <- subset(seuObjFilt, idents = "L_UNDET", invert = TRUE)

plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE, colors_use = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue"))
ggsave(filename = "SEURAT_NA_P56_UMAP_CELLTYPE_JAN2023.pdf", plot = plotCluUMAP1, width = 8, height = 6, units = "in", dpi = 300)

plotCluUMAP2 <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 2, raster = TRUE, colors_use = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue"))
ggsave(filename = "SEURAT_NA_P56_UMAP_CELLTYPE_JAN2023_LABELLED.pdf", plot = plotCluUMAP2, width = 10, height = 8, units = "in", dpi = 300)

# "A_PROG" <- "gold2"
# "B_NEURONAL_PROG" <- "gold4"
# "C_dSPN" <- "blue2"
# "D_iSPN" <- "red3"
# "E_eSPN" <- "green3"
# "F_ASTRO" <- "peru"
# "G_OLIG" <- "palegreen"
# "H_MICRO" <- "khaki3"
# "I_ENDO_EPEN" <- "tomato2"
# "J_INTER" <- "violet"
# "K_CORTICAL" <- "wheat2"

vpl3 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = c("nCount_RNA", "nFeature_RNA", "pMito_RNA"), group.by = "CellType", pt.size = 0, x_lab_rotate = TRUE)
ggsave(filename = "SEURAT_NA_P56_VIOLINPLOT_CELLTYPE_QC.pdf", plot = vpl3 , width = 9, height = 9, units = "in", dpi = 150, useDingbats=FALSE)

# vpl4 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = c("nCount_RNA", "nFeature_RNA", "pMito_RNA"), group.by = "CellTypeCluster", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "SEURAT_NA_P09_VIOLINPLOT_CELLTYPE_CLUSTER_QC.pdf", plot = vpl4 , width = 12, height = 9, units = "in", dpi = 150, useDingbats=FALSE)


## BARPLOT SHOWING CELLS PER CELLTYPE PER GENOTYPE
cellsPerCluster <- as.data.frame.matrix(table(seuObjSel@meta.data$CellType, seuObjSel@meta.data$Genotype))
# cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
cellsPerCluster$CellType <- row.names(cellsPerCluster)
# cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
write.table(cellsPerCluster, "SEURAT_NA_P56_CELLS_BY_CELLTYPE_GENOTYPE.txt", row.names = F, col.names = T, quote = F, sep = "\t")

cellsPerCluster2 <- melt(cellsPerCluster)
colnames(cellsPerCluster2) <- c("CELLTYPE", "GENOTYPE", "CELLS")
p7 <- ggplot(cellsPerCluster2) +
        geom_bar(aes(x = CELLTYPE, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_NA_P56_CELLS_BY_CELLTYPE_GENOTYPE.pdf", plot = p7, width = 9, height = 4.5, units = "in", dpi = 150)









##-------------------------------------------------------
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


load("SEURAT_NA_P56_INTERGRATION_HARMONY_ANNOTATED.RData")
# seuObjFilt

Idents(seuObjFilt) <- "CellType"

seuObjSel <- subset(seuObjFilt, idents = "L_UNDET", invert = TRUE)

DefaultAssay(seuObjSel) <- "RNA"
seuObjSel <- NormalizeData(seuObjSel)


mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Slc17a7", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")

vpl1 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = mygenes, x_lab_rotate = TRUE, group.by = "CellType")
ggsave(filename = "SEURAT_NA_P56_VIOLINPLOT_GENES_CELLTYPE.pdf", plot = vpl1, width = 8, height = 12, units = "in", dpi = 300)

dpl1 <- DotPlot_scCustom(seurat_object = seuObjSel, features = mygenes, x_lab_rotate = TRUE, group.by = "CellType")
ggsave(filename = "SEURAT_NA_P56_DOTPLOT_GENES_CELLTYPE.pdf", plot = dpl1, width = 8, height = 4, units = "in", dpi = 300)





# ##-------------------------------------------------------
# ## LOAD SEURAT OBJECTS FOR EACH GENOTYPE
# print("P56 CTL")
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/12_CTL_woAmbient_woDoublets/NA_SEURAT_CTL_P56.RData")
# # seuObj
# ctl.seuObj <- seuObj
# rm(seuObj)
# dim(ctl.seuObj)
# # 21967 70907

# print("P56 P1-CKO")
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/12_P1CKO_woAmbient_woDoublets/NA_SEURAT_P1CKO_P56.RData")
# # seuObj
# p1cko.seuObj <- seuObj
# rm(seuObj)
# dim(p1cko.seuObj)
# # 21967 55192

# print("P56 P2CKO")
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/12_P2CKO_woAmbient_woDoublets_Again2/NA_SEURAT_P2CKO_P56.RData")
# # seuObj
# p2cko.seuObj <- seuObj
# rm(seuObj)
# dim(p2cko.seuObj)
# # 21967 44192

# print("P56 P1P2CKO")
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/12_P1P2CKO_woAmbient_woDoublets/NA_SEURAT_P1P2CKO_P56.RData")
# # seuObj
# p1p2cko.seuObj <- seuObj
# rm(seuObj)
# dim(p1p2cko.seuObj)
# # 21967 67734


# ##-------------------------------------------------------
# ## FILTER SEURAT OBJECTS
# nuhi <- 10000
# pmhi <- 1.5
# # Lower cut-off for genes =(0.2*ncol(count_mat))/100 

# ctl.seuObjFilt <- subset(ctl.seuObj, subset = nFeature_RNA > (0.2 * ncol(ctl.seuObj)) / 100 & nCount_RNA < nuhi & pMito_RNA < pmhi)
# print(dim(ctl.seuObjFilt))
# # 21967 65405

# p1cko.seuObjFilt <- subset(p1cko.seuObj, subset = nFeature_RNA > (0.2 * ncol(p1cko.seuObj)) / 100 & nCount_RNA < nuhi & pMito_RNA < pmhi)
# print(dim(p1cko.seuObjFilt))
# # 21967 53727

# p2cko.seuObjFilt <- subset(p2cko.seuObj, subset = nFeature_RNA > (0.2 * ncol(p2cko.seuObj)) / 100 & nCount_RNA < nuhi & pMito_RNA < pmhi)
# print(dim(p2cko.seuObjFilt))
# # 21967 42277

# p1p2cko.seuObjFilt <- subset(p1p2cko.seuObj, subset = nFeature_RNA > (0.2 * ncol(p1p2cko.seuObj)) / 100 & nCount_RNA < nuhi & pMito_RNA < pmhi)
# print(dim(p1p2cko.seuObjFilt))
# # 21967 47628


# ##-------------------------------------------------------
# ## MERGE SEURAT OBJECTS
# seuObjFilt <- merge(x = ctl.seuObjFilt, y = c(p1cko.seuObjFilt, p2cko.seuObjFilt, p1p2cko.seuObjFilt))

# seuObjFilt$Genotype <- seuObjFilt$Age
# seuObjFilt$Age <- "P56"

# table(seuObjFilt$Genotype)
# #     CTL   P1CKO P1P2CKO   P2CKO 
# #   65405   53727   47628   42277

# table(seuObjFilt$GenoAgeSample)
# #     P56_CTL_NA03     P56_CTL_NA05     P56_CTL_NA07   P56_P1CKO_NA01 
# #            23191            23687            18527            18972 
# #   P56_P1CKO_NA06   P56_P1CKO_NA10 P56_P1P2CKO_NA02 P56_P1P2CKO_NA04 
# #            25901             8854             6714            25077 
# # P56_P1P2CKO_NA45   P56_P2CKO_NA11   P56_P2CKO_NA25   P56_P2CKO_NA41 
# #            15837             1893            13744              611 
# #   P56_P2CKO_NA44 
# #            26029


# ##-------------------------------------------------------
# ## PREPROCESS
# seuObjFilt <- NormalizeData(seuObjFilt) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE, npcs = 100)


# ##-------------------------------------------------------
# ## JACKSTRAW
# ## Compute and Score Jackstraw
# seuObjFilt <- JackStraw(seuObjFilt, num.replicate = 100, dims = 100)
# seuObjFilt <- ScoreJackStraw(seuObjFilt, dims = 1:100)

# ## Plot PCA loadings and Jackstraw scores
# p6 <- ElbowPlot(seuObjFilt, ndims = 100)
# p7 <- JackStrawPlot(seuObjFilt, dims = 1:100)

# ggsave(filename = "SEURAT_NA_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
# ggsave(filename = "SEURAT_NA_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

# ## IDENTIFY PCS
# pcScores <- seuObjFilt@reductions$pca@jackstraw$overall.p.values
# pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
# selpcs <- min(pcScoresNoSign[,1]) - 1
# print(selpcs)
# # 

# ##-------------------------------------------------------
# ## HARMONY
# seuObjFilt <- RunHarmony(object = seuObjFilt, group.by.vars = c("LibPrep", "SeqBatch", "Sex"), reduction = "pca", assay.use = "RNA", max.iter.harmony = 100) #, project.dim = FALSE)


# ##-------------------------------------------------------
# ## CLUSTERING
# seuObjFilt <- FindNeighbors(seuObjFilt, dims = 1:selpcs, reduction = "harmony")
# seuObjFilt <- FindClusters(seuObjFilt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), reduction = "harmony")
# seuObjFilt <- RunUMAP(seuObjFilt, dims = 1:selpcs, reduction = "harmony")


# ##-------------------------------------------------------
# ## SAVE RDATA
# save(seuObjFilt, file = "SEURAT_NA_P56_INTERGRATION_HARMONY.RData")


# ##-------------------------------------------------------
# ## PLOTS
# seutree <- clustree(seuObjFilt, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
# ggsave(filename = "SEURAT_NA_P56_INTERGRATION_HARMONY_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


# ## UMAP FOR ALL RESOLUTIONS
# ResolutionList <- grep("RNA_snn_res", colnames(seuObjFilt@meta.data), value = TRUE)

# DefaultAssay(seuObjFilt) <- "RNA"

# for (Resolution in ResolutionList)
#     {
#     print(paste("====> RESOLUTION ", Resolution, sep = ""))

#     pdf(paste0("SEURAT_NK_P56_INTEGRATE_UMAP_RES_", Resolution, ".pdf"), width = 7, height = 6)
#     g <- DimPlot(object = seuObjFilt, label = TRUE, reduction = "umap", group.by = Resolution)
#     print(g)
#     dev.off()

#     pdf(paste0("SEURAT_NK_P56_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width = 9, height = 3)
#     v <- VlnPlot(object = seuObjFilt, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
#     print(v)
#     dev.off()
#     }


# # mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Slc17a7", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")
# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cspg4", "Mag", "Cx3cr1", "Flt1", "Slc17a7", "Chat", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b", "Oprm1", "Isl1", "Pdyn", "Lypd1", "Nnat", "Ebf1", "Epha4", "Mef2c")


# fpl1 <- FeaturePlot_scCustom(seurat_object = seuObjFilt, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 16, height = 24, units = "in", dpi = 150)

# fpl2 <- FeaturePlot_scCustom(seurat_object = seuObjFilt, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 16, height = 24, units = "in", dpi = 150)



# plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjFilt, group.by = "Genotype", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_GENOTYPE.pdf", plot = plotCluUMAP1, width = 10, height = 8, units = "in", dpi = 300)

# plotCluUMAP1b <- DimPlot_scCustom(seurat_object = seuObjFilt, group.by = "Genotype", split.by = "Genotype", num_columns = 4, pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_GENOTYPE_FACET.pdf", plot = plotCluUMAP1b, width = 32, height = 8, units = "in", dpi = 300)

# plotCluUMAP1c <- DimPlot_scCustom(seurat_object = seuObjFilt, group.by = "GenoAgeSample", split.by = "GenoAgeSample", num_columns = 4, pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_GENOAGESAMPLE_FACET.pdf", plot = plotCluUMAP1c, width = 25, height = 20, units = "in", dpi = 300)

# plotCluUMAP1d <- DimPlot_scCustom(seurat_object = seuObjFilt, group.by = "GenoAgeSample", pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_GENOAGESAMPLE.pdf", plot = plotCluUMAP1d, width = 8, height = 6, units = "in", dpi = 300)


# ## RES 1.2
# ## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
# cellsPerCluster <- as.data.frame.matrix(table(seuObjFilt@meta.data$RNA_snn_res.1.2, seuObjFilt@meta.data$Genotype))
# cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
# cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
# write.table(cellsPerCluster, paste("NA_SEURAT_GENOTYPE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

# cellsPerCluster2 <- melt(cellsPerCluster)
# colnames(cellsPerCluster2) <- c("CLUSTER", "GENOTYPE", "CELLS")
# p7 <- ggplot(cellsPerCluster2) +
#         geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") + 
#         scale_y_continuous(labels = scales::percent_format()) + 
#         theme_classic() + 
#         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#         NULL
# ggsave(filename = paste("NA_SEURAT_GENOTYPE_PER_CLUSTER_1.2.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)



# ## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
# cellsPerCluster <- as.data.frame.matrix(table(seuObjFilt@meta.data$RNA_snn_res.1.2, seuObjFilt@meta.data$GenoAgeSample))
# cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
# cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
# write.table(cellsPerCluster, paste("NA_SEURAT_GENOAGESAMPLE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

# cellsPerCluster2 <- melt(cellsPerCluster)
# colnames(cellsPerCluster2) <- c("CLUSTER", "GENOAGESAMPLE", "CELLS")
# p7 <- ggplot(cellsPerCluster2) +
#         geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOAGESAMPLE), stat = "identity", position = "fill") + 
#         scale_y_continuous(labels = scales::percent_format()) + 
#         theme_classic() + 
#         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#         NULL
# ggsave(filename = paste("NA_SEURAT_GENOAGESAMPLE_PER_CLUSTER_1.2.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)



# my_levels <- seq(0, max(as.numeric(seuObjFilt$RNA_snn_res.1.2)) - 1)
# seuObjFilt$RNA_snn_res.1.2 <- factor(x = seuObjFilt$RNA_snn_res.1.2, levels = my_levels)

# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cspg4", "Mag", "Cx3cr1", "Flt1", "Slc17a7", "Chat", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b", "Oprm1", "Isl1", "Pdyn", "Lypd1", "Nnat", "Ebf1", "Epha4", "Mef2c")

# vpl2 <- Stacked_VlnPlot(seurat_object = seuObjFilt, features = mygenes, group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_1.2_B.pdf", plot = vpl2 , width = 12, height = 18, units = "in", dpi = 150, useDingbats=FALSE)

# vpl3 <- Stacked_VlnPlot(seurat_object = seuObjFilt, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_1.2_C.pdf", plot = vpl3 , width = 12, height = 4, units = "in", dpi = 150, useDingbats=FALSE)

# dpl2 <- DotPlot_scCustom(seurat_object = seuObjFilt, features = rev(mygenes), colors_use = viridis_plasma_dark_high, group.by = "RNA_snn_res.1.2", x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_DotPlot_1.2.pdf", plot = dpl2, width = 12, height = 16, units = "in", dpi = 150, useDingbats=FALSE)


# Idents(seuObjFilt) <- "RNA_snn_res.1.2"

# str.markers <- FindAllMarkers(seuObjFilt)  #, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

# write.table(str.markers, "SEURAT_NA_P56_INTERGRATION_HARMONY_CLU_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# save(str.markers, file = "SEURAT_NA_P56_INTERGRATION_HARMONY_CLU_MARKERS.RData")



#-------------------------------------------------------
# END
#-------------------------------------------------------








# ##-------------------------------------------------------
# ## SEURAT ANALYSIS | INTEGRATE DATA | PLOTS
# ##-------------------------------------------------------

# rm(list = ls())
# library(Seurat)
# library(dplyr)
# library(Matrix)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
# library(ggrepel)
# library(reticulate)
# library(clustree)

# load("SEURAT_NA_P09_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# # org.integrated, org.anchors


# plotCluUMAP1 <- DimPlot(object = org.integrated, group.by = "Source", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_SOURCE.pdf", plot = plotCluUMAP1, width = 10, height = 8, units = "in", dpi = 300)

# plotCluUMAP1b <- DimPlot(object = org.integrated, group.by = "Source", split.by = "Source", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_SOURCE_FACET.pdf", plot = plotCluUMAP1b, width = 32, height = 8, units = "in", dpi = 300)

# plotCluUMAP2 <- DimPlot(object = org.integrated, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_CELLTYPES.pdf", plot = plotCluUMAP2, width = 15, height = 12, units = "in", dpi = 300)

# plotCluUMAP3 <- DimPlot(object = org.integrated, group.by = "integrated_snn_res.0.8", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_CLUSTERS_RES_0.8.pdf", plot = plotCluUMAP3, width = 6, height = 6, units = "in", dpi = 300, useDingbats = FALSE)


# clusterByCellType <- as.data.frame.matrix(table(org.integrated@meta.data$integrated_snn_res.0.8, org.integrated@meta.data$CellType))
# write.table(clusterByCellType, "SEURAT_ORGANOIDS_CLUSTERS_BY_CELLTYPE.TSV", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# clusterBySource <- as.data.frame.matrix(table(org.integrated@meta.data$integrated_snn_res.0.8, org.integrated@meta.data$Source))
# write.table(clusterBySource, "SEURAT_ORGANOIDS_CLUSTERS_BY_SOURCE.TSV", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# clusterInfo <- merge(clusterByCellType, clusterBySource, by = "row.names")
# write.table(clusterInfo, "SEURAT_ORGANOIDS_CLUSTERS_BY_CELLTYPE_SOURCE.TSV", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


# DefaultAssay(org.integrated) <- "RNA"
# Idents(org.integrated) <- "integrated_snn_res.0.8"

# mygenes4 <- c("FOXP1", "MKI67", "SOX2", "ASCL1", "SOX4", "SOX11", "PAX6", "HES5", "PARD3", "EOMES", "NEUROG2", "BTG2", "NEUROD2", "NEUROD6", "TUBB3", "RBFOX3", "NCAM2", "BCL11B")

# my_levels <- seq(0, max(as.numeric(org.integrated$integrated_snn_res.0.8)) - 1)
# org.integrated@active.ident <- factor(x = org.integrated@active.ident, levels = my_levels)

# vpl1 <- VlnPlot(org.integrated, features = mygenes4, ncol = 3, pt.size = 0, assay = "RNA", combine = TRUE, sort = FALSE)
# ggsave(filename = "SEURAT_ORG_DATA_ViolinPlot_NoLog_0.8_3.pdf", plot = vpl1 , width = 32, height = 12, units = "in", dpi = 150, useDingbats=FALSE)

# dpl1 <- DotPlot(object = org.integrated, features = rev(mygenes4), dot.scale = 12, dot.min = 0, cols = c("white", "blue"), col.min = 0, col.max = 3, assay = "RNA") + theme(axis.text.x = element_text(angle = 90))
# ggsave(filename = "SEURAT_ORG_DATA_DotPlot_0.8_3.pdf", plot = dpl1, width = 12, height = 18, units = "in", dpi = 150, useDingbats=FALSE)


# # mygenes5 <- c("SOX2", "PAX6", "HES5", "EOMES", "NEUROG2", "BTG2", "NEUROD2", "TUBB3", "NEUROD6")

# # vpl2 <- VlnPlot(org.integrated, features = mygenes5, group.by = "integrated_snn_res.1.8", ncol = 3, pt.size = 0, assay = "RNA", combine = TRUE, sort = TRUE)
# # ggsave(filename = "SEURAT_ORG_DATA_ViolinPlot_NoLog_1.8_4.pdf", plot = vpl2 , width = 32, height = 6, units = "in", dpi = 150, useDingbats=FALSE)

# # dpl2 <- DotPlot(object = org.integrated, features = rev(mygenes5), dot.scale = 10, dot.min = 0, cols = c("white", "blue"), col.min = 0, col.max = 3, group.by = "integrated_snn_res.1.8", assay = "RNA") + theme(axis.text.x = element_text(angle = 90))
# # ggsave(filename = "SEURAT_ORG_DATA_DotPlot_1.8_4.pdf", plot = dpl2, width = 8, height = 24, units = "in", dpi = 150, useDingbats=FALSE)



##-------------------------------------------------------
## SEURAT ANALYSIS | INTEGRATE DATA | END
##-------------------------------------------------------



# celltypeByCluster <- as.data.frame.matrix(table(org.integrated@meta.data$CellType, org.integrated@meta.data$integrated_snn_res.0.8))
# write.table(celltypeByCluster, "SEURAT_ORGANOIDS_DATA_FILT_NORM_PCA_CLU_INTEGRATED_CELLTYPE_BY_CLUSTERS.TSV", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# celltypeBySource <- as.data.frame.matrix(table(org.integrated@meta.data$CellType, org.integrated@meta.data$Source))
# write.table(celltypeBySource, "SEURAT_ORGANOIDS_DATA_FILT_NORM_PCA_CLU_INTEGRATED_CELLTYPE_BY_SOURCE.TSV", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")




# metaDataTemp <- as.data.frame(org.integrated@meta.data)
# umapDataTemp <- as.data.frame(Embeddings(object = org.integrated, reduction = "umap"))
# metaData <- merge(umapDataTemp, metaDataTemp, by = "row.names")
# row.names(metaData) <- metaData$Row.names
# metaData$Row.names <- NULL


# emptyCells <- row.names(metaData[metaData$CellType == "",])
# org.integrated2 <- subset(org.integrated, cells = emptyCells, invert = TRUE)

# plotCluUMAP2 <- DimPlot(object = org.integrated2, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_CELLTYPES.pdf", plot = plotCluUMAP2, width = 16, height = 12, units = "in", dpi = 300)

# plotCluUMAP3 <- DimPlot(object = org.integrated, group.by = "integrated_snn_res.0.8", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_UMAP_CLUSTERS.pdf", plot = plotCluUMAP3, width = 16, height = 12, units = "in", dpi = 300)



# DefaultAssay(org.integrated) <- "RNA"

# mygenes4 <- c("FOXP1", "MKI67", "SOX2", "ASCL1", "SOX4", "SOX11", "PAX6")

# fpl1 <- FeaturePlot(object = org.integrated, features = mygenes4, cols = c("gray90", "blue"), reduction = "umap", ncol = 3, pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_FeaturePlot_orderT_3.pdf", plot = fpl1, width = 14, height = 12, units = "in", dpi = 300, useDingbats=FALSE)

# fpl2 <- FeaturePlot(object = org.integrated, features = mygenes4, cols = c("gray90", "blue"), reduction = "umap", ncol = 3, pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE)
# ggsave(filename = "SEURAT_ORG_DATA_FeaturePlot_orderF_3.pdf", plot = fpl2, width = 14, height = 12, units = "in", dpi = 300, useDingbats=FALSE)

# vpl1 <- VlnPlot(org.integrated, features = mygenes4, group.by = "integrated_snn_res.0.8", ncol = 3, pt.size = 0, assay = "RNA", combine = TRUE)
# ggsave(filename = "SEURAT_ORG_DATA_ViolinPlot_NoLog_0.8_3.pdf", plot = vpl1 , width = 20, height = 6, units = "in", dpi = 150, useDingbats=FALSE)

# dpl1 <- DotPlot(object = org.integrated, features = rev(mygenes4), dot.scale = 10, dot.min = 0, cols = c("white", "blue"), col.min = 0, col.max = 3, group.by = "integrated_snn_res.0.8", assay = "RNA") + theme(axis.text.x = element_text(angle = 90))
# ggsave(filename = "SEURAT_ORG_DATA_DotPlot_0.8_3.pdf", plot = dpl1, width = 6, height = 12, units = "in", dpi = 150, useDingbats=FALSE)






# meta_oRG <- metaData[grepl("gliogenic", metaData$cl_FullLineage),]

# table(meta_oRG$integrated_snn_res.0.8)
# sort(table(meta_oRG$integrated_snn_res.0.8), decreasing = TRUE)
#   12   24    0    7    8    5   15   16    3   13   22   23   26    1   10   11 
# 1497  229  117   63   12    8    2    2    2    1    1    1    1    0    0    0 
#   14   17   18   19    2   20   21   25   27   28   29   30   31   32   33   34 
#    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
#   35    4    6    9 
#    0    0    0    0







# clusterByCellType <- as.data.frame.matrix(table(org.integrated2@meta.data$integrated_snn_res.0.8, org.integrated2@meta.data$CellType))
# sort(clusterByCellType[row.names(clusterByCellType) == 6,], decreasing = TRUE)

# metaData$NewAnnotation <- metaData$CellType
# metaData$NewAnnotation[metaData$integrated_snn_res.0.8 %in% c(0, 5, 6)] <- "RG"
# metaData$NewAnnotation[metaData$integrated_snn_res.0.8 %in% c(2)] <- "IPC"
# metaData$NewAnnotation[metaData$integrated_snn_res.0.8 %in% c(1, 3, 4)] <- "NE"

##-------------------------------------------------------
## SEURAT ANALYSIS | INTEGRATE DATA | END
##-------------------------------------------------------









# ## marker genes
# markerGenes <- c("Aqp4", "Olig1", "Cx3cr1", "Flt2", "Slc17a7", "Chat", "Npy", "Sox4", "Mki67", "Ascl1", "Dlx2", "Sox11", "Sp9", "Ppp1r1b", "Drd1", "Drd2", "Tac1", "Penk")
# # goi <- intersect(row.names(org.integrated@assays$RNA), markerGenes)

# fpl1 <- FeaturePlot(object = org.integrated, features = markerGenes, cols = c("gray90", "darkblue"), reduction = "umap", ncol = 4, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = TRUE)
# ggsave(filename = paste("SEURAT_STR_DATA_UMAP_GENES_1.pdf", sep = ""), plot = fpl1, width = 14, height = 15, units = "in", dpi = 150)

# fpl2 <- FeaturePlot(object = org.integrated, features = markerGenes, cols = c("gray90", "darkblue"), reduction = "umap", ncol = 4, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = FALSE)
# ggsave(filename = paste("SEURAT_STR_DATA_UMAP_GENES_2.pdf", sep = ""), plot = fpl2, width = 14, height = 15, units = "in", dpi = 150)


# save(org.integrated, org.anchors, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")




# # # load("SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# # # org.integrated, org.anchors

# ## UMAP PLOT COLORED BY CLUSTERS & SOURCE
# # plotCluUMAP1 <- DimPlot(object = org.integrated, pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# # ggsave(filename = "SEURAT_STR_DATA_UMAP_CLUSTERS.pdf", plot = plotCluUMAP1, width = 13, height = 12, units = "in", dpi = 300)

# # plotCluUMAP2 <- DimPlot(object = org.integrated, pt.size = 0.1, reduction = "umap", label = FALSE, label.size = 3, repel = TRUE) + facet_wrap(~ident, nrow = 7)
# # ggsave(filename = "SEURAT_STR_DATA_UMAP_CLUSTERS_FACET.pdf", plot = plotCluUMAP2, width = 13, height = 12, units = "in", dpi = 300)

# plotCluUMAP3 <- DimPlot(object = org.integrated, pt.size = 0.1, group.by = "orig.ident", reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_STR_DATA_UMAP_SOURCES.pdf", plot = plotCluUMAP3, width = 13, height = 12, units = "in", dpi = 300)

# plotCluUMAP4 <- DimPlot(object = org.integrated, pt.size = 0.1, group.by = "orig.ident", reduction = "umap", label = FALSE, label.size = 3, repel = TRUE) + facet_wrap(~orig.ident, nrow = 4)
# ggsave(filename = "SEURAT_STR_DATA_UMAP_SOURCES_FACET.pdf", plot = plotCluUMAP4, width = 13, height = 12, units = "in", dpi = 300)







# # # load("SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# # # org.integrated, org.anchors

# org.integrated <- FindNeighbors(object = org.integrated, dims = 1:30) %>% FindClusters(resolution = 0.8)

# # UMAP PLOT COLORED BY CLUSTERS & SOURCE
# plotCluUMAP1 <- DimPlot(object = org.integrated, pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# ggsave(filename = "SEURAT_STR_DATA_UMAP_CLUSTERS.pdf", plot = plotCluUMAP1, width = 13, height = 12, units = "in", dpi = 300)

# plotCluUMAP2 <- DimPlot(object = org.integrated, pt.size = 0.1, reduction = "umap", label = FALSE, label.size = 3, repel = TRUE) + facet_wrap(~ident, nrow = 7)
# ggsave(filename = "SEURAT_STR_DATA_UMAP_CLUSTERS_FACET.pdf", plot = plotCluUMAP2, width = 12, height = 12, units = "in", dpi = 300)

# save(org.integrated, org.anchors, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED_2.RData")







# load("SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED_2.RData")

# seuMarkers <- FindAllMarkers(object = org.integrated)

# write.table(seuMarkers, "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED_2_DEG.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# save(seuMarkers, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED_2_DEG.RData")




# load("SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED_2.RData")

# ## marker genes
# markerGenes <- c("Foxp1", "Foxp2")
# goi <- intersect(row.names(org.integrated@assays$RNA), markerGenes)
# fpl <- FeaturePlot(object = org.integrated, features = goi, cols = c("gray90", "darkblue"), reduction = "umap", ncol = 2, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = TRUE)
# ggsave(filename = paste("SEURAT_STR_DATA_UMAP_GENES_2.pdf", sep = ""), plot = fpl, width = 8, height = 4, units = "in", dpi = 150)



# ## transfer data to classify the query cells based on reference data
# ## query dataset = "SAUNDERS"
# ## ref datasets = "ASHLEY", "CHEN", "MUNOZ1", "MUNOZ2", "QUAKE", "QUAKE2019"
# org.query <- org.list[["SAUNDERS"]]
# org.anchors <- FindTransferAnchors(reference = org.integrated, query = org.query, dims = 1:30)
# predictions <- TransferData(anchorset = org.anchors, refdata = org.integrated$orig.ident, dims = 1:30)
# org.query <- AddMetaData(org.query, metadata = predictions)


# ## match prediction
# org.query$prediction.match <- org.query$predicted.id == org.query$orig.ident
# table(org.query$prediction.match)

# save(org.integrated, org.anchors, org.list, org.query, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")


# seuObjAll <- CreateSeuratObject(seuObjAll.newdata) %>% SCTransform(vars.to.regress = c("nCount_RNA", "orig.ident")) %>% RunPCA(npcs = 100) %>% JackStraw(dims = 100) %>% ScoreJackStraw(dims = 1:100)
# # seuObjAll <- CreateSeuratObject(seuObjAll.newdata) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("nCount_RNA", "orig.ident")) %>% RunPCA(npcs = 100) %>% JackStraw(dims = 100) %>% ScoreJackStraw(dims = 1:100)

# pdf("SEURAT_STR_DATA_ElbowPlot.pdf", width = 8, height = 6)
# ElbowPlot(object = seuObjAll, ndims = 100)
# dev.off()

# pdf("SEURAT_STR_DATA_JackstrawPlot.pdf", width = 12, height = 6)
# JackStrawPlot(seuObjAll, dims = 1:100)
# dev.off()

# save(seuObjAll, file = "SEURAT_STR_DATA_FILT_NORM_PCA.RData")


# ##-------------------------------------------------------
# ## SEURAT ANALYSIS | PART 2 CLUSTERING
# ##-------------------------------------------------------
# load("SEURAT_STR_DATA_FILT.RData")
# ## seuObjAll.newdata
# load("SEURAT_STR_DATA_FILT_NORM_PCA.RData")
# # seuObjAll

# # seuObjAll <- CreateSeuratObject(seuObjAll.newdata) %>% SCTransform(vars.to.regress = "nCount_RNA") %>% RunPCA() %>% FindNeighbors(dims = 1:13) %>% RunUMAP(dims = 1:13) %>% FindClusters()
# seuObjAll <- FindNeighbors(object = seuObjAll, dims = 1:12) %>% RunUMAP(dims = 1:12) %>% FindClusters(resolution = 0.8)

# save(seuObjAll, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU.RData")

# pdf("SEURAT_STR_DATA_UMAP.pdf", width = 7, height = 6)
# DimPlot(seuObjAll, label = TRUE) + NoLegend()
# dev.off()


# ##-------------------------------------------------------
# ## SEURAT ANALYSIS | PART 3 CLUSTERING
# ##-------------------------------------------------------
# # load("SEURAT_STR_DATA_FILT_NORM_PCA_CLU.RData")
# seuMarkers <- FindAllMarkers(object = seuObjAll)
# save(seuMarkers, file = "SEURAT_STR_DATA_FILT_NORM_PCA_CLU_DEG.RData")


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------

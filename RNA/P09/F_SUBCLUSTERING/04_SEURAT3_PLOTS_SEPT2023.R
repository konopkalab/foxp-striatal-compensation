##-------------------------------------------------------
## SEURAT ANALYSIS | MODULES & LIBRARIES
##-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

## SCROLL DOWN FOR SEPT 2023 SECTION

# ##-------------------------------------------------------
# ## LOAD LIBRARIES
# rm(list = ls())
# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(ggplot2)
# library(reshape2)
# library(clustree)
# library(harmony)
# library(scCustomize)
# set.seed(10)


# ## https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html



# # ##-------------------------------------------------------
# # ## SUBCLUSTERING
# # load("./../SEURAT_NA_P09_INTERGRATION_HARMONY.RData")
# # # seuObjFilt

# # Idents(seuObjFilt) <- "RNA_snn_res.1.2"
# # DefaultAssay(seuObjFilt) <- "RNA"
# # seuObjFilt <- NormalizeData(seuObjFilt)

# # clusters2keep <- sort(c(12, 7, 33, 10, 1, 2, 4, 9, 0, 5, 8, 19, 26, 11))

# # seuObjSel <- subset(seuObjFilt, idents = clusters2keep, slot = "counts")

# # table(seuObjSel$GenoAgeSample)
# # #     CTL_P09_NA13     CTL_P09_NA14     CTL_P09_NA29   P1CKO_P09_NA32 
# # #             8720            12645            12734            11490 
# # #   P1CKO_P09_NA33   P1CKO_P09_NA42 P1P2CKO_P09_NA15 P1P2CKO_P09_NA20 
# # #             7266            11067            10034             7883 
# # # P1P2CKO_P09_NA24   P2CKO_P09_NA16   P2CKO_P09_NA22   P2CKO_P09_NA43 
# # #            11425             4252             7716            19462

# # table(seuObjSel$Genotype)
# # #     CTL   P1CKO P1P2CKO   P2CKO 
# # #   34099   29823   29342   31430

# # table(seuObjSel$Batch)
# # #     1     2     3     4     5 
# # # 35651 27024 12734 18756 30529

# # table(seuObjSel$Sex)
# # #     F     M 
# # # 87641 37053



# # ##-------------------------------------------------------
# # ## PREPROCESS
# # seuObjSel <- NormalizeData(seuObjSel) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE, npcs = 30)

# # ##-------------------------------------------------------
# # ## HARMONY
# # seuObjSel <- RunHarmony(seuObjSel, group.by.vars = c("Batch", "Sex"))

# # ##-------------------------------------------------------
# # ## CLUSTERING
# # seuObjSel <- RunUMAP(seuObjSel, reduction = "harmony", dims = 1:30)
# # myres <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
# # seuObjSel <- FindNeighbors(seuObjSel, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = myres)

# # ##-------------------------------------------------------
# # ## SAVE RDATA
# # save(seuObjSel, file = "SEURAT_NA_P09_INTERGRATION_HARMONY_SPN.RData")


# # ##-------------------------------------------------------
# # ## PLOTS
# # seutree <- clustree(seuObjSel, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
# # ggsave(filename = "SEURAT_NA_P09_INTERGRATION_HARMONY_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


# # ## UMAP FOR ALL RESOLUTIONS
# # ResolutionList <- grep("RNA_snn_res", colnames(seuObjSel@meta.data), value = TRUE)

# # DefaultAssay(seuObjSel) <- "RNA"

# # for (Resolution in ResolutionList)
# #     {
# #     print(paste("====> RESOLUTION ", Resolution, sep = ""))

# #     pdf(paste0("SEURAT_NK_P09_INTEGRATE_UMAP_RES_", Resolution, ".pdf"), width = 7, height = 6)
# #     g <- DimPlot(object = seuObjSel, label = TRUE, reduction = "umap", group.by = Resolution)
# #     print(g)
# #     dev.off()

# #     pdf(paste0("SEURAT_NK_P09_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width = 9, height = 3)
# #     v <- VlnPlot(object = seuObjSel, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
# #     print(v)
# #     dev.off()
# #     }





# # my_levels <- seq(0, max(as.numeric(seuObjSel$RNA_snn_res.0.8)) - 1)
# # seuObjSel$RNA_snn_res.0.8 <- factor(x = seuObjSel$RNA_snn_res.0.8, levels = my_levels)

# # mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cspg4", "Mag", "Cx3cr1", "Flt1", "Slc17a7", "Chat", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b", "Oprm1", "Isl1", "Pdyn", "Lypd1", "Nnat", "Ebf1", "Epha4", "Mef2c")

# # vpl2 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = mygenes, group.by = "RNA_snn_res.0.8", pt.size = 0, x_lab_rotate = TRUE)
# # ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_0.8_B.pdf", plot = vpl2 , width = 9, height = 18, units = "in", dpi = 150, useDingbats=FALSE)

# # vpl3 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.0.8", pt.size = 0, x_lab_rotate = TRUE)
# # ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_0.8_C.pdf", plot = vpl3 , width = 9, height = 5, units = "in", dpi = 150, useDingbats=FALSE)

# # dpl2 <- DotPlot_scCustom(seurat_object = seuObjSel, features = rev(mygenes), colors_use = viridis_plasma_dark_high, group.by = "RNA_snn_res.0.8", x_lab_rotate = TRUE)
# # ggsave(filename = "NA_SEURAT_PLOT_DotPlot_0.8.pdf", plot = dpl2, width = 12, height = 16, units = "in", dpi = 150, useDingbats=FALSE)


# # fpl1 <- FeaturePlot_scCustom(seurat_object = seuObjSel, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
# # ggsave(filename = "NA_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 24, height = 24, units = "in", dpi = 150)

# # fpl2 <- FeaturePlot_scCustom(seurat_object = seuObjSel, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
# # ggsave(filename = "NA_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 24, height = 24, units = "in", dpi = 150)

# # vpl1 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = mygenes, x_lab_rotate = TRUE, group.by = "RNA_snn_res.0.8")
# # ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_Res1.2.pdf", plot = vpl1, width = 24, height = 24, units = "in", dpi = 150)


# # plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "Genotype", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE)
# # ggsave(filename = "NA_SEURAT_UMAP_GENOTYPE.pdf", plot = plotCluUMAP1, width = 10, height = 8, units = "in", dpi = 300)

# # plotCluUMAP1b <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "Genotype", split.by = "Genotype", num_columns = 4, pt.size = 0.1, reduction = "umap", raster = TRUE)
# # ggsave(filename = "NA_SEURAT_UMAP_GENOTYPE_FACET.pdf", plot = plotCluUMAP1b, width = 32, height = 8, units = "in", dpi = 300)

# # plotCluUMAP1c <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "GenoAgeSample", split.by = "GenoAgeSample", num_columns = 4, pt.size = 0.1, reduction = "umap", raster = TRUE)
# # ggsave(filename = "NA_SEURAT_UMAP_GENOAGESAMPLE_FACET.pdf", plot = plotCluUMAP1c, width = 25, height = 20, units = "in", dpi = 300)

# # plotCluUMAP1d <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "GenoAgeSample", pt.size = 0.1, reduction = "umap", raster = TRUE)
# # ggsave(filename = "NA_SEURAT_UMAP_GENOAGESAMPLE.pdf", plot = plotCluUMAP1d, width = 8, height = 6, units = "in", dpi = 300)



# # ## RES 0.8
# # ## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 0.8
# # cellsPerCluster <- as.data.frame.matrix(table(seuObjSel@meta.data$RNA_snn_res.0.8, seuObjSel@meta.data$Genotype))
# # cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
# # cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
# # write.table(cellsPerCluster, paste("NA_SEURAT_GENOTYPE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

# # cellsPerCluster2 <- melt(cellsPerCluster)
# # colnames(cellsPerCluster2) <- c("CLUSTER", "GENOTYPE", "CELLS")
# # p7 <- ggplot(cellsPerCluster2) +
# #         geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") + 
# #         scale_y_continuous(labels = scales::percent_format()) + 
# #         theme_classic() + 
# #         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
# #         NULL
# # ggsave(filename = paste("NA_SEURAT_GENOTYPE_PER_CLUSTER_0.8.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)




# # ## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
# # cellsPerCluster <- as.data.frame.matrix(table(seuObjSel@meta.data$RNA_snn_res.0.8, seuObjSel@meta.data$GenoAgeSample))
# # cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
# # cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
# # write.table(cellsPerCluster, paste("NA_SEURAT_GENOAGESAMPLE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

# # cellsPerCluster2 <- melt(cellsPerCluster)
# # colnames(cellsPerCluster2) <- c("CLUSTER", "GENOAGESAMPLE", "CELLS")
# # p7 <- ggplot(cellsPerCluster2) +
# #         geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOAGESAMPLE), stat = "identity", position = "fill") + 
# #         scale_y_continuous(labels = scales::percent_format()) + 
# #         theme_classic() + 
# #         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
# #         NULL
# # ggsave(filename = paste("NA_SEURAT_GENOAGESAMPLE_PER_CLUSTER_0.8.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)



# ##-------------------------------------------------------
# ## UPDATE ANNOTATION
# rm(list = ls())
# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(ggplot2)
# library(reshape2)
# library(clustree)
# library(harmony)
# library(scCustomize)
# set.seed(10)


# load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN.RData")
# # seuObjSel

# Idents(seuObjSel) <- "RNA_snn_res.0.8"

# table(seuObjSel@active.ident)
# #     0     1    10    11    12    13    14    15    16    17     2     3     4 
# # 16111 15045  3895  3347  3188  2560  1685  1046   807   643 12364 11461 10855 
# #     5     6     7     8     9 
# #  9748  9715  9011  8032  5181

# seuObjSel$CellType <- gsub("^9$", "A_PROG", seuObjSel$RNA_snn_res.0.8)
# seuObjSel$CellType <- gsub("^6$|^11$", "B_NEURONAL_PROG", seuObjSel$CellType)
# seuObjSel$CellType <- gsub("^0$|^1$|^5$|^8$|^12$", "C_dSPN", seuObjSel$CellType)
# seuObjSel$CellType <- gsub("^2$|^3$|^4$|^7$", "D_iSPN", seuObjSel$CellType)
# seuObjSel$CellType <- gsub("^10$|^13$|^14$", "E_eSPN", seuObjSel$CellType)
# seuObjSel$CellType <- gsub("^15$|^16$|^17$", "F_UNDET", seuObjSel$CellType)

# seuObjSel$CellTypeCluster <- paste(seuObjSel$CellType, seuObjSel$RNA_snn_res.0.8, sep = "_")

# dimp1 <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "CellType", pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE.pdf", plot = dimp1, width = 8, height = 6, units = "in", dpi = 300)

# dimp1b <- DimPlot_scCustom(seurat_object = seuObjSel, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER.pdf", plot = dimp1b, width = 12, height = 6, units = "in", dpi = 300)

# vpl3 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = c("nCount_RNA", "nFeature_RNA"), group.by = "CellType", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_0.8.pdf", plot = vpl3 , width = 9, height = 5, units = "in", dpi = 150, useDingbats=FALSE)

# Idents(seuObjSel) <- "CellType"

# seuObjFilt <- subset(seuObjSel, idents = "F_UNDET", invert = TRUE)

# dimp2 <- DimPlot(seuObjFilt, group.by = "CellType", pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLEANED.pdf", plot = dimp2, width = 8, height = 6, units = "in", dpi = 300)

# dimp2b <- DimPlot(seuObjFilt, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER_CLEANED.pdf", plot = dimp2b, width = 8, height = 6, units = "in", dpi = 300)

# dimp3 <- DimPlot(seuObjFilt, group.by = "CellType", pt.size = 0.1, reduction = "umap", raster = FALSE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLEANED2.pdf", plot = dimp2, width = 8, height = 6, units = "in", dpi = 300)

# dimp3b <- DimPlot(seuObjFilt, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", raster = FALSE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER_CLEANED2.pdf", plot = dimp2b, width = 8, height = 6, units = "in", dpi = 300)

# vpl4 <- Stacked_VlnPlot(seurat_object = seuObjFilt, features = c("nCount_RNA", "nFeature_RNA"), group.by = "CellType", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_0.8_CLEANED.pdf", plot = vpl4 , width = 6, height = 5, units = "in", dpi = 150, useDingbats=FALSE)

# vpl5 <- Stacked_VlnPlot(seurat_object = seuObjFilt, features = c("nCount_RNA", "nFeature_RNA"), group.by = "CellTypeCluster", pt.size = 0, x_lab_rotate = TRUE)
# ggsave(filename = "NA_SEURAT_PLOT_ViolinPlot_0.8_CLEANED_Cluster.pdf", plot = vpl5 , width = 12, height = 5, units = "in", dpi = 150, useDingbats=FALSE)

# save(seuObjSel, file = "SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_ANNOTATED.RData")
# save(seuObjFilt, file = "SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")



# rm(list = ls())
# load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt

# Idents(seuObjFilt) <- "CellType"

# p09_spn_celltype_sample <- as.data.frame.matrix(table(seuObjFilt$CellType, seuObjFilt$GenoAgeSample))
# write.table(p09_spn_celltype_sample, "P09_SPN_CELLTYPE_BY_SAMPLE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# p09_spn_celltypecluster_sample <- as.data.frame.matrix(table(seuObjFilt$CellTypeCluster, seuObjFilt$GenoAgeSample))
# write.table(p09_spn_celltypecluster_sample, "P09_SPN_CELLTYPE_CLUSTER_BY_SAMPLE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


# plotCluUMAP1 <- DimPlot(seuObjFilt, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE.pdf", plot = plotCluUMAP1, width = 6.5, height = 6, units = "in", dpi = 300)

# plotCluUMAP2 <- DimPlot(seuObjFilt, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER_LABELLED.pdf", plot = plotCluUMAP2, width = 10, height = 8, units = "in", dpi = 300)

# plotCluUMAP2b <- DimPlot(seuObjFilt, group.by = "CellTypeCluster", pt.size = 0.1, reduction = "umap", label = FALSE, label.size = 3, raster = TRUE)
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_CLUSTER.pdf", plot = plotCluUMAP2b, width = 10, height = 8, units = "in", dpi = 300)



# ##-------------------------------------------------------
# ## LOAD LIBRARIES
# rm(list = ls())
# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(ggplot2)
# library(reshape2)
# library(clustree)
# library(harmony)
# library(scCustomize)
# set.seed(10)


# rm(list = ls())
# load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt

# Idents(seuObjFilt) <- "CellType"

# plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObjFilt, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE, colors_use = c("gold2", "gold4", "blue2", "red3", "green3"))
# ggsave(filename = "NA_SEURAT_UMAP_CELLTYPE_11072022.pdf", plot = plotCluUMAP1, width = 7.2, height = 6, units = "in", dpi = 300)







# ##-------------------------------------------------------
# ## FETCH THE LIST OF GENES EXPRESSED ACROSS ALL SPNs
# ## without P and NP
# rm(list = ls())

# load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt

# seuObjFiltSel <- subset(seuObjFilt, ident = c("C_dSPN", "D_iSPN", "E_eSPN"))

# gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

# gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

# expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

# exp_genes_sorted <- names(expressed.genes)

# write.table(exp_genes_sorted, "NA_P09_ALL_SPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
# write.table(gene.expression.sorted2, "NA_P09_ALL_SPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")




# ##-------------------------------------------------------
# ## FETCH THE LIST OF GENES EXPRESSED ACROSS dSPNs ONLY
# rm(list = ls())

# load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# # seuObjFilt

# seuObjFiltSel <- subset(seuObjFilt, ident = c("C_dSPN"))

# gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

# gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

# expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

# exp_genes_sorted <- names(expressed.genes)

# write.table(exp_genes_sorted, "NA_P09_dSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


# gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
# write.table(gene.expression.sorted2, "NA_P09_dSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")







##-------------------------------------------------------
## SEPT 2023
##-------------------------------------------------------
## FETCH THE LIST OF GENES EXPRESSED ACROSS iSPNs ONLY
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



load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

seuObjFiltSel <- subset(seuObjFilt, ident = c("D_iSPN"))

gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

exp_genes_sorted <- names(expressed.genes)

write.table(exp_genes_sorted, "NA_P09_iSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
write.table(gene.expression.sorted2, "NA_P09_iSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





rm(list = ls())
load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

seuObjFiltSel <- subset(seuObjFilt, ident = c("E_eSPN"))

gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

exp_genes_sorted <- names(expressed.genes)

write.table(exp_genes_sorted, "NA_P09_eSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
write.table(gene.expression.sorted2, "NA_P09_eSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





rm(list = ls())
load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

seuObjFiltSel <- subset(seuObjFilt, ident = c("A_PROG"))

gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

exp_genes_sorted <- names(expressed.genes)

write.table(exp_genes_sorted, "NA_P09_PROG_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
write.table(gene.expression.sorted2, "NA_P09_PROG_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





rm(list = ls())
load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

seuObjFiltSel <- subset(seuObjFilt, ident = c("B_NEURONAL_PROG"))

gene.expression <- apply(seuObjFiltSel@assays$RNA@data, 1, mean)

gene.expression.sorted <- sort(gene.expression, decreasing = TRUE)

expressed.genes <- gene.expression.sorted[gene.expression.sorted > 0]

exp_genes_sorted <- names(expressed.genes)

write.table(exp_genes_sorted, "NA_P09_NEURO_PROG_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


gene.expression.sorted2 <- as.data.frame(gene.expression.sorted)
write.table(gene.expression.sorted2, "NA_P09_NEURO_PROG_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





##-------------------------------------------------------
## SEPT 2023 | CTL
##-------------------------------------------------------
## FETCH THE LIST OF GENES EXPRESSED ACROSS iSPNs ONLY
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



load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

table(seuObjFilt$GenoAgeSample)
#     CTL_P09_NA13     CTL_P09_NA14     CTL_P09_NA29   P1CKO_P09_NA32 
#             8647            12579            12663            10972 
#   P1CKO_P09_NA33   P1CKO_P09_NA42 P1P2CKO_P09_NA15 P1P2CKO_P09_NA20 
#             7006            11012             9760             7846 
# P1P2CKO_P09_NA24   P2CKO_P09_NA16   P2CKO_P09_NA22   P2CKO_P09_NA43 
#            11238             4212             7690            18573

table(seuObjFilt$Genotype)
#     CTL   P1CKO P1P2CKO   P2CKO 
#   33889   28990   28844   30475

seuObjFiltCTL <- subset(seuObjFilt, subset = Genotype == c("CTL"))

table(seuObjFiltCTL$Genotype)
#   CTL 
# 33889 

table(seuObjFiltCTL$GenoAgeSample)
# CTL_P09_NA13 CTL_P09_NA14 CTL_P09_NA29 
#         8647        12579        12663

seuObjFiltCTL_prog <- subset(seuObjFiltCTL, ident = c("A_PROG"))
seuObjFiltCTL_neuroprog <- subset(seuObjFiltCTL, ident = c("B_NEURONAL_PROG"))
seuObjFiltCTL_dSPN <- subset(seuObjFiltCTL, ident = c("C_dSPN"))
seuObjFiltCTL_iSPN <- subset(seuObjFiltCTL, ident = c("D_iSPN"))
seuObjFiltCTL_eSPN <- subset(seuObjFiltCTL, ident = c("E_eSPN"))

gene.expression.prog <- apply(seuObjFiltCTL_prog@assays$RNA@data, 1, mean)
gene.expression.neuroprog <- apply(seuObjFiltCTL_neuroprog@assays$RNA@data, 1, mean)
gene.expression.dSPN <- apply(seuObjFiltCTL_dSPN@assays$RNA@data, 1, mean)
gene.expression.iSPN <- apply(seuObjFiltCTL_iSPN@assays$RNA@data, 1, mean)
gene.expression.eSPN <- apply(seuObjFiltCTL_eSPN@assays$RNA@data, 1, mean)

gene.expression.prog.sorted <- sort(gene.expression.prog, decreasing = TRUE)
gene.expression.neuroprog.sorted <- sort(gene.expression.neuroprog, decreasing = TRUE)
gene.expression.dSPN.sorted <- sort(gene.expression.dSPN, decreasing = TRUE)
gene.expression.iSPN.sorted <- sort(gene.expression.iSPN, decreasing = TRUE)
gene.expression.eSPN.sorted <- sort(gene.expression.eSPN, decreasing = TRUE)

expressed.genes.prog <- gene.expression.prog.sorted[gene.expression.prog.sorted > 0]
expressed.genes.neuroprog <- gene.expression.neuroprog.sorted[gene.expression.neuroprog.sorted > 0]
expressed.genes.dSPN <- gene.expression.dSPN.sorted[gene.expression.dSPN.sorted > 0]
expressed.genes.iSPN <- gene.expression.iSPN.sorted[gene.expression.iSPN.sorted > 0]
expressed.genes.eSPN <- gene.expression.eSPN.sorted[gene.expression.eSPN.sorted > 0]

exp_genes_sorted_prog <- names(expressed.genes.prog)
exp_genes_sorted_neuroprog <- names(expressed.genes.neuroprog)
exp_genes_sorted_dSPN <- names(expressed.genes.dSPN)
exp_genes_sorted_iSPN <- names(expressed.genes.iSPN)
exp_genes_sorted_eSPN <- names(expressed.genes.eSPN)

gene.expression.prog.sorted2 <- as.data.frame(gene.expression.prog.sorted)
gene.expression.neuroprog.sorted2 <- as.data.frame(gene.expression.neuroprog.sorted)
gene.expression.dSPN.sorted2 <- as.data.frame(gene.expression.dSPN.sorted)
gene.expression.iSPN.sorted2 <- as.data.frame(gene.expression.iSPN.sorted)
gene.expression.eSPN.sorted2 <- as.data.frame(gene.expression.eSPN.sorted)


write.table(exp_genes_sorted_prog, "NA_P09_CTL_PROG_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exp_genes_sorted_neuroprog, "NA_P09_CTL_NEURO_PROG_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exp_genes_sorted_dSPN, "NA_P09_CTL_dSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exp_genes_sorted_iSPN, "NA_P09_CTL_iSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exp_genes_sorted_eSPN, "NA_P09_CTL_eSPNs_EXPRESSED_GENES.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(gene.expression.prog.sorted2, "NA_P09_CTL_PROG_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene.expression.neuroprog.sorted2, "NA_P09_CTL_NEURO_PROG_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene.expression.dSPN.sorted2, "NA_P09_CTL_dSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene.expression.iSPN.sorted2, "NA_P09_CTL_iSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene.expression.eSPN.sorted2, "NA_P09_CTL_eSPNs_ALL_GENES.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")



#-------------------------------------------------------
## Combine into a single table
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


# prog <- read.table("NA_P09_CTL_PROG_ALL_GENES.txt", sep = "\t", header = FALSE)
# neuroprog <- read.table("NA_P09_CTL_NEURO_PROG_ALL_GENES.txt", sep = "\t", header = FALSE)
dSPN <- read.table("NA_P09_CTL_dSPNs_ALL_GENES.txt", sep = "\t", header = FALSE)
iSPN <- read.table("NA_P09_CTL_iSPNs_ALL_GENES.txt", sep = "\t", header = FALSE)
eSPN <- read.table("NA_P09_CTL_eSPNs_ALL_GENES.txt", sep = "\t", header = FALSE)

colnames(dSPN) <- c("Genes", "dSPN")
colnames(iSPN) <- c("Genes", "iSPN")
colnames(eSPN) <- c("Genes", "eSPN")

dSPN_sorted <- dSPN[order(dSPN$Genes),]
iSPN_sorted <- iSPN[order(iSPN$Genes),]
eSPN_sorted <- eSPN[order(eSPN$Genes),]

## combine individual tables into a giant data frame
dataCombined <- list("dSPN" = dSPN_sorted, "iSPN" = iSPN_sorted, "eSPN" = eSPN_sorted)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
write.table(combinedData, "SEPT2023_NA_P09_CTL_ALL_SPNs_ALL_GENES.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")






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


load("SEURAT_NA_P09_INTERGRATION_HARMONY_SPN_CLEANED.RData")
# seuObjFilt

table(seuObjFilt@active.ident)
#  D_iSPN          C_dSPN          E_eSPN          A_PROG B_NEURONAL_PROG 
#   43691           52124            8140            5181           13062

table(seuObjFilt$GenoAgeSample)
#     CTL_P09_NA13     CTL_P09_NA14     CTL_P09_NA29   P1CKO_P09_NA32 
#             8647            12579            12663            10972 
#   P1CKO_P09_NA33   P1CKO_P09_NA42 P1P2CKO_P09_NA15 P1P2CKO_P09_NA20 
#             7006            11012             9760             7846 
# P1P2CKO_P09_NA24   P2CKO_P09_NA16   P2CKO_P09_NA22   P2CKO_P09_NA43 
#            11238             4212             7690            18573

table(seuObjFilt$Genotype)
#     CTL   P1CKO P1P2CKO   P2CKO 
#   33889   28990   28844   30475

seuObjFiltCTL <- subset(seuObjFilt, subset = Genotype == c("CTL"))

table(seuObjFiltCTL$Genotype)
#   CTL 
# 33889 

table(seuObjFiltCTL$GenoAgeSample)
# CTL_P09_NA13 CTL_P09_NA14 CTL_P09_NA29 
#         8647        12579        12663

seuObjFiltCTL_prog <- subset(seuObjFiltCTL, ident = c("A_PROG"))
seuObjFiltCTL_neuroprog <- subset(seuObjFiltCTL, ident = c("B_NEURONAL_PROG"))
seuObjFiltCTL_dspn <- subset(seuObjFiltCTL, ident = c("C_dSPN"))
seuObjFiltCTL_ispn <- subset(seuObjFiltCTL, ident = c("D_iSPN"))
seuObjFiltCTL_espn <- subset(seuObjFiltCTL, ident = c("E_eSPN"))

gene.expression.prog <- as.data.frame(apply(seuObjFiltCTL_prog@assays$RNA@counts, 1, sum))
gene.expression.neuroprog <- as.data.frame(apply(seuObjFiltCTL_neuroprog@assays$RNA@counts, 1, sum))
gene.expression.dspn <- as.data.frame(apply(seuObjFiltCTL_dspn@assays$RNA@counts, 1, sum))
gene.expression.ispn <- as.data.frame(apply(seuObjFiltCTL_ispn@assays$RNA@counts, 1, sum))
gene.expression.espn <- as.data.frame(apply(seuObjFiltCTL_espn@assays$RNA@counts, 1, sum))

colnames(gene.expression.prog) <- c("Prog")
colnames(gene.expression.neuroprog) <- c("NeuroProg")
colnames(gene.expression.dspn) <- c("dSPN")
colnames(gene.expression.ispn) <- c("iSPN")
colnames(gene.expression.espn) <- c("eSPN")

gene.expression.prog$Genes <- row.names(gene.expression.prog)
gene.expression.neuroprog$Genes <- row.names(gene.expression.neuroprog)
gene.expression.dspn$Genes <- row.names(gene.expression.dspn)
gene.expression.ispn$Genes <- row.names(gene.expression.ispn)
gene.expression.espn$Genes <- row.names(gene.expression.espn)


gene.expression.prog.sorted <- gene.expression.prog[order(gene.expression.prog$Genes),]
gene.expression.neuroprog.sorted <- gene.expression.neuroprog[order(gene.expression.neuroprog$Genes),]
gene.expression.dspn.sorted <- gene.expression.dspn[order(gene.expression.dspn$Genes),]
gene.expression.ispn.sorted <- gene.expression.ispn[order(gene.expression.ispn$Genes),]
gene.expression.espn.sorted <- gene.expression.espn[order(gene.expression.espn$Genes),]


## combine individual tables into a giant data frame
dataCombined <- list("prog" = gene.expression.prog.sorted,
                     "neuroprog" = gene.expression.neuroprog.sorted,
                     "dspn" = gene.expression.dspn.sorted,
                     "ispn" = gene.expression.ispn.sorted,
                     "espn" = gene.expression.espn.sorted)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
write.table(combinedData, "SEPT2023_NA_P09_CTL_ALL_SPNs_ALL_GENES_COUNTS.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


#-------------------------------------------------------
# END
#-------------------------------------------------------



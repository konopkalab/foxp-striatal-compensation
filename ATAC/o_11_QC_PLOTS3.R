##---------------------------------------------------------------------------##
## Modules
# module purge && module load shared slurm python/3.7.x-anaconda
# module load hdf5_18/1.8.17
# module load gcc/8.3.0
# module load htslib
# module load gsl/2.4
# module load macs/2.1.2
# module load R/4.1.1-gccmkl

##---------------------------------------------------------------------------##
## R Libraries
rm(list = ls())
library(devtools)
library(dplyr)
library(scran)
library(Seurat)
library(Signac)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(chromVAR)
library(motifmatchr)
# library(ArchR)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(ggplot2)
library(scCustomize)
# addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

# # register(MulticoreParam(24, progressbar = TRUE))
# register(SerialParam())

##---------------------------------------------------------------------------##
## Load Data
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/I_LABELTRANSFER/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")
# naatac

DefaultAssay(naatac) <- "peaks"

metaatac <- as.data.frame(naatac@meta.data)
metaatac$Sample <- gsub("_CTL$|_P1CKO$|_P2CKO$|_P1P2CKO$", "", metaatac$Dataset)
metaatac$GenoSample <- paste(metaatac$Genotype, metaatac$Sample, sep = "_")

##---------------------------------------------------------------------------##
## Number of Nuclei per Genotype
nuclei.by.genosample <- as.data.frame(table(metaatac$predicted.id))
colnames(nuclei.by.genosample) <- c("CellType", "Nuclei")
row.names(nuclei.by.genosample) <- nuclei.by.genosample$CellType

#           CellType Nuclei
# 1           A_PROG   2723   gold2
# 2  B_NEURONAL_PROG   7694   gold4
# 3           C_dSPN  18519   blue2
# 4           D_iSPN  17308   red3
# 5           E_eSPN   2458   green3
# 6          F_ASTRO  11788   peru
# 7           G_OLIG   6230   palegreen
# 8          H_MICRO   1337   khaki3
# 9      I_ENDO_EPEN   1813   tomato2
# 10         J_INTER   3237   violet
# 11      K_CORTICAL   1928   steelblue
# 12         L_UNDET    900   grey70

## Barplot for Number of Nuclei
p1 <- ggplot(nuclei.by.genosample, aes(x = CellType, y = Nuclei, fill = CellType)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70")) +
      labs(x = "Genotypes", y = "# Nuclei", fill = "Genotype") +
      ggtitle("Number of Nuclei") + 
      NULL
ggsave(filename = "NA_ATAC_NUCLEI_PER_CELLTYPE_BARPLOT.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 300)


##---------------------------------------------------------------------------##
## Boxplot for Number of Peaks
p2 <- ggplot(metaatac, aes(x = predicted.id, y = nCount_peaks, color = predicted.id)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70")) +
      labs(x = "CellType", y = "# Peaks", fill = "CellType") +
      ggtitle("Number of Peaks") +
      NULL
ggsave(filename = "NA_ATAC_PEAKS_PER_CELLTYPE_BOXPLOT.pdf", plot = p2, width = 8, height = 4, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Boxplot for Number Reads in Peaks
p3 <- ggplot(metaatac, aes(x = predicted.id, y = peak_region_fragments, color = predicted.id)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70")) +
      labs(x = "CellType", y = "# Peaks", fill = "CellType") +
      ggtitle("Number of Reads in Peaks") +
      NULL
ggsave(filename = "NA_ATAC_READS_IN_PEAKS_PER_CELLTYPE_BOXPLOT.pdf", plot = p3, width = 8, height = 4, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Boxplot for Percent Reads in Peaks
p4 <- ggplot(metaatac, aes(x = predicted.id, y = pct_reads_in_peaks, color = predicted.id)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70")) +
      labs(x = "CellType", y = "# Peaks", fill = "CellType") +
      ggtitle("Number of Reads in Peaks") +
      ylim(0, 100) + 
      NULL
ggsave(filename = "NA_ATAC_PERCENT_READS_IN_PEAKS_PER_CELLTYPE_BOXPLOT.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 300)



# ##---------------------------------------------------------------------------##
# ## Scatterplot for Percent Reads in Peaks against Number of Reads
# p5 <- ggplot(metaatac, aes(x = passed_filters, y = pct_reads_in_peaks, color = Genotype)) +
#       geom_point(size = 0.2, alpha = 0.2) +
#       theme_classic() +
#       # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       scale_color_manual(values = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10")) +
#       labs(x = "# Total Reads", y = "% Reads in Peaks", fill = "Genotype") +
#       ylim(0, 100) + 
#       # ggtitle("Number of Reads in Peaks") +
#       NULL
# ggsave(filename = "NA_ATAC_PERCENT_READS_IN_PEAKS_VS_TOTAL_READS_SCATTERPLOT.pdf", plot = p5, width = 7, height = 6, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## UMAP colored by Cluster
pDim1 <- DimPlot_scCustom(naatac, group.by = "peaks_snn_res.1.2", pt.size = 2, raster = TRUE, raster.dpi = c(2048, 2048), label = TRUE, ggplot_default_colors = TRUE)
ggsave("NA_ATAC_UMAP_CLUSTERS_RES_1.2.pdf", plot = pDim1, width = 7, height = 6, units = "in", dpi = 300)

## UMAP colored by CellType
pDim2 <- DimPlot_scCustom(naatac, group.by = "predicted.id", pt.size = 2, raster = TRUE, raster.dpi = c(2048, 2048), label = TRUE, label.size = 2, repel = TRUE, colors_use = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70"))
ggsave("NA_ATAC_UMAP_CELLTYPE_wLABELS.pdf", plot = pDim2, width = 7.5, height = 6, units = "in", dpi = 300)

pDim2b <- DimPlot_scCustom(naatac, group.by = "predicted.id", pt.size = 2, raster = TRUE, raster.dpi = c(2048, 2048), label = FALSE, colors_use = c("gold2", "gold4", "blue2", "red3", "green3", "peru", "palegreen", "khaki3", "tomato2", "violet", "steelblue", "grey70"))
ggsave("NA_ATAC_UMAP_CELLTYPE_woLABELS.pdf", plot = pDim2b, width = 7.5, height = 6, units = "in", dpi = 300)

## UMAP colored by Genotype
pDim3 <- DimPlot_scCustom(naatac, group.by = "Genotype", pt.size = 2, raster = TRUE, raster.dpi = c(2048, 2048), label = FALSE, colors_use = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10"))
ggsave("NA_ATAC_UMAP_GENOTYPE_woLABELS.pdf", plot = pDim3, width = 6.5, height = 6, units = "in", dpi = 300)

## UMAP colored by GenoSample
pDim4 <- DimPlot_scCustom(naatac, group.by = "Dataset", pt.size = 2, raster = TRUE, raster.dpi = c(2048, 2048), label = FALSE, ggplot_default_colors = TRUE)
ggsave("NA_ATAC_UMAP_SAMPLE_woLABELS.pdf", plot = pDim4, width = 7, height = 6, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Number of Nuclei per CellType by GenoSample

celltype.by_genomsample <- as.data.frame.matrix(table(metaatac$predicted.id, metaatac$GenoSample))
write.table(celltype.by_genomsample, file = "NA_ATAC_CELLTYPE_BY_GENOSAMPLE.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

##---------------------------------------------------------------------------##
## END
##---------------------------------------------------------------------------##

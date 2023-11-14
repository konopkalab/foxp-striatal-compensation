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
nuclei.by.genosample <- as.data.frame(table(metaatac$Genotype))
colnames(nuclei.by.genosample) <- c("Geno", "Nuclei")
nuclei.by.genosample$Genotype <- gsub("^CTL$", "A_CTL", nuclei.by.genosample$Geno)
nuclei.by.genosample$Genotype <- gsub("^P1CKO$", "B_P1CKO", nuclei.by.genosample$Genotype)
nuclei.by.genosample$Genotype <- gsub("^P2CKO$", "C_P2CKO", nuclei.by.genosample$Genotype)
nuclei.by.genosample$Genotype <- gsub("^P1P2CKO$", "D_P1P2CKO", nuclei.by.genosample$Genotype)

row.names(nuclei.by.genosample) <- nuclei.by.genosample$Genotype

## Barplot for Number of Nuclei
p1 <- ggplot(nuclei.by.genosample, aes(x = Genotype, y = Nuclei, fill = Genotype)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10")) +
      labs(x = "Genotypes", y = "# Nuclei", fill = "Genotype") +
      ggtitle("Number of Nuclei") + 
      NULL
ggsave(filename = "NA_ATAC_NUCLEI_PER_GENOTYPE_BARPLOT.pdf", plot = p1, width = 4, height = 4, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Boxplot for Number of Peaks
p2 <- ggplot(metaatac, aes(x = Genotype, y = nCount_peaks, color = Genotype)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10")) +
      labs(x = "Genotype", y = "# Peaks", fill = "Genotype") +
      ggtitle("Number of Peaks") +
      NULL
ggsave(filename = "NA_ATAC_PEAKS_PER_GENOTYPE_BOXPLOT.pdf", plot = p2, width = 4, height = 4, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Boxplot for Number Reads in Peaks
p3 <- ggplot(metaatac, aes(x = Genotype, y = peak_region_fragments, color = Genotype)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10")) +
      labs(x = "Genotype", y = "# Peaks", fill = "Genotype") +
      ggtitle("Number of Reads in Peaks") +
      NULL
ggsave(filename = "NA_ATAC_READS_IN_PEAKS_PER_GENOTYPE_BOXPLOT.pdf", plot = p3, width = 4, height = 4, units = "in", dpi = 300)



##---------------------------------------------------------------------------##
## Boxplot for Percent Reads in Peaks
p4 <- ggplot(metaatac, aes(x = Genotype, y = pct_reads_in_peaks, color = Genotype)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10")) +
      labs(x = "Genotype", y = "# Peaks", fill = "Genotype") +
      ggtitle("Number of Reads in Peaks") +
      ylim(0, 100) + 
      NULL
ggsave(filename = "NA_ATAC_PERCENT_READS_IN_PEAKS_PER_GENOTYPE_BOXPLOT.pdf", plot = p4, width = 4, height = 4, units = "in", dpi = 300)



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



# ##---------------------------------------------------------------------------##
# ## UMAP colored by Cluster
# pDim1 <- DimPlot_scCustom(naatac, group.by = "peaks_snn_res.1.2", pt.size = 0.1, raster = TRUE, label = TRUE, ggplot_default_colors = TRUE)
# ggsave("NA_ATAC_UMAP_CLUSTERS_RES_1.2.pdf", plot = pDim1, width = 7, height = 6, units = "in", dpi = 300)

# ## UMAP colored by CellType
# pDim2 <- DimPlot_scCustom(naatac, group.by = "predicted.id", pt.size = 0.1, raster = TRUE, label = TRUE, label.size = 2, repel = TRUE, ggplot_default_colors = TRUE)
# ggsave("NA_ATAC_UMAP_CELLTYPE_wLABELS.pdf", plot = pDim2, width = 7.5, height = 6, units = "in", dpi = 300)

# pDim2b <- DimPlot_scCustom(naatac, group.by = "predicted.id", pt.size = 0.1, raster = TRUE, label = FALSE, ggplot_default_colors = TRUE)
# ggsave("NA_ATAC_UMAP_CELLTYPE_woLABELS.pdf", plot = pDim2b, width = 7.5, height = 6, units = "in", dpi = 300)

# ## UMAP colored by Genotype
# pDim3 <- DimPlot_scCustom(naatac, group.by = "Genotype", pt.size = 0.1, raster = TRUE, label = FALSE, colors_use = c("#191919", "#BB0BF8", "#5DCA9F", "#EEBC10"))
# ggsave("NA_ATAC_UMAP_GENOTYPE_woLABELS.pdf", plot = pDim3, width = 6.5, height = 6, units = "in", dpi = 300)

# ## UMAP colored by GenoSample
# pDim4 <- DimPlot_scCustom(naatac, group.by = "Dataset", pt.size = 0.1, raster = TRUE, label = FALSE, ggplot_default_colors = TRUE)
# ggsave("NA_ATAC_UMAP_SAMPLE_woLABELS.pdf", plot = pDim4, width = 7, height = 6, units = "in", dpi = 300)


##---------------------------------------------------------------------------##
## END
##---------------------------------------------------------------------------##

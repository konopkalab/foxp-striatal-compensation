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
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(harmony)
library(gridExtra)
library(rtracklayer)
library(scCustomize)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)



##-------------------------------------------------------
## Integrating with scRNA-seq data
## Load the pre-processed snRNA-seq data
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SEURAT_NA_P09_INTERGRATION_HARMONY_ANNOTATED.RData")
narna <- seuObjFilt
rm(seuObjFilt)

DefaultAssay(narna) <- "RNA"
narna <- NormalizeData(narna)



## Load snATAC-seq data
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/H_INTEGRATION/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED.RData")
naatac <- INTEGRATED
rm(INTEGRATED)

DefaultAssay(naatac) <- "RNA"


## Find Variable Features for snRNA-seq data
narna <- FindVariableFeatures(object = narna, nfeatures = 5000)

## Identify transfer anchors
transfer.anchors <- FindTransferAnchors(reference = narna, query = naatac, reduction = 'cca', dims = 1:50)

## Extract predicted labels
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = narna$CellType, weight.reduction = naatac[['integrated_lsi']], dims = 2:50)

## Add predicted labels to snATAC-seq object
naatac <- AddMetaData(object = naatac, metadata = predicted.labels)


## Save RData
save(naatac, transfer.anchors, predicted.labels, narna, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER.RData")


## DimPlots
plot1 <- DimPlot(narna, group.by = 'CellType', label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle('snRNA-seq')
plot2 <- DimPlot(naatac, group.by = 'predicted.id', label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle('snATAC-seq')
pInt <- grid.arrange(plot1, plot2, ncol = 2)
ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_UMAP.pdf", plot = pInt, width = 12, height = 6, units = "in", dpi = 300)


DefaultAssay(naatac) <- "RNA"
mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Tac1", "Drd2", "Penk", "Casz1", "Aqp4", "Cspg4", "Olig1", "Olig2", "Mag", "Cx3cr1", "Flt1")

plot3 <- FeaturePlot_scCustom(seurat_object = naatac, features = mygenes, pt.size = 0.1, raster = TRUE, min.cutoff = 0, max.cutoff = 2)
ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_GENES.pdf", plot = plot3, width = 16, height = 12, units = "in", dpi = 300)

pvln1 <- Stacked_VlnPlot(seurat_object = naatac, features = mygenes, group.by = "peaks_snn_res.0.8")
ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_GENES2.pdf", plot = pvln1, width = 16, height = 12, units = "in", dpi = 300)

pDim1 <- DimPlot(naatac, group.by = "peaks_snn_res.0.8", pt.size = 0.1, raster = TRUE, label = TRUE)
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_UMAP2.pdf", plot = pDim1, width = 7, height = 6, units = "in", dpi = 300)


## Additional resolutions for clustering
DefaultAssay(naatac) <- "peaks"
naatac <- FindClusters(naatac, verbose = FALSE, algorithm = 3, resolution = c(0.2, 0.4, 0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

pDim1 <- DimPlot(naatac, group.by = "peaks_snn_res.1.2", pt.size = 0.1, raster = TRUE, label = TRUE)
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_UMAP2.pdf", plot = pDim1, width = 7, height = 6, units = "in", dpi = 300)

## Save RData
save(naatac, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")

##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------

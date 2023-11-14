##-------------------------------------------------------
## SEURAT ANALYSIS | MODULES & LIBRARIES
##-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

##-------------------------------------------------------
## Load libraries
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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


## REFERENCE ANNOTATION
gtf <- rtracklayer::import("/work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_CELLRANGER_ATACSEQ/ForSignac/genes.gtf")
# fa <- Rsamtools::FaFile("/work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_CELLRANGER_ATACSEQ/ForSignac/genome.fa.gz")
mm10_annotation <- gtf[gtf$gene_type == "protein_coding"]
genome(mm10_annotation) <- "mm10"
seqlevelsStyle(mm10_annotation) <- "UCSC"
mm10_annotation$gene_biotype <- mm10_annotation$gene_type
mm10_annotation$tx_id <- mm10_annotation$transcript_id
mm10_annotation <- keepStandardChromosomes(mm10_annotation, pruning.mode = "coarse")



load("NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")
# naatac_dspn, daPeaks_ctl_p1cko, daPeaks_ctl_p2cko, daPeaks_ctl_p1p2cko


daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko
daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko
daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko



## dSPN | CTL vs. P1CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC > 0, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC < 0, ]
open_dspn_ctl <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]
closest_dspn_ctl <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl))
closest_dspn_p1cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1cko))
bg_ctl_p1cko <- unique(sort(c(unique(sort(closest_dspn_ctl$gene_name)), unique(sort(closest_dspn_p1cko$gene_name)))))
length(bg_ctl_p1cko)
# 1814

## dSPN | CTL vs. P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC > 0, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC < 0, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]
closest_dspn_ctl2 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl2))
closest_dspn_p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p2cko))
bg_ctl_p2cko <- unique(sort(c(unique(sort(closest_dspn_ctl2$gene_name)), unique(sort(closest_dspn_p2cko$gene_name)))))
length(bg_ctl_p2cko)
# 12652

## dSPN | CTL vs. P1P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC > 0, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC < 0, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]
closest_dspn_ctl3 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl3))
closest_dspn_p1p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1p2cko))
bg_ctl_p1p2cko <- unique(sort(c(unique(sort(closest_dspn_ctl3$gene_name)), unique(sort(closest_dspn_p1p2cko$gene_name)))))
length(bg_ctl_p1p2cko)
# 3045

length(unique(sort(c(bg_ctl_p1cko, bg_ctl_p2cko, bg_ctl_p1p2cko))))
# 12859












## LABEL TRANSFERRED snATAC-seq OBJECT
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/I_LABELTRANSFER/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")
# naatac

DefaultAssay(naatac) <- "peaks"

Annotation(naatac) <- mm10_annotation

table(naatac$Dataset)
#     NA104_CTL   NA105_P1CKO   NA106_P2CKO NA107_P1P2CKO     NA108_CTL 
#          2968          2656          2614           587         11060 
#   NA109_P1CKO   NA110_P2CKO NA111_P1P2CKO     NA112_CTL   NA113_P1CKO 
#         15557          4765          8311          8690          7095 
#   NA114_P2CKO NA115_P1P2CKO 
#          6286          5346

# naatac$Genotype <- gsub("NA104_CTL|NA108_CTL|NA112_CTL", "CTL", naatac$Dataset)
# naatac$Genotype <- gsub("NA105_P1CKO|NA109_P1CKO|NA113_P1CKO", "P1CKO", naatac$Genotype)
# naatac$Genotype <- gsub("NA106_P2CKO|NA110_P2CKO|NA114_P2CKO", "P2CKO", naatac$Genotype)
# naatac$Genotype <- gsub("NA107_P1P2CKO|NA111_P1P2CKO|NA115_P1P2CKO", "P1P2CKO", naatac$Genotype)

naatac$NewIdent <- paste(naatac$predicted.id, naatac$Genotype, sep = "_")

labeltransferbycluster <- as.data.frame.matrix(table(naatac$peaks_snn_res.1.2, naatac$predicted.id))
write.table(labeltransferbycluster, "NA_ATAC_LABELTRANSFER_BY_CLUSTER.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

pDim1 <- DimPlot(naatac, group.by = "peaks_snn_res.1.2", pt.size = 0.1, raster = TRUE, label = TRUE) + NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_RES_1.2.pdf", plot = pDim1, width = 6, height = 6, units = "in", dpi = 300)

pDim2 <- DimPlot(naatac, group.by = "Genotype", pt.size = 0.1, raster = TRUE, label = FALSE) #+ NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_GENOTYPE.pdf", plot = pDim2, width = 7, height = 6, units = "in", dpi = 300)

pDim2b <- DimPlot(naatac, group.by = "Genotype", split.by = "Genotype", pt.size = 0.1, raster = TRUE, label = FALSE) + NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_GENOTYPE_FACET.pdf", plot = pDim2b, width = 24, height = 6, units = "in", dpi = 300)

pDim3 <- DimPlot(naatac, group.by = "Dataset", pt.size = 0.1, raster = TRUE, label = FALSE) #+ NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_SAMPLE.pdf", plot = pDim3, width = 7.5, height = 6, units = "in", dpi = 300)

pDim3b <- DimPlot(naatac, group.by = "Dataset", split.by = "Dataset", pt.size = 0.1, raster = TRUE, label = FALSE, ncol = 4) + NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_SAMPLE_FACET.pdf", plot = pDim3b, width = 24, height = 18, units = "in", dpi = 300)

pDim4 <- DimPlot(naatac, group.by = "predicted.id", pt.size = 0.1, raster = TRUE, label = TRUE, label.size = 2) #+ NoLegend()
ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP_CELLTYPE.pdf", plot = pDim4, width = 7.5, height = 6, units = "in", dpi = 300)




DefaultAssay(naatac) <- "RNA"
mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Tac1", "Drd2", "Penk", "Casz1", "Aqp4", "Cspg4", "Olig1", "Olig2", "Mag", "Cx3cr1", "Flt1")
pvln1 <- Stacked_VlnPlot(seurat_object = naatac, features = mygenes, group.by = "peaks_snn_res.1.2")
ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_GENES_RES_1.2.pdf", plot = pvln1, width = 20, height = 12, units = "in", dpi = 300)

DefaultAssay(naatac) <- "peaks"

Idents(naatac) <- "peaks_snn_res.1.2"

table(naatac@active.ident)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7734 6522 3877 3633 3538 3471 3158 2901 2827 2710 2701 2672 2554 2172 2095 2090 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 1830 1805 1683 1458 1438 1424 1421 1307 1253 1211 1073  869  861  853  851  781 
#   32   33   34   35 
#  607  384  114   57

## subset for dSPNs
naatac_dspn <- subset(naatac, idents = c(1, 3, 6, 8))

table(naatac_dspn@active.ident)
#    1    3    6    8 
# 6522 3633 3158 2827

DefaultAssay(naatac_dspn) <- "peaks"
Idents(naatac_dspn) <- "Genotype"

## differential peaks 
daPeaks_ctl_p1cko <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P1CKO", min.pct = 0.05, logfc.threshold = 0.1, test.use = 'LR', latent.vars = 'peak_region_fragments')
daPeaks_ctl_p2cko <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P2CKO", min.pct = 0.05, logfc.threshold = 0.1, test.use = 'LR', latent.vars = 'peak_region_fragments')
daPeaks_ctl_p1p2cko <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P1P2CKO", min.pct = 0.05, logfc.threshold = 0.1, test.use = 'LR', latent.vars = 'peak_region_fragments')

save(naatac_dspn, daPeaks_ctl_p1cko, daPeaks_ctl_p2cko, daPeaks_ctl_p1p2cko, file = "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")


write.table(daPeaks_ctl_p1cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(daPeaks_ctl_p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(daPeaks_ctl_p1p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$p_val_adj <= 0.05,]
daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko[daPeaks_ctl_p2cko$p_val_adj <= 0.05,]
daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko[daPeaks_ctl_p1p2cko$p_val_adj <= 0.05,]


## dSPN | CTL vs. P1CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]
closest_dspn_ctl <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl))
closest_dspn_p1cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1cko))
write.table(closest_dspn_ctl, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Down_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(closest_dspn_p1cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Up_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


## dSPN | CTL vs. P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]
closest_dspn_ctl2 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl2))
closest_dspn_p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p2cko))
write.table(closest_dspn_ctl2, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Down_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(closest_dspn_p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Up_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")



## dSPN | CTL vs. P1P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]
closest_dspn_ctl3 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl3))
closest_dspn_p1p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1p2cko))
write.table(closest_dspn_ctl3, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Down_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(closest_dspn_p1p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Up_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")





# head(daPeaks_ctl_p1cko)

# pda1 <- Stacked_VlnPlot(seurat_object = naatac, features = rownames(daPeaks_ctl_p1cko)[1], pt.size = 0.1, idents = c("C_dSPN_CTL", "C_dSPN_P1CKO"))
# ggsave(filename = "DA_CTL_P1CKO_VIOLIN_PLOT_01.pdf", plot = pda1, width = 6, height = 4, units = "in", dpi = 300)

# open_dspn_ctl <- rownames(daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$avg_log2FC > 0.25, ])
# open_dspn_p1cko <- rownames(daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$avg_log2FC < -0.25, ])
# closest_dspn_ctl <- ClosestFeature(naatac, open_dspn_ctl)
# closest_dspn_p1cko <- ClosestFeature(naatac, open_dspn_p1cko)

## set plotting order
pcov1 <- CoveragePlot(object = naatac_dspn, region = "Foxp1", extend.upstream = 10000, extend.downstream = 10000, ncol = 1)
ggsave(filename = "COVERAGEPLOT_dSPN_Foxp1_2.pdf", plot = pcov1, width = 12, height = 10, units = "in", dpi = 300)

pcov2 <- CoveragePlot(object = naatac_dspn, region = "Drd", extend.upstream = 10000, extend.downstream = 10000, ncol = 1)
ggsave(filename = "COVERAGEPLOT_dSPN_Drd_2.pdf", plot = pcov2, width = 12, height = 10, units = "in", dpi = 300)

pcov3 <- CoveragePlot(object = naatac_dspn, region = "Foxp2", extend.upstream = 10000, extend.downstream = 10000, ncol = 1)
ggsave(filename = "COVERAGEPLOT_dSPN_Foxp2_2.pdf", plot = pcov3, width = 12, height = 10, units = "in", dpi = 300)


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------
load("NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")

daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$p_val_adj <= 0.05,]
daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko[daPeaks_ctl_p2cko$p_val_adj <= 0.05,]
daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko[daPeaks_ctl_p1p2cko$p_val_adj <= 0.05,]


## absolute(L2FC) >= 0.25

## dSPN | CTL vs. P1CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.25, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.25, ]
open_dspn_ctl <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]
# closest_dspn_ctl <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl))
# closest_dspn_p1cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1cko))
# write.table(closest_dspn_ctl, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Down_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p1cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Up_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


## dSPN | CTL vs. P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.25, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.25, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]
# closest_dspn_ctl2 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl2))
# closest_dspn_p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p2cko))
# write.table(closest_dspn_ctl2, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Down_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Up_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")



## dSPN | CTL vs. P1P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.25, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.25, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]
# closest_dspn_ctl3 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl3))
# closest_dspn_p1p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1p2cko))
# write.table(closest_dspn_ctl3, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Down_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p1p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Up_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")





## absolute(L2FC) >= 0.1375

## dSPN | CTL vs. P1CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]
# closest_dspn_ctl <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl))
# closest_dspn_p1cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1cko))
# write.table(closest_dspn_ctl, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Down_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p1cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1CKO_Closest_Up_in_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


## dSPN | CTL vs. P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]
# closest_dspn_ctl2 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl2))
# closest_dspn_p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p2cko))
# write.table(closest_dspn_ctl2, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Down_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P2CKO_Closest_Up_in_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")



## dSPN | CTL vs. P1P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]
# closest_dspn_ctl3 <- ClosestFeature(naatac_dspn, row.names(open_dspn_ctl3))
# closest_dspn_p1p2cko <- ClosestFeature(naatac_dspn, row.names(open_dspn_p1p2cko))
# write.table(closest_dspn_ctl3, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Down_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(closest_dspn_p1p2cko, "NA_ATAC_DifferentialPeaks_dSPN_CTL_vs_P1P2CKO_Closest_Up_in_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")








# ##-------------------------------------------------------
# ## Integrating with scRNA-seq data
# ## Load the pre-processed snRNA-seq data
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/20_INTEGRATION_wHarmony/SEURAT_NA_P09_INTERGRATION_HARMONY_ANNOTATED.RData")
# narna <- seuObjFilt
# rm(seuObjFilt)

# DefaultAssay(narna) <- "RNA"
# narna <- NormalizeData(narna)

# ## Load snATAC-seq data
# load("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY.RData")
# naatac <- HARMONISED
# rm(HARMONISED)

# DefaultAssay(naatac) <- "RNA"


# ## Find Variable Features for snRNA-seq data
# narna <- FindVariableFeatures(object = narna, nfeatures = 5000)

# ## Identify transfer anchors
# transfer.anchors <- FindTransferAnchors(reference = narna, query = naatac, reduction = 'cca', dims = 1:50)

# ## Extract predicted labels
# predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = narna$CellType, weight.reduction = naatac[['lsi']], dims = 2:50)

# ## Add predicted labels to snATAC-seq object
# naatac <- AddMetaData(object = naatac, metadata = predicted.labels)


# ## Save RData
# save(naatac, transfer.anchors, predicted.labels, narna, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER.RData")


# ## DimPlots
# plot1 <- DimPlot(narna, group.by = 'CellType', label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle('snRNA-seq')
# plot2 <- DimPlot(naatac, group.by = 'predicted.id', label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle('snATAC-seq')
# pInt <- grid.arrange(plot1, plot2, ncol = 2)
# ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_INT_UMAP2.pdf", plot = pInt, width = 12, height = 6, units = "in", dpi = 300)


# DefaultAssay(naatac) <- "RNA"
# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Tac1", "Drd2", "Penk", "Casz1", "Aqp4", "Cspg4", "Olig1", "Olig2", "Mag", "Cx3cr1", "Flt1")

# plot3 <- FeaturePlot_scCustom(seurat_object = naatac, features = mygenes, pt.size = 0.1, raster = TRUE, min.cutoff = 0, max.cutoff = 2)
# ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_GENES.pdf", plot = plot3, width = 16, height = 12, units = "in", dpi = 300)

# pvln1 <- Stacked_VlnPlot(seurat_object = naatac, features = mygenes, group.by = "peaks_snn_res.0.8")
# ggsave(filename = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_GENES2.pdf", plot = pvln1, width = 20, height = 12, units = "in", dpi = 300)

# pDim1 <- DimPlot(naatac, group.by = "peaks_snn_res.0.8", pt.size = 0.1, raster = TRUE, label = TRUE)
# ggsave("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_GA_HARMONY_LABTRANSFER_UMAP2.pdf", plot = pDim1, width = 7, height = 6, units = "in", dpi = 300)


# ##-------------------------------------------------------
# ## SEURAT ANALYSIS | END
# ##-------------------------------------------------------

## Modules
# module purge && module load shared slurm python/3.7.x-anaconda
# module load hdf5_18/1.8.17
# module load gcc/8.3.0
# module load htslib
# module load gsl/2.4
# module load macs/2.1.2
# module load R/4.1.1-gccmkl


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
# addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

# # register(MulticoreParam(24, progressbar = TRUE))
# register(SerialParam())

## FOX Motifs selected by NA
na_motifs <- read.table("NA_ALL_FOX_MOTIFS.txt", sep = "\t", header = FALSE)

## Tables for differential peaks
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/J_DIFF_PEAKS/BATCH_SEX/NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")
# daPeaks_ctl_p1cko
# daPeaks_ctl_p2cko
# daPeaks_ctl_p1p2cko
# naatac_dspn
rm(naatac_dspn)

## Differential peaks & filtering differential peaks
daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$p_val_adj <= 0.05,]
daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko[daPeaks_ctl_p2cko$p_val_adj <= 0.05,]
daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko[daPeaks_ctl_p1p2cko$p_val_adj <= 0.05,]

## dSPN | CTL vs. P1CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl1 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]

## dSPN | CTL vs. P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]

## dSPN | CTL vs. P1P2CKO
open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]


## Seurat/Signac object post motif analysis and tables for enriched motifs
# rm(list = ls())
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/K_MOTIF_ANALYSIS/BATCH_SEX/NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_2.RData")
# "enrichedMotifsCTLP1CKOup"  
# "enrichedMotifsCTLP1CKOdn"
# "enrichedMotifsCTLP2CKOup"  
# "enrichedMotifsCTLP2CKOdn"
# "enrichedMotifsCTLP1P2CKOup"
# "enrichedMotifsCTLP1P2CKOdn"
# "naatac_dspn"

## Fetch motif matrix
naatac_motif_matrix <- GetMotifData(naatac_dspn, assay = "peaks", slot = "data")
write.table(naatac_motif_matrix, "NA_ATAC_ALLPEAKS_DSPN_MOTIF_MATRIX.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")



naatac_mm_temp <- as.data.frame(naatac_motif_matrix[,colnames(naatac_motif_matrix) %in% na_motifs$V1])
naatac_mm_temp2 <- naatac_mm_temp[rowSums(naatac_mm_temp[]) > 0,]

naatac_mm_open_p1p2cko <- naatac_mm_temp2[row.names(naatac_mm_temp2) %in% row.names(open_dspn_p1p2cko),]
write.table(naatac_mm_open_p1p2cko, "DARs_FOR_FOX_MOTIFS_OPEN_P1P2CKO.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")









# df1[rowSums(df1[])>0,]




# # FOXP1 == MA0481.3
# # FOXP2 == MA0593.1
# # FOXP3 == MA0850.1

# naatac_mm_p1_temp <- as.data.frame(naatac_motif_matrix[,c("MA0481.3")])
# colnames(naatac_mm_p1_temp) <- "FOXP1"


# naatac_mm_p2_temp <- naatac_motif_matrix[,c("MA0593.1")]
# naatac_mm_p3_temp <- naatac_motif_matrix[,c("MA0850.1")]


# naatac_mm_p1 <- naatac_mm_p1_temp[naatac_mm_p1_temp == 1]
# naatac_mm_p2 <- naatac_mm_p2_temp[naatac_mm_p2_temp == 1]
# naatac_mm_p3 <- naatac_mm_p3_temp[naatac_mm_p3_temp == 1]



# ## Load NA motif lists
# namotifs <- read.table("chromVAR.txt", header = TRUE, sep = "\t")




# ## Load ChromVAR data
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/K_MOTIF_ANALYSIS/BATCH_SEX/NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_2.RData")
# # naatac_dspn, diffActCTLP1CKO, diffActCTLP2CKO, diffActCTLP1P2CKO

# motifConvertPlot <- function(motifids)
#     {
#     prefx <- gsub(".txt", "", motifids)
#     motifTable <- read.table(motifids, header = FALSE, sep = "\t")

#     pMotif <- MotifPlot(object = naatac_dspn, motifs = motifTable[,1], assay = 'peaks')
#     ggsave(filename = paste(prefx, "_Motifs.pdf", sep = ""), plot = pMotif, width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)

#     motifs_genes <- getMatrixByID(JASPAR2020, ID = motifTable[,1])
#     motifs_genes2 <- t(as.data.frame(lapply(motifs_genes, function(x) { x@name })))
#     write.table(motifs_genes2, paste(prefx, "_MotifGenes.txt", sep = ""), row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
#     }


# motifConvertPlot("NA_CHROMVAR_P1_UP.txt")
# motifConvertPlot("NA_CHROMVAR_P1_DN.txt")

# motifConvertPlot("NA_CHROMVAR_P2_UP.txt")
# motifConvertPlot("NA_CHROMVAR_P2_DN.txt")

# motifConvertPlot("NA_CHROMVAR_P1P2_UP.txt")
# motifConvertPlot("NA_CHROMVAR_P1P2_DN.txt")


# # ctlp1cko_motifs <- read.table("fromNA_motifs_P1CKO.txt", header = FALSE, sep = "\t")

# # pMotif_ctlp1cko <- MotifPlot(object = naatac_dspn, motifs = ctlp1cko_motifs[,1], assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_ChromVarMotifs_NA.pdf", plot = pMotif_ctlp1cko, width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)

# # ctlp1cko_motifs_genes <- getMatrixByID(JASPAR2020, ID = ctlp1cko_motifs[,1])
# # ctlp1cko_motifs_genes2 <- t(as.data.frame(lapply(ctlp1cko_motifs_genes, function(x) { x@name })))
# # write.table(ctlp1cko_motifs_genes2, "fromNA_motifs_P1CKO_MotifNames.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





# # ctlp2cko_motifs <- read.table("fromNA_motifs_P2CKO.txt", header = FALSE, sep = "\t")

# # pMotif_ctlp2cko <- MotifPlot(object = naatac_dspn, motifs = ctlp2cko_motifs[,1], assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_ChromVarMotifs_NA.pdf", plot = pMotif_ctlp2cko, width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)

# # ctlp2cko_motifs_genes <- getMatrixByID(JASPAR2020, ID = ctlp2cko_motifs[,1])
# # ctlp2cko_motifs_genes2 <- t(as.data.frame(lapply(ctlp2cko_motifs_genes, function(x) { x@name })))
# # write.table(ctlp2cko_motifs_genes2, "fromNA_motifs_P2CKO_MotifNames.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")






# # ctlp1p2cko_motifs <- read.table("fromNA_motifs_P1P2CKO.txt", header = FALSE, sep = "\t")

# # pMotif_ctlp1p2cko <- MotifPlot(object = naatac_dspn, motifs = ctlp1p2cko_motifs[,1], assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_ChromVarMotifs_NA.pdf", plot = pMotif_ctlp1p2cko, width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)

# # ctlp1p2cko_motifs_genes <- getMatrixByID(JASPAR2020, ID = ctlp1p2cko_motifs[,1])
# # ctlp1p2cko_motifs_genes2 <- t(as.data.frame(lapply(ctlp1p2cko_motifs_genes, function(x) { x@name })))
# # write.table(ctlp1p2cko_motifs_genes2, "fromNA_motifs_P1P2CKO_MotifNames.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")





# # ## Load snATAC-seq seurat/signac/macs2 data
# # # load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/97_SIGNAC_COMMON_PEAKS/NA_ATAC_COMMONPEAKS_SEURAT_COMBINED_NS_TSSE_FILT_SVD_GA_HARMONY_LABTRANSFER.RData")
# # # # naatac
# # # DefaultAssay(naatac) <- "ATAC"

# # ## Load snATAC-seq seurat/signac/macs2 and differential peaks data for dSPNs
# # load("NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS.RData")
# # # naatac_dspn, daPeaks_ctl_p1cko, daPeaks_ctl_p2cko, daPeaks_ctl_p1p2cko
# # DefaultAssay(naatac_dspn) <- "peaks"

# # ## Fetch matrix info from JASPAR
# # pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))


# # ## Add motif information
# # naatac_dspn <- AddMotifs(object = naatac_dspn, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)


# # ## Differential peaks & filtering differential peaks
# # daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$p_val_adj <= 0.05,]
# # daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko[daPeaks_ctl_p2cko$p_val_adj <= 0.05,]
# # daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko[daPeaks_ctl_p1p2cko$p_val_adj <= 0.05,]

# # ## dSPN | CTL vs. P1CKO
# # open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.1375, ]
# # open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.1375, ]
# # open_dspn_ctl1 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
# # open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]

# # ## dSPN | CTL vs. P2CKO
# # open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.1375, ]
# # open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.1375, ]
# # open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
# # open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]

# # ## dSPN | CTL vs. P1P2CKO
# # open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.1375, ]
# # open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.1375, ]
# # open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
# # open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]


# # ## Find peaks open in dSPNs (to be used as background peaks)
# # openPeaksCTLP1CKO <- AccessiblePeaks(naatac_dspn, idents = c("C_dSPN_CTL", "C_dSPN_P1CKO"))
# # openPeaksCTLP2CKO <- AccessiblePeaks(naatac_dspn, idents = c("C_dSPN_CTL", "C_dSPN_P2CKO"))
# # openPeaksCTLP1P2CKO <- AccessiblePeaks(naatac_dspn, idents = c("C_dSPN_CTL", "C_dSPN_P1P2CKO"))

# # ## Fetch meta features for all peaks
# # metaFeature <- GetAssayData(naatac_dspn, assay = "peaks", slot = "meta.features")

# # ## Match the overall GC content in the peak set
# # peaksMatchedCTLP1CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl1), ], n = 50000)
# # peaksMatchedCTLP1CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1CKO, ], query.feature = metaFeature[row.names(open_dspn_p1cko), ], n = 50000)

# # peaksMatchedCTLP2CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP2CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl2), ], n = 50000)
# # peaksMatchedCTLP2CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP2CKO, ], query.feature = metaFeature[row.names(open_dspn_p2cko), ], n = 50000)

# # peaksMatchedCTLP1P2CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1P2CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl3), ], n = 50000)
# # peaksMatchedCTLP1P2CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1P2CKO, ], query.feature = metaFeature[row.names(open_dspn_p1p2cko), ], n = 50000)


# # ## Finding overrepresented motifs 
# # ## use matched peaks in earlier step as background peaks
# # enrichedMotifsCTLP1CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl1), background = peaksMatchedCTLP1CKOdn)
# # enrichedMotifsCTLP1CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p1cko), background = peaksMatchedCTLP1CKOup)

# # enrichedMotifsCTLP2CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl2), background = peaksMatchedCTLP2CKOdn)
# # enrichedMotifsCTLP2CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p2cko), background = peaksMatchedCTLP2CKOup)

# # enrichedMotifsCTLP1P2CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl3), background = peaksMatchedCTLP1P2CKOdn)
# # enrichedMotifsCTLP1P2CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p1p2cko), background = peaksMatchedCTLP1P2CKOup)


# # ## Filter overrepresented motifs
# # ## use 'p.adjust' and 'fold.enrichment'
# # enrichedMotifsCTLP1CKOupFilt <- enrichedMotifsCTLP1CKOup[enrichedMotifsCTLP1CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP1CKOup$fold.enrichment >= 2,]
# # enrichedMotifsCTLP1CKOdnFilt <- enrichedMotifsCTLP1CKOdn[enrichedMotifsCTLP1CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP1CKOdn$fold.enrichment >= 2,]
# # enrichedMotifsCTLP1CKOupFiltSorted <- enrichedMotifsCTLP1CKOupFilt[order(enrichedMotifsCTLP1CKOupFilt$fold.enrichment, decreasing = TRUE),]
# # enrichedMotifsCTLP1CKOdnFiltSorted <- enrichedMotifsCTLP1CKOdnFilt[order(enrichedMotifsCTLP1CKOdnFilt$fold.enrichment, decreasing = TRUE),]

# # pMotif_ctlp1cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1CKOupFiltSorted)))
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_EnrichedMotifs_UP.pdf", plot = pMotif_ctlp1cko_up, width = 10, height = 6, units = "in", dpi = 300)

# # # pMotif_ctlp1cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1CKOdnFiltSorted)))
# # # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_EnrichedMotifs_DN.pdf", plot = pMotif_ctlp1cko_dn, width = 10, height = 6, units = "in", dpi = 300)



# # enrichedMotifsCTLP2CKOupFilt <- enrichedMotifsCTLP2CKOup[enrichedMotifsCTLP2CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP2CKOup$fold.enrichment >= 2,]
# # enrichedMotifsCTLP2CKOdnFilt <- enrichedMotifsCTLP2CKOdn[enrichedMotifsCTLP2CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP2CKOdn$fold.enrichment >= 2,]
# # enrichedMotifsCTLP2CKOupFiltSorted <- enrichedMotifsCTLP2CKOupFilt[order(enrichedMotifsCTLP2CKOupFilt$fold.enrichment, decreasing = TRUE),]
# # enrichedMotifsCTLP2CKOdnFiltSorted <- enrichedMotifsCTLP2CKOdnFilt[order(enrichedMotifsCTLP2CKOdnFilt$fold.enrichment, decreasing = TRUE),]

# # pMotif_ctlp2cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP2CKOupFiltSorted)))
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_EnrichedMotifs_UP.pdf", plot = pMotif_ctlp2cko_up, width = 10, height = 6, units = "in", dpi = 300)

# # # pMotif_ctlp2cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP2CKOdnFiltSorted)))
# # # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_EnrichedMotifs_DN.pdf", plot = pMotif_ctlp2cko_dn, width = 10, height = 6, units = "in", dpi = 300)




# # enrichedMotifsCTLP1P2CKOupFilt <- enrichedMotifsCTLP1P2CKOup[enrichedMotifsCTLP1P2CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP1P2CKOup$fold.enrichment >= 2,]
# # enrichedMotifsCTLP1P2CKOdnFilt <- enrichedMotifsCTLP1P2CKOdn[enrichedMotifsCTLP1P2CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP1P2CKOdn$fold.enrichment >= 2,]
# # enrichedMotifsCTLP1P2CKOupFiltSorted <- enrichedMotifsCTLP1P2CKOupFilt[order(enrichedMotifsCTLP1P2CKOupFilt$fold.enrichment, decreasing = TRUE),]
# # enrichedMotifsCTLP1P2CKOdnFiltSorted <- enrichedMotifsCTLP1P2CKOdnFilt[order(enrichedMotifsCTLP1P2CKOdnFilt$fold.enrichment, decreasing = TRUE),]

# # pMotif_ctlp1p2cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1P2CKOupFiltSorted)))
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_EnrichedMotifs_UP.pdf", plot = pMotif_ctlp1p2cko_up, width = 10, height = 6, units = "in", dpi = 300)

# # pMotif_ctlp1p2cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1P2CKOdnFiltSorted)))
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_EnrichedMotifs_DN.pdf", plot = pMotif_ctlp1p2cko_dn, width = 10, height = 6, units = "in", dpi = 300)




# # save(naatac_dspn, 
# #     enrichedMotifsCTLP1CKOup, enrichedMotifsCTLP1CKOdn, 
# #     enrichedMotifsCTLP2CKOup, enrichedMotifsCTLP2CKOdn, 
# #     enrichedMotifsCTLP1P2CKOup, enrichedMotifsCTLP1P2CKOdn, 
# #     file = "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS.RData")


# # write.table(enrichedMotifsCTLP1CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1CKOup.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# # write.table(enrichedMotifsCTLP1CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1CKOdn.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# # write.table(enrichedMotifsCTLP2CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P2CKOup.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# # write.table(enrichedMotifsCTLP2CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P2CKOdn.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# # write.table(enrichedMotifsCTLP1P2CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1P2CKOup.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# # write.table(enrichedMotifsCTLP1P2CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1P2CKOdn.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")






# # ## Computing motif activities
# # ## differentially-active motifs
# # ## ChromVAR identifies motifs associated with variability in chromatin accessibility between cells
# # naatac_dspn <- RunChromVAR(object = naatac_dspn, genome = BSgenome.Mmusculus.UCSC.mm10)

# # DefaultAssay(naatac_dspn) <- 'chromvar'

# # ## Differential activity scores
# # diffActCTLP1CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1CKO", only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# # diffActCTLP2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P2CKO", only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# # diffActCTLP1P2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1P2CKO", only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)


# # diffActCTLP1CKOFilt <- diffActCTLP1CKO[diffActCTLP1CKO$p_val_adj <= 0.05,]
# # diffActCTLP1CKOFiltSort <- diffActCTLP1CKOFilt[order(diffActCTLP1CKOFilt$avg_diff),]

# # diffActCTLP2CKOFilt <- diffActCTLP2CKO[diffActCTLP2CKO$p_val_adj <= 0.05,]
# # diffActCTLP2CKOFiltSort <- diffActCTLP2CKO[order(diffActCTLP2CKO$avg_diff),]

# # diffActCTLP1P2CKOFilt <- diffActCTLP1P2CKO[diffActCTLP1P2CKO$p_val_adj <= 0.05,]
# # diffActCTLP1P2CKOFiltSort <- diffActCTLP1P2CKO[order(diffActCTLP1P2CKO$avg_diff),]

# # pMotif_ctlp1cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP1CKOFiltSort)), assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_ChromVarMotifs.pdf", plot = pMotif_ctlp1cko, width = 10, height = 6, units = "in", dpi = 300)

# # pMotif_ctlp2cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP2CKOFiltSort)), assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_ChromVarMotifs.pdf", plot = pMotif_ctlp2cko, width = 10, height = 6, units = "in", dpi = 300)

# # pMotif_ctlp1p2cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP1P2CKOFiltSort)), assay = 'peaks')
# # ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_ChromVarMotifs.pdf", plot = pMotif_ctlp1p2cko, width = 10, height = 6, units = "in", dpi = 300)



# # save(naatac_dspn, 
# #      diffActCTLP1CKO,
# #      diffActCTLP2CKO,
# #      diffActCTLP1P2CKO,
# #      file = "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR.RData")



# # write.table(diffActCTLP1CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P1CKO.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# # write.table(diffActCTLP2CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P2CKO.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
# # write.table(diffActCTLP1P2CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P1P2CKO.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")







# # # Look at the activity of Foxp1 and Foxp2
# # p1 <- DimPlot(object = naatac_dspn, label = TRUE, pt.size = 0.1)
# # p2 <- FeaturePlot(object = naatac_dspn, features = "MA0481.3", min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)
# # p3 <- FeaturePlot(object = naatac_dspn, features = "MA0593.1", min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)
# # ggsave("NA_ATAC_ALLPEAKS_DSPN_ChromVAR_Foxp1_Foxp2.pdf", plot = grid.arrange(p1, p2, p3, ncol = 3), width = 18, height = 5, units = "in", dpi = 300)

# # ## END






# # DefaultAssay(naatac_dspn) <- 'chromvar'

# # ## Differential activity scores
# # diffActCTLP1CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# # diffActCTLP2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# # diffActCTLP1P2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)


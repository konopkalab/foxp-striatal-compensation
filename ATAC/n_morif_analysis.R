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
library(ArchR)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

# register(MulticoreParam(24, progressbar = TRUE))
register(SerialParam())

## Load snATAC-seq seurat/signac/macs2 data
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/I_LABELTRANSFER/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")
# # naatac
# DefaultAssay(naatac) <- "ATAC"

## Load snATAC-seq seurat/signac/macs2 and differential peaks data for dSPNs
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/J_DIFF_PEAKS/BATCH_SEX/NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")
# naatac_dspn, daPeaks_ctl_p1cko, daPeaks_ctl_p2cko, daPeaks_ctl_p1p2cko
DefaultAssay(naatac_dspn) <- "peaks"

## Fetch matrix info from JASPAR
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))


## Add motif information
naatac_dspn <- AddMotifs(object = naatac_dspn, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)


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


## Find peaks open in dSPNs (to be used as background peaks)
openPeaksCTLP1CKO <- AccessiblePeaks(naatac_dspn, idents = c("CTL", "P1CKO"))
openPeaksCTLP2CKO <- AccessiblePeaks(naatac_dspn, idents = c("CTL", "P2CKO"))
openPeaksCTLP1P2CKO <- AccessiblePeaks(naatac_dspn, idents = c("CTL", "P1P2CKO"))

## Fetch meta features for all peaks
metaFeature <- GetAssayData(naatac_dspn, assay = "peaks", slot = "meta.features")

## Match the overall GC content in the peak set
peaksMatchedCTLP1CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl1), ], n = 50000)
peaksMatchedCTLP1CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1CKO, ], query.feature = metaFeature[row.names(open_dspn_p1cko), ], n = 50000)

peaksMatchedCTLP2CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP2CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl2), ], n = 50000)
peaksMatchedCTLP2CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP2CKO, ], query.feature = metaFeature[row.names(open_dspn_p2cko), ], n = 50000)

peaksMatchedCTLP1P2CKOdn <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1P2CKO, ], query.feature = metaFeature[row.names(open_dspn_ctl3), ], n = 50000)
peaksMatchedCTLP1P2CKOup <- MatchRegionStats(meta.feature = metaFeature[openPeaksCTLP1P2CKO, ], query.feature = metaFeature[row.names(open_dspn_p1p2cko), ], n = 50000)


## Finding overrepresented motifs 
## use matched peaks in earlier step as background peaks
enrichedMotifsCTLP1CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl1), background = peaksMatchedCTLP1CKOdn)
enrichedMotifsCTLP1CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p1cko), background = peaksMatchedCTLP1CKOup)

enrichedMotifsCTLP2CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl2), background = peaksMatchedCTLP2CKOdn)
enrichedMotifsCTLP2CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p2cko), background = peaksMatchedCTLP2CKOup)

enrichedMotifsCTLP1P2CKOdn <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_ctl3), background = peaksMatchedCTLP1P2CKOdn)
enrichedMotifsCTLP1P2CKOup <- FindMotifs(object = naatac_dspn, features = row.names(open_dspn_p1p2cko), background = peaksMatchedCTLP1P2CKOup)


## Filter overrepresented motifs
## use 'p.adjust' and 'fold.enrichment'
enrichedMotifsCTLP1CKOupFilt <- enrichedMotifsCTLP1CKOup[enrichedMotifsCTLP1CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP1CKOup$fold.enrichment >= 2,]
enrichedMotifsCTLP1CKOdnFilt <- enrichedMotifsCTLP1CKOdn[enrichedMotifsCTLP1CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP1CKOdn$fold.enrichment >= 2,]
enrichedMotifsCTLP1CKOupFiltSorted <- enrichedMotifsCTLP1CKOupFilt[order(enrichedMotifsCTLP1CKOupFilt$fold.enrichment, decreasing = TRUE),]
enrichedMotifsCTLP1CKOdnFiltSorted <- enrichedMotifsCTLP1CKOdnFilt[order(enrichedMotifsCTLP1CKOdnFilt$fold.enrichment, decreasing = TRUE),]

pMotif_ctlp1cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1CKOupFiltSorted)))
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_EnrichedMotifs_UP_2.pdf", plot = pMotif_ctlp1cko_up, width = 10, height = 6, units = "in", dpi = 300)

# pMotif_ctlp1cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1CKOdnFiltSorted)))
# ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_EnrichedMotifs_DN_2.pdf", plot = pMotif_ctlp1cko_dn, width = 10, height = 6, units = "in", dpi = 300)



enrichedMotifsCTLP2CKOupFilt <- enrichedMotifsCTLP2CKOup[enrichedMotifsCTLP2CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP2CKOup$fold.enrichment >= 2,]
enrichedMotifsCTLP2CKOdnFilt <- enrichedMotifsCTLP2CKOdn[enrichedMotifsCTLP2CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP2CKOdn$fold.enrichment >= 2,]
enrichedMotifsCTLP2CKOupFiltSorted <- enrichedMotifsCTLP2CKOupFilt[order(enrichedMotifsCTLP2CKOupFilt$fold.enrichment, decreasing = TRUE),]
enrichedMotifsCTLP2CKOdnFiltSorted <- enrichedMotifsCTLP2CKOdnFilt[order(enrichedMotifsCTLP2CKOdnFilt$fold.enrichment, decreasing = TRUE),]

pMotif_ctlp2cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP2CKOupFiltSorted)))
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_EnrichedMotifs_UP_2.pdf", plot = pMotif_ctlp2cko_up, width = 10, height = 6, units = "in", dpi = 300)

# pMotif_ctlp2cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP2CKOdnFiltSorted)))
# ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_EnrichedMotifs_DN.pdf", plot = pMotif_ctlp2cko_dn, width = 10, height = 6, units = "in", dpi = 300)




enrichedMotifsCTLP1P2CKOupFilt <- enrichedMotifsCTLP1P2CKOup[enrichedMotifsCTLP1P2CKOup$p.adjust <= 0.05 & enrichedMotifsCTLP1P2CKOup$fold.enrichment >= 2,]
enrichedMotifsCTLP1P2CKOdnFilt <- enrichedMotifsCTLP1P2CKOdn[enrichedMotifsCTLP1P2CKOdn$p.adjust <= 0.05 & enrichedMotifsCTLP1P2CKOdn$fold.enrichment >= 2,]
enrichedMotifsCTLP1P2CKOupFiltSorted <- enrichedMotifsCTLP1P2CKOupFilt[order(enrichedMotifsCTLP1P2CKOupFilt$fold.enrichment, decreasing = TRUE),]
enrichedMotifsCTLP1P2CKOdnFiltSorted <- enrichedMotifsCTLP1P2CKOdnFilt[order(enrichedMotifsCTLP1P2CKOdnFilt$fold.enrichment, decreasing = TRUE),]

pMotif_ctlp1p2cko_up <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1P2CKOupFiltSorted)))
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_EnrichedMotifs_UP_2.pdf", plot = pMotif_ctlp1p2cko_up, width = 10, height = 6, units = "in", dpi = 300)

pMotif_ctlp1p2cko_dn <- MotifPlot(object = naatac_dspn, motifs = head(rownames(enrichedMotifsCTLP1P2CKOdnFiltSorted)))
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_EnrichedMotifs_DN_2.pdf", plot = pMotif_ctlp1p2cko_dn, width = 10, height = 6, units = "in", dpi = 300)




save(naatac_dspn, 
    enrichedMotifsCTLP1CKOup, enrichedMotifsCTLP1CKOdn, 
    enrichedMotifsCTLP2CKOup, enrichedMotifsCTLP2CKOdn, 
    enrichedMotifsCTLP1P2CKOup, enrichedMotifsCTLP1P2CKOdn, 
    file = "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_2.RData")


write.table(enrichedMotifsCTLP1CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1CKOup_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(enrichedMotifsCTLP1CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1CKOdn_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(enrichedMotifsCTLP2CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P2CKOup_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(enrichedMotifsCTLP2CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P2CKOdn_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(enrichedMotifsCTLP1P2CKOup, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1P2CKOup_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(enrichedMotifsCTLP1P2CKOdn, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CTL_P1P2CKOdn_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")






## Computing motif activities
## differentially-active motifs
## ChromVAR identifies motifs associated with variability in chromatin accessibility between cells
naatac_dspn <- RunChromVAR(object = naatac_dspn, genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(naatac_dspn) <- 'chromvar'

## Differential activity scores
diffActCTLP1CKO <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P1CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
diffActCTLP2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
diffActCTLP1P2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "CTL", ident.2 = "P1P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)


diffActCTLP1CKOFilt <- diffActCTLP1CKO[diffActCTLP1CKO$p_val_adj <= 0.05,]
diffActCTLP1CKOFiltSort <- diffActCTLP1CKOFilt[order(diffActCTLP1CKOFilt$avg_diff),]

diffActCTLP2CKOFilt <- diffActCTLP2CKO[diffActCTLP2CKO$p_val_adj <= 0.05,]
diffActCTLP2CKOFiltSort <- diffActCTLP2CKO[order(diffActCTLP2CKO$avg_diff),]

diffActCTLP1P2CKOFilt <- diffActCTLP1P2CKO[diffActCTLP1P2CKO$p_val_adj <= 0.05,]
diffActCTLP1P2CKOFiltSort <- diffActCTLP1P2CKO[order(diffActCTLP1P2CKO$avg_diff),]

pMotif_ctlp1cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP1CKOFiltSort)), assay = 'peaks')
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1CKO_ChromVarMotifs_2.pdf", plot = pMotif_ctlp1cko, width = 10, height = 6, units = "in", dpi = 300)

pMotif_ctlp2cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP2CKOFiltSort)), assay = 'peaks')
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P2CKO_ChromVarMotifs_2.pdf", plot = pMotif_ctlp2cko, width = 10, height = 6, units = "in", dpi = 300)

pMotif_ctlp1p2cko <- MotifPlot(object = naatac_dspn, motifs = head(rownames(diffActCTLP1P2CKOFiltSort)), assay = 'peaks')
ggsave(filename = "NA_ATAC_ALL_PEAKS_dSPN_CTL_P1P2CKO_ChromVarMotifs_2.pdf", plot = pMotif_ctlp1p2cko, width = 10, height = 6, units = "in", dpi = 300)



save(naatac_dspn, 
     diffActCTLP1CKO,
     diffActCTLP2CKO,
     diffActCTLP1P2CKO,
     file = "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_2.RData")



write.table(diffActCTLP1CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P1CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(diffActCTLP2CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(diffActCTLP1P2CKO, "NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_MOTIFS_CHROMVAR_CTL_P1P2CKO_2.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")







# Look at the activity of Foxp1 and Foxp2
p1 <- DimPlot(object = naatac_dspn, label = TRUE, pt.size = 0.1)
p2 <- FeaturePlot(object = naatac_dspn, features = "MA0481.3", min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)
p3 <- FeaturePlot(object = naatac_dspn, features = "MA0593.1", min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)
ggsave("NA_ATAC_ALLPEAKS_DSPN_ChromVAR_Foxp1_Foxp2_2.pdf", plot = grid.arrange(p1, p2, p3, ncol = 3), width = 18, height = 5, units = "in", dpi = 300)

## END






# DefaultAssay(naatac_dspn) <- 'chromvar'

# ## Differential activity scores
# diffActCTLP1CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# diffActCTLP2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)
# diffActCTLP1P2CKO <- FindMarkers(object = naatac_dspn, ident.1 = "C_dSPN_CTL", ident.2 = "C_dSPN_P1P2CKO", only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.1, min.pct = 0.05)


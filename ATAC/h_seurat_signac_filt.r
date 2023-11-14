##-------------------------------------------------------
## SEURAT ANALYSIS | MODULES & LIBRARIES
##-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

# NA104, NA108, NA112 --> CTL
# NA105b, NA109, NA113 --> P1CKO
# NA106, NA110, NA114 --> P1CKO
# NA107, NA111, NA115 --> P1P2CKO


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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## LOAD ARCHR DOUBLETS META
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/02_SEURAT_SIGNAC_DOUBLETS_ARCHR1/NA_ATAC_ARCHR_DOUBLETS_META.RData")
# naarchrmeta, nafiltarchrmeta
# Filtering 13433 cells from ArchRProject!
#         NA105_P1CKO : 159 of 3997 (4%)
#         NA106_P2CKO : 317 of 5636 (5.6%)
#         NA111_P1P2CKO : 2032 of 14255 (14.3%)
#         NA109_P1CKO : 2941 of 23685 (12.4%)
#         NA114_P2CKO : 638 of 7990 (8%)
#         NA104_CTL : 178 of 4222 (4.2%)
#         NA113_P1CKO : 847 of 9205 (9.2%)
#         NA115_P1P2CKO : 464 of 6819 (6.8%)
#         NA107_P1P2CKO : 10 of 1044 (1%)
#         NA108_CTL : 2955 of 23401 (12.6%)
#         NA110_P2CKO : 452 of 6727 (6.7%)
#         NA112_CTL : 2440 of 15623 (15.6%)

##-------------------------------------------------------
## Custom function to filter cells based on ArchR doublets
filterARCHR <- function(sam, sig)
    {
    # sam <- "NA104_CTL"
    # load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/A_CTL/NA104_CTL_MACS2_SIGNAC.RData")
    load(sig)
    Idents(MACS2_SEURAT) <- 'Dataset'

    metaArchRTemp <- nafiltarchrmeta[nafiltarchrmeta$Sample == sam,]
    temp1 <- as.data.frame(matrix(unlist(strsplit(row.names(metaArchRTemp), "#")), ncol = 2, byrow = TRUE))
    row.names(temp1) <- row.names(metaArchRTemp)
    metaArchR <- merge(metaArchRTemp, temp1, by = "row.names")
    row.names(metaArchR) <- metaArchR$V2
    seumeta <- as.data.frame(MACS2_SEURAT@meta.data)
    seumetaFilt <- merge(seumeta, metaArchR, by = "row.names")
    row.names(seumetaFilt) <- seumetaFilt$Row.names
    seumetaFilt$Row.names <- NULL
    cells2keep <- row.names(seumetaFilt)

    MACS2_SUB <- subset(MACS2_SEURAT, cells = cells2keep, slot = "counts")

    arTSSEnrichment <- seumetaFilt$TSSEnrichment
    names(arTSSEnrichment) <- row.names(seumetaFilt)
    arReadsInTSS <- seumetaFilt$ReadsInTSS
    names(arReadsInTSS) <- row.names(seumetaFilt)
    arReadsInPromoter <- seumetaFilt$ReadsInPromoter
    names(arReadsInPromoter) <- row.names(seumetaFilt)
    arReadsInBlacklist <- seumetaFilt$ReadsInBlacklist
    names(arReadsInBlacklist) <- row.names(seumetaFilt)
    arPromoterRatio <- seumetaFilt$PromoterRatio
    names(arPromoterRatio) <- row.names(seumetaFilt)
    arPassQC <- seumetaFilt$PassQC
    names(arPassQC) <- row.names(seumetaFilt)
    arNucleosomeRatio <- seumetaFilt$NucleosomeRatio
    names(arNucleosomeRatio) <- row.names(seumetaFilt)
    arnMultiFrags <- seumetaFilt$nMultiFrags
    names(arnMultiFrags) <- row.names(seumetaFilt)
    arnMonoFrags <- seumetaFilt$nMonoFrags
    names(arnMonoFrags) <- row.names(seumetaFilt)
    arnFrags <- seumetaFilt$nFrags
    names(arnFrags) <- row.names(seumetaFilt)
    arnDiFrags <- seumetaFilt$nDiFrags
    names(arnDiFrags) <- row.names(seumetaFilt)
    arDoubletScore <- seumetaFilt$DoubletScore
    names(arDoubletScore) <- row.names(seumetaFilt)
    arDoubletEnrichment <- seumetaFilt$DoubletEnrichment
    names(arDoubletEnrichment) <- row.names(seumetaFilt)
    arBlacklistRatio <- seumetaFilt$BlacklistRatio
    names(arBlacklistRatio) <- row.names(seumetaFilt)

    MACS2_SUB$arTSSEnrichment <- arTSSEnrichment
    MACS2_SUB$arReadsInTSS <- arReadsInTSS
    MACS2_SUB$arReadsInPromoter <- arReadsInPromoter
    MACS2_SUB$arReadsInBlacklist <- arReadsInBlacklist
    MACS2_SUB$arPromoterRatio <- arPromoterRatio
    MACS2_SUB$arPassQC <- arPassQC
    MACS2_SUB$arNucleosomeRatio <- arNucleosomeRatio
    MACS2_SUB$arnMultiFrags <- arnMultiFrags
    MACS2_SUB$arnMonoFrags <- arnMonoFrags
    MACS2_SUB$arnFrags <- arnFrags
    MACS2_SUB$arnDiFrags <- arnDiFrags
    MACS2_SUB$arDoubletScore <- arDoubletScore
    MACS2_SUB$arDoubletEnrichment <- arDoubletEnrichment
    MACS2_SUB$arBlacklistRatio <- arBlacklistRatio

    MACS2_FILT <- subset(x = MACS2_SUB, subset = peak_region_fragments > 1500 &
                                                 peak_region_fragments < 100000 &
                                                 pct_reads_in_peaks > 10 &
                                                 arBlacklistRatio < 0.03 &
                                                 nucleosome_signal < 2 &
                                                 TSS.enrichment > 2
                        )

    print(dim(MACS2_FILT))

    pqc1 <- VlnPlot(object = MACS2_SUB, features = c('pct_reads_in_peaks', 'peak_region_fragments', 'arReadsInPromoter', 'arReadsInTSS', 'TSS.enrichment', 'arTSSEnrichment', 'arBlacklistRatio', 'nucleosome_signal', 'arNucleosomeRatio', 'arPromoterRatio'), pt.size = 0.1, ncol = 5)
    ggsave(filename = paste(sam, "_FILT_QC_PLOT_2.pdf", sep = ""), plot = pqc1, width = 15, height = 10, units = "in", dpi = 300)

    save(MACS2_FILT, file = paste(sam, "_MACS2_SIGNAC_FILT.RData", sep = ""))
    return(MACS2_FILT)
    }


NA104_CTL_FILT <- filterARCHR("NA104_CTL", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/A_CTL/NA104_CTL_MACS2_SIGNAC.RData")
NA108_CTL_FILT <- filterARCHR("NA108_CTL", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/A_CTL/NA108_CTL_MACS2_SIGNAC.RData")
NA112_CTL_FILT <- filterARCHR("NA112_CTL", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/A_CTL/NA112_CTL_MACS2_SIGNAC.RData")

NA105_P1CKO_FILT <- filterARCHR("NA105_P1CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/B_P1CKO/NA105_P1CKO_MACS2_SIGNAC.RData")
NA109_P1CKO_FILT <- filterARCHR("NA109_P1CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/B_P1CKO/NA109_P1CKO_MACS2_SIGNAC.RData")
NA113_P1CKO_FILT <- filterARCHR("NA113_P1CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/B_P1CKO/NA113_P1CKO_MACS2_SIGNAC.RData")

NA106_P2CKO_FILT <- filterARCHR("NA106_P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/C_P2CKO/NA106_P2CKO_MACS2_SIGNAC.RData")
NA110_P2CKO_FILT <- filterARCHR("NA110_P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/C_P2CKO/NA110_P2CKO_MACS2_SIGNAC.RData")
NA114_P2CKO_FILT <- filterARCHR("NA114_P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/C_P2CKO/NA114_P2CKO_MACS2_SIGNAC.RData")

NA107_P1P2CKO_FILT <- filterARCHR("NA107_P1P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/D_P1P2CKO/NA107_P1P2CKO_MACS2_SIGNAC.RData")
NA111_P1P2CKO_FILT <- filterARCHR("NA111_P1P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/D_P1P2CKO/NA111_P1P2CKO_MACS2_SIGNAC.RData")
NA115_P1P2CKO_FILT <- filterARCHR("NA115_P1P2CKO", "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/D_P1P2CKO/NA115_P1P2CKO_MACS2_SIGNAC.RData")

save(NA104_CTL_FILT, NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT,
     NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT,
     NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT,
     file = "NA_ATAC_SIGNAC_MACS2_FILT.RData"
    )


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------

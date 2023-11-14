##--------------------------------------------
## Modules
# module purge && module load shared slurm python/3.7.x-anaconda
# module load hdf5_18/1.8.17
# module load gcc/8.3.0
# module load htslib
# module load gsl/2.4
# module load macs/2.1.2
# module load R/4.1.1-gccmkl

##--------------------------------------------
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
library(ggbio)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##--------------------------------------------
## snRNA-seq bigwig
bw_p09_ctl_for <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_CTLTagDir/P09_CTL.for.bw"
bw_p09_ctl_rev <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_CTLTagDir/P09_CTL.rev.bw"

bw_p09_p1cko_for <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P1CKOTagDir/P09_P1CKO.for.bw"
bw_p09_p1cko_rev <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P1CKOTagDir/P09_P1CKO.rev.bw"

bw_p09_p2cko_for <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P2CKOTagDir/P09_P2CKO.for.bw"
bw_p09_p2cko_rev <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P2CKOTagDir/P09_P2CKO.rev.bw"

bw_p09_p1p2cko_for <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P1P2CKOTagDir/P09_P1P2CKO.for.bw"
bw_p09_p1p2cko_rev <- "/project/Neuroinformatics_Core/Konopka_lab/akulk1/NA_RNA_BAM/P09_P1P2CKOTagDir/P09_P1P2CKO.rev.bw"


##--------------------------------------------
## snATAC-seq signac/seurat data
# ## Load snATAC-seq seurat/signac/macs2 data
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/I_LABELTRANSFER/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")
# # naatac

## Load snATAC-seq seurat/signac/macs2 and differential peaks data for dSPNs
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/J_DIFF_PEAKS/BATCH_SEX/NA_ATAC_ALLPEAKS_DSPN_DIFFERENTIAL_PEAKS_2.RData")
# naatac_dspn, daPeaks_ctl_p1cko, daPeaks_ctl_p2cko, daPeaks_ctl_p1p2cko

DefaultAssay(naatac_dspn) <- "peaks"

##--------------------------------------------
## snATAC-seq differential peaks
daPeaks_ctl_p1cko_sig <- daPeaks_ctl_p1cko[daPeaks_ctl_p1cko$p_val_adj <= 0.05,]
daPeaks_ctl_p2cko_sig <- daPeaks_ctl_p2cko[daPeaks_ctl_p2cko$p_val_adj <= 0.05,]
daPeaks_ctl_p1p2cko_sig <- daPeaks_ctl_p1p2cko[daPeaks_ctl_p1p2cko$p_val_adj <= 0.05,]

open_dspn_ctlTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1ckoTemp <- daPeaks_ctl_p1cko_sig[daPeaks_ctl_p1cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1cko <- open_dspn_p1ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1ckoTemp)),]

open_dspn_ctlTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p2ckoTemp <- daPeaks_ctl_p2cko_sig[daPeaks_ctl_p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl2 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p2cko <- open_dspn_p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p2ckoTemp)),]

open_dspn_ctlTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC >= 0.1375, ]
open_dspn_p1p2ckoTemp <- daPeaks_ctl_p1p2cko_sig[daPeaks_ctl_p1p2cko_sig$avg_log2FC <= -0.1375, ]
open_dspn_ctl3 <- open_dspn_ctlTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_ctlTemp)),]
open_dspn_p1p2cko <- open_dspn_p1p2ckoTemp[!grepl("^chrM|^chrX|^chrY", row.names(open_dspn_p1p2ckoTemp)),]

na_open_dko <- read.table("dKO_UP.bed", header = TRUE, sep = "\t")
row.names(na_open_dko) <- paste(na_open_dko$chr, na_open_dko$start, na_open_dko$end, sep = "-")


##--------------------------------------------
## regions or genes of interst
# regPde1c <- "chr6-56039160-56669191"
# regPde1c <- "chr6-56039160-56400028"
regPde1c <- "chr6-56577503-56578000"

##--------------------------------------------
## Pde1c (regPde1c)
ranges_open_ctl1 <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_ctl)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_p1cko <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_p1cko)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_ctl2 <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_ctl2)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_p2cko <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_p2cko)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_ctl3 <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_ctl3)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_p1p2cko <- GenomicRanges::intersect(StringToGRanges(row.names(open_dspn_p1p2cko)), StringToGRanges(regPde1c), ignore.strand = TRUE)
ranges_open_dko <- GenomicRanges::intersect(StringToGRanges(row.names(na_open_dko)), StringToGRanges(regPde1c), ignore.strand = TRUE)


##--------------------------------------------
## bigwig list
bwlist <- list("RNA_CTL_FOR" = bw_p09_ctl_for, "RNA_CTL_REV" = bw_p09_ctl_rev, 
               "RNA_P1CKO_FOR" = bw_p09_p1cko_for, "RNA_P1CKO_REV" = bw_p09_p1cko_rev, 
               "RNA_P2CKO_FOR" = bw_p09_p2cko_for, "RNA_P2CKO_REV" = bw_p09_p2cko_rev, 
               "RNA_P1P2CKO_FOR" = bw_p09_p1p2cko_for, "RNA_P1P2CKO_REV" = bw_p09_p1p2cko_rev)

##--------------------------------------------
## range list
rangeslist <- GRangesList("Open_CTL1" = ranges_open_ctl1,
                        "Open_P1cKO" = ranges_open_p1cko,
                        "Open_CTL2" = ranges_open_ctl2,
                        "Open_P2cKO" = ranges_open_p2cko,
                        "Open_CTL3" = ranges_open_ctl3,
                        "Open_P1P2cKO" = ranges_open_p1p2cko,
                        "Open_dKO" = ranges_open_dko)

print(unlist(rangeslist))
# GRanges object with 2 ranges and 0 metadata columns:
#                seqnames            ranges strand
#                   <Rle>         <IRanges>  <Rle>
#   Open_P1P2cKO     chr6 56577503-56577991      *
#       Open_dKO     chr6 56577503-56577991      *


##--------------------------------------------
## Coverage Plot for selected region of a Gene
cov_plot <- CoveragePlot(object = naatac_dspn, 
                region = regPde1c, 
                annotation = TRUE, 
                peaks = TRUE, 
                bigwig = bwlist, 
                bigwig.type = "coverage", 
                ranges = unlist(rangeslist),
                extend.upstream = 10000,
                extend.downstream = 2500
                )
ggsave(filename = "COVERAGEPLOT_dSPN_Pde1c_1.pdf", plot = cov_plot, width = 12, height = 12, units = "in", dpi = 300)


cov_plot <- CoveragePlot(object = naatac_dspn, 
                region = regPde1c, 
                annotation = TRUE, 
                peaks = TRUE, 
                # bigwig = bwlist, 
                bigwig.type = "coverage", 
                ranges = unlist(rangeslist),
                extend.upstream = 10000,
                extend.downstream = 2500
                )
ggsave(filename = "COVERAGEPLOT_dSPN_Pde1c_1_woRNA.pdf", plot = cov_plot, width = 12, height = 12, units = "in", dpi = 300)




##--------------------------------------------
## Coverage Plot for Complete Gene
cov_plot <- CoveragePlot(object = naatac_dspn, 
                region = "Pde1c", 
                annotation = TRUE, 
                peaks = TRUE, 
                bigwig = bwlist, 
                bigwig.type = "coverage", 
                ranges = unlist(rangeslist),
                extend.upstream = 1000,
                extend.downstream = 10000
                )
ggsave(filename = "COVERAGEPLOT_dSPN_Pde1c_2.pdf", plot = cov_plot, width = 16, height = 12, units = "in", dpi = 300)








## Coverage Plot for Individual Peaks/Ranges
myranges <- as.data.frame(ranges(unlist(rangeslist)))
myranges$new <- paste("chr6", myranges$start, myranges$end, sep = "-")
mypeaks <- unique(myranges$new)

print(mypeaks)

for(i in 1:length(mypeaks))
    {
    print(mypeaks[[i]])
    cov_plot2 <- CoveragePlot(
                object = naatac_dspn,
                region = StringToGRanges(mypeaks[[i]]),
                annotation = TRUE,
                peaks = TRUE,
                tile = FALSE,
                bigwig = bwlist, 
                bigwig.type = "coverage", 
                extend.upstream = 2500, 
                extend.downstream = 2500,
                ranges = StringToGRanges(mypeaks[[i]])
                )
    ggsave(filename = paste("COVERAGEPLOT_dSPN_Foxp1_", mypeaks[[i]], ".pdf", sep = ""), plot = cov_plot2, width = 8, height = 12, units = "in", dpi = 300)
    }


##--------------------------------------------
## END
##--------------------------------------------







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
library(SeuratDisk)
set.seed(1234)


# # register(MulticoreParam(24, progressbar = TRUE))
# register(SerialParam())

##---------------------------------------------------------------------------##
## Load Data
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/I_LABELTRANSFER/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED_LABTRANSFER_2.RData")
# naatac

DefaultAssay(naatac) <- "peaks"

## Extract meta table
metaatac <- as.data.frame(naatac@meta.data)
metaatac$Sample <- gsub("_CTL$|_P1CKO$|_P2CKO$|_P1P2CKO$", "", metaatac$Dataset)
metaatac$GenoSample <- paste(metaatac$Genotype, metaatac$Sample, sep = "_")

## Extract counts table
countsatac <- GetAssayData(naatac, slot = "counts")
mygenes <- as.data.frame(row.names(countsatac))
colnames(mygenes) <- "GeneID"
mygenes$GeneName <- mygenes$GeneID
mygenes$Category <- rep("Gene Expression", nrow(mygenes))
mycellbarcodes <- colnames(countsatac)

write.table(metaatac, "meta.tsv", col.names = T, sep = '\t', quote = F, row.names = T)
write.table(mygenes, "features.tsv", col.names = F, sep = '\t', quote = F, row.names = F)
write.table(mycellbarcodes, "barcodes.tsv", col.names = F, sep = '\t', quote = F, row.names = F)
writeMM(countsatac, "matrix.mtx")


##---------------------------------------------------------------------------##
## END
##---------------------------------------------------------------------------##

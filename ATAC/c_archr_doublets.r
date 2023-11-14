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
addArchRGenome("mm10")
set.seed(1234)



## Input files
na104_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA104_CTL/fragments.tsv.gz'
na105_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA105_P1CKO/fragments.tsv.gz'
na106_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA106_P2CKO/fragments.tsv.gz'
na107_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA107_P1P2CKO/fragments.tsv.gz'

na108_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA108_CTL/fragments.tsv.gz'
na109_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA109_P1CKO/fragments.tsv.gz'
na110_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA110_P2CKO/fragments.tsv.gz'
na111_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA111_P1P2CKO/fragments.tsv.gz'

na112_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA112_CTL_F/fragments.tsv.gz'
na113_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA113_P1CKO_M/fragments.tsv.gz'
na114_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA114_P2CKO_M/fragments.tsv.gz'
na115_frag <- '/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA115_P1P2CKO_M/fragments.tsv.gz'


inputFiles <- c(na104_frag, na105_frag, na106_frag, na107_frag,
                na108_frag, na109_frag, na110_frag, na111_frag,
                na112_frag, na113_frag, na114_frag, na115_frag)
names(inputFiles) <- c("NA104_CTL", "NA105_P1CKO", "NA106_P2CKO", "NA107_P1P2CKO", "NA108_CTL", "NA109_P1CKO", "NA110_P2CKO", "NA111_P1P2CKO", "NA112_CTL", "NA113_P1CKO", "NA114_P2CKO", "NA115_P1P2CKO")


## Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 0, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading=FALSE,
  force=TRUE
  )


## Inferring Doublets
### Change n_threads=1 to oversee parallelization issue on ArchR for addDoubletScores() function
addArchRThreads(threads = 1)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE
  )

## Creating an ArchRProject
projNA <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "NA_ARCHR",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
  )


## Filtering Doublets
projNAfilt <- filterDoublets(ArchRProj = projNA)
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

save(inputFiles, ArrowFiles, projNA, projNAfilt, file = "NA_ATAC_ARCHR_DOUBLETS.RData")

naarchrmeta <- as.data.frame(projNA@cellColData)  
nafiltarchrmeta <- as.data.frame(projNAfilt@cellColData)  

save(naarchrmeta, nafiltarchrmeta, file = "NA_ATAC_ARCHR_DOUBLETS_META.RData")


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
## Load filtered seurat/signac objects
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/E_FILT/NA_ATAC_SIGNAC_MACS2_FILT.RData")
# NA104_CTL_FILT, NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT,
# NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT,
# NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT,


##-------------------------------------------------------
## Default assay to peaks
DefaultAssay(NA104_CTL_FILT) <- 'peaks'
DefaultAssay(NA108_CTL_FILT) <- 'peaks'
DefaultAssay(NA112_CTL_FILT) <- 'peaks'

# DefaultAssay(NA105_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA109_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA113_P1CKO_FILT) <- 'peaks'

# DefaultAssay(NA106_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA110_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA114_P2CKO_FILT) <- 'peaks'

# DefaultAssay(NA107_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA111_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA115_P1P2CKO_FILT) <- 'peaks'


##-------------------------------------------------------
## Merge all 12 samples from four genotypes
# MERGED <- merge(x = NA104_CTL_FILT, y = list(NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT, NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT, NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT), add.cell.ids = c("NA104_CTL", "NA105_P1CKO", "NA106_P2CKO", "NA107_P1P2CKO", "NA108_CTL", "NA109_P1CKO", "NA110_P2CKO", "NA111_P1P2CKO", "NA112_CTL", "NA113_P1CKO", "NA114_P2CKO", "NA115_P1P2CKO"))
MERGED.CTL <- merge(x = NA104_CTL_FILT, y = list(NA108_CTL_FILT, NA112_CTL_FILT), add.cell.ids = c("NA104_CTL", "NA108_CTL", "NA112_CTL"))


print(table(MERGED.CTL$Dataset))
# NA104_CTL NA108_CTL NA112_CTL 
#      2968     11060      8690

save(MERGED.CTL, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_CTL_TEMP.RData")

##-------------------------------------------------------
## Process merged dataset
MERGED.CTL <- RunTFIDF(MERGED.CTL)
MERGED.CTL <- FindTopFeatures(MERGED.CTL, min.cutoff = 20)
MERGED.CTL <- RunSVD(MERGED.CTL)
MERGED.CTL <- RunUMAP(MERGED.CTL, dims = 2:50, reduction = 'lsi')
MERGED.CTL <- FindNeighbors(MERGED.CTL, dims = 2:50, reduction = 'lsi')
MERGED.CTL <- FindClusters(MERGED.CTL, verbose = FALSE, algorithm = 3)

MERGED_UMAP <- DimPlot(MERGED.CTL, group.by = c("Dataset", "seurat_clusters"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_CTL.pdf", plot = MERGED_UMAP, width = 15, height = 6, units = "in", dpi = 300)

MERGED_UMAP_FACET <- DimPlot(MERGED.CTL, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_CTL_FACET.pdf", plot = MERGED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)

save(MERGED.CTL, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_CTL.RData")


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------




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
## Load filtered seurat/signac objects
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/E_FILT/NA_ATAC_SIGNAC_MACS2_FILT.RData")
# NA104_CTL_FILT, NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT,
# NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT,
# NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT,


##-------------------------------------------------------
## Default assay to peaks
# DefaultAssay(NA104_CTL_FILT) <- 'peaks'
# DefaultAssay(NA108_CTL_FILT) <- 'peaks'
# DefaultAssay(NA112_CTL_FILT) <- 'peaks'

DefaultAssay(NA105_P1CKO_FILT) <- 'peaks'
DefaultAssay(NA109_P1CKO_FILT) <- 'peaks'
DefaultAssay(NA113_P1CKO_FILT) <- 'peaks'

# DefaultAssay(NA106_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA110_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA114_P2CKO_FILT) <- 'peaks'

# DefaultAssay(NA107_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA111_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA115_P1P2CKO_FILT) <- 'peaks'


##-------------------------------------------------------
## Merge all 12 samples from four genotypes
# MERGED <- merge(x = NA104_CTL_FILT, y = list(NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT, NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT, NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT), add.cell.ids = c("NA104_CTL", "NA105_P1CKO", "NA106_P2CKO", "NA107_P1P2CKO", "NA108_CTL", "NA109_P1CKO", "NA110_P2CKO", "NA111_P1P2CKO", "NA112_CTL", "NA113_P1CKO", "NA114_P2CKO", "NA115_P1P2CKO"))
# MERGED.CTL <- merge(x = NA104_CTL_FILT, y = list(NA108_CTL_FILT, NA112_CTL_FILT), add.cell.ids = c("NA104_CTL", "NA108_CTL", "NA112_CTL"))
MERGED.P1CKO <- merge(x = NA105_P1CKO_FILT, y = list(NA109_P1CKO_FILT, NA113_P1CKO_FILT), add.cell.ids = c("NA105_P1CKO", "NA109_P1CKO", "NA113_P1CKO"))

print(table(MERGED.P1CKO$Dataset))
# NA105_P1CKO NA109_P1CKO NA113_P1CKO 
#        2656       15557        7095

save(MERGED.P1CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1CKO_TEMP.RData")

##-------------------------------------------------------
## Process merged dataset
MERGED.P1CKO <- RunTFIDF(MERGED.P1CKO)
MERGED.P1CKO <- FindTopFeatures(MERGED.P1CKO, min.cutoff = 20)
MERGED.P1CKO <- RunSVD(MERGED.P1CKO)
MERGED.P1CKO <- RunUMAP(MERGED.P1CKO, dims = 2:50, reduction = 'lsi')
MERGED.P1CKO <- FindNeighbors(MERGED.P1CKO, dims = 2:50, reduction = 'lsi')
MERGED.P1CKO <- FindClusters(MERGED.P1CKO, verbose = FALSE, algorithm = 3)

MERGED_UMAP <- DimPlot(MERGED.P1CKO, group.by = c("Dataset", "seurat_clusters"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1CKO.pdf", plot = MERGED_UMAP, width = 15, height = 6, units = "in", dpi = 300)

MERGED_UMAP_FACET <- DimPlot(MERGED.P1CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1CKO_FACET.pdf", plot = MERGED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)

save(MERGED.P1CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1CKO.RData")


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------





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
## Load filtered seurat/signac objects
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/E_FILT/NA_ATAC_SIGNAC_MACS2_FILT.RData")
# NA104_CTL_FILT, NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT,
# NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT,
# NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT,


##-------------------------------------------------------
## Default assay to peaks
# DefaultAssay(NA104_CTL_FILT) <- 'peaks'
# DefaultAssay(NA108_CTL_FILT) <- 'peaks'
# DefaultAssay(NA112_CTL_FILT) <- 'peaks'

# DefaultAssay(NA105_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA109_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA113_P1CKO_FILT) <- 'peaks'

DefaultAssay(NA106_P2CKO_FILT) <- 'peaks'
DefaultAssay(NA110_P2CKO_FILT) <- 'peaks'
DefaultAssay(NA114_P2CKO_FILT) <- 'peaks'

# DefaultAssay(NA107_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA111_P1P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA115_P1P2CKO_FILT) <- 'peaks'


##-------------------------------------------------------
## Merge all 12 samples from four genotypes
# MERGED <- merge(x = NA104_CTL_FILT, y = list(NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT, NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT, NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT), add.cell.ids = c("NA104_CTL", "NA105_P1CKO", "NA106_P2CKO", "NA107_P1P2CKO", "NA108_CTL", "NA109_P1CKO", "NA110_P2CKO", "NA111_P1P2CKO", "NA112_CTL", "NA113_P1CKO", "NA114_P2CKO", "NA115_P1P2CKO"))
# MERGED.CTL <- merge(x = NA104_CTL_FILT, y = list(NA108_CTL_FILT, NA112_CTL_FILT), add.cell.ids = c("NA104_CTL", "NA108_CTL", "NA112_CTL"))
# MERGED.P1CKO <- merge(x = NA105_P1CKO_FILT, y = list(NA109_P1CKO_FILT, NA113_P1CKO_FILT), add.cell.ids = c("NA105_P1CKO", "NA109_P1CKO", "NA113_P1CKO"))
MERGED.P2CKO <- merge(x = NA106_P2CKO_FILT, y = list(NA110_P2CKO_FILT, NA114_P2CKO_FILT), add.cell.ids = c("NA106_P2CKO", "NA110_P2CKO", "NA114_P2CKO"))

print(table(MERGED.P2CKO$Dataset))
# NA106_P2CKO NA110_P2CKO NA114_P2CKO 
#        2614        4765        6286

save(MERGED.P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P2CKO_TEMP.RData")

##-------------------------------------------------------
## Process merged dataset
MERGED.P2CKO <- RunTFIDF(MERGED.P2CKO)
MERGED.P2CKO <- FindTopFeatures(MERGED.P2CKO, min.cutoff = 20)
MERGED.P2CKO <- RunSVD(MERGED.P2CKO)
MERGED.P2CKO <- RunUMAP(MERGED.P2CKO, dims = 2:50, reduction = 'lsi')
MERGED.P2CKO <- FindNeighbors(MERGED.P2CKO, dims = 2:50, reduction = 'lsi')
MERGED.P2CKO <- FindClusters(MERGED.P2CKO, verbose = FALSE, algorithm = 3)

MERGED_UMAP <- DimPlot(MERGED.P2CKO, group.by = c("Dataset", "seurat_clusters"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P2CKO.pdf", plot = MERGED_UMAP, width = 15, height = 6, units = "in", dpi = 300)

MERGED_UMAP_FACET <- DimPlot(MERGED.P2CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P2CKO_FACET.pdf", plot = MERGED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)

save(MERGED.P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P2CKO.RData")


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------





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
## Load filtered seurat/signac objects
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/E_FILT/NA_ATAC_SIGNAC_MACS2_FILT.RData")
# NA104_CTL_FILT, NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT,
# NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT,
# NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT,


##-------------------------------------------------------
## Default assay to peaks
# DefaultAssay(NA104_CTL_FILT) <- 'peaks'
# DefaultAssay(NA108_CTL_FILT) <- 'peaks'
# DefaultAssay(NA112_CTL_FILT) <- 'peaks'

# DefaultAssay(NA105_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA109_P1CKO_FILT) <- 'peaks'
# DefaultAssay(NA113_P1CKO_FILT) <- 'peaks'

# DefaultAssay(NA106_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA110_P2CKO_FILT) <- 'peaks'
# DefaultAssay(NA114_P2CKO_FILT) <- 'peaks'

DefaultAssay(NA107_P1P2CKO_FILT) <- 'peaks'
DefaultAssay(NA111_P1P2CKO_FILT) <- 'peaks'
DefaultAssay(NA115_P1P2CKO_FILT) <- 'peaks'


##-------------------------------------------------------
## Merge all 12 samples from four genotypes
# MERGED <- merge(x = NA104_CTL_FILT, y = list(NA105_P1CKO_FILT, NA106_P2CKO_FILT, NA107_P1P2CKO_FILT, NA108_CTL_FILT, NA109_P1CKO_FILT, NA110_P2CKO_FILT, NA111_P1P2CKO_FILT, NA112_CTL_FILT, NA113_P1CKO_FILT, NA114_P2CKO_FILT, NA115_P1P2CKO_FILT), add.cell.ids = c("NA104_CTL", "NA105_P1CKO", "NA106_P2CKO", "NA107_P1P2CKO", "NA108_CTL", "NA109_P1CKO", "NA110_P2CKO", "NA111_P1P2CKO", "NA112_CTL", "NA113_P1CKO", "NA114_P2CKO", "NA115_P1P2CKO"))
# MERGED.CTL <- merge(x = NA104_CTL_FILT, y = list(NA108_CTL_FILT, NA112_CTL_FILT), add.cell.ids = c("NA104_CTL", "NA108_CTL", "NA112_CTL"))
# MERGED.P1CKO <- merge(x = NA105_P1CKO_FILT, y = list(NA109_P1CKO_FILT, NA113_P1CKO_FILT), add.cell.ids = c("NA105_P1CKO", "NA109_P1CKO", "NA113_P1CKO"))
# MERGED.P2CKO <- merge(x = NA106_P2CKO_FILT, y = list(NA110_P2CKO_FILT, NA114_P2CKO_FILT), add.cell.ids = c("NA106_P2CKO", "NA110_P2CKO", "NA114_P2CKO"))
MERGED.P1P2CKO <- merge(x = NA107_P1P2CKO_FILT, y = list(NA111_P1P2CKO_FILT, NA115_P1P2CKO_FILT), add.cell.ids = c("NA107_P1P2CKO", "NA111_P1P2CKO", "NA115_P1P2CKO"))

print(table(MERGED.P1P2CKO$Dataset))
# NA107_P1P2CKO NA111_P1P2CKO NA115_P1P2CKO 
#           587          8311          5346

save(MERGED.P1P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1P2CKO_TEMP.RData")

##-------------------------------------------------------
## Process merged dataset
MERGED.P1P2CKO <- RunTFIDF(MERGED.P1P2CKO)
MERGED.P1P2CKO <- FindTopFeatures(MERGED.P1P2CKO, min.cutoff = 20)
MERGED.P1P2CKO <- RunSVD(MERGED.P1P2CKO)
MERGED.P1P2CKO <- RunUMAP(MERGED.P1P2CKO, dims = 2:50, reduction = 'lsi')
MERGED.P1P2CKO <- FindNeighbors(MERGED.P1P2CKO, dims = 2:50, reduction = 'lsi')
MERGED.P1P2CKO <- FindClusters(MERGED.P1P2CKO, verbose = FALSE, algorithm = 3)

MERGED_UMAP <- DimPlot(MERGED.P1P2CKO, group.by = c("Dataset", "seurat_clusters"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1P2CKO.pdf", plot = MERGED_UMAP, width = 15, height = 6, units = "in", dpi = 300)

MERGED_UMAP_FACET <- DimPlot(MERGED.P1P2CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1P2CKO_FACET.pdf", plot = MERGED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)

save(MERGED.P1P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1P2CKO.RData")


##-------------------------------------------------------
## SEURAT ANALYSIS | END
##-------------------------------------------------------


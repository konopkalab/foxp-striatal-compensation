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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## Load merged signac/seurat object for genotype
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/F_MERGE/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_CTL.RData")
# MERGED.CTL


##-------------------------------------------------------
## Batch correction using harmony
DefaultAssay(MERGED.CTL) <- "peaks"
HARMONISED.CTL <- RunHarmony(object = MERGED.CTL, group.by.vars = c("Batch", "Sex"), reduction = "lsi", assay.use = "peaks", max.iter.harmony = 100, project.dim = F)

HARMONISED.CTL <- RunUMAP(HARMONISED.CTL, dims = 2:50, reduction = "harmony")
HARMONISED.CTL <- FindNeighbors(HARMONISED.CTL, dims = 2:50, reduction = 'lsi')
HARMONISED.CTL <- FindClusters(HARMONISED.CTL, verbose = FALSE, algorithm = 3)

HARMONISED_UMAP <- DimPlot(HARMONISED.CTL, group.by = c("Dataset", "seurat_clusters", "Batch"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_CTL.pdf", plot = HARMONISED_UMAP, width = 21, height = 6, units = "in", dpi = 300)

HARMONISED_UMAP_FACET <- DimPlot(HARMONISED.CTL, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_CTL_FACET.pdf", plot = HARMONISED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)


save(HARMONISED.CTL, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_CTL.RData")


##-------------------------------------------------------
## END
##-------------------------------------------------------




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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## Load merged signac/seurat object for genotype
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/F_MERGE/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1CKO.RData")
# MERGED.P1CKO


##-------------------------------------------------------
## Batch correction using harmony
DefaultAssay(MERGED.P1CKO) <- "peaks"
HARMONISED.P1CKO <- RunHarmony(object = MERGED.P1CKO, group.by.vars = c("Batch", "Sex"), reduction = "lsi", assay.use = "peaks", max.iter.harmony = 100, project.dim = F)

HARMONISED.P1CKO <- RunUMAP(HARMONISED.P1CKO, dims = 2:50, reduction = "harmony")
HARMONISED.P1CKO <- FindNeighbors(HARMONISED.P1CKO, dims = 2:50, reduction = 'lsi')
HARMONISED.P1CKO <- FindClusters(HARMONISED.P1CKO, verbose = FALSE, algorithm = 3)

HARMONISED_UMAP <- DimPlot(HARMONISED.P1CKO, group.by = c("Dataset", "seurat_clusters", "Batch"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1CKO.pdf", plot = HARMONISED_UMAP, width = 21, height = 6, units = "in", dpi = 300)

HARMONISED_UMAP_FACET <- DimPlot(HARMONISED.P1CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1CKO_FACET.pdf", plot = HARMONISED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)


save(HARMONISED.P1CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P1CKO.RData")


##-------------------------------------------------------
## END
##-------------------------------------------------------




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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## Load merged signac/seurat object for genotype
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/F_MERGE/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P2CKO.RData")
# MERGED.P2CKO


##-------------------------------------------------------
## Batch correction using harmony
DefaultAssay(MERGED.P2CKO) <- "peaks"
HARMONISED.P2CKO <- RunHarmony(object = MERGED.P2CKO, group.by.vars = c("Batch", "Sex"), reduction = "lsi", assay.use = "peaks", max.iter.harmony = 100, project.dim = F)

HARMONISED.P2CKO <- RunUMAP(HARMONISED.P2CKO, dims = 2:50, reduction = "harmony")
HARMONISED.P2CKO <- FindNeighbors(HARMONISED.P2CKO, dims = 2:50, reduction = 'lsi')
HARMONISED.P2CKO <- FindClusters(HARMONISED.P2CKO, verbose = FALSE, algorithm = 3)

HARMONISED_UMAP <- DimPlot(HARMONISED.P2CKO, group.by = c("Dataset", "seurat_clusters", "Batch"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P2CKO.pdf", plot = HARMONISED_UMAP, width = 21, height = 6, units = "in", dpi = 300)

HARMONISED_UMAP_FACET <- DimPlot(HARMONISED.P2CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P2CKO_FACET.pdf", plot = HARMONISED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)


save(HARMONISED.P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P2CKO.RData")


##-------------------------------------------------------
## END
##-------------------------------------------------------





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
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## Load merged signac/seurat object for genotype
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/F_MERGE/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_P1P2CKO.RData")
# MERGED.P1P2CKO


##-------------------------------------------------------
## Batch correction using harmony
DefaultAssay(MERGED.P1P2CKO) <- "peaks"
HARMONISED.P1P2CKO <- RunHarmony(object = MERGED.P1P2CKO, group.by.vars = c("Batch", "Sex"), reduction = "lsi", assay.use = "peaks", max.iter.harmony = 100, project.dim = F)

HARMONISED.P1P2CKO <- RunUMAP(HARMONISED.P1P2CKO, dims = 2:50, reduction = "harmony")
HARMONISED.P1P2CKO <- FindNeighbors(HARMONISED.P1P2CKO, dims = 2:50, reduction = 'lsi')
HARMONISED.P1P2CKO <- FindClusters(HARMONISED.P1P2CKO, verbose = FALSE, algorithm = 3)

HARMONISED_UMAP <- DimPlot(HARMONISED.P1P2CKO, group.by = c("Dataset", "seurat_clusters", "Batch"), pt.size = 0.1)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1P2CKO.pdf", plot = HARMONISED_UMAP, width = 21, height = 6, units = "in", dpi = 300)

HARMONISED_UMAP_FACET <- DimPlot(HARMONISED.P1P2CKO, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_MERGED_P1P2CKO_FACET.pdf", plot = HARMONISED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)


save(HARMONISED.P1P2CKO, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P1P2CKO.RData")


##-------------------------------------------------------
## END
##-------------------------------------------------------



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


# ##-------------------------------------------------------
# ## Load harmonized signac/seurat object for all genotypes
# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/G_HARMONY/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_CTL.RData")
# # HARMONISED.CTL

# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/G_HARMONY/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P1CKO.RData")
# # HARMONISED.P1CKO

# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/G_HARMONY/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P2CKO.RData")
# # HARMONISED.P2CKO

# load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/06_MACS2_SIGNAC_HARMONY_JAN2023/G_HARMONY/BATCH_SEX/NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_P1P2CKO.RData")
# # HARMONISED.P1P2CKO

# DefaultAssay(HARMONISED.CTL) <- "peaks"
# DefaultAssay(HARMONISED.P1CKO) <- "peaks"
# DefaultAssay(HARMONISED.P2CKO) <- "peaks"
# DefaultAssay(HARMONISED.P1P2CKO) <- "peaks"


# ##-------------------------------------------------------
# ## Merge
# MERGED <- merge(x = HARMONISED.CTL, y = list(HARMONISED.P1CKO, HARMONISED.P2CKO, HARMONISED.P1P2CKO))


# MERGED <- RunTFIDF(MERGED)
# MERGED <- FindTopFeatures(MERGED, min.cutoff = 20)
# MERGED <- RunSVD(MERGED)
# MERGED <- RunUMAP(MERGED, dims = 2:50, reduction = 'lsi')
# MERGED <- FindNeighbors(MERGED, dims = 2:50, reduction = 'lsi')
# MERGED <- FindClusters(MERGED, verbose = FALSE, algorithm = 3)

# MERGED_UMAP <- DimPlot(MERGED, group.by = c("Dataset", "seurat_clusters"), pt.size = 0.1)
# ggsave("NA_ATACSEQ_UMAP_MERGED_ALL.pdf", plot = MERGED_UMAP, width = 15, height = 6, units = "in", dpi = 300)

# MERGED_UMAP_FACET <- DimPlot(MERGED, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
# ggsave("NA_ATACSEQ_UMAP_MERGED_ALL_FACET.pdf", plot = MERGED_UMAP_FACET, width = 20, height = 8, units = "in", dpi = 300)

# ##-------------------------------------------------------
# ## save merged seurat objects
# save(MERGED, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_MERGED.RData")






##-------------------------------------------------------
## Load merged data
rm(list = ls())
load("NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_MERGED.RData")
# MERGED

MERGED$Genotype <- gsub("^NA104_CTL$|^NA108_CTL$|^NA112_CTL$", "CTL", MERGED$Dataset)
MERGED$Genotype <- gsub("^NA105_P1CKO$|^NA109_P1CKO$|^NA113_P1CKO$", "P1CKO", MERGED$Genotype)
MERGED$Genotype <- gsub("^NA106_P2CKO$|^NA110_P2CKO$|^NA114_P2CKO$", "P2CKO", MERGED$Genotype)
MERGED$Genotype <- gsub("^NA107_P1P2CKO$|^NA111_P1P2CKO$|^NA115_P1P2CKO$", "P1P2CKO", MERGED$Genotype)


##-------------------------------------------------------
## Integration
NA_ATAC_LIST <- SplitObject(MERGED, split.by = "Genotype")


# anchfeat <- unique(sort(c(row.names(HARMONISED.CTL), row.names(HARMONISED.P1CKO), row.names(HARMONISED.P2CKO), row.names(HARMONISED.P1P2CKO))))

## Find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = NA_ATAC_LIST,
  anchor.features = row.names(NA_ATAC_LIST[["CTL"]]),
  reduction = "rlsi",
  dims = 2:30
)


# Integrate LSI embeddings
INTEGRATED <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = MERGED[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
INTEGRATED <- RunUMAP(INTEGRATED, dims = 2:30, reduction = 'integrated_lsi')
INTEGRATED <- FindNeighbors(INTEGRATED, dims = 2:30, reduction = 'integrated_lsi')
INTEGRATED <- FindClusters(INTEGRATED, verbose = FALSE, algorithm = 3)

##-------------------------------------------------------
## save integrated seurat objects
save(INTEGRATED, file = "NA_ATAC_SIGNAC_MACS2_FILT_MERGED_HARMONY_INTEGRATED.RData")


##-------------------------------------------------------
## UMAP Plots
INT_UMAP <- DimPlot(INTEGRATED, group.by = c("Dataset", "Batch", "Sex", "Genotype"), pt.size = 0.1, ncol = 2)
ggsave("NA_ATACSEQ_UMAP_INTEGRATED_ALL.pdf", plot = INT_UMAP, width = 15, height = 12, units = "in", dpi = 300)

INT_UMAP_FACET <- DimPlot(INTEGRATED, group.by = "Dataset", pt.size = 0.1, split.by = "Dataset", ncol = 4)
ggsave("NA_ATACSEQ_UMAP_INTEGRATED_ALL_FACET.pdf", plot = INT_UMAP_FACET, width = 21, height = 15, units = "in", dpi = 300)


##-------------------------------------------------------
## END
##-------------------------------------------------------

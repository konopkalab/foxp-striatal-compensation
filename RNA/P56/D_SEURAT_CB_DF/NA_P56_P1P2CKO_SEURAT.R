##------------------------------------
## NA | WT | MAR 2021
##------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl

##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(harmony)
library(scCustomize)
set.seed(10)


##------------------------------------
## reference gene symbol list
ref <- read.table("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P09/00_REF/gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## load CellBender + DoubletFinder dataset
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/11_DOUBLETFINDER_P1P2CKO/SEURAT_NA_P56_P1P2CKO_CB_DF.RData")
# seuObjDF

table(seuObjDF$DoubletFinder)
# Discarded  Retained 
#      9976     89788

seuObjDFCleaned <- subset(seuObjDF, subset = DoubletFinder == "Retained")

table(seuObjDFCleaned$GenoAgeSample)
#  P1P2CKO_P56_NA02 P1P2CKO_P56_NA04A P1P2CKO_P56_NA04B  P1P2CKO_P56_NA09 
#              6762             25423             22021             13623 
#  P1P2CKO_P56_NA45 
#             21959

## cells retained
cells2keepTemp <- row.names(seuObjDFCleaned@meta.data[seuObjDFCleaned@meta.data$DoubletFinder == "Retained",])
cells2keepNA02 <- cells2keepTemp[grepl("^P56_P1P2CKO_NA02_", cells2keepTemp)]
cells2keepNA04A <- cells2keepTemp[grepl("^P56_P1P2CKO_NA04A_", cells2keepTemp)]
cells2keepNA04B <- cells2keepTemp[grepl("^P56_P1P2CKO_NA04B_", cells2keepTemp)]
cells2keepNA09 <- cells2keepTemp[grepl("^P56_P1P2CKO_NA09_", cells2keepTemp)]
cells2keepNA45 <- cells2keepTemp[grepl("^P56_P1P2CKO_NA45_", cells2keepTemp)]

cells2keepNA04A <- gsub("P56_P1P2CKO_NA04A_", "P56_P1P2CKO_NA04_", cells2keepNA04A)
cells2keepNA04B <- gsub("P56_P1P2CKO_NA04B_", "P56_P1P2CKO_NA04_", cells2keepNA04B)

cells2keep <- c(cells2keepNA02, unique(sort(c(cells2keepNA04A, cells2keepNA04B))), cells2keepNA09, cells2keepNA45)


## load dataset from 10X run
p1p1p2cko.na02 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/02_CELLBENDER/P1P2CKO_NA02/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
print(dim(p1p1p2cko.na02)) # 53715 8153

p1p1p2cko.na04a <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/02_CELLBENDER/P1P2CKO_NA04_01/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
print(dim(p1p1p2cko.na04a)) # 53715 27437

p1p1p2cko.na04b <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/02_CELLBENDER/P1P2CKO_NA04_02/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
print(dim(p1p1p2cko.na04b)) # 53715 28536

p1p1p2cko.na09 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/P56/02_CELLBENDER/P1P2CKO_NA09/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
print(dim(p1p1p2cko.na09)) # 53715 13623

p1p1p2cko.na45 <- Read10X_h5("/project/Neuroinformatics_Core/Konopka_lab/akulk1/NEW_SEP2022_NK_NA_EO_CO/03_CELLBENDER_NA/NA45_CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
print(dim(p1p1p2cko.na45)) # 53715 22015

## adding sample & condition info to column names
colnames(p1p1p2cko.na02) <- paste("P56_P1P2CKO_NA02", colnames(p1p1p2cko.na02), sep = "_")
colnames(p1p1p2cko.na04a) <- paste("P56_P1P2CKO_NA04", colnames(p1p1p2cko.na04a), sep = "_")
colnames(p1p1p2cko.na04b) <- paste("P56_P1P2CKO_NA04", colnames(p1p1p2cko.na04b), sep = "_")
colnames(p1p1p2cko.na09) <- paste("P56_P1P2CKO_NA09", colnames(p1p1p2cko.na09), sep = "_")
colnames(p1p1p2cko.na45) <- paste("P56_P1P2CKO_NA45", colnames(p1p1p2cko.na45), sep = "_")


## fetch genes or rows corresponding to gencode protein coding gene symbols
p1p1p2cko.na02.temp <- p1p1p2cko.na02[row.names(p1p1p2cko.na02) %in% ref$GeneSymbol,]
p1p1p2cko.na04a.temp <- p1p1p2cko.na04a[row.names(p1p1p2cko.na04a) %in% ref$GeneSymbol,]
p1p1p2cko.na04b.temp <- p1p1p2cko.na04b[row.names(p1p1p2cko.na04b) %in% ref$GeneSymbol,]
p1p1p2cko.na09.temp <- p1p1p2cko.na09[row.names(p1p1p2cko.na09) %in% ref$GeneSymbol,]
p1p1p2cko.na45.temp <- p1p1p2cko.na45[row.names(p1p1p2cko.na45) %in% ref$GeneSymbol,]


## make a data from from matrix
P56_P1P2CKO_NA02 <- as.data.frame(as.matrix(p1p1p2cko.na02.temp))
P56_P1P2CKO_NA04_A <- as.data.frame(as.matrix(p1p1p2cko.na04a.temp))
P56_P1P2CKO_NA04_B <- as.data.frame(as.matrix(p1p1p2cko.na04b.temp))
P56_P1P2CKO_NA09 <- as.data.frame(as.matrix(p1p1p2cko.na09.temp))
P56_P1P2CKO_NA45 <- as.data.frame(as.matrix(p1p1p2cko.na45.temp))


## add gene symbol as a column
P56_P1P2CKO_NA02$Genes <- row.names(P56_P1P2CKO_NA02)
P56_P1P2CKO_NA04_A$Genes <- row.names(P56_P1P2CKO_NA04_A)
P56_P1P2CKO_NA04_B$Genes <- row.names(P56_P1P2CKO_NA04_B)
P56_P1P2CKO_NA09$Genes <- row.names(P56_P1P2CKO_NA09)
P56_P1P2CKO_NA45$Genes <- row.names(P56_P1P2CKO_NA45)



## collapse the UMI counts for common cells by adding the counts for NA04
commonNA04 <- intersect(colnames(P56_P1P2CKO_NA04_A), colnames(P56_P1P2CKO_NA04_B))
P56_P1P2CKO_NA04_A_COMMON <- P56_P1P2CKO_NA04_A[,colnames(P56_P1P2CKO_NA04_A) %in% commonNA04]
P56_P1P2CKO_NA04_B_COMMON <- P56_P1P2CKO_NA04_B[,colnames(P56_P1P2CKO_NA04_B) %in% commonNA04]

dataCommonNA04 <- list("P56_P1P2CKO_NA04_A" = P56_P1P2CKO_NA04_A_COMMON, "P56_P1P2CKO_NA04_B" = P56_P1P2CKO_NA04_B_COMMON)

combinedDataNA04 <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCommonNA04)
row.names(combinedDataNA04) <- combinedDataNA04$Genes
combinedDataNA04$Genes <- NULL
# colnames(combinedDataNA04) <- substr(colnames(combinedDataNA04), 1, nchar(colnames(combinedDataNA04)) - 2)

combinedMetaNA04 <- as.data.frame(matrix(unlist(strsplit(colnames(combinedDataNA04), "_")), ncol = 4, byrow = TRUE))
row.names(combinedMetaNA04) <- colnames(combinedDataNA04)
colnames(combinedMetaNA04) <- c("Age", "Genotype", "Library", "CellBarcode")
combinedMetaNA04$CellBarcode <- substr(combinedMetaNA04$CellBarcode, 1, nchar(combinedMetaNA04$CellBarcode) - 2)

newDataTempNA04 <- as.data.frame(t(combinedDataNA04))
newDataNA04 <- merge(newDataTempNA04, combinedMetaNA04, by = "row.names")
row.names(newDataNA04) <- newDataNA04$Row.names
newDataNA04$Row.names <- NULL
cellBC <- list(newDataNA04$CellBarcode)
newDataNA04$Genotype <- NULL
newDataNA04$Library <- NULL
newDataNA04$CellBarcode <- NULL
newDataNA04$Age <- NULL

newDataNA04AggrTemp <- aggregate(newDataNA04, by = cellBC, FUN = "sum")
row.names(newDataNA04AggrTemp) <- newDataNA04AggrTemp$Group.1
newDataNA04AggrTemp$Group.1 <- NULL

P56_P1P2CKO_NA04 <- as.data.frame(t(newDataNA04AggrTemp))
colnames(P56_P1P2CKO_NA04) <- paste("P56_P1P2CKO_NA04", colnames(P56_P1P2CKO_NA04), sep = "_")
P56_P1P2CKO_NA04$Genes <- row.names(P56_P1P2CKO_NA04)


## combine individual tables into a giant data frame
dataCombined <- list("P56_P1P2CKO_NA02" = P56_P1P2CKO_NA02, "P56_P1P2CKO_NA04" = P56_P1P2CKO_NA04, "P56_P1P2CKO_NA09" = P56_P1P2CKO_NA09, "P56_P1P2CKO_NA45" = P56_P1P2CKO_NA45)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

na.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

na.data[is.na(na.data)] <- 0
print(dim(na.data)) ## 21967 71188


p56.p1p1p2cko.data <- na.data[,colnames(na.data) %in% cells2keep]
print(dim(p56.p1p1p2cko.data))
# 21967 67734


save(p56.p1p1p2cko.data, file = "NA_SEURAT_P1P2CKO_P56_CB_DF_COUNTS.RData")




##------------------------------------
rm(list = ls())
load("NA_SEURAT_P1P2CKO_P56_CB_DF_COUNTS.RData")
# p56.p1p1p2cko.data

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.p1p1p2cko.data, project = "NA_P1P2CKO_P56")
print(dim(seuObj)) ## 21967 67734


##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 4, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Genotype", "Age", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

seuObj$Genotype <- seuObj$orig.ident
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")


seuObj$LibPrep <- gsub("P56_P1P2CKO_NA02", 1, seuObj$GenoAgeSample)
seuObj$LibPrep <- gsub("P56_P1P2CKO_NA04", 2, seuObj$LibPrep)
seuObj$LibPrep <- gsub("P56_P1P2CKO_NA09", 4, seuObj$LibPrep)
seuObj$LibPrep <- gsub("P56_P1P2CKO_NA45", 11, seuObj$LibPrep)

seuObj$SeqBatch <- gsub("P56_P1P2CKO_NA02", 1, seuObj$GenoAgeSample)
seuObj$SeqBatch <- gsub("P56_P1P2CKO_NA04", 1, seuObj$SeqBatch)
seuObj$SeqBatch <- gsub("P56_P1P2CKO_NA09", 3, seuObj$SeqBatch)
seuObj$SeqBatch <- gsub("P56_P1P2CKO_NA45", 10, seuObj$SeqBatch)

seuObj$Sex <- gsub("P56_P1P2CKO_NA02", "M", seuObj$GenoAgeSample)
seuObj$Sex <- gsub("P56_P1P2CKO_NA04", "F", seuObj$Sex)
seuObj$Sex <- gsub("P56_P1P2CKO_NA09", "M", seuObj$Sex)
seuObj$Sex <- gsub("P56_P1P2CKO_NA45", "F", seuObj$Sex)

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

print(table(seuObj@active.ident))


## Visualize Data QC
p1 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_NA_QC_1.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NA_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NA_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

# ##------------------------------------
# ## Data Filtering
# seuObj <- subset(seuObj, subset = nCount_RNA < 3000 & pMito_RNA < 0.1)
# print(dim(seuObj)) ## 21967  13943

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)

ggsave(filename = "SEURAT_NA_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NA_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Scaling
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_NA_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NA_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)


save(seuObj, file = "NA_SEURAT_P1P2CKO_P56_TEMP.RData")

##------------------------------------
## Harmony
seuObj <- RunHarmony(object = seuObj, group.by.vars = c("LibPrep", "SeqBatch", "Sex"), reduction = "pca", assay.use = "RNA") #, project.dim = FALSE)


##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)


seuObj <- FindNeighbors(seuObj, dims = 1:selpcs,reduction = "harmony")
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), reduction = "harmony")
seuObj <- RunUMAP(seuObj, dims = 1:selpcs, reduction = "harmony")


##------------------------------------
## Save RData
save(seuObj, file = "NA_SEURAT_P1P2CKO_P56.RData")



HARMONISED_UMAP <- DimPlot(seuObj, group.by = c("GenoAgeSample", "LibPrep", "SeqBatch", "Sex", "RNA_snn_res.0.8"), pt.size = 0.01)
ggsave("NA_SEURAT_P1P2CKO_P56_UMAP.pdf", plot = HARMONISED_UMAP, width = 20, height = 12, units = "in", dpi = 300)

HARMONISED_UMAP_FACET <- DimPlot(seuObj, group.by = "GenoAgeSample", pt.size = 0.01, split.by = "GenoAgeSample", ncol = 3)
ggsave("NA_SEURAT_P1P2CKO_P56_UMAP_FACET.pdf", plot = HARMONISED_UMAP_FACET, width = 20, height = 6, units = "in", dpi = 300)









# ##------------------------------------
# ## Plots

# seutree <- clustree(seuObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
# ggsave(filename = "NA_SEURAT_PLOT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


# ResolutionList <- grep("RNA_snn_res", colnames(seuObj@meta.data), value = TRUE)

# for (Resolution in ResolutionList)
#     {
#     print(paste("====> RESOLUTION ", Resolution, sep = ""))

#     pdf(paste0("SEURAT_NA_INTEGRATE_UMAP_RES_", Resolution, ".pdf"), width=7, height=6)
#     g <- DimPlot(object = seuObj, label = TRUE, reduction = "umap", group.by = Resolution)
#     print(g)
#     dev.off()

#     pdf(paste0("SEURAT_NA_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width=7, height=3)
#     v <- VlnPlot(object = seuObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
#     print(v)
#     dev.off()
#     }


# ##------------------------------------
# ## Save RData
# save(seuObj, file = "NA_SEURAT_P1P2CKO_P56.RData")














# Idents(seuObj) <- "RNA_snn_res.0.8"

# myres <- 0.8
# dir.cbeate(paste("RES_", myres, sep = ""))
# setwd(paste("RES_", myres, sep = ""))



# pumap1a <- DimPlot(seuObj, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE, label.size = 4, repel = TRUE)
# ggsave("NA_SEURAT_PLOT_UMAP_CLUSTERS.PDF", plot = pumap1a, width = 7, height = 6, units = "in", dpi = 150)

# pumap1b <- DimPlot(seuObj, reduction = "umap", group.by = "RNA_snn_res.0.8") + facet_wrap(~RNA_snn_res.0.8, nrow = 5)
# ggsave("NA_SEURAT_PLOT_UMAP_CLUSTERS_FACET.PDF", plot = pumap1b, width = 24, height = 24, units = "in", dpi = 150)

# pumap5a <- DimPlot(seuObj, reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 4, repel = TRUE)
# ggsave("NA_SEURAT_PLOT_UMAP_GENOTYPE.PDF", plot = pumap5a, width = 7, height = 6, units = "in", dpi = 150)

# pumap5b <- DimPlot(seuObj, reduction = "umap", group.by = "orig.ident") + facet_wrap(~orig.ident, nrow = 1)
# ggsave("NA_SEURAT_PLOT_UMAP_GENOTYPE_FACET.PDF", plot = pumap5b, width = 12, height = 4, units = "in", dpi = 150)

# pumap6 <- DimPlot(seuObj, reduction = "umap", group.by = "Sample") + facet_wrap(~Sample, nrow = 3)
# ggsave("NA_SEURAT_PLOT_UMAP_SAMPLE_FACET.PDF", plot = pumap6, width = 12, height = 12, units = "in", dpi = 150)



# ## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 0.8
# cellsPerCluster <- as.data.frame.matrix(table(seuObj@meta.data$RNA_snn_res.0.8, seuObj@meta.data$GenoAgeSample))
# cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
# write.table(cellsPerCluster, paste(seuObj@project.name, "GENOTYPE_PER_CLUSTER_0.8.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

# cellsPerCluster2 <- melt(cellsPerCluster)
# colnames(cellsPerCluster2) <- c("CLUSTER", "GENOTYPE", "CELLS")
# p7 <- ggplot(cellsPerCluster2) +
#         geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") + 
#         scale_y_continuous(labels = scales::percent_format()) + 
#         theme_classic() + 
#         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#         NULL
# ggsave(filename = paste(seuObj@project.name, "GENOTYPE_PER_CLUSTER_0.8.pdf", sep = "_"), plot = p7, width = 8, height = 4, units = "in", dpi = 150)


# ## VIOLINPLOT AND DOTPLOT | RES 0.8
# vpl1 <- VlnPlot(object = seuObj, features = mygenes, ncol = 4, pt.size = 0, group.by = "RNA_snn_res.0.8", assay = "RNA", slot = "data", log = TRUE)
# ggsave(filename = "NA_SEURAT_ViolinPlot_Log_0.8.pdf", plot = vpl1, width = 24, height = 16, units = "in", dpi = 150)

# dpl1 <- DotPlot(object = seuObj, features = rev(mygenes), dot.scale = 10, group.by = "RNA_snn_res.0.8", assay = "RNA") + theme(axis.text.x = element_text(angle = 90))
# ggsave(filename = "NA_SEURAT_DotPlot_0.8.pdf", plot = dpl1, width = 8, height = 16, units = "in", dpi = 150)


# vpl2 <- VlnPlot(object = seuObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = "RNA_snn_res.0.8")
# ggsave(filename = "NA_SEURAT_ViolinPlot_nUMI_0.8.pdf", plot = vpl2, width = 16, height = 4, units = "in", dpi = 150)

# setwd("./../")



# ##------------------------------------
# ## Find Cluster Markers
# rm(list = ls())
# load("NA_SEURAT_P1P2CKO_P1.RData")
# # seuObj

# Idents(seuObj) <- "RNA_snn_res.0.8"

# seuObj.markers <- FindAllMarkers(seuObj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
# write.table(seuObj.markers, "NA_SEURAT_RES_0.8_CLUSTER_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")



# ##------------------------------------
# ## Genes of interest
# rm(list = ls())
# load("NA_SEURAT_P1P2CKO_P1.RData")
# # seuObj

# myres <- 0.8
# setwd(paste("RES_", myres, sep = ""))

# Idents(seuObj) <- "RNA_snn_res.0.8"

# table(seuObj@active.ident)
# #    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# # 2960 2946 1516 1134 1030  973  796  749  731  678  670  644  569  538  335  311 
# #   16   17   18   19   20   21   22   23   24 
# #  288  167  141  134  117  112   83   62   61

# cellspercluster <- as.data.frame(table(seuObj@active.ident))
# colnames(cellspercluster) <- c("Cluster", "Cells")
# cellspercluster <- cellspercluster[order(cellspercluster$Cells, decreasing = TRUE),]
# write.table(cellspercluster, paste("SEURAT_NA_CELLS_PER_CLUSTER_", myres, ".tsv", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# # Aqp4 for astrocytes, 
# # Olig1 for OPCs, 
# # Cx3cr1 for microglia, 
# # Flt2 for endothelial, 
# # Slc17a7 for glutamatergic cortical neurons, 
# # interneuron populations (Chat, Npy), 
# # neurogenic and neural differentiation marker (Sox4, Sox11). 
# # proliferating cells (Mki67), 
# # neural progenitors (Ascl1), 
# # neural progenitors derived from the lateral ganglionic eminence (Dlx2). 
# # genes important for iSPN specification (Sp9), 
# # mature SPN marker (Ppp1r1b), 
# # major SPN subtypes (Drd1, Drd2)

# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Slc17a7", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")

# pvln <- VlnPlot(seuObj, features = mygenes, pt.size = 0, ncol = 1, group.by = "RNA_snn_res.0.8")
# ggsave(filename = paste("SEURAT_NA_GENES_VIOLINPLOT_", myres, "_1.pdf", sep = ""), plot = pvln, width = 6, height = 36, units = "in", dpi = 150)

# pftr1 <- FeaturePlot(seuObj, features = mygenes, order = FALSE, pt.size = 0.25, min.cutoff = 0, max.cutoff = 3)
# ggsave(filename = paste("SEURAT_NA_GENES_FEATUREPLOT_", myres, "_1.pdf", sep = ""), plot = pftr1, width = 12, height = 12, units = "in", dpi = 150)

# pftr2 <- FeaturePlot(seuObj, features = mygenes, order = FALSE, pt.size = 0.25, min.cutoff = 0, max.cutoff = 3)
# ggsave(filename = paste("SEURAT_NA_GENES_FEATUREPLOT_", myres, "_2.pdf", sep = ""), plot = pftr2, width = 12, height = 12, units = "in", dpi = 150)

# pdot1 <- DotPlot(seuObj, features = rev(mygenes), cols = c("white", "blue"), dot.scale = 8, col.min = 0, col.max = 3, group.by = "RNA_snn_res.0.8") + RotatedAxis()
# ggsave(filename = paste("SEURAT_NA_GENES_DOTPLOT_", myres, "_1.pdf", sep = ""), plot = pdot1, width = 8, height = 8, units = "in", dpi = 150)

# setwd("./../")




# ##------------------------------------
# ## Annotate Clusters

# rm(list = ls())
# library(Seurat)
# library(dplyr)
# library(Matrix)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
# library(ggrepel)
# library(reticulate)
# library(WGCNA)


# # load("NA_SEURAT.RData")
# # # seuObj

# cluMarkers <- read.table("NA_SEURAT_RES_0.8_CLUSTER_MARKERS.txt", sep = "\t", header = TRUE)
# tab <- cluMarkers
# tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_logFC >= 0.25 & tab$pct.1 >= 0.3,]
# # tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_logFC >= 0.3,]
# # tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.3,]
# tab <- tab[c(7,6)]
# tab$cluster <- as.factor(paste("Cluster_", sprintf("%02d", as.numeric(as.character(tab$cluster))), sep=""))
# colnames(tab)=c("Gene","DEFINITION")
# Genes=as.data.frame(table(tab$DEFINITION))


# load("./../Saunders2018DEG.RData")
# # degSaunders
# GeneSets <- degSaunders


# for(i in 1:length(GeneSets)){
# 	colnames(GeneSets[[i]])[1] <- "Gene"
# }

# ln=length(GeneSets)
# cl=length(Genes$Var1)
# TEMP=list()
# INT=list()
# for (i in 1:ln)
# {
# TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
# INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
# }
# names(INT)=names(GeneSets)
# names(TEMP)=names(GeneSets)

# # Create the matrix for each GeneSet
# NROWS <- sapply(GeneSets,nrow)

# #
# #               GeneSet
# #              +       -
# #          +-------+-------+
# #        + |   a   |   b   |
# #  Module  +-------+-------+
# #        - |   c   |   d   |
# #          +-------+-------+

# for (i in 1:length(INT))
# {
# INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
# INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
# INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585 #13517 #nrow(tab) #8321 # nrow(cluMarkers)(genes in seurat object)
# }

# # sum(Genes$Freq)
# RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
# 	f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
# 	return(c(row,
# 			P_val = f$p.value,
# 			LogP = -log10(f$p.value), 
# 			OR = f$estimate[[1]],
# 			OR_Low = f$conf.int[1],
# 			OR_Up = f$conf.int[2]))
# }

# # run
# FisherMat=list()
# for (i in 1:length(INT))
# {
# FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
# rownames(FisherMat[[i]]) <- INT[[i]]$Var1
# FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
# }
# names(FisherMat)<-names(INT)
# save(FisherMat, TEMP, file= "NA_SEURAT_RES_0.8_FisherOutput_Saunders_Enrich.RData")


# # Create matrices of Pval
# tmp<-list()
# FisherP<-matrix()
# rowNames <- rownames(FisherMat[[1]])
# colNames <- names(FisherMat)
# for (i in 1:length(INT))
# {
# tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
# FisherP <- do.call(cbind,tmp)
# }
# rownames(FisherP) <- rowNames
# colnames(FisherP) <- colNames

# # Create matrices of OR
# tmp<-list()
# FisherOR<-matrix()
# rowNames <- rownames(FisherMat[[1]])
# colNames <- names(FisherMat)
# for (i in 1:length(INT))
# {
# tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
# FisherOR <- do.call(cbind,tmp)
# }
# rownames(FisherOR) <- rowNames
# colnames(FisherOR) <- colNames

# # Pval Adjusted
# library(magrittr)
# FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
# rownames(FisherAdj) <- rowNames
# colnames(FisherAdj) <- colNames

# FisherAdj[FisherAdj>0.05]=1
# FisherOR[FisherOR < 1]=0


# pdf("NA_SEURAT_RES_0.8_FisherOutput_Saunders_Enrich.pdf", width=12, height=12, pointsize=12)
# par(mar=c(15, 7, 2, 2))
# df=-log10(FisherAdj)
# LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
# LabelMat[LabelMat == "0\n(1)"] <- NA
# labeledHeatmap(Matrix = df, 
# xLabels = colnames(df), 
# yLabels = rownames(df), 
# colorLabels =FALSE,
# colors=colorRampPalette(c("white", "red"))(50),
# textMatrix=LabelMat, 
# setStdMargins = FALSE, 
# cex.text = 0.5,
# xLabelsAngle = 90)
# dev.off()



# FisherORt <- as.data.frame(t(FisherOR))
# colnames(FisherORt) <- paste(colnames(FisherORt), "OR", sep = "_")

# FisherAdjt <- as.data.frame(t(FisherAdj))
# colnames(FisherAdjt) <- paste(colnames(FisherAdjt), "Pval", sep = "_")

# FisherData <- merge(FisherORt, FisherAdjt, by = "row.names")
# row.names(FisherData) <- FisherData$Row.names
# FisherData$Row.names <- NULL

# FisherData2 <- FisherData[,order(colnames(FisherData))]
# write.table(FisherData2, "NA_SEURAT_RES_0.8_FisherOutput_Saunders_PlotData.txt", row.names = T, col.names = T, quote = F, sep = "\t")


# sink("NA_SEURAT_RES_0.8_AnnotationData.txt")
# aa <- FisherData2
# for(i in 1:ncol(aa))
#  {	 
#  if(i %% 2 != 0)
#   {
#   cluor <- i
#   clup <- i + 1
#   print("------------------")
#   bb <- aa[,c(cluor, clup)]
#   cc <- bb[order(bb[,1], decreasing = T),]
#   dd <- bb[order(bb[,2], decreasing = F),]
#   print(gsub("_OR", "", colnames(cc)[1]))
#   print(cc[1:2,])
#   print(dd[1:2,])
#   print("------------------")
#   }
#  }

# sink()



# ##------------------------------------
# ## Update Seurat Object
# rm(list = ls())
# library(Seurat)
# library(dplyr)
# library(Matrix)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
# library(ggrepel)
# library(reticulate)
# library(WGCNA)

# load("NA_SEURAT_P1P2CKO_P1.RData")
# # seuObj

# Idents(seuObj) <- "RNA_snn_res.0.8"

# metaTemp <- as.data.frame(seuObj@meta.data)
# metaTemp$Cluster <- metaTemp$RNA_snn_res.0.8
# metaTemp$CellBarcode <- row.names(metaTemp)

# annotations <- read.table("RES_0.8/SEURAT_NA_ANNOTATION_PER_CLUSTER_0.8.tsv", header = TRUE, sep = "\t")
# # annotations$clusterid <- annotations$Cluster

# annotations2 <- merge(metaTemp, annotations, by = "Cluster")

# annotations3 <- annotations2$CellType
# names(annotations3) <- annotations2$CellBarcode

# seuObj[["CellType"]] <- annotations3
# seuObj[["CellType_Cluster"]] <- paste(seuObj@meta.data$CellType, seuObj@meta.data$RNA_snn_res.0.8, sep = "_")

# p1 <- DimPlot(seuObj, reduction = "umap", pt.size = 0.1, group.by = "CellType", label = TRUE, label.size = 3, repel = TRUE) 
# ggsave(filename = "RES_0.8/NA_SEURAT_RES_0.8_UMAP_CELLTYPE.pdf", plot = p1, width = 7.5, height = 6, units = "in", dpi = 300)

# p1 <- DimPlot(seuObj, reduction = "umap", pt.size = 0.1, group.by = "CellType", label = FALSE, label.size = 3, repel = TRUE) 
# ggsave(filename = "RES_0.8/NA_SEURAT_RES_0.8_UMAP_CELLTYPE2.pdf", plot = p1, width = 7.5, height = 6, units = "in", dpi = 300)

# p2 <- DimPlot(seuObj, reduction = "umap", pt.size = 0.1, group.by = "CellType_Cluster", label = TRUE, label.size = 2, repel = TRUE) #+ facet_wrap(~CellType_Cluster, nrow = 4)
# ggsave(filename = "RES_0.8/NA_SEURAT_RES_0.8_UMAP_CELLTYPE_CLUSTER.pdf", plot = p2, width = 9, height = 6, units = "in", dpi = 150)

# save(seuObj, file = "NA_SEURAT_P1P2CKO_P1_UPDATED.RData")





# ## ------------------------------------
# ## ------------------------------------
# ## NEW WPRE
# ## Mark P1-CTL-01 and P1-CTL-02 cellbarcodes for WPRE 
# rm(list = ls())
# library(Seurat)
# library(dplyr)
# library(Matrix)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
# library(ggrepel)
# library(reticulate)
# library(WGCNA)


# ## Load Seurat clustering data
# load("NA_SEURAT_P1P2CKO_P1_UPDATED.RData")
# # seuObj

# seuMetaTemp <- as.data.frame(seuObj@meta.data)

# ## Load P1-CTL-01 and P1-CTL-02 from 10X run## load dataset from 10X run
# p1.p1p2cko.1 <- Read10X(data.dir = "/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE-5/WT03WPRE/outs/filtered_feature_bc_matrix")
# print(dim(p1.p1p2cko.1)) # 53716 23368

# p1.p1p2cko.2 <- Read10X(data.dir = "/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE-5/WT04WPRE/outs/filtered_feature_bc_matrix")
# print(dim(p1.p1p2cko.2)) # 53716 14099

# ## adding sample & condition info to column names
# colnames(p1.p1p2cko.1) <- paste("CTL_P1_01", colnames(p1.p1p2cko.1), sep = "_")
# colnames(p1.p1p2cko.2) <- paste("CTL_P1_02", colnames(p1.p1p2cko.2), sep = "_")

# p1.ctlw.1 <- as.data.frame(p1.p1p2cko.1["WPRE",])
# p1.ctlw.2 <- as.data.frame(p1.p1p2cko.2["WPRE",])
# colnames(p1.ctlw.1) <- "WPRE"
# colnames(p1.ctlw.2) <- "WPRE"


# # cellbarcodes.ctl1 <- scan("/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE/WT03WPRE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", what = "", sep = "\n") # 21387
# # cellbarcodes.ctl2 <- scan("/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE/WT04WPRE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", what = "", sep = "\n") # 11855

# # cb.ctl1 <- as.data.frame(cellbarcodes.ctl1)
# # colnames(cb.ctl1) <- "CellBarcode1"
# # cb.ctl1$Geno1 <- rep("CTL", nrow(cb.ctl1))
# # cb.ctl1$Age1 <- rep("P1", nrow(cb.ctl1))
# # cb.ctl1$Lib1 <- rep("01", nrow(cb.ctl1))
# # row.names(cb.ctl1) <- paste("CTL_P1_01", cb.ctl1$CellBarcode1, sep = "_" )

# # cb.ctl2 <- as.data.frame(cellbarcodes.ctl2)
# # colnames(cb.ctl2) <- "CellBarcode1"
# # cb.ctl2$Geno1 <- rep("CTL", nrow(cb.ctl2))
# # cb.ctl2$Age1 <- rep("P1", nrow(cb.ctl2))
# # cb.ctl2$Lib1 <- rep("02", nrow(cb.ctl2))
# # row.names(cb.ctl2) <- paste("CTL_P1_02", cb.ctl2$CellBarcode1, sep = "_" )


# ## Combine CTL1 and CTL2 tables
# ctl.data <- rbind(p1.ctlw.1, p1.ctlw.2)
# ctl.data$CELLBARCODE <- row.names(ctl.data)
# common.cb <- intersect(ctl.data$CELLBARCODE, row.names(seuMetaTemp))
# common.data <- as.data.frame(ctl.data[ctl.data$CELLBARCODE %in% common.cb,]) # as.data.frame(ctl.data[row.names(ctl.data) %in% common.cb,])
# common.data$WPRELIB <- rep(1, nrow(common.data))

# seuMetaUpdated <- merge(seuMetaTemp, common.data, by = "row.names", all.x = TRUE)
# row.names(seuMetaUpdated) <- seuMetaUpdated$Row.names
# seuMetaUpdated$Row.names <- NULL
# seuMetaUpdated$WPRELIB[is.na(seuMetaUpdated$WPRELIB)] <- 0
# seuMetaUpdated$CELLBARCODE[is.na(seuMetaUpdated$CELLBARCODE)] <- "MISSING"
# seuMetaUpdated$WPRE[is.na(seuMetaUpdated$WPRE)] <- 999

# seuwpre <- factor(seuMetaUpdated$WPRE)
# names(seuwpre) <- row.names(seuMetaUpdated)

# seuObj$WPRE <- seuwpre

# seuwpre <- factor(seuMetaUpdated$WPRELIB)
# names(seuwpre) <- row.names(seuMetaUpdated)

# seuObj$WPRELIB <- seuwpre


# save(seuObj, file = "NA_SEURAT_P1P2CKO_P1_UPDATED_WPRE_NEW.RData")




# myres <- 0.8

# pwpre <- DimPlot(seuObj, reduction = "umap", group.by = "WPRELIB", cols = c("grey90", "blue"), pt.size = 0.1, order = c(1, 0))
# ggsave(filename = "NEW_SEURAT_NA_CLUST_0.8_WPRELIB_1.pdf", plot = pwpre, width = 6.5, height = 6, units = "in", dpi = 150)

# pwpre <- DimPlot(seuObj, reduction = "umap", group.by = "WPRE", pt.size = 0.1)
# ggsave(filename = "NEW_SEURAT_NA_CLUST_0.8_WPRE_1.pdf", plot = pwpre, width = 6.5, height = 6, units = "in", dpi = 150)

# wprepercluster <- as.data.frame.matrix(table(seuObj@meta.data$CellType_Cluster, seuObj@meta.data$WPRE))
# write.table(wprepercluster, "NEW_SEURAT_NA_CLUST_0.8_WPRE.txt", row.names = T, col.names = T, quote = F, sep = "\t")



# metaDataNew <- as.data.frame(seuObj@meta.data)
# metaDataNew2 <- metaDataNew[!metaDataNew$WPRE %in% c(0, 999),]
# metaDataNew2$WPRE_RATIO <- 100 * (as.numeric(metaDataNew2$WPRE) / as.numeric(metaDataNew2$nCount_RNA))




# tags <- c("(0_0.05]", "(0.05_0.1]", "(0.1_0.2]", "(0.2_0.4]", "(0.4_0.6]", "(0.6_0.8]", "(0.8_1.0]", "(1.0_1.2]", "(1.2_1.4]", "(1.4_1.6]", "(1.6_1.8]", "(1.8_2.0]", "(2.0_2.2]", "(2.2_3.4]")

# v <- metaDataNew2 %>% select(CellType_Cluster, WPRE_RATIO) #pick the variable 
# vgroup <- as_tibble(v) %>% 
#   mutate(tag = case_when(
#     WPRE_RATIO <= 0.05 ~ tags[1],
#     WPRE_RATIO > 0.05 & WPRE_RATIO <= 0.1 ~ tags[2],
#     WPRE_RATIO > 0.1 & WPRE_RATIO <= 0.2 ~ tags[3],
#     WPRE_RATIO > 0.2 & WPRE_RATIO <= 0.4 ~ tags[4],
#     WPRE_RATIO > 0.4 & WPRE_RATIO <= 0.6 ~ tags[5],
#     WPRE_RATIO > 0.6 & WPRE_RATIO <= 0.8 ~ tags[6],
#     WPRE_RATIO > 0.8 & WPRE_RATIO <= 1.0 ~ tags[7],
#     WPRE_RATIO > 1.0 & WPRE_RATIO <= 1.2 ~ tags[8],
#     WPRE_RATIO > 1.2 & WPRE_RATIO <= 1.4 ~ tags[9],
#     WPRE_RATIO > 1.4 & WPRE_RATIO <= 1.6 ~ tags[10],
#     WPRE_RATIO > 1.6 & WPRE_RATIO <= 1.8 ~ tags[11],
#     WPRE_RATIO > 1.8 & WPRE_RATIO <= 2.0 ~ tags[12],
#     WPRE_RATIO > 2.0 & WPRE_RATIO <= 2.2 ~ tags[13],
#     WPRE_RATIO > 2.2 & WPRE_RATIO <= 3.4 ~ tags[14]
#     ))
# summary(vgroup)

# vgroup$tag <- factor(vgroup$tag,
#                        levels = tags,
#                        ordered = FALSE)
# summary(vgroup$tag)

# #   (0_0.05] (0.05_0.1]  (0.1_0.2]  (0.2_0.4]  (0.4_0.6]  (0.6_0.8]  (0.8_1.0] 
# #       3682       2716       1135        282         54         30         14 
# #  (1.0_1.2]  (1.2_1.4]  (1.4_1.6]  (1.6_1.8]  (1.8_2.0]  (2.0_2.2]  (2.2_3.4] 
# #         12          5          1          2          1          2          3

# vgroup2 <- as.data.frame(vgroup)

# vgroup3 <- vgroup2 %>%
#                 group_by(CellType_Cluster, tag) %>%
#                 summarise(WPRE2 = n())        

# vgroup4 <- dcast(vgroup3, CellType_Cluster ~ tag)
# write.table(vgroup4, "NEW_SEURAT_NA_CLUST_0.8_WPRE_RATIO.txt", row.names = F, col.names = T, quote = F, sep = "\t")





# # library(dplyr)
# # metaData <- as.data.frame(seuObj@meta.data)
# # metaData2 <- metaData[,c("CellType_Cluster", "GenoAgeSample", "WPRE")]
# # metaData21 <- metaData2[metaData2$WPRE != 0,]
# # metaData22 <- metaData21[metaData21$WPRE != 999,]
# # metaData3 <- metaData22 %>% 
# #         group_by(CellType_Cluster, GenoAgeSample) %>%
# #         summarise(WPRE2 = n())
# # metaData3 <- as.data.frame(metaData3)
# # metaData4 <- dcast(metaData3, CellType_Cluster ~ Genotype)

# # write.table(metaData4, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPECLUSTER_GENOTYPE.txt", row.names = F, col.names = T, quote = F, sep = "\t")
     

# # metaData32 <- metaData2 %>% 
# #         group_by(CellType_Cluster, Genotype) %>%
# #         summarise(NUCLEI = n())
# # metaData32 <- as.data.frame(metaData32)
# # metaData42 <- dcast(metaData32, CellType_Cluster ~ Genotype)
# # colnames(metaData42) <- c("CellType_Cluster", "CTL_TOTAL")

# # metaData5 <- merge(metaData4, metaData42, by = "CellType_Cluster")

# # write.table(metaData5, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPECLUSTER_GENOTYPE2.txt", row.names = F, col.names = T, quote = F, sep = "\t")





# # metaData <- as.data.frame(seuObj@meta.data)
# # metaData2 <- metaData[,c("CellType", "Genotype", "WPRELIB")]
# # metaData21 <- metaData2[metaData2$WPRELIB != 0,]
# # metaData3 <- metaData21 %>% 
# #         group_by(CellType, Genotype) %>%
# #         summarise(WPRE2 = n())
# # metaData3 <- as.data.frame(metaData3)
# # metaData4 <- dcast(metaData3, CellType ~ Genotype)

# # write.table(metaData4, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE.txt", row.names = F, col.names = T, quote = F, sep = "\t")
     

# # metaData32 <- metaData2 %>% 
# #         group_by(CellType, Genotype) %>%
# #         summarise(NUCLEI = n())
# # metaData32 <- as.data.frame(metaData32)
# # metaData42 <- dcast(metaData32, CellType ~ Genotype)
# # colnames(metaData42) <- c("CellType", "CTL_TOTAL")

# # metaData5 <- merge(metaData4, metaData42, by = "CellType")

# # write.table(metaData5, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE2.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# # metaData5$PCT_CTL <- 100*(metaData5$CTL/metaData5$CTL_TOTAL)

# # p2 <- ggplot(metaData5, aes(x = CellType, y = PCT_CTL, fill = CellType)) + 
# #       geom_bar(position = "fill") +
# #       scale_y_continuous(labels = scales::percent) +
# #       NULL
# # ggsave(filename = "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE2.pdf", plot = p2, width = 8, height = 4, units = "in", dpi = 150)            



# # ## ------------------------------------
# # ## ------------------------------------
# # ## OLD WPRE
# # ## Mark P1-CTL-01 and P1-CTL-02 cellbarcodes for WPRE 
# # rm(list = ls())
# # library(Seurat)
# # library(dplyr)
# # library(Matrix)
# # library(ggplot2)
# # library(gridExtra)
# # library(reshape2)
# # library(ggrepel)
# # library(reticulate)
# # library(WGCNA)


# # ## Load Seurat clustering data
# # load("NA_SEURAT_P1P2CKO_P1_UPDATED.RData")
# # # seuObj

# # seuMetaTemp <- as.data.frame(seuObj@meta.data)

# # ## Load P1-CTL-01 and P1-CTL-02 from 10X run
# # cellbarcodes.ctl1 <- scan("/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE/WT03WPRE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", what = "", sep = "\n") # 21387
# # cellbarcodes.ctl2 <- scan("/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/NOVASEQ_JUN2021/03_COUNT_WT-WPRE/WT04WPRE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", what = "", sep = "\n") # 11855

# # cb.ctl1 <- as.data.frame(cellbarcodes.ctl1)
# # colnames(cb.ctl1) <- "CellBarcode1"
# # cb.ctl1$Geno1 <- rep("CTL", nrow(cb.ctl1))
# # cb.ctl1$Age1 <- rep("P1", nrow(cb.ctl1))
# # cb.ctl1$Lib1 <- rep("01", nrow(cb.ctl1))
# # row.names(cb.ctl1) <- paste("CTL_P1_01", cb.ctl1$CellBarcode1, sep = "_" )

# # cb.ctl2 <- as.data.frame(cellbarcodes.ctl2)
# # colnames(cb.ctl2) <- "CellBarcode1"
# # cb.ctl2$Geno1 <- rep("CTL", nrow(cb.ctl2))
# # cb.ctl2$Age1 <- rep("P1", nrow(cb.ctl2))
# # cb.ctl2$Lib1 <- rep("02", nrow(cb.ctl2))
# # row.names(cb.ctl2) <- paste("CTL_P1_02", cb.ctl2$CellBarcode1, sep = "_" )


# # ## Combine CTL2 and cKO tables
# # ctl.data <- rbind(cb.ctl1, cb.ctl2)
# # common.cb <- intersect(row.names(ctl.data), row.names(seuMetaTemp))
# # common.data <- ctl.data[row.names(ctl.data) %in% common.cb,]
# # common.data$WPRELIB <- rep(1, nrow(common.data))

# # seuMetaUpdated <- merge(seuMetaTemp, common.data, by = "row.names", all.x = TRUE)
# # row.names(seuMetaUpdated) <- seuMetaUpdated$Row.names
# # seuMetaUpdated$Row.names <- NULL
# # seuMetaUpdated$WPRELIB[is.na(seuMetaUpdated$WPRELIB)] <- 0


# # seuwpre <- factor(seuMetaUpdated$WPRELIB)
# # names(seuwpre) <- row.names(seuMetaUpdated)

# # seuObj$WPRELIB <- seuwpre

# # save(seuObj, file = "NA_SEURAT_P1P2CKO_P1_UPDATED_WPRE.RData")




# # myres <- 0.8

# # pwpre <- DimPlot(seuObj, reduction = "umap", group.by = "WPRELIB", cols = c("grey90", "blue"), pt.size = 0.1, order = c(1, 0))
# # ggsave(filename = "SEURAT_NA_CLUST_0.8_WPRELIB_1.pdf", plot = pwpre, width = 6.5, height = 6, units = "in", dpi = 150)

# # wprepercluster <- as.data.frame.matrix(table(seuObj@meta.data$CellType_Cluster, seuObj@meta.data$WPRELIB))
# # write.table(wprepercluster, "SEURAT_NA_CLUST_0.8_WPRELIB.txt", row.names = T, col.names = T, quote = F, sep = "\t")




# # library(dplyr)
# # metaData <- as.data.frame(seuObj@meta.data)
# # metaData2 <- metaData[,c("CellType_Cluster", "Genotype", "WPRELIB")]
# # metaData21 <- metaData2[metaData2$WPRELIB != 0,]
# # metaData3 <- metaData21 %>% 
# #         group_by(CellType_Cluster, Genotype) %>%
# #         summarise(WPRE2 = n())
# # metaData3 <- as.data.frame(metaData3)
# # metaData4 <- dcast(metaData3, CellType_Cluster ~ Genotype)

# # write.table(metaData4, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPECLUSTER_GENOTYPE.txt", row.names = F, col.names = T, quote = F, sep = "\t")
     

# # metaData32 <- metaData2 %>% 
# #         group_by(CellType_Cluster, Genotype) %>%
# #         summarise(NUCLEI = n())
# # metaData32 <- as.data.frame(metaData32)
# # metaData42 <- dcast(metaData32, CellType_Cluster ~ Genotype)
# # colnames(metaData42) <- c("CellType_Cluster", "CTL_TOTAL")

# # metaData5 <- merge(metaData4, metaData42, by = "CellType_Cluster")

# # write.table(metaData5, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPECLUSTER_GENOTYPE2.txt", row.names = F, col.names = T, quote = F, sep = "\t")





# # metaData <- as.data.frame(seuObj@meta.data)
# # metaData2 <- metaData[,c("CellType", "Genotype", "WPRELIB")]
# # metaData21 <- metaData2[metaData2$WPRELIB != 0,]
# # metaData3 <- metaData21 %>% 
# #         group_by(CellType, Genotype) %>%
# #         summarise(WPRE2 = n())
# # metaData3 <- as.data.frame(metaData3)
# # metaData4 <- dcast(metaData3, CellType ~ Genotype)

# # write.table(metaData4, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE.txt", row.names = F, col.names = T, quote = F, sep = "\t")
     

# # metaData32 <- metaData2 %>% 
# #         group_by(CellType, Genotype) %>%
# #         summarise(NUCLEI = n())
# # metaData32 <- as.data.frame(metaData32)
# # metaData42 <- dcast(metaData32, CellType ~ Genotype)
# # colnames(metaData42) <- c("CellType", "CTL_TOTAL")

# # metaData5 <- merge(metaData4, metaData42, by = "CellType")

# # write.table(metaData5, "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE2.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# # # metaData5$PCT_CTL <- 100*(metaData5$CTL/metaData5$CTL_TOTAL)

# # # p2 <- ggplot(metaData5, aes(x = CellType, y = PCT_CTL, fill = CellType)) + 
# # #       geom_bar(position = "fill") +
# # #       scale_y_continuous(labels = scales::percent) +
# # #       NULL
# # # ggsave(filename = "SEURAT_NA_CLUST_0.8_WPRELIB_CELLTYPE_GENOTYPE2.pdf", plot = p2, width = 8, height = 4, units = "in", dpi = 150)            



# ------------------------------------
# End
# ------------------------------------

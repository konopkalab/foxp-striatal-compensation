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
library(rtracklayer)
# library(EnsDb.Mmusculus.v79)
# library(BSgenome.Mmusculus.UCSC.GRCm38p6)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


##-------------------------------------------------------
## REFERENCE GENOME
# ## extract gene annotations from EnsDb
# # annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) #, standard.chromosomes = TRUE, biotypes = "protein_coding")
annotations <- GRangesForUCSCGenome(genome = "mm10")

# ## change to UCSC style since the data was mapped to mm10
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "mm10"


##-------------------------------------------------------
## REFERENCE ANNOTATION
gtf <- rtracklayer::import("/work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_CELLRANGER_ATACSEQ/ForSignac/genes.gtf")
# fa <- Rsamtools::FaFile("/work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_CELLRANGER_ATACSEQ/ForSignac/genome.fa.gz")
mm10_annotation <- gtf[gtf$gene_type == "protein_coding"]
genome(mm10_annotation) <- "mm10"
seqlevelsStyle(mm10_annotation) <- "UCSC"
mm10_annotation$gene_biotype <- mm10_annotation$gene_type
mm10_annotation <- keepStandardChromosomes(mm10_annotation, pruning.mode = "coarse")



##-------------------------------------------------------
## Custom function to convert MACS2 output table into a GRanges object
macs2GRanges <-function(peaks) {
  # generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
    strand="*",
    count=peaks$pileup,
    score=peaks$X.log10.pvalue.,
    FE=peaks$fold_enrichment,
    fdr=peaks$X.log10.qvalue.,
    maxpos=peaks$abs_summit
    )
  return(myrange)
  }


##-------------------------------------------------------
## Custom function to Read MACS2 output and make a GR
grMACS2 <- function(macs2PeaksPath)
 {
 MACS2_OUT <- read.table(macs2PeaksPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
 MACS2_PEAKS <- macs2GRanges(MACS2_OUT)
 return(MACS2_PEAKS)
 }


##-------------------------------------------------------
## Identify common peaks across samples of same genotype
na107macs2 <- grMACS2("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA107_P1P2CKO_MACS2/NA107_P1P2CKO_peaks.xls")
na111macs2 <- grMACS2("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA111_P1P2CKO_MACS2/NA111_P1P2CKO_peaks.xls")
na115macs2 <- grMACS2("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA115_P1P2CKO_MACS2/NA115_P1P2CKO_peaks.xls")

common.peaks.p1p2cko <- reduce(x = c(na107macs2, na111macs2, na115macs2))


##-------------------------------------------------------
## Function for Signac using MACS2 peaks
processSignac <- function(crCountsPath, crMetaPath, crFragPath, macs2PeaksPath, archrMetaPath, prefx, covSex, covBatch)
    {
    # crCountsPath <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA104_CTL/filtered_peak_bc_matrix.h5"
    # crMetaPath <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA104_CTL/singlecell.csv"
    # crFragPath <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA104_CTL/fragments.tsv.gz"
    # macs2PeaksPath <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA104_CTL_MACS2/NA104_CTL_peaks.xls"
    # archrMetaPath <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/02_SEURAT_SIGNAC_DOUBLETS_ARCHR1/NA_ATAC_ARCHR_DOUBLETS_META.RData"
    # prefx <- "NA104_CTL"
    # covSex <- "M"
    # covBatch <- 1

    ## Read CellRanger output
    CR_COUNTS <- Read10X_h5(crCountsPath)
    CR_META <- read.csv(crMetaPath, header = TRUE, row.names = 1)

    ## Read MACS2 output
    MACS2_OUT <- read.table(macs2PeaksPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    MACS2_PEAKS <- macs2GRanges(MACS2_OUT)

    ## ArchR Doublets Meta
    load(archrMetaPath)
    # naarchrmeta, nafiltarchrmeta
    AR_META <- naarchrmeta[naarchrmeta$Sample %in% prefx,]
    row.names(AR_META) <- gsub(paste(prefx, "#", sep = ""), "", row.names(AR_META))
    colnames(AR_META) <- paste("AR", colnames(AR_META), sep = "_")

    ## Create chromatin assay
    CR_CHRASSAY <- CreateChromatinAssay(counts = CR_COUNTS, sep=c(":","-"), genome = "mm10", fragments = crFragPath, annotation = mm10_annotation)

    ## Create Seurat object
    CR_SEURAT <- CreateSeuratObject(CR_CHRASSAY, assay = "peaks", project = prefx, meta.data = CR_META)

    ## Create counts table using MACS2 output
    # MACS2_COUNTS <- FeatureMatrix(fragments = Fragments(CR_SEURAT), features = MACS2_PEAKS, cells = colnames(CR_SEURAT))
    MACS2_COUNTS <- FeatureMatrix(fragments = Fragments(CR_SEURAT), features = common.peaks.p1p2cko, cells = colnames(CR_SEURAT))

    ## Create chromatin assay using MACS2 counts
    MACS2_CHRASSAY <- CreateChromatinAssay(counts = MACS2_COUNTS, genome = "mm10", fragments = Fragments(CR_SEURAT), annotation = mm10_annotation)

    ## Create Seurat object using MACS2-based chromatin assay
    MACS2_SEURAT <- CreateSeuratObject(MACS2_CHRASSAY, assay = "peaks", meta.data = CR_META)

    ## Update annotation
    Annotation(MACS2_SEURAT) <- mm10_annotation
    
    ## Update meta
    MACS2_SEURAT$Dataset <- prefx
    MACS2_SEURAT$Sex <- covSex
    MACS2_SEURAT$Batch <- covBatch

    ## Nucleosome signal
    MACS2_SEURAT <- NucleosomeSignal(object = MACS2_SEURAT)

    ## TSS Enrichment
    MACS2_SEURAT <- TSSEnrichment(MACS2_SEURAT, fast = FALSE)

    ## Plot TTS enrichment
    pTSSE <- TSSPlot(MACS2_SEURAT) + NoLegend()
    ggsave(paste(prefx, "TSS_Enrichment_1.pdf", sep = "_"), plot = pTSSE, width = 10, height = 5, units = "in", dpi = 300)

    ## Add Meta
    MACS2_SEURAT$pct_reads_in_peaks <- MACS2_SEURAT$peak_region_fragments / MACS2_SEURAT$passed_filters * 100
    MACS2_SEURAT$blacklist_ratio <- MACS2_SEURAT$blacklist_region_fragments / MACS2_SEURAT$peak_region_fragments

    ## Plot QC
    pQC1 <- VlnPlot(object = MACS2_SEURAT, features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1, ncol = 5)
    ggsave(paste(prefx, "QC_PLOT_1.pdf", sep = "_"), plot = pQC1, width = 20, height = 5, units = "in", dpi = 300)

    ## Compute the term-frequency inverse-document-frequency normalization on a matrix used in LSA (Latent Semantic Analysis)
    MACS2_SEURAT <- RunTFIDF(MACS2_SEURAT)

    ## Find most frequently observed features.
    MACS2_SEURAT <- FindTopFeatures(MACS2_SEURAT, min.cutoff = 'q0')

    ## Run singular value decomposition, finds a few approximate largest (or, optionally, smallest) singular values and corresponding singular vectors of a sparse or dense matrix.
    MACS2_SEURAT <- RunSVD(object = MACS2_SEURAT)

    ## Plot depth cor
    pDC <- DepthCor(MACS2_SEURAT)
    ggsave(paste(prefx, "Depth_Cor.pdf", sep = "_"), plot = pDC, width = 6, height = 4, units = "in", dpi = 300)

    ## Run UMAP, Find Neighbours, Find Clusters
    MACS2_SEURAT <- RunUMAP(MACS2_SEURAT, dims = 2:50, reduction = 'lsi') %>% FindNeighbors(dims = 2:50, reduction = 'lsi') %>% FindClusters(algorithm = 3, resolution = 0.8, dims = 2:50, reduction = 'lsi')

    ## Plot UMAP
    pUMAP <- DimPlot(MACS2_SEURAT, group.by = 'peaks_snn_res.0.8', pt.size = 0.1, label = TRUE) + NoLegend()
    ggsave(paste(prefx, "UMAP_CLUSTERS.pdf", sep = "_"), plot = pUMAP, width = 6, height = 6, units = "in", dpi = 300)

    ## Gene activity matrix
    MACS2_GA <- GeneActivity(MACS2_SEURAT)

    ## Add the gene activity matrix to the Seurat object as a new assay
    MACS2_SEURAT[['RNA']] <- CreateAssayObject(counts = MACS2_GA)
    MACS2_SEURAT <- NormalizeData(object = MACS2_SEURAT, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(MACS2_SEURAT$nCount_RNA))

    ## Feature plot using RNA as default assay
    DefaultAssay(MACS2_SEURAT) <- 'RNA'
    pFeature <- FeaturePlot(object = MACS2_SEURAT, features = c("Foxp1", "Foxp2", "Drd", "Tac1", "Drd2", "Penk", "Casz1", "Aqp4", "Mag", "Cx3cr1", "Flt1", "Ppp1r1b"), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
    ggsave(paste(prefx, "UMAP_GENES.pdf", sep = "_"), plot = pFeature, width = 15, height = 16, units = "in", dpi = 300)

    save(MACS2_PEAKS, CR_CHRASSAY, CR_SEURAT, MACS2_COUNTS, MACS2_CHRASSAY, MACS2_SEURAT, AR_META, common.peaks.p1p2cko, file = paste(prefx, "MACS2_SIGNAC.RData", sep = "_"))
    }


print("=====> PROCESSING NA107_P1P2CKO")
# crCountsPath, crMetaPath, crFragPath, macs2PeaksPath, archrMetaPath, prefx, covSex, covBatch
processSignac("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA107_P1P2CKO/filtered_peak_bc_matrix.h5",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA107_P1P2CKO/singlecell.csv",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_COMBINED_A_B/NA107_P1P2CKO/fragments.tsv.gz",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA107_P1P2CKO_MACS2/NA107_P1P2CKO_peaks.xls",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/02_SEURAT_SIGNAC_DOUBLETS_ARCHR1/NA_ATAC_ARCHR_DOUBLETS_META.RData",
              "NA107_P1P2CKO",
              "M",
              1)


print("=====> PROCESSING NA111_P1P2CKO")
processSignac("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA111_P1P2CKO/filtered_peak_bc_matrix.h5",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA111_P1P2CKO/singlecell.csv",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_C/NA111_P1P2CKO/fragments.tsv.gz",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA111_P1P2CKO_MACS2/NA111_P1P2CKO_peaks.xls",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/02_SEURAT_SIGNAC_DOUBLETS_ARCHR1/NA_ATAC_ARCHR_DOUBLETS_META.RData",
              "NA111_P1P2CKO",
              "M",
              2)

print("=====> PROCESSING NA115_P1P2CKO")
processSignac("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA115_P1P2CKO_M/filtered_peak_bc_matrix.h5",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA115_P1P2CKO_M/singlecell.csv",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/00_CELLRANGER_D/NA115_P1P2CKO_M/fragments.tsv.gz",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/04_MACS2/NA115_P1P2CKO_MACS2/NA115_P1P2CKO_peaks.xls",
              "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NA_DATA/ATACSEQ/02_SEURAT_SIGNAC_DOUBLETS_ARCHR1/NA_ATAC_ARCHR_DOUBLETS_META.RData",
              "NA115_P1P2CKO",
              "M",
              3)


##-------------------------------------------------------
## END
##-------------------------------------------------------

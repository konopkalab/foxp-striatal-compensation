##------------------------------------
## NK | WT | MAR 2021
##------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl



##------------------------------------
## Libraries
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)
library(WGCNA)
options(future.globals.maxSize= 50000 * 1024^2)



cluMarkers <- read.table("SEURAT_NA_P09_INTERGRATION_HARMONY_CLU_MARKERS.txt", sep = "\t", header = TRUE)
tab <- cluMarkers
tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.5 & tab$pct.1 >= 0.25,]
# tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.5,]
# tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.25,]
tab <- tab[c(7,6)]
tab$cluster <- as.factor(paste("Cluster_", sprintf("%02d", as.numeric(as.character(tab$cluster))), sep=""))
colnames(tab)=c("Gene","DEFINITION")
Genes=as.data.frame(table(tab$DEFINITION))


load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/05_SEURAT/MISC/Saunders2018DEG.RData")
# degSaunders
GeneSets <- degSaunders


for(i in 1:length(GeneSets)){
	colnames(GeneSets[[i]])[1] <- "Gene"
}

ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- length(unique(sort(cluMarkers$gene)))-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585 #13517 #nrow(tab) #8321 # length(unique(sort(cluMarkers$gene))) 4540   # nrow(cluMarkers)
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
	f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
	return(c(row,
			P_val = f$p.value,
			LogP = -log10(f$p.value), 
			OR = f$estimate[[1]],
			OR_Low = f$conf.int[1],
			OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat, TEMP, file= "NA_ANNOTATION_FisherOutput_Saunders_Enrich.RData")


# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0


pdf("NA_ANNOTATION_FisherOutput_Saunders_Enrich.pdf", width=12, height=16, pointsize=12)
par(mar=c(15, 7, 2, 2))
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5,
xLabelsAngle = 90)
dev.off()



FisherORt <- as.data.frame(t(FisherOR))
colnames(FisherORt) <- paste(colnames(FisherORt), "OR", sep = "_")

FisherAdjt <- as.data.frame(t(FisherAdj))
colnames(FisherAdjt) <- paste(colnames(FisherAdjt), "Pval", sep = "_")

FisherData <- merge(FisherORt, FisherAdjt, by = "row.names")
row.names(FisherData) <- FisherData$Row.names
FisherData$Row.names <- NULL

FisherData2 <- FisherData[,order(colnames(FisherData))]
write.table(FisherData2, "NA_ANNOTATION_FisherOutput_Saunders_Enrich_PlotData.txt", row.names = T, col.names = T, quote = F, sep = "\t")


sink("NA_ANNOTATION_FisherOutput_Saunders_Enrich_AnnotationData.txt")
aa <- FisherData2
for(i in 1:ncol(aa))
 {	 
 if(i %% 2 != 0)
  {
  cluor <- i
  clup <- i + 1
  print("------------------")
  bb <- aa[,c(cluor, clup)]
  cc <- bb[order(bb[,1], decreasing = T),]
  dd <- bb[order(bb[,2], decreasing = F),]
  print(gsub("_OR", "", colnames(cc)[1]))
  print(cc[1:2,])
  print(dd[1:2,])
  print("------------------")
  }
 }

sink()


##------------------------------------
## End
##------------------------------------
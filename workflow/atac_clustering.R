if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

library(reshape2)
library(dplyr)
library(Seurat)
library(Signac)
library(patchwork)
library(ggplot2)
library(tidyr)
library(stringr)
#library(janitor)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDbData)
library(data.table)

setwd("results/")
pdf(file=paste0("plot/atac_clustering.pdf"))
gex<-readRDS("gex/features/gex_seurat.rds")

###ATAC clustering###

#Find gex markers
gex.markers <- FindAllMarkers(gex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

for(i in levels(gex@active.ident)){
  assign(paste("cluster",i,".markers", sep=""),(FindMarkers(gex, ident.1 = i, min.pct = 0.25)))
}

#read in scATAC data

counts<-readRDS("atac/filtered_matrix/MACS2/FILTER/matrix.rds")

chrom_assay <- CreateChromatinAssay(
  counts = counts
)

atac <- CreateSeuratObject(
  counts = chrom_assay,assay="peaks", min.cells = "atac/filtered_matrix/MACS2/FILTER/barcodes.txt"
)

Fragments(atac) <- CreateFragmentObject(
  path = "atac/summary/atac.fragments.tsv.gz",
  cells = colnames(atac)
)

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(atac) <- annotations

##Reductions
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


p1 <- DimPlot(gex, label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2

# quantify gene activity
gene.activities <- GeneActivity(atac, features = VariableFeatures(gex), gene.id = T)

# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = gex, query = atac, features = VariableFeatures(object = gex),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = gex$seurat_clusters,
                                     weight.reduction = atac[["lsi"]], dims = 2:30)

atac <- AddMetaData(atac, metadata = celltype.predictions)

p1 <- DimPlot(atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p1 

cluster_pred<-as.data.frame(atac$predicted.id)
cluster_pred$Barcode<- rownames(cluster_pred)
colnames(cluster_pred)<-c("Cluster","Barcode")
write.table(cluster_pred[,2:1], "atac/features/fromgex_cluster_atac.tsv", row.names = FALSE, quote=FALSE, col.names = TRUE, sep="\t")

dev.off()
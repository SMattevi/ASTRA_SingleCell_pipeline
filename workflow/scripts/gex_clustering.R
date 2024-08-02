# install.packages("remotes")
# remotes::install_github("trevorld/r-optparse")
# library("optparse")
library("reshape2")
library("dplyr")
library("Seurat")
library("Signac")
library("patchwork")
library("ggplot2")
library("tidyr")
library("stringr")
library("biomaRt")
library("data.table")

# option_list = list(
#   make_option(c("-s", "--sample"), type="character", default=NULL,
#               help="sample id", metavar="character")
# );

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# print(opt)

# setwd(paste0("results_",opt$sample,"/"))

counts_gex<-read.table("gex/features/counts.tsv.gz",header=T)

countst<-reshape2::dcast(counts_gex,gene~cell, fill="0")
countst <- data.frame(countst, row.names = 1)
gex<-CreateSeuratObject(counts =countst)

gex <- NormalizeData(gex, normalization.method = "LogNormalize", scale.factor = 10000)
gex <- FindVariableFeatures(gex, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(gex), 10)

all.genes <- rownames(gex)
gex <- ScaleData(gex, features = all.genes)
gex <- RunPCA(gex, features = VariableFeatures(object = gex) )

VizDimLoadings(gex, dims = 1:2, reduction = "pca")
DimPlot(gex, reduction = "pca")

gex <- JackStraw(gex, num.replicate = 100)
gex <- ScoreJackStraw(gex, dims = 1:20)

JackStrawPlot(gex, dims = 1:15)
ElbowPlot(gex)

gex <- FindNeighbors(gex, dims = 1:10)
gex <- FindClusters(gex, resolution = 0.5) #parameter

gex <- RunUMAP(gex, dims = 1:10)
gex <- RunTSNE(gex, dims = 1:10,check_duplicates = FALSE)
gex<-RunSVD(gex)

cl<-as.data.frame(gex$seurat_clusters)
cl$gene<-row.names(cl)

colnames(cl)<-c("Cluster","Barcode") 

write.table(cl[,c("Barcode","Cluster")],"gex/features/cluster_gex.tsv",sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

saveRDS(gex,"gex/features/gex_seurat.rds")

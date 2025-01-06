# https://prabhakarlab.github.io/Banksy/articles/multi-sample.html
# Use banksy for multi-sample clustering with Harmony integration 
suppressPackageStartupMessages({
library(spatialLIBD)
library(ExperimentHub)
library(harmony)
library(Banksy)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
})
setwd("/hamit/")

# Load individual seura objec and store as a list
sampleName <- c("sample1", "sample2", "sample3")
seu_list <-list()
seu_list <- lapply(sampleName, function(x) {
  readRDS(paste0("./spatial1/02_QC/spatialObjFiltered_", x, ".rds"))
})

# Perform basic normalization using seura function 
seu_list <- lapply(seu_list, function(x) {
  NormalizeData(x, scale.factor = 3000, normalization.method = "SCT")
})

# Compute HVGs for each dataset and take the union
hvgs <- lapply(seu_list, function(x) {
  VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)

# Subset to HVGs
seu_list <- lapply(seu_list, function(x) x[hvgs,])
seu <- Reduce(merge, seu_list)
locs <- do.call(rbind.data.frame, lapply(spe_list, spatialCoords))
seu@meta.data <- cbind(seu@meta.data, locs)

# Grouping variable
head(seu@meta.data)
table(seu$sample_id)

# Set spatial coordinates
sdimx <- 'pxl_col_in_fullres'
sdimy <- 'pxl_row_in_fullres'

#----------------------- Banksy ----------------------- #
seu <- RunBanksy(seu, 
    lambda = 0.2, 
    assay = 'originalexp', 
    slot = 'data',
    dimx = sdimx, 
    dimy = sdimy, 
    features = 'all',
    group = 'sample_id', 
    split.scale = FALSE, 
    k_geom = 6) #
seu <- RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 15)
seu <- RunHarmony(seu, group.by.vars='sample_id')

#----------------------- Harmony ----------------------- #
seu <- RunUMAP(seu, dims = 1:15, reduction = 'harmony')
seu <- FindNeighbors(seu, dims = 1:15, reduction = 'harmony')
seu <- FindClusters(seu, resolution = 0.4)
save(seu, file = "./spatial1/03_Clustering/SeuratClustered.rds")

pdf("./spatial1/03_Clustering/Clustering.pdf")
DimPlot(seu, pt.size = 0.25, label = TRUE, label.size = 3, cols = 3)
FeatureScatter(seu, 'staggered_sdimx', 'staggered_sdimy', cols = 3, pt.size = 0.75)
dev.off()


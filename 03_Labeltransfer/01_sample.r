# Apply log normalization method on both snRNAseq and visium data to project 11 clusters from snRNAseq to visium slide
# snRNAseq seurat object is merged data from 14 samples and visium is per sample
suppressPackageStartupMessages({
library(Seurat)
library(SeuratData)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggplot2)
})
setwd("/hamit/")

# Import sample name
args <- commandArgs(TRUE)
sampleName <- args[1]
print(paste0(sampleName, " is being processed"))

# Load single cell seurat data
load(file = "./scrna1/04_Clustering/ClusteredBroad.rds") # seuratObj

# Load spatial data from specific sample
load(file = paste0("./spatial1/02_QC/spatialObjFiltered_", sampleName, ".rds"))

# Set pdf
namePdf <- paste0("./spatial1/03_LabelTransfer/Labeled_", as.character(sampleName), ".pdf")
pdf(file = namePdf)

#---------------------- Normalization ----------------------#
options(future.globals.maxSize = 3e+010)
print("Normalizatio spatial data")
spatialObj <- NormalizeData(spatialObj, normalization.method = "LogNormalize", scale.factor = 10000)
spatialObj <- FindVariableFeatures(spatialObj, selection.method = "vst", nfeatures = 2000)
spatialObj <- ScaleData(spatialObj, vars.to.regress = c("percent.mt"))
spatialObj <- RunPCA(spatialObj)

print("Normalization scRNA data")
seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
seuratObj <- ScaleData(seuratObj, vars.to.regress = c("mitoPercent"))
seuratObj <- RunPCA(seuratObj)
seuratObj <- RunUMAP(seuratObj, dims = 1:15)

# Check clustering, should be same as the snRNAseq
DimPlot(seuratObj, group.by = "idents", label = TRUE)

#---------------------- Label transfer ----------------------#
# Note: selection of reduction to use could depend on type of data transfering https://github.com/satijalab/seurat/issues/8489
# In this case, using pca as a reduction result in cluster with zero count in visium data, cca is used instead
anchors <- FindTransferAnchors(reference = seuratObj, 
    query = spatialObj, 
    normalization.method = "LogNormalize", 
    reduction = "cca")

# Transfer the annotation label
predictions.assay <- TransferData(anchorset = anchors, 
    refdata = seuratObj$idents, 
    prediction.assay = TRUE,
    weight.reduction = "cca", 
    dims = 1:15)

spatialObj[["predictions"]] <- predictions.assay
DefaultAssay(spatialObj) <- "predictions" 

# Save normalized object
save(spatialObj, file = paste0("./spatial1/03_labeltransfer/Labeled_", sampleName, ".rds"))

# ---------------------- Plot ---------------------- #
# Plot broad cluster expression on slide
plots <- SpatialFeaturePlot(spatialObj, 
    features = c(paste0("Exc", 1:4), "Inh", "Oligo", "OPC", paste0("Astro", 1:2), "Micro", "Endo"),
    pt.size.factor = 2, image.alpha = 0.3, ncol = 3, crop = TRUE)
plots <- lapply(plots, function(p) {
    p + theme(
        legend.key.size = unit(0.4, "cm"), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        axis.text = element_text(size = 6), 
        axis.title = element_text(size = 6))
})
plot_grid(plotlist = plots, ncol = 4, nrow = 3, scale = 1.2)

# Predict cell types on spatially restricted area and plot top 4 clusters
# https://github.com/satijalab/seurat-object/issues/25
spatialObj <- FindSpatiallyVariableFeatures(spatialObj,
    assay = "predictions", 
    selection.method = "moransi",
    features = rownames(spatialObj), 
    r.metric = 5, 
    slot = "data")
# Check if "moransi.spatially.variable.rank" is in column
colnames(spatialObj[["predictions"]]@meta.features)
# Check which cluser is the most spatially restricted
top.clusters <- rownames(dplyr::slice_min(spatialObj[["predictions"]]@meta.features, 
    order_by = spatialObj[["predictions"]]@meta.features$moransi.spatially.variable.rank, 
    n = 4))
SpatialPlot(object = spatialObj, features = top.clusters, ncol = 2, 
    pt.size.factor = 2, image.alpha = 0.3, ncol = 3, crop = TRUE)

dev.off()
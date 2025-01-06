# Loop through each sample and build spatial seurat object
# Basic gene filtering is done 
# output: spatial seurat object stored as .rds for each subject
suppressPackageStartupMessages({
library(Seurat)
library(tidyverse)
})

setwd("/hamit/")

# List sample names
sampleName <- c("sample1", "sample2", "sample3")

# Create seurat object
for (x in sampleName) {
    samplePath <- file.path("./data/rawadata", x, "outs")
    spatialObj <- Load10X_Spatial(
        samplePath,
        filename = "filtered_feature_bc_matrix.h5",
        assay = "Spatial",
        slice = "slice1")

    # --------------- Filter low quality spots ---------------- #
    # Check total spots and genes
    dim(spatialObj)

    #01 Filter spots with no gene expression
    noExpr <- which(spatialObj$nFeature_Spatial == 0)
    length(noExpr) 
    spatialObj <- spatialObj[-noExpr, ]

    #02 Filter spots with zero counts
    noCount <- which(spatialObj$nCount_Spatial == 0)
    length(noCount)
    spatialObj <- spatialObj[-noCount, ]

    #03 Drop spots outside the tissue
    if ("in_tissue" %in% colnames(spatialObj@meta.data)) {
        spatialObj <- subset(spatialObj, subset = in_tissue == 1)}

    # Check total spots and genes
    dim(spatialObj)

    # Save the object 
    save(spatialObj, file = paste0("./spatial1/01_buildSeurat/spatialObj", x, ".rds"))
}

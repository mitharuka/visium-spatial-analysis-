# Loop through each sample and build spatial seurat object
# Basic gene filtering is done 
# output: spatial seurat object stored as .rds for each subject
suppressPackageStartupMessages({
library(Seurat)
library(tidyverse)
library(scran)
})

setwd("/hamit/")

# List sample names
sampleName <- c("sample1", "sample2", "sample3")

# Create seurat object
for (x in sampleName) {
    load(paste0("./spatial1/01_buildSeurat/spatialObj", x, ".rds"))

    # Calculate QC metrics and add to metadata
    spatialObj[["percent.mt"]] <- PercentageFeatureSet(spatialObj, pattern = "^MT-") 
    spatialObj$detected <- spatialObj$nFeature_Spatial  # Number of detected features
    spatialObj$sum <- spatialObj$nCount_Spatial  # Total counts/reads
    spatialObj$percent.mt <- spatialObj[["percent.mt"]]  # Percentage of mitochondrial genes

    # Convert to  spe object to use isOUtlier function
    spe <- as.SingleCellExperiment(spatialObj)

    # Check mitochondria percentage in each spot
    spe$high_mito_id <- isOutlier(spe$percent.mt, nmads = 2, type = "higher")
    table(spe@colData$high_mito_id)

    # Check number of counter detecter per spot
    spe$low_sum_id <- isOutlier(spe$sum, log = TRUE, nmads = 3, type = "lower")
    table(spe$low_sum_id)

    # Check number of genes detected per spot
    spe$low_detected_id <- isOutlier(spe$detected, log = TRUE, nmads = 3, type = "lower")
    table(spe$low_detected_id)

    # Get the outlier spots from the SPE object and add to meta data
    spe$discard_auto_id <- spe$high_mito_id | spe$low_sum_id | spe$low_detected_id
    table(spe$discard_auto_id)

    print("Percentage of spots to be drop")
    print(100 * sum(spe$discard_auto_id) / ncol(spe))

    # Get metadata from spe object and combine o seurat
    spatialObj$feature_outlier <- spe$discard_auto_id
    save(spatialObj, file = paste0("./spatial1/02_QC/spatialObjNonFiltered_", x, ".rds"))

    # Remove outlier 
    spatialObj <- subset(spatialObj, subset = feature_outlier == FALSE)

    # Save the object 
    save(spatialObj, file = paste0("./spatial1/02_QC/spatialObjFiltered_", x, ".rds"))
}

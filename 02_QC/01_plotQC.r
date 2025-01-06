# Plot basic QC metrics for the spatial data
# Output: spatial1/02_QC/QC.pdf
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(SpatialExperiment)
})
setwd("/hamit/")

# Load spatial object from each subject and merge into one 
sampleName <- c("sample1", "sample2", "sample3")

# Create a list of Seurat objects
seuratList <- list()
for (x in sampleName) {
  load(paste0("spatial1/01_buildSeurat/spatialObj", x, ".rds"))
  spatialObj$orig.ident <- x
  seuratList[[x]] <- spatialObj
  rm(spatialObj)
}
spatialObj <- merge(seuratList[[1]], y = seuratList[-1], add.cell.ids = sampleName)
spatialObj <- JoinLayers(spatialObj)
spatialObj

# Calculate QC metrics
spatialObj[["percent.mt"]] <- PercentageFeatureSet(spatialObj, pattern = "^MT-") 
spatialObj$sum <- spatialObj$nCount_Spatial  # Total counts/reads
spatialObj$detected <- spatialObj$nFeature_Spatial  # Number of detected features
spatialObj$percent.mt <- spatialObj[["percent.mt"]]  # Percentage of mitochondrial genes

# Extract count 
data <- spatialObj@meta.data

# -------------------- Plot -------------------- # 
pdf("./spatial1/02_QC/QC.pdf")

# UMI counts per celll
ggplot(data, aes(color=orig.ident, x=sum, fill=orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        ylab("Cell density") +
        geom_vline(xintercept = 500) +
        ggtitle("Number of counts per cell")

# Number of genes detected per cell
ggplot(data, aes(color=orig.ident, x=detected, fill= orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        ylab("Number of genes") +
        geom_vline(xintercept = 1000) +
        ggtitle("Number of genes detected per cell")

# Number of mitochondria cells
ggplot(data, aes(color=orig.ident, x=percent.mt , fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2) +
    ggtitle("Mitochondrial ratio")

# Correlation between genes deceted and UMI counts
ggplot(data, aes(x=sum, y=detected, color=percent.mt )) + 
    	geom_point() + 
	    scale_colour_gradient(low = "gray90", high = "black") +
    	stat_smooth(method=lm) +
  	  scale_x_log10() + 
  	  scale_y_log10() + 
  	  theme_classic() +
  	  geom_vline(xintercept = 500) +
  	  geom_hline(yintercept = 250) +
  	  facet_wrap(~orig.ident) +
      ggtitle("Correlation between genes detected and UMI counts")

# Plot number of counts on the spatial slide for each sample 
plotList <- list()
for (x in sampleName) {
  tempPlotList <- list()
  for (i in c("nCount_Spatial", "nFeature_Spatial")) {
    seuratSubset <- subset(spatialObj, subset = orig.ident == x)
    tempPlotList[[i]] <- SpatialFeaturePlot(seuratSubset, features = i, alpha = c(0.1, 1)) +
                         theme(legend.position = "right")
  }
  plotList[[x]] <-  wrap_plots(tempPlotList, ncol = 1, nrow = 2)
}

# Pring imgages for each sample 
for (plot in plotList) {
  print(plot)
}

dev.off()


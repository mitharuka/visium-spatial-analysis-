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
sampleList <- c("sample1", "sample2", "sample3")

for (x in sampleList) {
    load(file = paste0("./spatial1/03_Labeltransfer/Labeled_", x, ".rds"))
    pdf(file = paste0("./spatial1/03_Labeltransfer/LabeledGeneMarker_", x, ".pdf"))

    # Plot image based on mitochondria
    plot1 <- SpatialFeaturePlot(spatialObj, 
        features = "percent.mt",
        pt.size.factor = 2, image.alpha = 0.3, crop = TRUE) +
            theme(legend.position = "right",
                legend.key.size = unit(0.2, "cm"),
                legend.text = element_text(size = 6), 
                legend.title = element_text(size = 6), 
                axis.text = element_text(size = 6), 
                axis.title = element_text(size = 6))

    # Excitatory clusters
    plot2 <- SpatialFeaturePlot(spatialObj, 
        features = paste0("Exc", 1:4),
        pt.size.factor = 2, image.alpha = 0.3, crop = TRUE)
    plot2 <- lapply(plot2, function(p) {
        p + theme(
            legend.position = "right",
            legend.key.size = unit(0.2, "cm"), 
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6), 
            axis.text = element_text(size = 6), 
            axis.title = element_text(size = 6))
    })

    # Plot expression of gene marker
    plot3 <- SpatialFeaturePlot(spatialObj, 
        features = c("MALM2", "PROX1", "SV2B", "SATB2", "TYRO3", "PFKP", "HS3ST4"), 
        pt.size.factor = 2, image.alpha = 0.3, crop = TRUE)
    plot3 <- lapply(plot3, function(p) {
        p + theme(
            legend.position = "right",
            legend.key.size = unit(0.2, "cm"), 
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6), 
            axis.text = element_text(size = 6), 
            axis.title = element_text(size = 6))
    })

    p <- plot_grid(plotlist = c(list(plot1), plot2, plot3),  ncol = 3, nrow = 4, scale = 1.0)
    print(p)

    dev.off()
}


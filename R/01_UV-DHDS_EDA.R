## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")
library(pheatmap)
library(ggpubr)
library(PoiClaClu)

##*********************************************************************************************************
## Exploratory data analysis

## load data
processed_data <- readRDS("Data/processed_data.rds")
dds <- processed_data$dds
se  <- processed_data$se
vsd_raw <- processed_data$vsd_unadjusted
vsd_adjusted <- processed_data$vsd_adj
rm(processed_data)

##******
## plot correlation matrix

## calculate Poisson distances and define parameters for graphical representation
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- colData(se)$sample
colnames(samplePoisDistMatrix) <- NULL
colorsPoisd <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

## plot heatmap
svg("Results/heatmap_poisson_distance.svg", width=15, height=13)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colorsPoisd)
dev.off()

##******
## check for outliers
svg("Results/cooks_distance.svg", width=15, height=6)
par(mar = c(8, 5, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
title("Cooks distance across samples")
dev.off()


##******
## plot PCA
p1 <- pca_plot(data = vsd_raw, pcsToUse = 1:2, intgroup = c("cell", "donor"), 
         title = "PCA with unadjusted VST data", subtitle = "PC1 vs. PC2")

## plot PCA data after removal of donor effect
p2 <- pca_plot(data = vsd_adjusted, pcsToUse = 1:2, intgroup = c("irradiation", "cell"), 
         title = "PCA with VST data adjusting for donor effect", subtitle = "PC1 vs. PC2")
p3 <- pca_plot(data = vsd_adjusted, pcsToUse = 2:3, intgroup = c("irradiation", "cell"), 
         title = "PCA with VST data adjusting for donor effect", subtitle = "PC2 vs. PC3")
p4 <- pca_plot(data = vsd_adjusted, pcsToUse = 2:3, intgroup = c("time", "irradiation"), 
         title = "PCA with VST data adjusting for donor effect", subtitle = "PC2 vs. PC3")

## save plot
svg("Results/PCA.svg", width=15, height=12)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()


##******
##* Plot dispersion
svg("Results/dispersion.svg", width= 7, height= 5)
plotDispEsts(dds)
dev.off()



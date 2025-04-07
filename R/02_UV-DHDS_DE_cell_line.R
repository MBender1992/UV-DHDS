## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(ashr)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(RColorBrewer)
library(ggprism)

##*********************************************************************************************************
## Differential expression analysis for differences in cell lines

## load data
processed_data <- readRDS("Data/processed_data.rds")
dds <- processed_data$dds
se  <- processed_data$se
anno <- processed_data$anno
vsd_adjusted <- processed_data$vsd_adj
pThres <- dds@thresholds$pvalue
lfcThres <- dds@thresholds$lfc
rm(processed_data)

## extract results
res <- results(dds, lfcThreshold = lfcThres, alpha = pThres, 
               contrast = contraster(dds,
                                     group1 = list(c("irradiation", "Control"), c("cell", "Melanocytes")),
                                     group2 = list(c("irradiation", "Control"), c("cell", "DSC"))))

##*********************************************************************************************************
## Ma plots

## shrink LFCs and annotate results object
resShrunk <- lfcShrink(dds, res = res, contrast = c("cell", "Melanocytes", "DSC"), type="ashr")
allRes <- cbind(ENSEMBL = rownames(res), res)
allRes <- if (dim(anno)[1] == dim(allRes)[1]) {
  left_join(as.data.frame(allRes), anno)
} else {
  stop("Dimensions of annotation object and result object are different.")
}
allRes <- allRes[!str_detect(allRes$ENSEMBL, "\\."), ]
sigRes <- allRes[!is.na(allRes$padj) & allRes$padj < pThres,]


# plot with original lfc
svg("Results/cell_line/MA_plot_original.svg", width=9, height=6)
plotMA(res, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with original log fold changes")
dev.off()

# shrink estimates and exlpore differences in MA plots
svg("Results/cell_line/MA_plot_shrunken.svg", width=9, height=6)
plotMA(resShrunk, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with shrunken log fold changes")
dev.off()


##*********************************************************************************************************
## Volcano plots

##******
##* Volcano plot with labelling of most strongly DE genes

## get results without preliminary filter
res_volcano <- results(dds)

# adding ENSEMBL gene ID as a column in significant differentially expression gene table.
dat_volcano <- cbind(ENSEMBL = rownames(res_volcano), res_volcano)
dat_volcano <- left_join(as.data.frame(dat_volcano), anno)
res_volcano <- NULL

## parameters for Volcano plot
pointSize <- 2
labSize   <- 4

# plot Volcano
# this only serves as diagnostic and representation tool. the actual fold change thresholds and pvalues change later due to the way the results function works
svg("Results/cell_line/volcano_topSignif.svg", width=12, height=9)
EnhancedVolcano(as.data.frame(dat_volcano), lab = dat_volcano$SYMBOL, selectLab = ifelse(dat_volcano$padj <= 1e-30 & abs(dat_volcano$log2FoldChange) >=8, dat_volcano$SYMBOL, NA),
                x = 'log2FoldChange', y = 'padj', gridlines.major = FALSE, gridlines.minor = FALSE,
                xlim = c(-21, 21), ylim = c(-5, 200), title = 'Melanocytes vs. DSCs', 
                pCutoff = pThres, FCcutoff = 2, pointSize = pointSize, labSize = labSize) 
dev.off()

##******
##* Volcano plot with labelling of skin marker genes

## put marker genes at the end of the dataframe for better plotting
dat_volcano_aux <- dat_volcano[!dat_volcano$SYMBOL %in% SkinMarkerGenes, ]
dat_volcano_markers <- rbind(dat_volcano_aux, dat_volcano[dat_volcano$SYMBOL %in% SkinMarkerGenes, ])

svg("Results/cell_line/volcano_markers.svg", width=12, height=9)
EnhancedVolcano(as.data.frame(dat_volcano_markers), lab = dat_volcano_markers$SYMBOL, selectLab = SkinMarkerGenes,
                x = 'log2FoldChange', y = 'padj', gridlines.major = FALSE, gridlines.minor = FALSE,
                xlim = c(-21, 21), ylim = c(-5, 200), title = 'Melanocytes vs. DSCs', 
                pCutoff = pThres, FCcutoff = 2, pointSize = pointSize, labSize = labSize)
dev.off()


##*********************************************************************************************************
## Heatmap

## specify fontsize
heatmapFontsize <- 18

## specify colors for heatmap
colors  <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_fun <- colorRamp2(c(-2, 0, 2), c(c(colors[1], colors[51], colors[100])))

## Heatmap of all DE genes
ind <- colData(dds)$irradiation == "Control"
matAll <- assay(vsd_adjusted)[sigRes$ENSEMBL, ind]

## scale matrix
datScaled <- matAll %>%  t() %>% scale() %>%  t()

## specify annotation bar
annoHeatmap <- as.data.frame(colData(vsd_adjusted)[ind, c("cell", "time")])
colnames(annoHeatmap) <- c("Cell_type",  "Time")

## specify annotation bar colors
annColors <- HeatmapAnnotation(df =annoHeatmap,
                               col = list(Cell_type = c("DSC" = ggsci::pal_npg("nrc")(4)[1], "Melanocytes" = ggsci::pal_npg("nrc")(4)[2]),
                                          Time      = c("6h" = ggsci::pal_npg("nrc")(4)[3], "24h" = ggsci::pal_npg("nrc")(4)[4])),
                               annotation_legend_param = list(Cell_type = list(nrow=1), 
                                                              Time = list(nrow = 1),
                                                              title_gp = gpar(fontsize = heatmapFontsize), 
                                                              labels_gp = gpar(fontsize = heatmapFontsize)),
                               annotation_name_gp = gpar(fontsize = heatmapFontsize, fontface = "bold"))

## plot Heatmap
Ht <- Heatmap(datScaled,col= col_fun,
              top_annotation = annColors,
              clustering_method_row = "average",
              clustering_method_columns = "average",
              clustering_distance_row = "pearson",
              clustering_distance_column = "euclidean",
              show_row_dend = FALSE,
              show_row_names =  FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 10),
              column_title = paste("Differentially Expressed genes (", nrow(sigRes), ")", sep =""),
              heatmap_legend_param = list(
                title = "row Z-score",
                at = seq(-2,2,by=1),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

svg("Results/cell_line/Heatmap.svg", width=12, height=15)
draw(Ht, merge_legend = TRUE)
dev.off()


##*********************************************************************************************************
## GO analysis

##******
## store upregulated genes in a list
genes_upregulated <- allRes[allRes$padj < pThres & allRes$log2FoldChange > 0 & !is.na(allRes$padj) & !is.na(allRes$ENTREZID), ]$ENTREZID
## plot enriched go terms
ls_go_upregulated <- plot_go(genes_upregulated, showCategory = 30, sym.colors = TRUhttp://127.0.0.1:16689/graphics/plot_zoom_png?width=1364&height=863E, path = "Results/cell_line/GOBP/", return.res = TRUE)

## search for certain go terms
## everything connected to melanocytes and pigmentation
search_go(enrich.res = ls_go_upregulated$enrich_go, all.res = allRes, path = "Results/cell_line/GOBP/", search.string = "[Mm]el|[Pp]ig")

##******
## store downregulated genes in a list
genes_downregulated <- allRes[allRes$padj < pThres & allRes$log2FoldChange < 0 & !is.na(allRes$padj) & !is.na(allRes$ENTREZID), ]$ENTREZID
## plot enriched pathways
ls_go_downregulated <- plot_go(genes_downregulated, showCategory = 30, sym.colors = TRUE, path = "Results/cell_line/GOBP/", return.res = TRUE)






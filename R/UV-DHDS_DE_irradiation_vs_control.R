## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(ComplexHeatmap)
library(RColorBrewer)
library(openxlsx)
library(gridExtra)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)

##*********************************************************************************************************
## Custom function only used within this script

## function to save multiple results from dds object into list
apply_contrasts <- function(trt, ctrl, annObj = NULL, shrink = FALSE){
  res <- results(dds, lfcThreshold = lfcThres, alpha = pThres, 
          contrast = contraster(dds,
                                group1 = list(c("condition", trt)),
                                group2 = list(c("condition", ctrl))))
  if(shrink == TRUE){
    message("Output contains shrunken log fold changes.")
    res <- lfcShrink(dds, res = res, contrast = c("condition", trt, ctrl), type="ashr")
  } else {
    message("Output contains original fold changes.")
  }
  
  if(!is.null(annObj)){
    allRes <- cbind(ENSEMBL = rownames(res), res)
    allRes <- if (dim(annObj)[1] == dim(allRes)[1]) {
      left_join(as.data.frame(allRes), annObj)
    } else {
      stop("Dimensions of annotation object and result object are different.")
    }
    allRes <- allRes[!str_detect(allRes$ENSEMBL, "\\."), ]
    
    if(shrink == TRUE){
      res_print <- as.data.frame(allRes)
      write.xlsx(res_print, paste0("Results/irradiation_vs_control/xlsx/",trt, "_vs_", ctrl, "_shrunken_LFC.xlsx"))
    }
    res <- allRes
  }
  res
}

## wrapper to plot volacno plot
plot_volcano <- function(results, title = ""){
  EnhancedVolcano(as.data.frame(results), lab = results$SYMBOL, 
                  x = 'log2FoldChange', y = 'padj', gridlines.major = FALSE, gridlines.minor = FALSE,
                  xlim = c(-8, 8), ylim = c(-0.5, 40),  title = title, 
                  pCutoff = pThres, FCcutoff = 1, pointSize = 3.5, labSize = 5.5) 
}

## function to list significantly up or downregulated genes (as ENTREZIDs) for cluster profiler
list_signif_genes <- function(x, data,  direction = c("greater", "lesser")){
  tmp <- data[[x]]
  if(direction == "greater"){
    res <- tmp[!is.na(tmp$padj) & tmp$padj < pThres & !is.na(tmp$ENTREZID) & tmp$log2FoldChange > 0, ]$ENTREZID
  }   else if(direction == "lesser"){
    res <- tmp[!is.na(tmp$padj) & tmp$padj < pThres & !is.na(tmp$ENTREZID) & tmp$log2FoldChange < 0, ]$ENTREZID
  } else {
    stop("Please specify a direction.")
  }
  res
} ## --> INCLUDE IN EKBSEQ

## function to do GO enrichment analysis on multiple lists of DE genes and plot results as dotplot
plot_go_clusters <- function(gene.list, path){
  names <- str_remove(deparse(substitute(gene.list)), "ls_")
  ## calculate clustered pathways
  ck <- compareCluster(geneCluster = gene.list, fun = enrichGO, ont = "BP", keyType = "ENTREZID", OrgDb = org.Hs.eg.db, pvalueCutoff = pThres)
  ck <- enrichplot::pairwise_termsim(ck) 
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  p <- dotplot(ck, font.size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  svg(paste0(path, names, "_GOBP.svg"), width=6, height=14)
  print(p)
  dev.off()
}

## direct contrasts between conditions (siehe Blog, mit diff2-diff1)

##*********************************************************************************************************
## Data wrangling

## extract unique conditions
conditions <- as.character(unique(colData$condition))

## define contrasts
input <- data.frame(treated  = c(conditions[seq(2,11,3)], conditions[seq(3,12,3)]),
                    controls = rep(conditions[seq(1, 10, 3)],2))
input$contrasts <- paste0(input[,1], "_vs_", input[,2])

## save multiple contrasts
res_shrunk <- mapply(apply_contrasts, trt = input[,1], ctrl = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = TRUE, annObj = anno))

## define names of contrasts
names(res_shrunk) <- input$contrasts

## extract log fold changes
lfcs <- lapply(1:length(res_shrunk), function(i){
  res_shrunk[[i]]$pflag <- ifelse(res_shrunk[[i]]$padj < pThres & !is.na(res_shrunk[[i]]$padj), TRUE, FALSE)
  data.frame(res_shrunk[[i]]$log2FoldChange, res_shrunk[[i]]$pflag)
})

## bind lfcs to data.frame
logFC <- as.matrix(do.call(cbind, lfcs))
rownames(logFC) <- rownames(res_shrunk[[1]])

## keep rows with all genes differentially express in at least one condition
logFC <- logFC[,seq(1,15,2)][which(rowSums(logFC[,seq(2,16,2)]) != 0),]
colnames(logFC) <- input$contrasts

##*********************************************************************************************************
## plot Heatmap of fold changes

## define fontsize for heatmap
heatmapFontsize <- 18

## define color function
colors  <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(c(colors[1], colors[51], colors[100])))

colnames(logFC) <- stringr::str_remove(colnames(logFC), "_vs.+")

Ht <- Heatmap(logFC, col= col_fun,
        #top_annotation = annColors,
        clustering_method_row = "average",
        clustering_method_columns = "average",
        clustering_distance_row = "pearson",
        clustering_distance_column = "pearson",
        show_row_dend = FALSE,
        show_row_names =  FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 9),
        heatmap_legend_param = list(
          title = "Log2 FC",
          at = seq(-2,2,by=1),
          color_bar="continuous",
          title_position ="topcenter",
          legend_direction = "horizontal",
          legend_width = unit(4, "cm")
        ))

svg("Results/irradiation_vs_control/heatmap_lfc.svg", width=12, height=22)
draw(Ht)
dev.off()


##*********************************************************************************************************
## Plot MA plots

res_MA <- mapply(apply_contrasts, trt = input[,1], ctrl = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE))

## loop over plotMA to get MA plots for the 8 conditions
ma_plots <- lapply(1:length(res_shrunk), function(i){
  svg(paste0("Results/irradiation_vs_control/MA_plots/MA_plot_",names(res_MA)[i],  ".svg"), width=5, height=4)
  plotMA(res_MA[[i]], ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = names(res_MA)[i])
  dev.off()
  })


##*********************************************************************************************************
## Plot volcano plots

## define new results object without shrinkage to plot as volcano plot
res_volcano <- mapply(apply_contrasts, trt = input[,1], ctrl = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE, annObj = anno))

## loop over plot_volcano to get volcano plots for the 8 conditions
volcanos <- lapply(1:length(res_volcano), function(i){
   plot_volcano(res_volcano[[i]], title = names(res_volcano)[i])
})

## save volcano plots
svg("Results/irradiation_vs_control/volcano_arranged.svg", width=11, height=20)
do.call("grid.arrange", c(volcanos, ncol=2))
dev.off()


##*********************************************************************************************************
## GO analysis

## store upregulated genes in a list
ls_upregulated <- lapply(1:length(res_shrunk), list_signif_genes, data = res_shrunk, direction = "greater")
names(ls_upregulated) <- names(res_shrunk)
## plot enriched pathways
plot_go_clusters(ls_upregulated, path = "Results/irradiation_vs_control/GOBP/")

## store downregulated genes in a list
ls_downregulated <- lapply(1:length(res_shrunk), list_signif_genes, data = res_shrunk, direction = "lesser")
names(ls_downregulated) <- names(res_shrunk)
## plot enriched pathways
plot_go_clusters(ls_downregulated, path = "Results/irradiation_vs_control/GOBP/")


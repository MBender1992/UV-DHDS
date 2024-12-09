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


##*********************************************************************************************************
## Custom function only used within this script

## function to save multiple results from dds object into list
apply_contrasts2 <- function(trt1, ctrl1, trt2, ctrl2, annObj = NULL, shrink = FALSE){
  
  # irradiated in group1
  diff1 <-  contraster(dds,
                       group1 = list(c("condition", trt1)), # in this context irradiated melanocytes
                       group2 = list(c("condition", ctrl1))) # in this context control melanocytes
  
  # irradiated in group 2
  diff2 <- contraster(dds,
                      group1 = list(c("condition", trt2)),
                      group2 = list(c("condition", ctrl2)))
  
  res <- results(dds, lfcThreshold = lfcThres, alpha = pThres,
                 contrast = diff1-diff2)
  
  if(shrink == TRUE){
    message("Output contains shrunken log fold changes.")
    res <- lfcShrink(dds, res = res, contrast = diff1-diff2, type="ashr")
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
      write.xlsx(res_print, paste0("Results/irradiation_DSC_vs_Melanocytes/xlsx/",trt1, "_vs_", trt2, "_shrunken_LFC.xlsx"))
    }
    res <- allRes
  }
  res
}

## wrapper to plot volacno plot
plot_volcano <- function(results, title = ""){
  EnhancedVolcano(as.data.frame(results), lab = results$SYMBOL, 
                  x = 'log2FoldChange', y = 'padj', gridlines.major = FALSE, gridlines.minor = FALSE,
                  xlim = c(-10, 10), ylim = c(-0.5, 25),  title = title, 
                  pCutoff = pThres, FCcutoff = 1, pointSize = 3.5, labSize = 5.5) 
}

##*********************************************************************************************************
## Data wrangling

## extract unique conditions
conditions <- as.character(unique(colData$condition))

## define contrasts
input <- data.frame(treated1  = conditions[c(8,9,11,12)],
                    controls1 = conditions[c(7,7,10,10)],
                    treated2  = conditions[c(2,3,5,6)],
                    controls2 = conditions[c(1,1,4,4)])
input$contrasts <- paste0(input[,3], "_vs_", input[,1])

## save multiple contrasts
res_shrunk <- mapply(apply_contrasts2, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = TRUE, annObj = anno))

## define names of contrasts
names(res_shrunk) <- input$contrasts

##*********************************************************************************************************
## Plot MA plots

res_MA <- mapply(apply_contrasts2, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE))

## loop over plotMA to get MA plots for the 8 conditions
ma_plots <- lapply(1:length(res_shrunk), function(i){
  svg(paste0("Results/irradiation_DSC_vs_Melanocytes/MA_plots/MA_plot_",names(res_MA)[i],  ".svg"), width=5, height=4)
  plotMA(res_MA[[i]], ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = names(res_MA)[i])
  dev.off()
  })

##*********************************************************************************************************
## Plot volcano plots

## define new results object without shrinkage to plot as volcano plot
res_volcano <- mapply(apply_contrasts2, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE, annObj = anno))

## loop over plot_volcano to get volcano plots for the 8 conditions
volcanos <- lapply(1:length(res_volcano), function(i){
   plot_volcano(res_volcano[[i]], title = names(res_volcano)[i])
})

## save volcano plots
svg("Results/irradiation_DSC_vs_Melanocytes/volcano_arranged.svg", width=12, height=16)
do.call("grid.arrange", c(volcanos, ncol=2))
dev.off()




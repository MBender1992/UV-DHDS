## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

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
library(limma)
library(ggprism)
library(UCell)
library(GO.db)

##*********************************************************************************************************
## Custom function only used within this script

## function to save multiple results from dds object into list
apply_contrasts <- function(trt1, ctrl1, trt2, ctrl2, annObj = NULL, shrink = FALSE){
  
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
                  pCutoff = pThres, FCcutoff = 1, pointSize = 2, labSize = 4) 
}

## blot gene expression data as boxlots from adjusted vsd object
plot_gene_expression <- function(vsd_obj, anno_obj, genes, ncol = 1){
  ind <- which(!is.na(anno_obj$SYMBOL) & anno_obj$SYMBOL %in% genes)
  vsd_counts <- assay(vsd_obj)
  df <- t(vsd_counts[ind,])
  colnames(df) <- anno$SYMBOL[match(colnames(df), anno$ENSEMBL)]
  df <- cbind(colData(vsd_obj), df) %>% 
    as.data.frame() %>%
    gather(gene, vsd, -colnames(.)[!colnames(.) %in% genes])
  
  df %>%
    ggplot(aes(condition, vsd, color = irradiation, shape = cell)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.7)) +
    theme_classic(base_size = 8)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = unit(c(1,1,1,1), "cm")) +
    guides(x = "prism_offset", y = "prism_offset") +
    scale_color_manual(values = c("grey50", "#2B7BBA", "#083674")) +
    facet_wrap(~gene, scales = "free_y", ncol = ncol)
}

## wrapper around complex heatmap
plot_heatmap <- function(vsd_obj, anno_obj, genes, row_split = 5, column_split = 3, return_data = FALSE){
  matAll <- assay(vsd_obj)[anno_obj$SYMBOL %in% genes & rowSums(assay(dds)) != 0 & !is.na(anno_obj$SYMBOL), ]
  rownames(matAll) <- anno_obj$SYMBOL[match(rownames(matAll), anno_obj$ENSEMBL)]
  datScaled <- matAll %>%  t() %>% scale() %>%  t()

  ht <- Heatmap(datScaled, col= col_fun,
          top_annotation = ha,
          clustering_method_row = "average",
          clustering_method_columns = "average",
          clustering_distance_row = "pearson",
          clustering_distance_column = "pearson",
          show_row_dend = TRUE,
          show_row_names =  TRUE,
          show_column_names = TRUE,
          column_split = column_split,
          split = row_split,
          border = TRUE,
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 9),
          heatmap_legend_param = list(
            title = "Log2 FC",
            at = seq(-2,2,by=1),
            color_bar="continuous",
            title_position ="topcenter",
            legend_direction = "horizontal",
            legend_width = unit(4, "cm")
          ))
  
  if(return_data == TRUE){
    return(datScaled)
  } else {
    return(ht)
  }
} 


##********************************************************************************************************
## Data wrangling

## load data
processed_data <- readRDS("Data/processed_data.rds")
dds <- processed_data$dds
se  <- processed_data$se
anno <- processed_data$anno
vsd_adjusted <- processed_data$vsd_adj
pThres <- dds@thresholds$pvalue
lfcThres <- dds@thresholds$lfc
rm(processed_data)

## extract unique conditions
conditions <- as.character(unique(colData$condition))

## define contrasts
input <- data.frame(treated1  = conditions[c(8,9,11,12)],
                    controls1 = conditions[c(7,7,10,10)],
                    treated2  = conditions[c(2,3,5,6)],
                    controls2 = conditions[c(1,1,4,4)])
input$contrasts <- paste0(input[,3], "_vs_", input[,1])

## save multiple contrasts
res_shrunk <- mapply(apply_contrasts, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = TRUE, annObj = anno))

## define names of contrasts
names(res_shrunk) <- input$contrasts

##*********************************************************************************************************
## Plot MA plots

res_MA <- mapply(apply_contrasts, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE))

## loop over plotMA to get MA plots for the 8 conditions
ma_plots <- lapply(1:length(res_shrunk), function(i){
  svg(paste0("Results/irradiation_DSC_vs_Melanocytes/MA_plots/MA_plot_",names(res_MA)[i],  ".svg"), width=5, height=4)
  plotMA(res_MA[[i]], ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = names(res_MA)[i])
  dev.off()
  })

##*********************************************************************************************************
## Plot volcano plots

## define new results object without shrinkage to plot as volcano plot
res_volcano <- mapply(apply_contrasts, trt1 = input[,3], ctrl1 = input[,4], trt2 = input[,1], ctrl2 = input[,2], SIMPLIFY = FALSE, MoreArgs = list(shrink = FALSE, annObj = anno))

## loop over plot_volcano to get volcano plots for the 8 conditions
volcanos <- lapply(1:length(res_volcano), function(i){
   plot_volcano(res_volcano[[i]], title = names(res_volcano)[i])
})

## save volcano plots
svg("Results/irradiation_DSC_vs_Melanocytes/volcano_arranged.svg", width=12, height=16)
do.call("grid.arrange", c(volcanos, ncol=2))
dev.off()


##*********************************************************************************************************
## GO analysis

## store upregulated genes in a list
ls_upregulated <- list_signif_genes(list = res_shrunk, p.threshold = pThres, direction = "greater")
names(ls_upregulated) <- names(res_shrunk)
## plot enriched pathways
plot_go_clusters(ls_upregulated, sym.colors = TRUE, path = "Results/irradiation_DSC_vs_Melanocytes/GOBP/")

## store downregulated genes in a list
ls_downregulated <- list_signif_genes(list = res_shrunk, p.threshold = pThres,  direction = "lesser")
names(ls_downregulated) <- names(res_shrunk)
## plot enriched pathways
plot_go_clusters(ls_downregulated, path = "Results/irradiation_DSC_vs_Melanocytes/GOBP/")


##*********************************************************************************************************
## targeted gene analysis

## get Kegg pathways
tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
                     column="SYMBOL", keytype="ENTREZID")

## extract names of interesting pathways
kegg_names <- getKEGGPathwayNames(species="hsa")
kegg_names %>%
  filter(str_detect(Description, "[Cc]ell [Cc]ycle|UV|Apoptosis|[Rr]epair|[Cc]hromatin"))

## retrieve GO terms
go_to_genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = keys(org.Hs.eg.db, keytype = "GO"),
                                     columns = c("GO", "ENTREZID"),
                                     keytype = "GO")
go_to_genes <- unique(go_to_genes)
go_to_genes$TERM <- AnnotationDbi::mapIds(GO.db, 
                                          keys = go_to_genes$GO,
                                          column = "TERM", 
                                          keytype = "GOID", 
                                          multiVals = "first")
go_to_genes$SYMBOL <- anno$SYMBOL[match(go_to_genes$ENTREZID, anno$ENTREZID)]

## extract cell cycle genes
cc_genes <- tab %>% 
  filter(PathwayID == "hsa04110") %>%
  .$Symbol

## extract apoptosis genes
apoptosis_genes <- tab %>% 
  filter(PathwayID == "hsa04210") %>%
  .$Symbol

## extract NER genes
ner_genes <- tab %>% 
  filter(PathwayID == "hsa03420") %>%
  .$Symbol

## extract chromatin remodeling genes
chr_remodeling_genes <- tab %>% 
  filter(PathwayID == "hsa03082") %>%
  .$Symbol

## cell differentiation genes
cell_differentiation_genes <- go_to_genes %>%
  filter(GO == "GO:0030154") %>%
  filter(!duplicated(SYMBOL) & !is.na(SYMBOL)) %>%
  .$SYMBOL

## EV genes
ev_genes <- go_to_genes %>%
  filter(GO == "GO:1903561") %>%
  filter(!duplicated(SYMBOL) & !is.na(SYMBOL)) %>%
  .$SYMBOL

## filter processes which contain GDF15
go_to_genes %>% 
  filter(SYMBOL == "GDF15") %>%
  filter(!duplicated(TERM))
 

## plot GDF and cell cycle genes
p <- plot_gene_expression(vsd_adjusted, anno, genes = cc_genes[str_detect(cc_genes, "CDK|E2F|RB[L|\\d]")], ncol = 4)
ggsave("Results/cell_cycle_genes.svg", p, width = 10, height = 15)
## plot apoptose genes
p <- plot_gene_expression(vsd_adjusted, anno, genes = apoptosis_genes[str_detect(apoptosis_genes, "CASP|BCL|CTS")], ncol = 4)
ggsave("Results/apoptosis_genes.svg", p, width = 10, height = 15)
## plot apoptose genes (which were elevated in 6h DSCs)
p <- plot_gene_expression(vsd_adjusted, anno, genes = c("TRADD", "PIK3R2", "TUBA4A", "EIF2S1", "BCL2L11", "GADD45A", "PMAIP1", "FAS", "TNFRSF10A", "BAK1", "PIDD1", "ERN1"), ncol = 3)
ggsave("Results/apoptosis_DSC_6h_high_genes.svg", p, width = 10, height = 15)
## plot repair genes
p <- plot_gene_expression(vsd_adjusted, anno, genes = c("ERCC5","XPA","ERCC3","POLR2J3","UVSSA","ERCC6","POLD3","ERCC8","GTF2H2C","POLR2B","POLD4",
                                                        "ERCC1","CETN2","GTF2H4","ERCC4","PCNA","DDB2","POLR2D","XPC","POLR2A"), ncol = 4)
ggsave("Results/repair_DSC_6h_high_genes.svg", p, width = 10, height = 15)

##******
## Heatmap preparation

## specify fontsize
heatmapFontsize <- 18

## specify colors for heatmap
colors  <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c(c(colors[1], colors[51], colors[100])))

## specify annotation bar
annoHeatmap <- as.data.frame(colData(vsd_adjusted)[, c("cell", "donor", "time", "irradiation")])
colnames(annoHeatmap) <- c("Cell_type", "Donor",  "Time", "Irradiation")
donors <- colorRampPalette(rev(brewer.pal(n = 7, name = "Set1")))(6)
names(donors) <- unique(colData(vsd_adjusted)$donor)
ha <- HeatmapAnnotation(
  Cell_type = annoHeatmap$Cell_type,
  Donor = annoHeatmap$Donor, # Spalte extrahieren
  Time = annoHeatmap$Time,
  Irradiation = annoHeatmap$Irradiation,
  col = list(
    Cell_type = c("DSC" = "#cab2d6", "Melanocytes" = "#33a02c"),
    Donor = donors,
    Time = c("6h" = "#fb6a4a", "24h" = "#67000d"),
    Irradiation = c("Control" = "grey70", "1x 300J/m² UVB" = "#2B7BBA",  "3x 300J/m² UVB" = "#083674")
  )
)

##******
## Plot heatmaps

## cc genes
svg("Results/Heatmap_cc_genes.svg", width = 10, height = 15)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = cc_genes, row_split = 4, column_split = 3), column_title = "Cell cycle genes (KEGG:hsa04110)")
dev.off()

## apoptosis genes
svg("Results/Heatmap_apoptosis_genes.svg", width = 10, height = 15)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = apoptosis_genes), column_title = "Apoptosis genes (KEGG:hsa04210)")
dev.off()

## NER genes
svg("Results/Heatmap_NER_genes.svg", width = 10, height = 15)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = ner_genes, row_split = 4, column_split = 3), column_title = "NER genes (KEGG:hsa03420)")
dev.off()

## NER genes
svg("Results/Heatmap_Chromatin_remodeling_genes.svg", width = 10, height = 15)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = c(chr_remodeling_genes, "EZH2", "DNMT1", "DNMT3A", "DNMT3B", "DROSHA", "DICER1")), 
         column_title = "Chromatin remodeling genes (KEGG:hsa03082) \n and DNMTs, EZH2, DROSHA & DICER")
dev.off()

## differentiation genes
svg("Results/Heatmap_cell_differentiation_genes.svg", width = 10, height = 35)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = cell_differentiation_genes), column_title = "Cell differentiation (GO:0030154)")
dev.off()

svg("Results/Heatmap_EV_cargo.svg", width = 10, height = 15)
draw(plot_heatmap(vsd_obj = vsd_adjusted, anno_obj = anno, genes = ev_genes, column_split = 2, row_split = 4), column_title = "EV cargo (GO:1903561)")
dev.off()



##*********************************************************************************************************
## Calculate gene set scores

## make list of gene sets
gene.sets <- list(CC_signature = cc_genes,
                  Apoptosis_signature = apoptosis_genes,
                  NER_signature = ner_genes,
                  Chr_remodeling_signature = chr_remodeling_genes,
                  Differentiation_signature = cell_differentiation_genes)

## score gene sets                
mat <- assay(vsd_adjusted)
rownames(mat) <- anno$SYMBOL[match(rownames(mat), anno$ENSEMBL)]
scores <- ScoreSignatures_UCell(mat, features=gene.sets, ncores = 10)
metadata <- colData(vsd_adjusted)
scores <- cbind(metadata, scores) %>%
  as.data.frame() %>%
  select(ID, sample, donor, condition, cell, irradiation, time, CC_signature_UCell, Apoptosis_signature_UCell, NER_signature_UCell,  Chr_remodeling_signature_UCell, Differentiation_signature_UCell) 

## plot scores
p <- scores %>%
  gather(signature, score, -c(ID, sample, donor, condition, cell, irradiation, time)) %>%
  group_by(signature) %>%
  mutate(zscore = scale(score)) %>%
  ggplot(aes(condition, zscore, color = irradiation, shape = cell)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_point(position = position_dodge(0.7)) +
  theme_classic(base_size = 8)+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  guides(x = "prism_offset", y = "prism_offset") +
  scale_color_manual(values = c("grey50", "#2B7BBA", "#083674")) +
  facet_wrap(~signature, scales = "free_y", ncol = 1)

## save plot
ggsave("Results/gene_signatures.svg", p, width = 5, height = 9)


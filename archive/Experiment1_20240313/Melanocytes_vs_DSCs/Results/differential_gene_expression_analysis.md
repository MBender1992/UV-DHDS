---
title: "Next-generation sequencing - analysis with DESeq2"
subtitle: "Comparison of gene expression in DSCs compared to melanocytes - differential expression analysis"
author: "Marc Bender"
date: '2024-06-04'
knit: (function(input, ...) {
    rmarkdown::render(
      input,
      output_dir = "../Results",
      output_format = "pdf_document"
    )
  })
output: 
  word_document:
    toc: true
    keep_md: yes
  pdf_document:
    toc: true
    keep_md: yes
  html_notebook:
    theme: united
    toc: true
  html_document:
    df_print: paged
    toc: true
    code_folding: "hide"
header-includes: 
- \usepackage[utf8]{inputenc} 
- \usepackage{pdflscape}
- \newcommand{\blandscape}{\begin{landscape}}
- \newcommand{\elandscape}{\end{landscape}}
---










\newpage 
# Resources and introduction

This document includes the differential expression analysis for the comparison of RNASeq data from dermal stem cells (DSCs) and Melanocytes. The exploratory data analysis and functional analysis (pathway analysis and GO terms) are included in separate files. Statistical methods are described in a separate .docx file. 

For additional info on analysis workflows check: 
DESEQ workflow: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 
https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html  
https://www.bigbioinformatics.org/r-and-rnaseq-analysis 
https://github.com/bigbioinformatics/r-programming-and-rnaseq-workshop 
ClusterProfiler: https://pubmed.ncbi.nlm.nih.gov/34557778/ 

\newpage
# Differential expression analysis

## Overview


```
## 
## out of 36904 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.58 (up)    : 3162, 8.6%
## LFC < -0.58 (down) : 3292, 8.9%
## outliers [1]       : 157, 0.43%
## low counts [2]     : 12443, 34%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

\newpage
## Volcano Plot (most differentially expressed genes labelled)

Volcano plots are commonly used to display the results of RNA-seq or other omics experiments. A volcano plot is a type of scatterplot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes with large fold changes that are also statistically significant. These may be the most biologically significant genes. In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left, and the most statistically significant genes are towards the top.

This plot shows the most significantly changed genes (<= 1e-30, logFC >= |8|). 

![Volcano plot showing differentially expressed genes between Melanocytes and DSCs.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/Volcano_plot_sig-1.pdf) 

\newpage
## Volcano Plot (skin cell marker genes labelled)

This plot shows 24 selected genes from a publication assessing differences between Melanocytes (12 genes) and DSCs (12 genes) (Zabierowski et al. 2012) as well as 12 genes characteristic for Fibroblasts. 

![Volcano plot showing differentially expressed genes between Melanocytes and DSCs with labelled marker genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/Volcano_plot_markers-1.pdf) 

 \newpage
## Volcano Plot (melanoma marker genes labelled)

This plot shows genes frequently mutated in melanoma [taken from Cherepakhin et. al (2022)] and differentially expressed genes [taken form Xie et. al (2022)]. In total 443 differentially regulated (or mutationally amplified/repressed) genes were either labelled as "upregulated in melanoma" or "downregulated in melanoma". Out of these 443 genes, 356 genes were expressed above the expression threshold. The following table shows the total number of genes which were identified to be up/downregulated in melanoma and the corresponding results of up/downregulation in our experiments: 



|                 |Melanoma_up |Melanoma_down |
|:----------------|:-----------|:-------------|
|Melanocytes_up   |NA (NA%)    |NA (NA%)      |
|Melanocytes_down |NA (NA%)    |NA (NA%)      |



\newpage
![Volcano plot showing differentially expressed genes between Melanocytes and DSCs with tumour marker genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/Volcano_plot_melanoma-1.pdf) 

\newpage
## MA plot

An MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated
coefficients in the model, e.g. the comparisons of interest, across all genes. On the y-axis, the “M”
stands for “minus” – subtraction of log values is equivalent to the log of the ratio – and on the x-axis,
the “A” stands for “average”. You may hear this plot also referred to as a mean-difference plot, or a
Bland-Altman plot.

Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison
of DSCs vs. Melanocytes. There are three types of shrinkage estimators in DESeq2, which are covered
in the DESeq2 vignette. Here we specify the apeglm method for shrinking coefficients, which is good for shrinking
the noisy LFC estimates while giving low bias LFC estimates for true large differences (Zhu, Ibrahim, and Love 2018).
To use apeglm we specify a coefficient from the model to shrink, either by name or number as the coefficient appears
in resultsNames(dds).

![MA plot with original log fold changes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/MA_plot-1.pdf) 

<br> 

![MA plot with shrunken log fold changes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/MA_plot_shrunk-1.pdf) 

\newpage
## Comparison of marker genes

## Expression of melanocyte markers

Marker genes were labelled according to the publication of Zabierowski 2012, lab practice and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4278762/. 

![Expression of melanocyte markers.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/mc_marker_genes-1.pdf) 

\newpage
## Expression of dsc markers

![Expression of DSC markers.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/dsc_marker_genes-1.pdf) 

\newpage
## Expression of fibroblast markers

![Expression of fibroblast markers.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/fb_marker_genes-1.pdf) 



\newpage
## Heatmap

If using either of Euclidean distance or Pearson correlation, data should follow a Gaussian /
normal (parametric) distribution. So, if coming from a microarray, anything from RMA normalisation is fine,
whereas, if coming from RNA-seq, any data deriving from a transformed normalised count metric should be fine,
such as variance-stabilised, regularised log, or log CPM expression levels.

If you are performing clustering on non-normal data, like 'normalised' [non-transformed] RNA-seq counts,
FPKM expression units, etc., then use Spearman correlation (non-parametric).



\newpage
### All differentially expressed (DE) genes
![Heatmap of all differentially expressed genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/Heatmap-1.pdf) 

\newpage
### Top 30 DE genes
![Heatmap of the Top30 differentially expressed genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/HeatmapTop30-1.pdf) 


\newpage
### Comparison of genes identified with the DRAGEN pipeline

This plot includes the same 30 genes which are depicted in the Heatmap included in the output generated from the DRAGEN analysis software. 

![Heatmap of genes identified as Top30 differentially expressed in the DRAGEN analysis pipeline.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/DSC_NGS/Results/differential_gene_expression_analysis_files/figure-latex/HeatmapDragen-1.pdf) 

\newpage
# Appenix

## Session Info

```
## R version 4.4.0 (2024-04-24 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 10 x64 (build 19041)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
## [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
## [5] LC_TIME=German_Germany.utf8    
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] circlize_0.4.16             RColorBrewer_1.1-3         
##  [3] ggprism_1.0.5               ComplexHeatmap_2.20.0      
##  [5] pheatmap_1.0.12             EnhancedVolcano_1.22.0     
##  [7] ggrepel_0.9.5               org.Hs.eg.db_3.19.1        
##  [9] AnnotationDbi_1.66.0        lubridate_1.9.3            
## [11] forcats_1.0.0               stringr_1.5.1              
## [13] dplyr_1.1.4                 purrr_1.0.2                
## [15] readr_2.1.5                 tidyr_1.3.1                
## [17] tibble_3.2.1                ggplot2_3.5.1              
## [19] tidyverse_2.0.0             DESeq2_1.44.0              
## [21] SummarizedExperiment_1.34.0 Biobase_2.64.0             
## [23] MatrixGenerics_1.16.0       matrixStats_1.3.0          
## [25] GenomicRanges_1.56.0        GenomeInfoDb_1.40.0        
## [27] IRanges_2.38.0              S4Vectors_0.42.0           
## [29] BiocGenerics_0.50.0        
## 
## loaded via a namespace (and not attached):
##   [1] rstudioapi_0.16.0       jsonlite_1.8.8          shape_1.4.6.1          
##   [4] magrittr_2.0.3          magick_2.8.3            farver_2.1.1           
##   [7] rmarkdown_2.26          GlobalOptions_0.1.2     zlibbioc_1.50.0        
##  [10] ragg_1.3.1              vctrs_0.6.5             memoise_2.0.1          
##  [13] askpass_1.2.0           tinytex_0.51            htmltools_0.5.8.1      
##  [16] S4Arrays_1.4.0          curl_5.2.1              SparseArray_1.4.1      
##  [19] cachem_1.0.8            uuid_1.2-0              mime_0.12              
##  [22] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
##  [25] Matrix_1.7-0            R6_2.5.1                fastmap_1.1.1          
##  [28] GenomeInfoDbData_1.2.12 shiny_1.8.1.1           clue_0.3-65            
##  [31] digest_0.6.35           colorspace_2.1-0        textshaping_0.3.7      
##  [34] RSQLite_2.3.6           labeling_0.4.3          fansi_1.0.6            
##  [37] timechange_0.3.0        httr_1.4.7              abind_1.4-5            
##  [40] compiler_4.4.0          bit64_4.0.5             fontquiver_0.2.1       
##  [43] withr_3.0.0             doParallel_1.0.17       BiocParallel_1.38.0    
##  [46] DBI_1.2.2               highr_0.10              openssl_2.1.2          
##  [49] DelayedArray_0.30.0     rjson_0.2.21            ggsci_3.0.3            
##  [52] gfonts_0.2.0            tools_4.4.0             zip_2.3.1              
##  [55] httpuv_1.6.15           glue_1.7.0              promises_1.3.0         
##  [58] cluster_2.1.6           generics_0.1.3          gtable_0.3.5           
##  [61] tzdb_0.4.0              data.table_1.15.4       hms_1.1.3              
##  [64] xml2_1.3.6              utf8_1.2.4              XVector_0.44.0         
##  [67] foreach_1.5.2           pillar_1.9.0            later_1.3.2            
##  [70] lattice_0.22-6          bit_4.0.5               ekbSeq_0.0.2           
##  [73] tidyselect_1.2.1        fontLiberation_0.1.0    locfit_1.5-9.9         
##  [76] Biostrings_2.72.0       knitr_1.46              fontBitstreamVera_0.1.1
##  [79] crul_1.4.2              xfun_0.43               stringi_1.8.4          
##  [82] UCSC.utils_1.0.0        yaml_2.3.8              evaluate_0.23          
##  [85] codetools_0.2-20        httpcode_0.3.0          officer_0.6.6          
##  [88] gdtools_0.3.7           cli_3.6.2               xtable_1.8-4           
##  [91] systemfonts_1.0.6       munsell_0.5.1           Rcpp_1.0.12            
##  [94] png_0.1-8               parallel_4.4.0          blob_1.2.4             
##  [97] scales_1.3.0            crayon_1.5.2            flextable_0.9.6        
## [100] GetoptLong_1.0.5        rlang_1.1.3             KEGGREST_1.44.0
```

---
title: "Next-generation sequencing - analysis with DESeq2"
subtitle: "Comparison of gene expression in DSCs after irradiation - differential expression analysis"
author: "Marc Bender"
date: '2024-07-04'
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
## out of 35739 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.38 (up)    : 68, 0.19%
## LFC < -0.38 (down) : 441, 1.2%
## outliers [1]       : 0, 0%
## low counts [2]     : 15354, 43%
## (mean count < 3)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

\newpage
## Volcano Plot (most differentially expressed genes labelled)

Volcano plots are commonly used to display the results of RNA-seq or other omics experiments. A volcano plot is a type of scatterplot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes with large fold changes that are also statistically significant. These may be the most biologically significant genes. In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left, and the most statistically significant genes are towards the top.

This plot shows the most significantly changed genes (<= 1e-30, logFC >= |8|). 

![Volcano plot showing differentially expressed genes between trt and ctrl.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/Volcano_plot_sig-1.pdf) 


\newpage
## Volcano Plot (melanoma markers labelled)

![Volcano plot showing differentially expressed genes between trt and ctrl.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/Volcano_plot_melanoma-1.pdf) 


\newpage
## MA plot

An MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated
coefficients in the model, e.g. the comparisons of interest, across all genes. On the y-axis, the “M”
stands for “minus” – subtraction of log values is equivalent to the log of the ratio – and on the x-axis,
the “A” stands for “average”. You may hear this plot also referred to as a mean-difference plot, or a
Bland-Altman plot.

Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison
of trt vs. ctrl. There are three types of shrinkage estimators in DESeq2, which are covered
in the DESeq2 vignette. Here we specify the apeglm method for shrinking coefficients, which is good for shrinking
the noisy LFC estimates while giving low bias LFC estimates for true large differences (Zhu, Ibrahim, and Love 2018).
To use apeglm we specify a coefficient from the model to shrink, either by name or number as the coefficient appears
in resultsNames(dds).

![MA plot with original log fold changes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/MA_plot-1.pdf) 

<br> 

![MA plot with shrunken log fold changes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/MA_plot_shrunk-1.pdf) 

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
![Heatmap of all differentially expressed genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/Heatmap-1.pdf) 

\newpage
### Top 30 DE genes
![Heatmap of the Top30 differentially expressed genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/HeatmapTop30-1.pdf) 

\newpage
# Publication Figures

\newpage
## Figure1

![Figure 1.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/differential_gene_expression_files/figure-latex/Figure1-1.pdf) 


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
##  [1] ekbSeq_0.0.4                RColorBrewer_1.1-3         
##  [3] circlize_0.4.16             ggpubr_0.6.0               
##  [5] ggprism_1.0.5               ComplexHeatmap_2.20.0      
##  [7] pheatmap_1.0.12             EnhancedVolcano_1.22.0     
##  [9] ggrepel_0.9.5               org.Hs.eg.db_3.19.1        
## [11] AnnotationDbi_1.66.0        lubridate_1.9.3            
## [13] forcats_1.0.0               stringr_1.5.1              
## [15] dplyr_1.1.4                 purrr_1.0.2                
## [17] readr_2.1.5                 tidyr_1.3.1                
## [19] tibble_3.2.1                ggplot2_3.5.1              
## [21] tidyverse_2.0.0             DESeq2_1.44.0              
## [23] SummarizedExperiment_1.34.0 Biobase_2.64.0             
## [25] MatrixGenerics_1.16.0       matrixStats_1.3.0          
## [27] GenomicRanges_1.56.0        GenomeInfoDb_1.40.0        
## [29] IRanges_2.38.0              S4Vectors_0.42.0           
## [31] BiocGenerics_0.50.0        
## 
## loaded via a namespace (and not attached):
##   [1] fs_1.6.4                 bitops_1.0-7             enrichplot_1.24.0       
##   [4] HDO.db_0.99.1            httr_1.4.7               doParallel_1.0.17       
##   [7] ggsci_3.0.3              tools_4.4.0              backports_1.4.1         
##  [10] utf8_1.2.4               R6_2.5.1                 lazyeval_0.2.2          
##  [13] mgcv_1.9-1               GetoptLong_1.0.5         withr_3.0.0             
##  [16] prettyunits_1.2.0        gridExtra_2.3            cli_3.6.2               
##  [19] textshaping_0.3.7        scatterpie_0.2.2         officer_0.6.6           
##  [22] labeling_0.4.3           slam_0.1-50              tm_0.7-13               
##  [25] goseq_1.56.0             askpass_1.2.0            Rsamtools_2.20.0        
##  [28] systemfonts_1.0.6        yulab.utils_0.1.4        gson_0.1.0              
##  [31] txdbmaker_1.0.0          gfonts_0.2.0             DOSE_3.30.0             
##  [34] rstudioapi_0.16.0        RSQLite_2.3.6            httpcode_0.3.0          
##  [37] treemap_2.4-4            gridGraphics_0.5-1       generics_0.1.3          
##  [40] shape_1.4.6.1            BiocIO_1.14.0            car_3.1-2               
##  [43] zip_2.3.1                GO.db_3.19.1             Matrix_1.7-0            
##  [46] fansi_1.0.6              abind_1.4-5              lifecycle_1.0.4         
##  [49] yaml_2.3.8               carData_3.0-5            qvalue_2.36.0           
##  [52] SparseArray_1.4.1        BiocFileCache_2.12.0     blob_1.2.4              
##  [55] promises_1.3.0           crayon_1.5.2             lattice_0.22-6          
##  [58] cowplot_1.1.3            GenomicFeatures_1.56.0   KEGGREST_1.44.0         
##  [61] magick_2.8.3             pillar_1.9.0             knitr_1.46              
##  [64] fgsea_1.30.0             rjson_0.2.21             codetools_0.2-20        
##  [67] fastmatch_1.1-4          glue_1.7.0               ggfun_0.1.4             
##  [70] fontLiberation_0.1.0     data.table_1.15.4        treeio_1.28.0           
##  [73] vctrs_0.6.5              png_0.1-8                gtable_0.3.5            
##  [76] cachem_1.0.8             xfun_0.43                S4Arrays_1.4.0          
##  [79] mime_0.12                tidygraph_1.3.1          iterators_1.0.14        
##  [82] tinytex_0.51             nlme_3.1-164             ggtree_3.12.0           
##  [85] bit64_4.0.5              fontquiver_0.2.1         progress_1.2.3          
##  [88] filelock_1.0.3           colorspace_2.1-0         rrvgo_1.16.0            
##  [91] DBI_1.2.2                tidyselect_1.2.1         bit_4.0.5               
##  [94] compiler_4.4.0           curl_5.2.1               httr2_1.0.1             
##  [97] BiasedUrn_2.0.11         flextable_0.9.6          xml2_1.3.6              
## [100] NLP_0.2-1                fontBitstreamVera_0.1.1  DelayedArray_0.30.0     
## [103] shadowtext_0.1.3         rtracklayer_1.64.0       scales_1.3.0            
## [106] rappdirs_0.3.3           digest_0.6.35            rmarkdown_2.26          
## [109] XVector_0.44.0           htmltools_0.5.8.1        pkgconfig_2.0.3         
## [112] umap_0.2.10.0            highr_0.10               dbplyr_2.5.0            
## [115] fastmap_1.1.1            rlang_1.1.3              GlobalOptions_0.1.2     
## [118] UCSC.utils_1.0.0         shiny_1.8.1.1            farver_2.1.1            
## [121] jsonlite_1.8.8           BiocParallel_1.38.0      GOSemSim_2.30.0         
## [124] RCurl_1.98-1.14          magrittr_2.0.3           ggplotify_0.1.2         
## [127] GenomeInfoDbData_1.2.12  wordcloud_2.6            patchwork_1.2.0         
## [130] munsell_0.5.1            Rcpp_1.0.12              ape_5.8                 
## [133] viridis_0.6.5            gdtools_0.3.7            reticulate_1.37.0       
## [136] stringi_1.8.4            ggraph_2.2.1             zlibbioc_1.50.0         
## [139] MASS_7.3-60.2            plyr_1.8.9               parallel_4.4.0          
## [142] Biostrings_2.72.0        graphlayouts_1.1.1       splines_4.4.0           
## [145] hms_1.1.3                geneLenDataBase_1.40.1   locfit_1.5-9.9          
## [148] igraph_2.0.3             uuid_1.2-0               ggsignif_0.6.4          
## [151] reshape2_1.4.4           biomaRt_2.60.0           crul_1.4.2              
## [154] XML_3.99-0.16.1          evaluate_0.23            tzdb_0.4.0              
## [157] foreach_1.5.2            tweenr_2.0.3             httpuv_1.6.15           
## [160] openssl_2.1.2            polyclip_1.10-6          clue_0.3-65             
## [163] gridBase_0.4-7           ggforce_0.4.2            broom_1.0.5             
## [166] xtable_1.8-4             restfulr_0.0.15          tidytree_0.4.6          
## [169] RSpectra_0.16-1          rstatix_0.7.2            later_1.3.2             
## [172] viridisLite_0.4.2        ragg_1.3.1               aplot_0.2.2             
## [175] clusterProfiler_4.12.0   memoise_2.0.1            GenomicAlignments_1.40.0
## [178] cluster_2.1.6            timechange_0.3.0
```

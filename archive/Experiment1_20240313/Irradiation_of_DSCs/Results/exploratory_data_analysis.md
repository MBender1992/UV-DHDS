---
title: "Next-generation sequencing - analysis with DESeq2"
subtitle: "Comparison of gene expression in DSCs after irradiation - exploratory data analysis"
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
    fig_caption: yes
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

This document includes the exploratory data analysis for the comparison of RNASeq data from dermal stem cells (DSCs) and Melanocytes. The differential expression analysis and functional analysis (pathway analysis and GO terms) are included in separate files. Statistical methods are described in a separate .docx file. 

For additional info on analysis workflows check: 
DESEQ workflow: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 
https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html  
https://www.bigbioinformatics.org/r-and-rnaseq-analysis 
https://github.com/bigbioinformatics/r-programming-and-rnaseq-workshop 
ClusterProfiler: https://pubmed.ncbi.nlm.nih.gov/34557778/ 

\newpage
# Type of analysis

```
## [1] "This analysis displays results adjusted for the donor effect"
```

\newpage
# Thresholds, design formula and number of transcripts


```{=latex}
\global\setlength{\Oldarrayrulewidth}{\arrayrulewidth}

\global\setlength{\Oldtabcolsep}{\tabcolsep}

\setlength{\tabcolsep}{2pt}

\renewcommand*{\arraystretch}{1.5}



\providecommand{\ascline}[3]{\noalign{\global\arrayrulewidth #1}\arrayrulecolor[HTML]{#2}\cline{#3}}

\begin{longtable}[c]{|p{4.21in}|p{2.79in}}



\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endfirsthead 

\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endhead



\multicolumn{2}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 7in+2\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{\textsuperscript{†}}}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{All\ transcripts\ which\ are\ above\ the\ minimum\ read\ count\ threshold\ (10)\ in\ at}}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ least\ 5\ samples.}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endlastfoot



\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{\textbf{Thresholds}}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{P-value}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{0.05}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{Log-fold\ change}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{0.3785}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{Minimum\ read\ count}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{10}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{Smallest\ group\ size}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{5}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{\textbf{Design\ Formula}}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{\textasciitilde donor\ +\ treatment}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{\textbf{Transcript\ number}}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{Total\ number\ of\ transcripts}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{60609}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 4.21in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{Filtered\ number\ of\ transcripts}}\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{\textsuperscript{†}}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 2.79in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{4}\selectfont{60564}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\end{longtable}



\arrayrulecolor[HTML]{000000}

\global\setlength{\arrayrulewidth}{\Oldarrayrulewidth}

\global\setlength{\tabcolsep}{\Oldtabcolsep}

\renewcommand*{\arraystretch}{1}
```


```{=latex}
\global\setlength{\Oldarrayrulewidth}{\arrayrulewidth}

\global\setlength{\Oldtabcolsep}{\tabcolsep}

\setlength{\tabcolsep}{2pt}

\renewcommand*{\arraystretch}{1.5}



\providecommand{\ascline}[3]{\noalign{\global\arrayrulewidth #1}\arrayrulecolor[HTML]{#2}\cline{#3}}

\begin{longtable}[c]{|p{0.85in}|p{0.64in}|p{0.59in}|p{0.64in}|p{0.59in}|p{0.64in}|p{0.59in}|p{0.64in}|p{0.59in}|p{0.64in}|p{0.59in}}



\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.85in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{27\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{27\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{28\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{28\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{29\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{29\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{30\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{30\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{31\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{31\_23\_trt}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endfirsthead 

\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.85in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{27\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{27\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{28\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{28\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{29\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{29\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{30\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{30\_23\_trt}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{31\_23\_ctrl}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{DSC}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{31\_23\_trt}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endhead



\multicolumn{11}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 7in+20\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\textsuperscript{†}}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{All\ transcripts\ which\ are\ above\ the\ minimum\ read\ count\ threshold\ (10)\ in\ at}}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\linebreak }}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textbf{\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ least5\ samples.}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endlastfoot



\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.85in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{Number\ of\ reads}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\ (total)}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{3.98e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{5.34e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{4.52e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{4.79e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{4.35e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{3.50e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{3.94e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{4.84e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{4.21e+07}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{3.86e+07}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 0.85in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{Number\ of\ reads}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\ (filtered)}}\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{\textsuperscript{†}}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.02e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.30e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.21e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{1.99e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.11e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.71e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.43e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{1.78e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.64in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{2.45e+07}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.59in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{8}{8}\selectfont{1.93e+07}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\end{longtable}



\arrayrulecolor[HTML]{000000}

\global\setlength{\arrayrulewidth}{\Oldarrayrulewidth}

\global\setlength{\tabcolsep}{\Oldtabcolsep}

\renewcommand*{\arraystretch}{1}
```

\newpage
# Sample sheet data


```{=latex}
\global\setlength{\Oldarrayrulewidth}{\arrayrulewidth}

\global\setlength{\Oldtabcolsep}{\tabcolsep}

\setlength{\tabcolsep}{2pt}

\renewcommand*{\arraystretch}{1.5}



\providecommand{\ascline}[3]{\noalign{\global\arrayrulewidth #1}\arrayrulecolor[HTML]{#2}\cline{#3}}

\begin{longtable}[c]{|p{1.74in}|p{1.72in}|p{1.46in}|p{0.89in}|p{1.19in}}



\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{[Header]}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value1}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value2}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value3}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value4}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endfirsthead 

\hhline{>{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}->{\arrayrulecolor[HTML]{000000}\global\arrayrulewidth=0pt}-}

\multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{[Header]}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value1}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value2}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value3}}}} & \multicolumn{1}{>{\cellcolor[HTML]{CFCFCF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\textbf{Value4}}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}

\endhead



\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{FileFormatVersion}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{2}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{RunName}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{mRNA\_Test\_Buxtehude}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{InstrumentPlatform}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{NextSeq1k2k}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IndexOrientation}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Forward}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[Reads]}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Read1Cycles}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{76}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Read2Cycles}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{76}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Index1Cycles}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Index2Cycles}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[Sequencing\_Settings]}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{LibraryPrepKits}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStrandedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[BCLConvert\_Settings]}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SoftwareVersion}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{03.10.2012}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AdapterRead1}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTGTCTCTTATACACATCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AdapterRead2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTGTCTCTTATACACATCT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{FastqCompressionFormat}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{gzip}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[BCLConvert\_Data]}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Sample\_ID}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Index}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Index2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{1}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GAACTGAGCG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TCGTGGAGCG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AGGTCAGATA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTACAAGATA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{3}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CGTCTCATAT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TATAGTAGCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{4}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ATTCCATAAG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TGCCTGGTGG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{5}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GACGAGATTA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ACATTATCCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{6}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AACATCGCGC}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTCCACTTGT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{7}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTAGTGCTCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TGGAACAGTA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{8}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GATCAAGGCA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CCTTGTTAAT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{9}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GACTGAGTAG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTTGATAGTG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AGTCAGACGA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ACCAGCGACA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{11}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CCGTATGTTC}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CATACACTGT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{12}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GAGTCATAGG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTGTGGCGCT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{13}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTTGCCATTA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ATCACGAAGG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{14}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GAAGCGGCAC}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CGGCTCTACT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{15}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TCCATTGCCG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GAATGCACGA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{UHR}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CGGTTACGGC}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AAGACTATAG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[DragenRna\_Settings]}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SoftwareVersion}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{03.10.2012}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ReferenceGenomeDir}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{hg38\_alt\_aware}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{MapAlignOutFormat}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{bam}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{KeepFastq}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{true}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{DifferentialExpressionEnable}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{true}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[DragenRna\_Data]}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Sample\_ID}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Comparison1}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Comparison2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{1}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{3}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{4}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{5}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{6}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{7}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{8}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{9}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{control}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{11}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{12}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{13}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{14}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{15}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{comparison}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{UHR}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{na}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[Cloud\_Settings]}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GeneratedVersion}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{1.12.0.202401022340}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{[Cloud\_Data]}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{Sample\_ID}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ProjectName}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{LibraryName}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{LibraryPrep-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{KitName}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IndexAdapter-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{KitName}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{1}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{1\_GAACTGAGCG\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TCGTGGAGCG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{2}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{2\_AGGTCAGATA\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CTACAAGATA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{3}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{3\_CGTCTCATAT\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TATAGTAGCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{4}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{4\_ATTCCATAAG\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TGCCTGGTGG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{5}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{5\_GACGAGATTA\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ACATTATCCT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{6}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{6\_AACATCGCGC\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTCCACTTGT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{7}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{7\_CTAGTGCTCT\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{TGGAACAGTA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{8}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{8\_GATCAAGGCA\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CCTTGTTAAT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{9}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{9\_GACTGAGTAG\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTTGATAGTG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{10\_AGTCAGACGA\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ACCAGCGACA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{11}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{11\_CCGTATGTTC\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CATACACTGT}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{12}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{12\_GAGTCATAGG\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GTGTGGCGCT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{13}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{13\_CTTGCCATTA\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ATCACGAAGG}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{14}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{14\_GAAGCGGCAC\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{CGGCTCTACT}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{15}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{15\_TCCATTGCCG\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{GAATGCACGA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\cellcolor[HTML]{EFEFEF}\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\multicolumn{1}{>{\raggedright}m{\dimexpr 1.74in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{UHR}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.72in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.46in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{UHR\_CGGTTACGGC\_}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{AAGACTATAG}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 0.89in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{ILMNStran-}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{dedmRNA}}} & \multicolumn{1}{>{\raggedright}m{\dimexpr 1.19in+0\tabcolsep}}{\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{IDTIlmnNRNAUDI}}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{\linebreak }}\textcolor[HTML]{000000}{\fontsize{9}{9}\selectfont{SetALigation}}} \\

\noalign{\global\arrayrulewidth 0pt}\arrayrulecolor[HTML]{000000}





\end{longtable}



\arrayrulecolor[HTML]{000000}

\global\setlength{\arrayrulewidth}{\Oldarrayrulewidth}

\global\setlength{\tabcolsep}{\Oldtabcolsep}

\renewcommand*{\arraystretch}{1}
```

\newpage
# Exploratory data analysis

## Heatmap of Poisson distances between samples

A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples.
* IMPORTANT:
  + use Poisson distance for raw (non-normalized) count data
  + use Euclidean distance for data normalized by regularized-logarithm transformation (rlog) or variance stablization transfromation (vst)

![Correlation of poisson distance calculated on raw counts.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/poisd-1.pdf) 

\newpage
## Transform data

Which transformation to choose? The VST is much faster to compute and is less sensitive to high count outliers
than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST
when there is a wide range of sequencing depth across samples (an order of magnitude difference).
We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and
compare the meanSdPlot or PCA plots generated, as described below.

To show the effect of the transformation, in the figure below we plot the first sample against the second,
first simply using the log2 function (after adding 1, to avoid taking the log of zero), and then using the
VST and rlog-transformed values. For the log2 approach, we need to first estimate size factors to account for
sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the
vst and rlog.

![Representative depiction of count transformation.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/plot_transformations-1.pdf) 

<!-- \newpage -->
<!-- ## Effect of transformation on the variance -->

<!-- The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the  the regularized log transformation.  -->

<!-- Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions. While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions. -->

<!-- ```{r} -->
<!-- meanSdPlot(assay(vsd)) -->
<!-- ``` -->

\newpage
## Check for outliers

RNA-seq data sometimes contain isolated instances of very large counts that are apparently unrelated to the experimental or study design, and which may be considered outliers. There are many reasons why outliers can arise, including rare technical or experimental artifacts, read mapping problems in the case of genetically differing samples, and genuine, but rare biological events. In many cases, users appear primarily interested in genes that show a consistent behavior, and this is the reason why by default, genes that are affected by such outliers are set aside by DESeq2, or if there are sufficient samples, outlier counts are replaced for model fitting. These two behaviors are described below.

The DESeq function calculates, for every gene and for every sample, a diagnostic test for outliers called Cook’s distance. Cook’s distance is a measure of how much a single sample is influencing the fitted coefficients for a gene, and a large value of Cook’s distance is intended to indicate an outlier count. The Cook’s distances are stored as a matrix available in assays(dds)[["cooks"]].

The results function automatically flags genes which contain a Cook’s distance above a cutoff for samples which have 3 or more replicates. The p values and adjusted p values for these genes are set to NA. At least 3 replicates are required for flagging, as it is difficult to judge which sample might be an outlier with only 2 replicates. This filtering can be turned off with results(dds, cooksCutoff=FALSE).

With many degrees of freedom – i.,e., many more samples than number of parameters to be estimated – it is undesirable to remove entire genes from the analysis just because their data include a single count outlier. When there are 7 or more replicates for a given sample, the DESeq function will automatically replace counts with large Cook’s distance with the trimmed mean over all samples, scaled up by the size factor or normalization factor for that sample. This approach is conservative, it will not lead to false positives, as it replaces the outlier value with the value predicted by the null hypothesis. This outlier replacement only occurs when there are 7 or more replicates, and can be turned off with DESeq(dds, minReplicatesForReplace=Inf).

The default Cook’s distance cutoff for the two behaviors described above depends on the sample size and number of parameters to be estimated. The default is to use the 99% quantile of the F(p,m-p) distribution (with p the number of parameters including the intercept and m number of samples). The default for gene flagging can be modified using the cooksCutoff argument to the results function. For outlier replacement, DESeq preserves the original counts in counts(dds) saving the replacement counts as a matrix named replaceCounts in assays(dds). Note that with continuous variables in the design, outlier detection and replacement is not automatically performed, as our current methods involve a robust estimation of within-group variance which does not extend easily to continuous covariates. However, users can examine the Cook’s distances in assays(dds)[["cooks"]], in order to perform manual visualization and filtering if necessary.

Note on many outliers: if there are very many outliers (e.g. many hundreds or thousands) reported by summary(res), one might consider further exploration to see if a single sample or a few samples should be removed due to low quality. The automatic outlier filtering/replacement is most useful in situations which the number of outliers is limited. When there are thousands of reported outliers, it might make more sense to turn off the outlier filtering/replacement (DESeq with minReplicatesForReplace=Inf and results with cooksCutoff=FALSE) and perform manual inspection: First it would be advantageous to make a PCA plot as described above to spot individual sample outliers; Second, one can make a boxplot of the Cook’s distances to see if one sample is consistently higher than others (here this is not the case):

![Sample-wise cook's distances.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/cooks_distance-1.pdf) 

\newpage
## Principal component analysis (PCA)

Die Hauptkomponentenanalyse oder *principal component analysis* (PCA) ist ein Verfahren zur Reduktion der Datendimensionalität und findet Anwendung bei komplexen Datensätzen (aus n Proben mit p gemessenen Merkmalen), um Muster oder Strukturen in den Daten zu erkennen und die Interpretation der Daten zu vereinfachen. Dies erfolgt durch Identifikation neuer unkorrelierter Variablen, den sog. Hauptkomponenten (PC = *principal component*), welche sukzessive die erklärte Varianz innerhalb eines Datensets maximieren. Die Berechnung dieser Hauptkomponenten erfolgt dabei durch Zentrieren der Daten und einen nachfolgenden iterativen Prozess, bei dem eine zufällige Gerade durch den Mittelwert der Datenpunkte gelegt und so lange rotiert wird, bis der Abstand der auf die Gerade projizierten Datenpunkte zum Zentrum maximal ist.
Die Distanzen werden berechnet als:

\begin{align*}
SS(distances) = d^2_{1} + d^2_{2} + ... + d^2_{n}
\end{align*}


mit SS(distances) = Eigenwert (Quadratsumme der Distanzen), d = Distanz vom projizierten Datenpunkt zum Zentrum, n = Probengröße.

Die Gerade mit dem größten Eigenwert wird als PC1 bezeichnet und trägt am meisten zur Varianzaufklärung bei. Im Anschluss wird die nächste Gerade gesucht, die orthogonal zur ersten Gerade ist und ebenfalls durch den Mittelwert der Datenpunkte verläuft. Nach diesem Grundsatz können (bis zur n-ten Hauptkomponente) weitere Hauptkomponenten gefunden werden, wobei diese stets orthogonal zu den bisherigen Hauptkomponenten stehen müssen.

Die Reduktion der Dimensionalität beruht darauf, dass die jeweiligen Hauptkomponenten einen unterschiedlich großen (absteigend von PC1 zu PCn) Anteil der Datenvarianz erklären und es daher in den meisten Fällen ausreicht die ersten drei Hauptkomponenten zu analysieren (welche jeweils als 2-dimensionaler Plot gegeneinander aufgetragen werden), um Gemeinsamkeiten innerhalb der Proben zu identifizieren.
				Der Anteil der Varianz, der durch die jeweilige Hauptkomponente erklärt wird ($var_{explained}$), kann anhand folgender Formel  berechnet werden:

\begin{align*}
var(PC_{i}) = \frac{SS(distances (PC_{i}))}{n-1}
\end{align*}

\begin{align*}
var_{explained} = \frac{var(PC_{i})}{\sum^{n}_{j=1}var(PC_{j})}
\end{align*}


mit $PC_{i}$ = Hauptkomponente i. \newline

Die PCA wurde in dem Statistikprogramm R mithilfe der Pakete *factoextra* und *FactoMineR* durchgeführt. Die Ermittlung relevanter Hauptkomponenten erfolgte mittels Skree Plot. Dort werden entweder die Eigenwerte oder die erklärte Varianz gegen die Nummer der jeweiligen Hauptkomponente aufgetragen und mittels einer zunächst steil abfallenden Linie, die sich asymptotisch der Abzisse annähert, verbunden. An dem sog. *elbow* ist ein deutlicher Knick in der Linie zu erkennen. Faktoren, die links dieses Punktes liegen reichen größtenteils aus, um Muster in den Daten zu erklären. Die anderen Faktoren unterscheiden sich oft nicht deutlich von Zufallskorrelationen und sind daher nicht von Bedeutung.

![Principal component analysis with vst (left) and rlog (right).](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/PCA-1.pdf) 

\
Add a labelled variant to identify possible batch effects.

![Principal component analysis with vst (left) and rlog (right) with labelled samples.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/PCA_labelled-1.pdf) 

\newpage
## Dispersion plot

Plotting the dispersion estimates is a useful diagnostic. The dispersion plot below shows the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value, (this outlier detection is described in the manual page for estimateDispersionsMAP). The amount of shrinkage can vary, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.

![Dispersion plot.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/dispersion_plot-1.pdf) 

\
There was a minimal batch effect present in the PCA plot (sample 27). Therefore a more detailled analysis of the dispersion is useful to further analyze this effect. In the following all genes with a dispersion estimate >1 are shown for DSCs and Melanocytes, respectively. The first plot depicts rlog gene counts for genes with a dispersion estimate >1 in DSCs. The second plot depicts rlog gene counts for genes with a dispersion estimate >1 in Melanocytes. The batch is indicated by different shapes.

\newpage
\blandscape
![Rlog counts in ctrl samples for high dispersion genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/high_disp_dsc-1.pdf) 

\elandscape

\newpage
\blandscape
![Rlog counts in trt samples for high dispersion genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/high_disp_mc-1.pdf) 

\elandscape

\newpage
## Scatter plot matrices

These types of plots were suggested in https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2968-1.
They depict the correlation pairs of each sample with other samples. The spread of points should be closer to the identity line (red line) in technical replicates than in biological replicates because a larger spread indicates higher differences in expression values most often associated with differential expression. Divergence of the expected pattern (as described in the aforemented publication) should also warrant a closer inspection.

![Scatterplot matrix for rlog counts.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/scatter_matrix-1.pdf) 

\newpage
Scatter plot matrix indicating genes with high dispersion in red.

![Scatterplot matrix for rlog counts with colored high dispersion genes.](C:/MBender/Arbeit/Github/EKB/UV-DHDS/Irradiation_of_DSCs/Results/exploratory_data_analysis_files/figure-latex/scatter_matrix_disp-1.pdf) 

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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] kableExtra_1.4.0            ggprism_1.0.5              
##  [3] ekbSeq_0.0.4                GGally_2.2.1               
##  [5] pheatmap_1.0.12             flextable_0.9.6            
##  [7] PoiClaClu_1.0.2.1           lubridate_1.9.3            
##  [9] forcats_1.0.0               stringr_1.5.1              
## [11] dplyr_1.1.4                 purrr_1.0.2                
## [13] readr_2.1.5                 tidyr_1.3.1                
## [15] tibble_3.2.1                ggplot2_3.5.1              
## [17] tidyverse_2.0.0             DESeq2_1.44.0              
## [19] SummarizedExperiment_1.34.0 Biobase_2.64.0             
## [21] MatrixGenerics_1.16.0       matrixStats_1.3.0          
## [23] GenomicRanges_1.56.0        GenomeInfoDb_1.40.0        
## [25] IRanges_2.38.0              S4Vectors_0.42.0           
## [27] BiocGenerics_0.50.0        
## 
## loaded via a namespace (and not attached):
##   [1] progress_1.2.3           nnet_7.3-19              Biostrings_2.72.0       
##   [4] vctrs_0.6.5              digest_0.6.35            png_0.1-8               
##   [7] shape_1.4.6.1            ggrepel_0.9.5            rrvgo_1.16.0            
##  [10] httpcode_0.3.0           MASS_7.3-60.2            fontLiberation_0.1.0    
##  [13] reshape2_1.4.4           httpuv_1.6.15            foreach_1.5.2           
##  [16] qvalue_2.36.0            withr_3.0.0              xfun_0.43               
##  [19] ggfun_0.1.4              ggpubr_0.6.0             survival_3.6-4          
##  [22] memoise_2.0.1            crul_1.4.2               hexbin_1.28.3           
##  [25] clusterProfiler_4.12.0   gson_0.1.0               BiasedUrn_2.0.11        
##  [28] ggsci_3.0.3              systemfonts_1.0.6        ragg_1.3.1              
##  [31] tidytree_0.4.6           zoo_1.8-12               GlobalOptions_0.1.2     
##  [34] Formula_1.2-5            prettyunits_1.2.0        KEGGREST_1.44.0         
##  [37] promises_1.3.0           httr_1.4.7               rstatix_0.7.2           
##  [40] restfulr_0.0.15          rstudioapi_0.16.0        pan_1.9                 
##  [43] UCSC.utils_1.0.0         generics_0.1.3           DOSE_3.30.0             
##  [46] curl_5.2.1               zlibbioc_1.50.0          ggraph_2.2.1            
##  [49] polyclip_1.10-6          GenomeInfoDbData_1.2.12  SparseArray_1.4.1       
##  [52] xtable_1.8-4             doParallel_1.0.17        evaluate_0.23           
##  [55] S4Arrays_1.4.0           BiocFileCache_2.12.0     hms_1.1.3               
##  [58] glmnet_4.1-8             colorspace_2.1-0         filelock_1.0.3          
##  [61] NLP_0.2-1                reticulate_1.37.0        treemap_2.4-4           
##  [64] magrittr_2.0.3           later_1.3.2              viridis_0.6.5           
##  [67] ggtree_3.12.0            lattice_0.22-6           XML_3.99-0.16.1         
##  [70] shadowtext_0.1.3         cowplot_1.1.3            pillar_1.9.0            
##  [73] nlme_3.1-164             iterators_1.0.14         gridBase_0.4-7          
##  [76] compiler_4.4.0           plotROC_2.3.1            RSpectra_0.16-1         
##  [79] stringi_1.8.4            jomo_2.7-6               minqa_1.2.6             
##  [82] GenomicAlignments_1.40.0 plyr_1.8.9               crayon_1.5.2            
##  [85] abind_1.4-5              BiocIO_1.14.0            gridGraphics_0.5-1      
##  [88] locfit_1.5-9.9           graphlayouts_1.1.1       bit_4.0.5               
##  [91] fastmatch_1.1-4          codetools_0.2-20         textshaping_0.3.7       
##  [94] openssl_2.1.2            slam_0.1-50              GetoptLong_1.0.5        
##  [97] tm_0.7-13                mime_0.12                splines_4.4.0           
## [100] circlize_0.4.16          Rcpp_1.0.12              dbplyr_2.5.0            
## [103] HDO.db_0.99.1            knitr_1.46               blob_1.2.4              
## [106] utf8_1.2.4               clue_0.3-65              lme4_1.1-35.3           
## [109] fs_1.6.4                 checkmate_2.3.1          ggsignif_0.6.4          
## [112] ggplotify_0.1.2          Matrix_1.7-0             tzdb_0.4.0              
## [115] svglite_2.1.3            tweenr_2.0.3             pkgconfig_2.0.3         
## [118] tools_4.4.0              cachem_1.0.8             RSQLite_2.3.6           
## [121] viridisLite_0.4.2        DBI_1.2.2                fastmap_1.1.1           
## [124] rmarkdown_2.26           scales_1.3.0             grid_4.4.0              
## [127] Rsamtools_2.20.0         broom_1.0.5              officer_0.6.6           
## [130] patchwork_1.2.0          ggstats_0.6.0            carData_3.0-5           
## [133] rpart_4.1.23             farver_2.1.1             survminer_0.4.9         
## [136] tidygraph_1.3.1          scatterpie_0.2.2         mgcv_1.9-1              
## [139] yaml_2.3.8               rtracklayer_1.64.0       cli_3.6.2               
## [142] txdbmaker_1.0.0          table1_1.4.3             lifecycle_1.0.4         
## [145] askpass_1.2.0            backports_1.4.1          BiocParallel_1.38.0     
## [148] timechange_0.3.0         gtable_0.3.5             rjson_0.2.21            
## [151] umap_0.2.10.0            parallel_4.4.0           ape_5.8                 
## [154] jsonlite_1.8.8           mitml_0.4-5              bitops_1.0-7            
## [157] bit64_4.0.5              yulab.utils_0.1.4        zip_2.3.1               
## [160] geneLenDataBase_1.40.1   mice_3.16.0              highr_0.10              
## [163] GOSemSim_2.30.0          survMisc_0.5.6           lazyeval_0.2.2          
## [166] shiny_1.8.1.1            htmltools_0.5.8.1        enrichplot_1.24.0       
## [169] KMsurv_0.1-5             GO.db_3.19.1             rappdirs_0.3.3          
## [172] emR_0.6.3                tinytex_0.51             glue_1.7.0              
## [175] gfonts_0.2.0             httr2_1.0.1              XVector_0.44.0          
## [178] gdtools_0.3.7            RCurl_1.98-1.14          treeio_1.28.0           
## [181] gridExtra_2.3            boot_1.3-30              igraph_2.0.3            
## [184] R6_2.5.1                 km.ci_0.5-6              labeling_0.4.3          
## [187] GenomicFeatures_1.56.0   cluster_2.1.6            wordcloud_2.6           
## [190] aplot_0.2.2              nloptr_2.0.3             DelayedArray_0.30.0     
## [193] tidyselect_1.2.1         htmlTable_2.4.2          ggforce_0.4.2           
## [196] xml2_1.3.6               fontBitstreamVera_0.1.1  car_3.1-2               
## [199] AnnotationDbi_1.66.0     munsell_0.5.1            goseq_1.56.0            
## [202] fontquiver_0.2.1         data.table_1.15.4        htmlwidgets_1.6.4       
## [205] fgsea_1.30.0             ComplexHeatmap_2.20.0    RColorBrewer_1.1-3      
## [208] biomaRt_2.60.0           rlang_1.1.3              uuid_1.2-0              
## [211] fansi_1.0.6
```

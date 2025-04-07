## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")


##*********************************************************************************************************
## Load data

## load marker genes
markerGenes <- read_csv2("Data/marker_genes.csv")

## extract marker genes from csv
McMarkerGenes    <- markerGenes[markerGenes$TYPE == "Melanocytes",]$SYMBOL
DscMarkerGenes   <- markerGenes[markerGenes$TYPE == "DSCs",]$SYMBOL
FibroMarkerGenes <- markerGenes[markerGenes$TYPE == "Fibroblasts",]$SYMBOL
SkinMarkerGenes  <- c(McMarkerGenes, DscMarkerGenes, FibroMarkerGenes)

## load data
counts_20240924 <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20240924_lexogen_20241016.tsv', sep = '\t', header = TRUE)

## convert counts into matrix
counts <- column_to_rownames(counts_20240924[,-1], "gene_id") 

## assign coldata
colData <- read_xlsx("Data/metadata_v0.6_20241108.xlsx") %>% filter(run == "EKB_20240924") %>%
  mutate(time = factor(time, levels = c("6h", "24h")),
         irradiation = factor(irradiation, levels = c("control", "1x_UVB300", "3x_UVB300"), labels = c("Control", "1x 300J/m² UVB", "3x 300J/m² UVB")),
         condition = factor(condition, levels = c("DSC_control_6h", "DSC_1x_UVB300_6h",  "DSC_3x_UVB300_6h", "DSC_control_24h", "DSC_1x_UVB300_24h",  "DSC_3x_UVB300_24h",
                                                  "Melanocytes_control_6h", "Melanocytes_1x_UVB300_6h",  "Melanocytes_3x_UVB300_6h", "Melanocytes_control_24h", "Melanocytes_1x_UVB300_24h",  "Melanocytes_3x_UVB300_24h")))

## construct summarized experiment object
se <- SummarizedExperiment(assays= as.matrix(counts), colData=colData)

##*********************************************************************************************************
## Set thresholds and variable names and perform DE analysis

## always specify this filter criterion to take the desired difference into account in significance testing
pThres <- 0.01
lfcThres <- round(log2(1), 4)

## define contrast
dds <- DESeqDataSet(se, design = ~ donor + condition)

## run DESeq
dds <- DESeq(dds, parallel = TRUE)

## add metadata to DESeq object with attributes
attr(dds, "thresholds") <- list(
  pvalue = pThres, 
  lfc = lfcThres
)

##*********************************************************************************************************
## Data transformation

## Transform data
vsd_raw <- vst(dds, blind = FALSE)

## remove batch effect
vsd_adjusted <- vsd_raw
mat <- assay(vsd_adjusted)
mat <- limma::removeBatchEffect(mat, vsd_adjusted[["donor"]])
assay(vsd_adjusted) <- mat
mat <- NULL

## Construct annotation object
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(dds), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

## remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENSEMBL))

## store processed objects as list
processed_data <- list(dds = dds, 
                       se = se, 
                       vsd_unadjusted = vsd_raw, 
                       vsd_adj = vsd_adjusted, 
                       anno = anno,
                       SkinMarkerGenes = SkinMarkerGenes)

## processsed data
saveRDS(processed_data, "Data/processed_data.rds")



## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ashr)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")

## load data
counts <- read.table(file = 'TCGA/skcm_tcga_gdc_2025_01_16/data_mrna_seq_read_counts.txt', sep = '\t', header = TRUE)
patient <-  read.table(file = 'TCGA/skcm_tcga_gdc_2025_01_16/data_clinical_patient.txt', sep = '\t', header = TRUE)
sample <-  read.table(file = 'TCGA/skcm_tcga_gdc_2025_01_16/data_clinical_sample.txt', sep = '\t', header = TRUE)

## convert counts into matrix
counts <- column_to_rownames(counts, "Entrez_Gene_Id") 

## construct summarized experiment object
ind <- sample$SAMPLE_ID %in% str_replace_all(colnames(counts), "\\.", "-") # remove samples without count data
se <- SummarizedExperiment(assays= as.matrix(counts), colData=left_join(sample[ind,], patient)) # NOTE_ 2 patients have 2 samples in the dataset

##*********************************************************************************************************
## Set thresholds and variable names and perform DE analysis

## always specify this filter criterion to take the desired difference into account in significance testing
pThres <- 0.05
lfcThres <- round(log2(1), 4)

## define contrast
dds <- DESeqDataSet(se, design = ~ VITAL_STATUS)

## run DESeq
dds <- DESeq(dds, parallel = TRUE)

##*********************************************************************************************************
## Exploratory data analysis

##* Transform data
vsd_raw <- vst(dds, blind = FALSE)

pca_plot(data = vsd_raw, pcsToUse = 1:2, intgroup = c("VITAL_STATUS", "SEX"), 
         title = "PCA with unadjusted VST data", subtitle = "PC1 vs. PC2")

##* Construct annotation object
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(dds), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENTREZID")

## remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENTREZID))

res <- results(dds, alpha = pThres)

allRes <- cbind(ENTREZID = rownames(res), res)
allRes <- left_join(as.data.frame(allRes), anno)
allRes <- allRes[!str_detect(allRes$ENSEMBL, "\\."), ]
sigRes <- allRes[!is.na(allRes$padj) & allRes$padj < pThres,]

res[!is.na(res$padj) & res$padj < 0.05,]

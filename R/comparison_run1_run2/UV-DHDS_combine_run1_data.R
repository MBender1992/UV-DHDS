## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)

## load data
counts_EKB_20240313 <- read.table(file = 'Data/run1/dsc_dsc_ir_counts.tsv', sep = '\t', header = TRUE)[,1:11]
counts_EKB_20240313_MC <- read.csv("Data/run1/counts_first_sequencing_msc.csv")[,-1]
counts_TUD <- read.csv("Data/run1/dsc_msc_novogene.csv")[,-1]

tmp <- merge(counts_EKB_20240313, counts_EKB_20240313_MC, by.x = "Gene_id", by.y = "gene")
counts_run1 <- merge(tmp, counts_TUD, by.x = "Gene_id", by.y = "gene") 
colnames(counts_run1)[1] <- "gene_id"

write.table(counts_run1, "Data/UV-DHDS_count_matrix_EKB_20240313_and_TUD_Ishita.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

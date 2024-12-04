## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(readxl)

## load data
counts_lexogen <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20240924_lexogen_20241016.tsv', sep = '\t', header = TRUE)[,-c(1,22)]
counts_ishita  <- read.csv(file   = 'Data/UV-DHDS_count_matrix_EKB_20240924_Ishita_20241021.csv')

## rename colnames to match lexogen colnames
colnames(counts_ishita)[c(1:4, 34:38, 68)] <- c("gene_id", "M.a", "M.b", "M.c","D.a", "D.c", "D.d", "D.e", "D.f", "D.b")

## sort Ishita df to match lexogen df
counts_ishita <- counts_ishita[,match(colnames(counts_lexogen),colnames(counts_ishita))]

## combine counts
counts_combined <- merge(counts_lexogen, counts_ishita, by = "gene_id")
counts_combined <- column_to_rownames(counts_combined, "gene_id")

## assign coldata
colData <- read_xlsx("Data/metadata_v0.6_20241108.xlsx")[1:68,] %>% filter(ID != "D20") %>%
  mutate(time = factor(time, levels = c("6h", "24h")),
         irradiation = factor(irradiation, levels = c("control", "1x_UVB300", "3x_UVB300"), labels = c("Control", "1x 300J/m² UVB", "3x 300J/m² UVB")),
         condition = factor(condition, levels = c("DSC_control_6h", "DSC_1x_UVB300_6h",  "DSC_3x_UVB300_6h", "DSC_control_24h", "DSC_1x_UVB300_24h",  "DSC_3x_UVB300_24h",
                                                  "Melanocytes_control_6h", "Melanocytes_1x_UVB300_6h",  "Melanocytes_3x_UVB300_6h", "Melanocytes_control_24h", "Melanocytes_1x_UVB300_24h",  "Melanocytes_3x_UVB300_24h")))

colData <-  rbind(colData, colData)
colData$workflow <- c(rep("Lexogen", 67), rep("Ishita", 67))


## construct summarized experiment object
se <- SummarizedExperiment(assays= as.matrix(counts_combined), colData=colData)

## define contrast
dds_comparison <- DESeqDataSet(se, design = ~ donor + condition)

## run DESeq
dds_comparison <- DESeq(dds_comparison)

## Transform data
vsd_raw_comparison <- vst(dds_comparison, blind = FALSE)

svg("Results/comparison_workflow/PCA.svg", width=8, height=6)
pca_plot(data = vsd_raw_comparison, pcsToUse = 1:2, intgroup = c("workflow", "cell"), 
         title = "PCA with unadjusted VST data", subtitle = "PC1 vs. PC2")
dev.off()

##*********************************************************************************************************
## correlation analysis

se_Lexogen <- se[,se$workflow == "Lexogen"]
se_Ishita <- se[,se$workflow == "Ishita"]

## perform DESeq with lexogen workflow
dds_Lexogen <- DESeqDataSet(se_Lexogen, design = ~ donor + condition)
dds_Lexogen <- DESeq(dds_Lexogen)
res_Lexogen <- results(dds_Lexogen, lfcThreshold = lfcThres, alpha = pThres, 
                        contrast = contraster(dds_Lexogen,
                                              group1 = list(c("irradiation", "Control"), c("cell", "Melanocytes")),
                                              group2 = list(c("irradiation", "Control"), c("cell", "DSC"))))

## perform DESeq with ishitas workflow
dds_Ishita <- DESeqDataSet(se_Ishita, design = ~ donor + condition)
dds_Ishita <- DESeq(dds_Ishita)
res_Ishita <- results(dds_Ishita, lfcThreshold = lfcThres, alpha = pThres, 
                       contrast = contraster(dds_Ishita,
                                             group1 = list(c("irradiation", "Control"), c("cell", "Melanocytes")),
                                             group2 = list(c("irradiation", "Control"), c("cell", "DSC"))))

## plot results
ind <- which(!is.na(res_Lexogen$log2FoldChange) & !is.na(res_Ishita$log2FoldChange) & res_Lexogen$baseMean >= 10 & res_Ishita$baseMean >= 10)
plot(res_Lexogen$log2FoldChange[ind], res_Ishita$log2FoldChange[ind], xlab = "Log2 FC Lexogen", ylab = "Log2 FC Ishita", main = "Correlation of log2 FC for DSC vs Melanocytes")

## calculate correlation
cor(res_Lexogen$log2FoldChange[ind], res_Ishita$log2FoldChange[ind])


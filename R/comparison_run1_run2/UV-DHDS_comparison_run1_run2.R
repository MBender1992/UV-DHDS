## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(readxl)

## load data
counts_lexogen <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20240924_lexogen_20241016.tsv', sep = '\t', header = TRUE)[,-c(1,22)]
counts_20240313 <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20240313_and_TUD_Ishita.tsv', sep = '\t', header = TRUE)
counts_ishita <- read.csv(file = 'Data/UV-DHDS_count_matrix_EKB_20240924_Ishita_20241021.csv')
colnames(counts_ishita)[1] <- "gene_id"

## rename colnames to match lexogen colnames
colnames(counts_ishita)[c(1:4, 34:38, 68)] <- c("gene_id", "M.a", "M.b", "M.c","D.a", "D.c", "D.d", "D.e", "D.f", "D.b")

## sort Ishita df to match lexogen df
counts_ishita <- counts_ishita[,match(colnames(counts_lexogen),colnames(counts_ishita))]

## combine counts
counts_combined <- merge(counts_ishita, counts_20240313, by = "gene_id")
counts_combined <- column_to_rownames(counts_combined, "gene_id")

## assign coldata
colData <- read_xlsx("Data/metadata_v0.6_20241108.xlsx") %>% 
  filter(ID != "D20") %>%
  mutate(time = factor(time, levels = c("6h", "24h")),
         irradiation = factor(irradiation, levels = c("control", "1x_UVB300", "3x_UVB300"), labels = c("Control", "1x 300J/m² UVB", "3x 300J/m² UVB")),
         condition = factor(condition, levels = c("DSC_control_6h", "DSC_1x_UVB300_6h",  "DSC_3x_UVB300_6h", "DSC_control_24h", "DSC_1x_UVB300_24h",  "DSC_3x_UVB300_24h",
                                                  "Melanocytes_control_6h", "Melanocytes_1x_UVB300_6h",  "Melanocytes_3x_UVB300_6h", "Melanocytes_control_24h", "Melanocytes_1x_UVB300_24h",  "Melanocytes_3x_UVB300_24h"))) %>%
  mutate(ID = str_replace_all(.$ID, "-", "."))

## define index to include only data which was sequenced at EKB on the first run and the corresponding data of the new run for a comparison  between both libraries
ind <- which((colData$irradiation == "Control" & colData$time == "24h")| colData$condition == "DSC_1x_UVB300_24h")

## construct summarized experiment object
se <- SummarizedExperiment(assays= as.matrix(counts_combined[,ind]), colData=colData[ind,])

## define contrast
dds_comparison <- DESeqDataSet(se, design = ~ donor + condition)

## run DESeq
dds_comparison <- DESeq(dds_comparison)

## Transform data
vsd_raw_comparison <- vst(dds_comparison, blind = FALSE)

svg("Results/comparison_run1/PCA.svg", width=8, height=6)
pca_plot(data = vsd_raw_comparison, pcsToUse = 1:2, intgroup = c("run", "cell"), 
         title = "PCA with unadjusted VST data", subtitle = "PC1 vs. PC2")
dev.off()


## comparison of counts
ind <- which((colData$irradiation == "Control" & colData$time == "24h" & colData$run != "TUD" & colData$cell == "DSC" & colData$ID != "D.d"))
colData[ind,] %>% print(n="all")
counts_combined[,ind]

depth_lexogen <- rowSums(counts_combined[,ind][,1:5])
depth_illumina <- rowSums(counts_combined[,ind][,6:10])

depths <- data.frame(lexogen = depth_lexogen, illumina = depth_illumina)
depths$ratio <- log2((depths$illumina)/(depths$lexogen))

hist(depths$ratio[depths$ratio != 0], xlim = c(-20,20), xlab = "count ratio", breaks = 100, main = "Distribution of count ratio Illumina/Lexogen")

median(depths$ratio, na.rm = T)
mean(depths$ratio[!is.infinite(depths$ratio) & !is.nan(depths$ratio)])

## 
sum(depths$ratio < 0, na.rm = T)
sum(depths$ratio > 0, na.rm = T)

## 
sum(depths$lexogen == 0, na.rm = T)
sum(depths$illumina == 0, na.rm = T)

## 
sum(depths$lexogen > 1000, na.rm = T)
sum(depths$illumina > 1000, na.rm = T)


##*********************************************************************************************************
## compare fold changes of different sequencing approaches

se_20240313 <- se[,se$run == "EKB_20240313"]
se_lexogen <- se[,se$run == "EKB_20240924"]

## perform DESeq on old EKB data
dds_20240313 <- DESeqDataSet(se_20240313, design = ~ donor + condition)
dds_20240313 <- DESeq(dds_20240313)
res_20240313 <- results(dds_20240313, lfcThreshold = lfcThres, alpha = pThres, 
                        contrast = contraster(dds_20240313,
                                              group1 = list(c("condition", "DSC_1x_UVB300_24h")),
                                              group2 = list(c("condition", "DSC_control_24h"))))

## perform DESeq on new EKB data
dds_lexogen <- DESeqDataSet(se_lexogen, design = ~ donor + condition)
dds_lexogen <- DESeq(dds_lexogen)
res_lexogen <- results(dds_lexogen, lfcThreshold = lfcThres, alpha = pThres, 
                        contrast = contraster(dds_lexogen,
                                              group1 = list(c("condition", "DSC_1x_UVB300_24h")),
                                              group2 = list(c("condition", "DSC_control_24h"))))

## plot results
plot(res_20240313$log2FoldChange, res_lexogen$log2FoldChange, xlab = "Log2 FC Illumina", ylab = "Log2 FC Lexogen", main = "Correlation of log2 FC for irradiated vs. control in DSCs (24h)")
ind <- which(!is.na(res_20240313$log2FoldChange) & !is.na(res_lexogen$log2FoldChange))

## calculate correlation
cor(res_20240313$log2FoldChange[ind], res_lexogen$log2FoldChange[ind])


median(se_lexogen$rna_input_ng)

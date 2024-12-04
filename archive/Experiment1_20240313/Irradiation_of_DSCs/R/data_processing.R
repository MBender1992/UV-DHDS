## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************
## load packages
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(apeglm)
library(PoiClaClu)
library(ekbSeq)
library(pathfindR)

## load sample sheet
sampleSheet <- read_csv("Data/SampleSheet.csv")

## replace some values with \n for newlines
sampleSheet$Value2 <- str_replace_all(sampleSheet$Value2, "\\_", "\\_\n")
sampleSheet$Value2 <- str_replace(sampleSheet$Value2, "\\_\n", "\\_")
sampleSheet$Value3 <- str_replace(sampleSheet$Value3, "Stran", "Stran-\n")
sampleSheet$Value3 <- str_replace(sampleSheet$Value3, "Prep", "Prep-\n")
sampleSheet$Value4 <- str_replace(sampleSheet$Value4, "AUDI", "AUDI\n")
sampleSheet$Value4 <- str_replace(sampleSheet$Value4, "Adapter", "Adapter-\n")

## load marker genes
markerGenes <- read_csv2("../Melanocytes_vs_DSCs/Data/marker_genes.csv")

## extract marker genes from csv
MelanomaMarkerGenes <- markerGenes[str_detect(markerGenes$TYPE, "[Mm]elanoma"), ]

## load dataset from Darmstadt
datTUD <- read.table(file = 'Data/DSC_irradiated_counts_Darmstadt.tsv', sep = '\t', header = TRUE)
datEKB  <- datTUD[,c(1:11)]
datTUD <- datTUD[,c(1,12:16)]

## load dataset from buxtehude
colnames(datEKB) <- ifelse(colnames(datEKB) != "Gene_id", paste(colnames(datEKB), "_EKB", sep = ""), colnames(datEKB))

## rename columns from datasource 2
colnames(datTUD) <-  ifelse(colnames(datTUD) != "Gene_id", paste(colnames(datTUD), "_TUD", sep = ""), colnames(datTUD))

## join expression matrices
dat <- left_join(datEKB, datTUD, by = "Gene_id")

## assign count matrix
counts <- as.matrix(dat[,-1])
rownames(counts) <- NULL
colnames(counts) <- NULL

## define column data
colData <- data.frame(
  treatment = factor(c(rep("ctrl", 5), rep("trt", 5), rep("ctrl", 5))),
  donor = factor(rep(c("DSC_27_23", "DSC_28_23", "DSC_29_23", "DSC_30_23", "DSC_31_23"), 3)),
  run = c(rep("EKB", 10), rep("TUD", 5))
)

colData$sample <- paste(colData$donor, colData$treatment,sep = "_")

## construct summarized experiment object
se <- SummarizedExperiment(assays=list(counts = counts),
                           colData=colData)

## set rownames as ENSEMBLE identifieres
rownames(se) <- unlist(dat[,1])

## store results from DRAGEN analysis
DragenRes <- read_csv("Data/DE_results_DRAGEN.csv")
DragenRes <- head(DragenRes, 30)

##*********************************************************************************************************
## Set thresholds and variable names

## always specify this filter criterion to take the desired difference into account in significance testing
pThres <- 0.01
lfcThres <- round(log2(1), 4)
minimumCount <- 10 # Here we perform pre-filtering to keep only rows that have a count of at least 10 for a
# minimal number of samples. The count of 10 is a reasonable choice for bulk RNA-seq.
smallestGroupSize <- 5 # A recommendation for the minimal number of samples is to specify the smallest group
# size,e.g. here there are 5 samples in each group.
factor1 <- "treatment" # Variable of interests
factor2 <- "donor" # placeholder in case we want to adjust for a second factor. For more than 2 factors the workflow needs to be adjusted.

## define contrast
contrast <- ekbSeq::make_design_formula(fct1 = factor1, fct2 = factor2)
dds <- DESeqDataSet(se, design = contrast)

## collapse technical replicates
dds <- collapseReplicates(dds, dds$sample, dds$run)

## extract number of reads per sample
nread <- colSums(assay(dds))
names(nread) <- str_replace(names(nread), "\\_", "\n")

## run DESeq
dds <- DESeq(dds)

## remove duplicate entries
dds <- dds[!str_detect(rownames(counts(dds)), "\\."),]

##*********************************************************************************************************
## Exploratory data analysis

##******
##* Poisd

## calculate Poisson distances and define parameters for graphical representation
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- rownames(colData(se))
colnames(samplePoisDistMatrix) <- NULL
colorsPoisd <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

##******
##* Transform data

## variance stabilized transformation
vsd <- vst(dds, blind = FALSE)
## rlog transformation (regularized log transformation)
rld <- rlog(dds, blind = FALSE)
## estimate size factors
dds <- estimateSizeFactors(dds)

##******
##* Remove batch effects (in this case the donor effect) if applicable

if (!is.null(factor2)) {
  ## remove batch for vsd
  mat <- assay(vsd)
  mat <- limma::removeBatchEffect(mat, vsd[[factor2]])
  assay(vsd) <- mat

  ## remove batch for rlog
  mat <- assay(rld)
  mat <- limma::removeBatchEffect(mat, rld[[factor2]])
  assay(rld) <- mat
}

##*********************************************************************************************************
## Differential expression analysis

##******
##* Extract results and shrink lfcs
res <- results(dds, lfcThreshold = lfcThres, alpha = pThres)


## shrink LFCs
resShrunk <- lfcShrink(dds, res = res, coef = "treatment_trt_vs_ctrl", type = "apeglm")

##******
##* Annotate data and subset significant genes
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(res),
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
)

## remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENSEMBL))

## adding ENSEMBL gene ID as a column in significant differentially expression gene table.
allRes <- cbind(ENSEMBL = rownames(resShrunk), resShrunk)
allRes <- if (dim(anno)[1] == dim(allRes)[1]) {
  left_join(as.data.frame(allRes), anno)
} else {
  stop("Dimensions of annotation object and result object are different.")
}
allRes <- allRes[!str_detect(allRes$ENSEMBL, "\\."), ]

## subset significant results
sigRes <- subset(allRes, padj < pThres)


##******
## Extract genes with high dispersion
indDisp <- which(mcols(dds, use.names = TRUE)$dispersion >= 1 & rowSums(counts(dds) >= 10) >= smallestGroupSize)
rldHighDisp <- assay(rld)[indDisp, ]

highDispGenes <- allRes[allRes$ENSEMBL %in% rownames(rldHighDisp), ]$SYMBOL
highDispGenes <- highDispGenes[!is.na(highDispGenes)]

##******
## Setup pathfindR workflow

## define pathway input to use the appropriate gene list background
keep1 <- rowSums(counts(dds) >= 10) >= smallestGroupSize
keep2 <- !is.na(allRes$padj) & allRes$padj < pThres
keep  <- keep1 | keep2
allResPath <- allRes[keep, ] %>% as.data.frame()
rownames(allResPath) <- NULL
keep <- NULL

## load WikiPathway data
wpsConv <- readRDS("../wpsConv.rds")[[1]]
wpNames <- readRDS("../wpNames.rds")[[1]]
wpIds   <- readRDS("../wpIds.rds")[[1]]

## define pathfindR input and remove NAs
pathfindR_input <- allResPath %>%
  column_to_rownames("ENSEMBL") %>%
  dplyr::select(SYMBOL, log2FoldChange, padj) %>%
  na.omit()

## make analysis with WikiPathways
outputDf <- run_pathfindR(pathfindR_input,
                          gene_sets = "Custom",
                          custom_genes = wpsConv,
                          p_val_threshold = pThres,
                          custom_descriptions = wpNames)

## melanocytes are treated as cases throughout the analysis
cases <- colnames(dds)[str_detect(colnames(dds),"trt")]

## setup expression matrix as input for pathfindRs score_term function
expMat <-  assay(rld)[rownames(assay(rld)) %in% rownames(pathfindR_input),] %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(anno[,c("ENSEMBL", "SYMBOL")]) %>%
  dplyr::select(-ENSEMBL) %>%
  filter(!duplicated(SYMBOL) & !is.na(SYMBOL)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

## calculate Z-Score for pathway activation/repression
scoreMatrix <- score_terms(
  enrichment_table = outputDf,
  exp_mat = expMat,
  cases = cases,
  plot_hmap = FALSE)

## extract colnames for each factor level for later calculation
lvl1 <- colnames(dds)[colData(dds)[[factor1]] == "ctrl"]
lvl2 <- colnames(dds)[colData(dds)[[factor1]] == "trt"]

## convert Matrix to df and assign row means for each factor level
scoreMatrix <- as.data.frame(scoreMatrix)
scoreMatrix$meanCtrl <- rowMeans(scoreMatrix[, colnames(scoreMatrix) %in% lvl1])
scoreMatrix$meanTrt <- rowMeans(scoreMatrix[, colnames(scoreMatrix) %in% lvl2])

## assign upregulated and downregulated pathways as vector
indAct <- rownames(scoreMatrix)[scoreMatrix$meanTrt > scoreMatrix$meanCtrl]

## assign direction of pathway enrichment
outputDf$Direction <- ifelse(outputDf$ID %in% indAct, "Activated", "Repressed")

## define clusters in upregulated pathways
upClustered   <- cluster_enriched_terms(outputDf[outputDf$Direction == "Activated",], plot_dend = FALSE, plot_clusters_graph = FALSE)
## define clusters in downregulated pathways
downClustered <- cluster_enriched_terms(outputDf[outputDf$Direction == "Repressed",], plot_dend = FALSE, plot_clusters_graph = FALSE)

## bind clustered data
outputClustered <- rbind(upClustered, downClustered)

## extract representative pathways
representativeDf <- outputClustered[outputClustered$Status == "Representative",]
## calculate rank to include both fold change and pvalue as "importance" metric
representativeDf$rank <- sqrt(representativeDf$Fold_Enrichment * -log(representativeDf$lowest_p))

##******
## Save data if anything has been changed (default is false)
ekbSeq::save_res(
  save = FALSE,
  all.res = allRes,
  sig.res = sigRes,
  pathway.res = outputClustered,
  background.gene.list = allResPath,
  results.name = "Results/dsc_irradiation",
  rdata.name = "DSC_NGS_UV",
  rm.batch = factor2
)




# testString <- resEnrich$Description
# testString <- str_replace_all(testString, "cell ", "cell_")
# testString <- str_replace_all(testString, "positive ", " ")
# testString <- str_replace_all(testString, "negative ", " ")
# testString <- str_replace_all(testString, "regulation ", "regulation_")
# testString <- str_replace_all(testString, " process", "_process")
# testString <- str_replace_all(testString, "response ", "response_")
# testString <- str_replace_all(testString, "of ", "of_")
# testString <- str_replace_all(testString, "to ", "to_")
# testString <- str_replace_all(testString, "via ", "via_")
#
# gobpDownFreq <- data.frame(table(unlist(strsplit(unlist(testString), " "))))
# wordcloud2(gobpDownFreq, size = 0.5)


##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************
## load packages
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(apeglm)
library(PoiClaClu)
library(corto)
library(pathfindR)
library(goseq)
library(rrvgo)
library(ekbSeq)

## define custom functions (only used in this script)

## function to calculate cross table
calc_cross_tab <- function(type, lfc = "greater"){
  if(lfc == "greater"){
    n <- sum(allRes[allRes$SYMBOL %in% MelanomaMarkerGenes[MelanomaMarkerGenes$TYPE == type,]$SYMBOL, ]$log2FoldChange>0)
    p <- mean(allRes[allRes$SYMBOL %in% MelanomaMarkerGenes[MelanomaMarkerGenes$TYPE == type,]$SYMBOL, ]$log2FoldChange>0)
  } else if(lfc == "lesser") {
    n <- sum(allRes[allRes$SYMBOL %in% MelanomaMarkerGenes[MelanomaMarkerGenes$TYPE == type,]$SYMBOL, ]$log2FoldChange<0)
    p <- mean(allRes[allRes$SYMBOL %in% MelanomaMarkerGenes[MelanomaMarkerGenes$TYPE == type,]$SYMBOL, ]$log2FoldChange<0)
  } else {
    stop("Please specify whether lfc greater or lesser than 0 should be compared.")
  }
  paste(n, " (",round(p*100,1), "%)", sep = "")
}


##*********************************************************************************************************
##*load data

## load sample sheet
sampleSheet <- read_csv("Data/SampleSheet.csv")
sampleSheet <- format_sample_sheet(sampleSheet)

## load marker genes
markerGenes <- read_csv2("Data/marker_genes.csv")

## extract marker genes from csv
McMarkerGenes    <- markerGenes[markerGenes$TYPE == "Melanocytes",]$SYMBOL
DscMarkerGenes   <- markerGenes[markerGenes$TYPE == "DSCs",]$SYMBOL
FibroMarkerGenes <- markerGenes[markerGenes$TYPE == "Fibroblasts",]$SYMBOL
SkinMarkerGenes  <- c(McMarkerGenes, DscMarkerGenes, FibroMarkerGenes)
MelanomaMarkerGenes <- markerGenes[str_detect(markerGenes$TYPE, "[Mm]elanoma"),]

## load dataset from Darmstadt
datTUD <- read.table(file = 'Data/DSC_vs_MC_genes_counts_Ishita.tsv', sep = '\t', header = TRUE)

## load dataset from buxtehude
datEKB  <- read_csv("Data/DSC_vs_MC_genes_counts.csv")
datEKB[,1] <- str_remove(unlist(datEKB[,1]), "\\..+$")
colnames(datEKB) <- ifelse(colnames(datTUD) != "Gene_id", paste(colnames(datTUD), "_EKB", sep = ""), colnames(datTUD))

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
  cell_type = factor(rep(c(rep("DSC", 5), rep("Melanocytes", 5)),2)),
  sample = rep(c(
    "DSC_27_23", "DSC_28_23", "DSC_29_23", "DSC_30_23", "DSC_31_23",
    "MC_27_23", "MC_28_23", "MC_29_23", "MC_30_23", "MC_31_23"
  ),2), 
  run = c(rep("EKB", 10), rep("TUD", 10))
)

## construct summarized experiment object
se <- SummarizedExperiment(assays=list(counts = counts),
                     colData=colData)

## collapse technical replicates
seColl <- collapseReplicates(se, se$sample, se$run)

## set rownames as ENSEMBLE identifieres
rownames(seColl) <- unlist(dat[,1])

## store results from DRAGEN analysis
DragenRes <- read_csv("Data/DE_results_DRAGEN.csv")
DragenRes <- head(DragenRes, 30)

##*********************************************************************************************************
## Set thresholds and variable names

## always specify this filter criterion to take the desired difference into account in significance testing
pThres            <- 0.05
lfcThres          <- round(log2(1.5),4)
minimumCount      <- 10 # Here we perform pre-filtering to keep only rows that have a count of at least 10 for a 
# minimal number of samples. The count of 10 is a reasonable choice for bulk RNA-seq.
smallestGroupSize <- 5  # A recommendation for the minimal number of samples is to specify the smallest group     
# size,e.g. here there are 5 samples in each group.
factor1 <- "cell_type" # Variable of interests
factor2 <- NULL # placeholder in case we want to adjust for a second factor. For more than 2 factors the workflow needs to be adjusted.

## define contrast
contrast <- ekbSeq::make_design_formula(fct1 = factor1, fct2 = factor2)
dds <- DESeqDataSet(seColl, design = contrast)

## extract number of reads per sample
nread <- colSums(assay(dds))
names(nread) <- str_replace(names(nread), "\\_", "\n")

##*********************************************************************************************************
## Exploratory data analysis

##****** 
##*Poisd

## calculate Poisson distances and define parameters for graphical representation
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- rownames(colData(se))
colnames(samplePoisDistMatrix) <- NULL
colorsPoisd <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)

##******
##*Transform data

## variance stabilized transformation
vsd <- vst(dds, blind = FALSE)
## rlog transformation (regularized log transformation)
rld <- rlog(dds, blind = FALSE)
## estimate size factors
dds <- estimateSizeFactors(dds)


##*********************************************************************************************************
## Differential expression analysis

## run DESeq
dds <- DESeq(dds)

##******
##*Extract results and shrink lfcs
res <- results(dds, lfcThreshold = lfcThres, alpha = pThres)

## shrink LFCs
resShrunk <- lfcShrink(dds, res = res, coef="cell_type_Melanocytes_vs_DSC", type="apeglm")

##******
##*Annotate data and subset significant genes
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(res), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

## remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENSEMBL))

## adding ENSEMBL gene ID as a column in significant differentially expression gene table.
allRes <- cbind(ENSEMBL = rownames(resShrunk), resShrunk)
allRes <- if(dim(anno)[1] == dim(allRes)[1]){
  left_join(as.data.frame(allRes), anno)
} else {
  stop("Dimensions of annotation object and result object are different.")
}

## subset significant results
sigRes <- subset(allRes, padj < pThres)

##******
## Extract genes with high dispersion
highDispGenes <- high_dispersion_genes(dds.obj = dds, norm.counts = rld, res = allRes, disp.thres = 1)

quantile(S4Vectors::mcols(dds, use.names = TRUE)$dispersion[arg2], na.rm =T, probs = seq(0,1,0.05))

indDisp <- which(arg1 & arg2)
##******
## define corto object and master regulators
load(system.file("extdata","centroids.rda", package="corto", mustWork=TRUE))

## define expression matrix object which serves as input to the corto function
expMat <- assay(rld)
expMat <- expMat %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>%
  left_join(allRes[,c("ENSEMBL", "SYMBOL")]) %>% dplyr::select(-ENSEMBL) %>%
  filter(!duplicated(SYMBOL) & !is.na(SYMBOL)) %>% column_to_rownames("SYMBOL") %>%
  as.matrix()

## Master regulators without known melanocyte factors
regulon <- corto(expMat,centroids=centroids,p = 0.0001, nthreads=2,nbootstraps=10,verbose=TRUE)
masterRegNew <- mra(expMat[,6:10], expMat[,1:5], regulon = regulon)

## Master regulators with known melanocyte factors
centroids <- c(centroids, SkinMarkerGenes)
centroids[!duplicated(centroids)]
regulon <- corto(expMat,centroids=centroids,p = 0.0001, nthreads=2,nbootstraps=10,verbose=TRUE)
masterRegMarkers <- mra(expMat[,6:10], expMat[,1:5], regulon = regulon)


##*********************************************************************************************************
## calculate cross table for comparison of melanocyte and melanoma gene expression

## set up cross table
crossTab <- data.frame(matrix(nrow = 2, ncol = 2))
rownames(crossTab) <- c("Melanocytes_up", "Melanocytes_down")
colnames(crossTab) <- c("Melanoma_up", "Melanoma_down")

## populate cross table
crossTab[1,1] <- calc_cross_tab(type = "Melanoma_up", lfc = "greater")
crossTab[2,2] <- calc_cross_tab(type = "Melanoma_down", lfc = "lesser")
crossTab[1,2] <- calc_cross_tab(type = "Melanoma_down", lfc = "greater")
crossTab[2,1] <- calc_cross_tab(type = "Melanoma_up", lfc = "lesser")

##*********************************************************************************************************
## setup pathfindR input

## define pathway input to use the appropriate gene list background
keep1 <- rowSums(counts(dds) >= 10) >= smallestGroupSize
keep2 <- !is.na(allRes$padj) & allRes$padj < 0.05
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
                          custom_descriptions = wpNames)

## melanocytes are treated as cases throughout the analysis
cases <- colnames(dds)[str_detect(colnames(dds),"MC")]

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
lvl1 <- colnames(dds)[colData(dds)[[factor1]] == "DSC"]
lvl2 <- colnames(dds)[colData(dds)[[factor1]] == "Melanocytes"]

## convert Matrix to df and assign row means for each factor level
scoreMatrix <- as.data.frame(scoreMatrix)
scoreMatrix$meanDSC <- rowMeans(scoreMatrix[, colnames(scoreMatrix) %in% lvl1])
scoreMatrix$meanMC <- rowMeans(scoreMatrix[, colnames(scoreMatrix) %in% lvl2])

## assign upregulated and downregulated pathways as vector
indAct <- rownames(scoreMatrix)[scoreMatrix$meanMC > scoreMatrix$meanDSC]

## assign direction of pathway enrichment
outputDf$Direction <- ifelse(outputDf$ID %in% indAct, "Activated", "Repressed")

## define clusters in upregulated pathways
upClustered   <- cluster_enriched_terms(outputDf[outputDf$Direction == "Activated",], plot_dend = FALSE, plot_clusters_graph = FALSE)
largeUpClusters <- as.numeric(names(table(upClustered$Cluster)[table(upClustered$Cluster)>3]))

## define clusters in downregulated pathways
downClustered <- cluster_enriched_terms(outputDf[outputDf$Direction == "Repressed",], plot_dend = FALSE, plot_clusters_graph = FALSE)
largeDownClusters <- as.numeric(names(table(downClustered$Cluster)[table(downClustered$Cluster)>3]))

## bind clustered data
outputClustered <- rbind(upClustered, downClustered)

## extract representative pathways
representativeDf <- outputClustered[outputClustered$Status == "Representative",]

## calculate rank to include both fold change and pvalue as "importance" metric
representativeDf$rank <- sqrt(representativeDf$Fold_Enrichment * -log(representativeDf$lowest_p))

##*********************************************************************************************************
## calculate enriched pathways with goseq for high dispersion genes (indicating donor differneces)

## define gene List based on high dispersion genes
goseqInput <- allResPath[!is.na(allResPath$padj) | allResPath$SYMBOL %in% highDispGenes,]
geneList <- as.integer(goseqInput$SYMBOL %in% highDispGenes)
names(geneList) <- goseqInput$ENSEMBL

## extract reduced go terms from goseq analysis
reducedTerms <- ekbSeq::goseq_to_revigo(gene.list = geneList, genome = "hg19")

##******
## Save data if anything has been changed (default is false)
ekbSeq::save_res(
  save = FALSE, 
  all.res = allRes,
  sig.res = sigRes, 
  pathway.res = outputClustered,
  background.gene.list = allResPath,
  results.name = "Results/melanocytes_vs_dscs", 
  rdata.name = "DSC_NGS", 
  rm.batch = factor2
)


# 
# library(biomaRt)
# 
# # the mart, so-to-say the "universe", here H.Sapiens
# mymart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # list everything that this mart offers
# attrs  <- listAttributes(mart = mymart)
# 
# #/ two example genes:
# wanted <- anno[anno$SYMBOL %in% highDispGenes,]$ENSEMBL
# 
# # retrieve the data for these genes:
# genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
#                         mart = mymart, 
#                         filters = "ensembl_gene_id", 
#                         values = wanted)
# 
# table(genes$chromosome_name)



# Darmstadt ca 6000 DE genes
# Buxtehude ca. 5700 DE genes
# Overlap ca 5200 DE genes
# Die non-overlapping genes zum Großteil wenig signifikant



## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(UCell)

##******
## load data 

## load data from GSE178514; BM Mesenchymal stem cells
GSE178514 <- read.delim("mesenchymal_stem_cell_data/data/GSE178514_AllGeneFPKM_mRNA.txt", sep = "\t", header = TRUE) %>% as.data.frame() 
GSE178514$GeneID <- anno$SYMBOL[match(GSE178514$GeneID, anno$ENTREZID)]
GSE178514 <- GSE178514 %>% filter(!is.na(GeneID) & !duplicated(GeneID)) %>%
  column_to_rownames("GeneID")

## extract markers
GSE178514_markers <- rowMeans(GSE178514) %>%
  sort(decreasing = TRUE) %>%
  .[1:1000]

## load data from GSE95072
GSE95072 <- read.delim("mesenchymal_stem_cell_data/data/GSE95072_368-Lewis_DESeq2_NormCounts.txt", sep = "\t", header = TRUE) %>% as.data.frame() 
GSE95072 <- GSE95072 %>%
  column_to_rownames("X368.Lewis_DESeq2_NormCounts")

GSE95072_markers <- rowMeans(GSE95072) %>%
  sort(decreasing = TRUE) %>%
  .[1:1000]

##******
## specify gene sets based on literature research marker genes
gene.sets <- list(MSC_signature = c("THY1", "ENG", "NT5E", "CDH5","ALCAM"), # CD90 Thy1; CD105 ENG; CD73 NT5E
                  GSE178514_signature = names(GSE178514_markers),
                  GSE95072_signature = names(GSE95072_markers),
                  DSC_signature = c("NGFR", "S100A4", "MSX1", "CDH2"), 
                  Fibro_signature = c("VIM", "ELN", "FN1", "LAMA5", "LAMA2"),
                  Melanocyte_signature = McMarkerGenes)

## write signatures as csv
plyr::ldply(gene.sets, rbind) %>% 
  column_to_rownames(".id") %>%
  t() %>% 
  data.frame() %>%
  mutate_if(is.character, ~replace(., is.na(.), "")) %>%
  write.csv("results/gene_signatures.csv")

##******
## score GSE178514
mat <- GSE178514 %>% as.matrix()
GSE178514_scores <- ScoreSignatures_UCell(mat, features=gene.sets) %>%
  as.data.frame() %>%
  mutate(Passage = str_remove(rownames(.), "\\d{3}\\.FPKM")) 
write.csv(GSE178514_scores, "results/GSE178514_signature_scores.csv")

## score GSE95072
mat <- GSE95072 %>% as.matrix()
GSE95072_scores <- ScoreSignatures_UCell(mat, features=gene.sets) %>%
  as.data.frame() %>%
  mutate(Passage = str_remove(rownames(.), "\\d{3}\\.FPKM")) 
write.csv(GSE95072_scores, "results/GSE95072_signature_scores.csv")


##******
## extract normalized counts from our dataset and calculate scores based on the gene signatures
mat <- assay(vsd_adjusted)
rownames(mat) <- anno$SYMBOL[match(rownames(mat), anno$ENSEMBL)]
scores <- ScoreSignatures_UCell(mat, features=gene.sets)
metadata <- colData(vsd_adjusted)
UVDHDS_scores <- cbind(metadata, scores) %>%
  as.data.frame() %>%
  select(ID, sample, donor, MSC_signature_UCell, GSE178514_signature_UCell, GSE95072_signature_UCell,  DSC_signature_UCell, Fibro_signature_UCell, Melanocyte_signature_UCell) 
write.csv(UVDHDS_scores, "results/UVDHDS_signature_scores.csv")
  
p <-UVDHDS_scores %>%
  mutate(cell_type = str_extract(.$sample, "[^_]+")) %>%
  gather("signature", "score", -c(ID, sample, donor, cell_type)) %>%
  ggplot(aes(signature, score, color = donor)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(0.2)) +
  facet_wrap(~cell_type, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)) +
  ggsci::scale_color_npg()

ggsave("results/UVDHDS_scores.png", width = 300, height = 150, unit = "mm")


library(pathfindR)
library(clusterProfiler)

## function to plot networks
plot_network <- function(enrichedData, gene.list, search.string = NULL, showCategory = 5, widths = c(0.35,0.65)){
  cnp <- cnetplot(enrichedData, color.params = list(foldChange = gene.list), 
                  showCategory = showCategory,
                  cex.params = list(gene_label = 0.5, category_label = 0.5)) +
    ggplot2::scale_color_gradientn(name = "Log2 fold change", colours = c("#053061", "#053061", "#2166ac","#4393c3", "#4393c3",
                                                                          "#92c5de", "white", "#f4a582", "#d6604d", "#d6604d",
                                                                          "#b2182b","#67001f" ,"#67001f"),
                                   values = c(0, 0.075,0.2,0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.8, 0.925, 1), 
                                   limits = c(-5, 5)) +
   theme(plot.tag = element_text(face = "bold", size = 20))
  print(cnp)
}

## function to save GO terms as dotplot and network
go_custom <- function(direction, ont = "BP"){
  if(direction == "greater"){
    tmp <- sigRes[sigRes$log2FoldChange > 0,]
    path <- paste("Results/Exploratory/upregulated_GO",ont, sep = "")
  } else if(direction == "lesser"){
    tmp <- sigRes[sigRes$log2FoldChange < 0,]
    path <- paste("Results/Exploratory/downregulated_GO",ont, sep = "")
  } else{
    stop("Please specify direction as either greater or lesser to analyze up- or downregulated genes, respectively.")
  }
  geneList <- unlist(tmp$log2FoldChange)
  names(geneList) <- tmp$ENTREZID
  resEnrich <- enrichGO(names(geneList), OrgDb = "org.Hs.eg.db", ont = ont, readable = TRUE)
  ego2 <- simplify(resEnrich, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  png(paste(path, "_dotplot.png",sep = ""), units="in", width=9, height=12, res=600)
  print(dotplot(ego2, showCategory = 30))
  dev.off()
  
  png(paste(path, "_network.png",sep = ""), units="in", width=12, height=9, res=600)
  plot_network(ego2, gene.list = geneList)
  dev.off()
}


rdata <- "DSC_NGS_UV_adjusted.RData"
load(file = rdata)

largeUpClusters <- as.numeric(names(table(upClustered$Cluster)[table(upClustered$Cluster)>3]))
png("Results/upClustered_wikiPW.png", units="in", width=9, height=13, res=600)
enrichment_chart(upClustered[upClustered$Cluster %in% largeUpClusters,], plot_by_cluster = TRUE)
dev.off()

largeDownClusters <- as.numeric(names(table(downClustered$Cluster)[table(downClustered$Cluster)>3]))
png("Results/downClustered_wikiPW.png", units="in", width=9, height=13, res=600)
enrichment_chart(downClustered[downClustered$Cluster %in% largeDownClusters,], plot_by_cluster = TRUE)
dev.off()


## GOBP upregulated
go_custom(direction = "greater", ont = "BP")
## GOBP downregulated
go_custom(direction = "lesser", ont = "BP")

## GOBP upregulated
go_custom(direction = "greater", ont = "MF")
## GOBP downregulated
go_custom(direction = "lesser", ont = "MF")

## GOBP upregulated
go_custom(direction = "greater", ont = "CC")
## GOBP downregulated
go_custom(direction = "lesser", ont = "CC")
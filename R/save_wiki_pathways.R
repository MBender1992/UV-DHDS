# ################# NOT RUN
library(qusage)
# ## read gmt file containing wikipathways
wps <- read.gmt("../wikipathways-20240510-gmt-Homo_sapiens.gmt")
wpIds <- str_extract(names(wps), "WP\\d+")

## convert ENTREZ Ids to gene symbols
wpsConv <- lapply(1:length(wps), function(i){
  AnnotationDbi::select(org.Hs.eg.db, wps[[i]],
                        columns="SYMBOL",
                        keytype="ENTREZID")$SYMBOL
})


remove <- c("%WikiPathways_20240510%", "%Homo sapiens", "WP\\d+")
wpNames <- str_remove_all(names(wps), paste(remove, collapse = "|"))
names(wpNames) <- wpIds
names(wpsConv) <- wpIds

## save R objects
obj_save(wpsConv, wpNames, wpIds, folder = "..")









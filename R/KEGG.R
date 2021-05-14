if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
if (!require(ggplot2))BiocManager::install("ggplot2")
if (!require(clusterProfiler))BiocManager::install("clusterProfiler")
library("clusterProfiler")


#' Title
#'
#' @param genelist a list of gene
#' @param genetype a string
#' @param KEGGfile Default:KEGG.info
#'
#' @return
#' @export
#'
#' @examples
KEGG_enrich <- function(genelist, genetype, KEGGfile = "KEGG.info") {
  kegglist <- read.delim(file = KEGGfile, sep = "\t", header = FALSE)
  gene <- genelist
  gene1 <- t(gene) # 获取需要分析的基因列表

  # kegg富集
  kegg2gene <- data.frame(kegglist$V1, kegglist$V3) # keggid,Geneid
  kegg2name <- data.frame(kegglist$V1, kegglist$V2) # keggid,kegg描述
  gene1 <- t(gene$V1) # 获取需要分析的基因列表
  spo_KEGG <- enricher(gene1, TERM2GENE = kegg2gene, TERM2NAME = kegg2name)
  res <- length(spo_KEGG$Count)
  print(res)
  # print(spo_KEGG)

  # 判断是否有富集到结果，有，则出图。
  if (res > 1) {
    p2 <- dotplot(spo_KEGG)
    p2 # KEGG出图
    ggsave(p1, file = paste(genetype, ".KEGG.tiff", sep = ""))
    ggsave(p1, file = paste(genetype, ".KEGG.pdf", sep = ""))
    write.csv(spo_KEGG, file = paste(genetype, "_KEGG.csv", sep = ""))
  }
  return(spo_KEGG)
}

if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
if (!require(ggplot2))BiocManager::install("ggplot2")
if (!require(clusterProfiler))BiocManager::install("clusterProfiler")
library("clusterProfiler")
#' Title
#'
#' @param genefile
#' @param genetype
#' @param TFfile
#'
#' @return
#' @export
#'
#' @examples
TF_enrich <- function(genelist, genetype, TFfile = "TF.info") {
  #gene <- read.delim(file = genefile, sep = "\t", header = FALSE)
  gene <- genelist
  tflist <- read.delim(file = TFfile, sep = "\t", header = FALSE)

  spo2gene <- data.frame(tflist$V1, tflist$V3) # TF,Geneid
  spo2name <- data.frame(tflist$V1, tflist$V2) # TF,TF描述
  gene1 <- t(gene$V1) # 获取需要分析的基因列表
  # TF富集分析
  spo_TF <- enricher(gene1, TERM2GENE = spo2gene, TERM2NAME = spo2name)
  res <- length(spo_TF$Count)
  print(res)
  # 判断是否有富集结果
  if (res > 1) {
    p1 <- dotplot(spo_TF, showCategory = length(spo_TF$Count)) # showCategory设置出图，显示的行数
    p1 # TF出图
    cnetplot(spo_TF, circular = TRUE, colorEdge = TRUE)
    ggsave("ST1_TF_circle.pdf")
    ggsave(p1, file = paste(genetype, ".TF.tiff", sep = ""))
    ggsave(p1, file = paste(genetype, ".TF.pdf", sep = ""))
    write.csv(spo_TF, file = paste(genetype, "_TF.csv", sep = ""))
  }
  return(spo_TF)
}

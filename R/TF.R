#' TFenrich: in fact it can enrich everything,just prepare the right database.
#' @import clusterProfiler
#' @import ggplot2
#' @param genelist:a gene list vector.
#' @param genetype:a string.
#' @param TFfile: just like `GO.info`.
#'
#' @return return a S3 enrich result
#' @export TFenrich
#'
#' @examples
#' genelist <- read.delim("ST1.txt")
#' TFenrich(genelist,"TF1","data/TF.info")
TFenrich <- function(genelist, genetype, TFfile = "TF.info") {
  # gene <- read.delim(file = genefile, sep = "\t", header = FALSE)
  gene <- genelist
  tflist <- read.delim(file = TFfile, sep = "\t", header = FALSE)

  spo2gene <- data.frame(tflist$V1, tflist$V3) # TF,Geneid
  spo2name <- data.frame(tflist$V1, tflist$V2) # TF,TF描述
  gene1 <- t(gene) # 获取需要分析的基因列表
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

#' KEGGenrich: for enrich kegg
#'
#' @param genelist a list of gene
#' @param genetype a string
#' @param KEGGfile Default:KEGG.info
#'
#' @return KEGG the enrich results
#' @export KEGGenrich
#'
#' @examples
#' genelist <- read.delim("data/ST1.txt")
#' KEGGenrich(genelist,"ST1","data/KEGG.info")
#'
KEGGenrich <- function(genelist, genetype, KEGGfile = "KEGG.info") {
  kegglist <- read.delim(file = KEGGfile, sep = "\t", header = FALSE)
  ###删除KEGG的描述字符中，字符长度超过60的，逗号分割的第一个逗号后的字符
  term <- c("term")
  for (i in 1:length(kegglist[,2])){
    if(nchar(kegglist[,2][i])>60){
      term <- append(term,kegglist[,2][i])
      #print(kegglist[,2][i])
      kegglist[,2][i] <- strsplit(kegglist[,2][i],split= ",")[[1]][1]
    }
  }
  unique(term)


  gene <- genelist
  # kegg富集
  kegg2gene <- data.frame(kegglist[,1], kegglist[,3]) # keggid,Geneid
  kegg2name <- data.frame(kegglist[,1], kegglist[,2]) # keggid,kegg描述
  gene1 <- t(gene) # 获取需要分析的基因列表
  spo_KEGG <- enricher(gene1, TERM2GENE = kegg2gene, TERM2NAME = kegg2name)
  res <- length(spo_KEGG$Count)
  print(res)
  # print(spo_KEGG)

  # 判断是否有富集到结果，有，则出图。
  if (res > 1) {
    p1 <- dotplot(spo_KEGG)
    print(p1) # KEGG出图
    ggsave(p1, file = paste(genetype, ".KEGG.tiff", sep = ""))
    ggsave(p1, file = paste(genetype, ".KEGG.pdf", sep = ""))
    write.csv(spo_KEGG, file = paste(genetype, "_KEGG.csv", sep = ""))
  }
  return(spo_KEGG)
}
